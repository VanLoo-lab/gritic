import bz2
import pickle
import warnings
from collections.abc import Mapping
from numbers import Integral, Real
from pathlib import Path

import numpy as np
import pandas as pd

from gritic.dataloader import validate_pascal_snake_case_columns


NON_PARSIMONY_PENALTY_COEFFICIENT = 2.7

ROUTE_KEY_COLUMNS = ['Sample_ID', 'Segment_ID', 'Route']
SEGMENT_KEY_COLUMNS = ['Sample_ID', 'Segment_ID']
SEGMENT_METADATA_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Chromosome',
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
    'Total_CN',
    'N_Mutations',
    'Mutation_Rate',
    'WGD_Status',
]
ROUTE_TABLE_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Chromosome',
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
    'Total_CN',
    'N_Mutations',
    'Mutation_Rate',
    'Route',
    'Probability',
    'Average_N_Events',
    'Average_Pre_WGD_Losses',
    'Average_Post_WGD_Losses',
    'Time',
    'Density',
    'WGD_Status',
    'WGD_Timing',
    'WGD_Timing_CI_Low',
    'WGD_Timing_CI_High',
]
GAIN_TIMING_TABLE_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Route',
    'Node',
    'Node_Phasing',
    'Timing',
    'Timing_CI_Low',
    'Timing_CI_High',
]
GAIN_DRAW_COLUMNS = SEGMENT_METADATA_COLUMNS + [
    'Posterior_Sample_Index',
    'Route',
    'Node',
    'Node_Phasing',
    'Gain_Timing',
    'WGD_Timing',
    'Gain_Index',
]
ROUTE_DRAW_COLUMNS = SEGMENT_METADATA_COLUMNS + [
    'Posterior_Sample_Index',
    'Route',
    'WGD_Timing',
]
SUMMARY_VALUE_COLUMNS = [
    'Gain_Index',
    'Proportion',
    'Timing_Median',
    'Timing_Low_CI',
    'Timing_High_CI',
    'Pre_WGD_Probability',
    'Post_WGD_Probability',
    'WGD_Timing_Median',
    'WGD_Timing_Low_CI',
    'WGD_Timing_High_CI',
]
SUMMARY_COLUMNS = SEGMENT_METADATA_COLUMNS + SUMMARY_VALUE_COLUMNS

_SEGMENT_CONSTANT_COLUMNS = SEGMENT_METADATA_COLUMNS + [
    'WGD_Timing',
    'WGD_Timing_CI_Low',
    'WGD_Timing_CI_High',
]


def _missing_columns(table, required_columns):
    return [column for column in required_columns if column not in table.columns]


def _require_columns(table, required_columns, table_name):
    missing_columns = _missing_columns(table, required_columns)
    if missing_columns:
        raise ValueError(
            f"The {table_name} is missing columns: "
            f"{', '.join(missing_columns)}"
        )


def _invalid_identifier(value):
    return pd.isna(value) or not str(value).strip()


def _validate_identifiers(table, columns, table_name):
    for column in columns:
        invalid = table[column].map(_invalid_identifier)
        if invalid.any():
            raise ValueError(
                f'{column} must be present for every {table_name} row'
            )


def _format_keys(keys):
    return ', '.join(':'.join(map(str, key)) for key in keys)


def _validate_route_rows(route_table):
    _require_columns(route_table, ROUTE_KEY_COLUMNS, 'route table')
    _validate_identifiers(route_table, ROUTE_KEY_COLUMNS, 'route table')

    duplicate_routes = route_table.duplicated(
        subset=ROUTE_KEY_COLUMNS,
        keep=False,
    )
    if duplicate_routes.any():
        keys = (
            route_table.loc[duplicate_routes, ROUTE_KEY_COLUMNS]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )
        raise ValueError(
            'The route table must contain exactly one row per route; '
            f'duplicate routes: {_format_keys(keys)}'
        )


def apply_non_parsimony_penalty(route_table):
    """Return a route table with one post-hoc penalty per route."""
    required_columns = ROUTE_KEY_COLUMNS + [
        'Average_N_Events',
        'Probability',
    ]
    _require_columns(route_table, required_columns, 'route table')
    _validate_route_rows(route_table)

    penalized_table = route_table.copy()
    penalized_table['Probability'] = pd.to_numeric(
        penalized_table['Probability'],
        errors='coerce',
    )
    penalized_table['Average_N_Events'] = pd.to_numeric(
        penalized_table['Average_N_Events'],
        errors='coerce',
    )

    if (penalized_table['Probability'] < 0).any():
        raise ValueError('Route probabilities must be non-negative')
    if (penalized_table['Average_N_Events'] < 0).any():
        raise ValueError('Average_N_Events values must be non-negative')

    for _, segment_indexes in penalized_table.groupby(
        SEGMENT_KEY_COLUMNS,
        sort=False,
    ).groups.items():
        segment_indexes = list(segment_indexes)
        probabilities = penalized_table.loc[
            segment_indexes,
            'Probability',
        ].to_numpy(dtype=float)
        average_events = penalized_table.loc[
            segment_indexes,
            'Average_N_Events',
        ].to_numpy(dtype=float)

        valid_segment = (
            np.isfinite(probabilities).all()
            and np.isfinite(average_events).all()
            and probabilities.sum() > 0
        )
        if not valid_segment:
            penalized_table.loc[segment_indexes, 'Probability'] = np.nan
            continue

        with np.errstate(divide='ignore'):
            log_weights = (
                np.log(probabilities)
                - NON_PARSIMONY_PENALTY_COEFFICIENT * average_events
            )
        finite_weights = np.isfinite(log_weights)
        if not finite_weights.any():
            penalized_table.loc[segment_indexes, 'Probability'] = np.nan
            continue

        log_weights -= np.max(log_weights[finite_weights])
        weights = np.exp(log_weights)
        normalizer = weights.sum()
        if not np.isfinite(normalizer) or normalizer <= 0:
            penalized_table.loc[segment_indexes, 'Probability'] = np.nan
            continue
        penalized_table.loc[segment_indexes, 'Probability'] = (
            weights / normalizer
        )

    return penalized_table


def _read_table(path, table_name):
    try:
        table = pd.read_csv(
            path,
            sep='\t',
            dtype={
                'Sample_ID': str,
                'Segment_ID': str,
                'Route': str,
                'Chromosome': str,
            },
            keep_default_na=False,
        )
    except pd.errors.EmptyDataError as error:
        raise ValueError(f'The {table_name} has no header') from error
    validate_pascal_snake_case_columns(
        table.columns,
        table_name.capitalize(),
    )
    return table


def _validate_route_table(route_table, sample_id):
    _require_columns(route_table, ROUTE_TABLE_COLUMNS, 'route table')
    _validate_route_rows(route_table)

    unexpected_sample_ids = sorted(
        set(route_table['Sample_ID']) - {str(sample_id)}
    )
    if unexpected_sample_ids:
        raise ValueError(
            f'The route table for sample {sample_id} contains other '
            f"Sample_ID values: {', '.join(unexpected_sample_ids)}"
        )

    route_table['Probability'] = pd.to_numeric(
        route_table['Probability'],
        errors='coerce',
    )
    route_table['Average_N_Events'] = pd.to_numeric(
        route_table['Average_N_Events'],
        errors='coerce',
    )
    if (route_table['Probability'] < 0).any():
        raise ValueError('Route probabilities must be non-negative')
    if (route_table['Average_N_Events'] < 0).any():
        raise ValueError('Average_N_Events values must be non-negative')

    for segment_key, segment_table in route_table.groupby(
        SEGMENT_KEY_COLUMNS,
        sort=False,
    ):
        inconsistent_columns = [
            column
            for column in _SEGMENT_CONSTANT_COLUMNS
            if segment_table[column].nunique(dropna=False) != 1
        ]
        if inconsistent_columns:
            raise ValueError(
                'Route rows for segment '
                f'{segment_key[1]} have inconsistent segment metadata: '
                f"{', '.join(inconsistent_columns)}"
            )


def _validate_gain_timing_table(gain_timing_table, route_table, sample_id):
    _require_columns(
        gain_timing_table,
        GAIN_TIMING_TABLE_COLUMNS,
        'gain timing table',
    )
    if gain_timing_table.empty:
        return

    _validate_identifiers(
        gain_timing_table,
        ROUTE_KEY_COLUMNS + ['Node'],
        'gain timing table',
    )
    unexpected_sample_ids = sorted(
        set(gain_timing_table['Sample_ID']) - {str(sample_id)}
    )
    if unexpected_sample_ids:
        raise ValueError(
            f'The gain timing table for sample {sample_id} contains other '
            f"Sample_ID values: {', '.join(unexpected_sample_ids)}"
        )

    duplicate_nodes = gain_timing_table.duplicated(
        subset=ROUTE_KEY_COLUMNS + ['Node'],
        keep=False,
    )
    if duplicate_nodes.any():
        duplicate_table = gain_timing_table.loc[
            duplicate_nodes,
            ROUTE_KEY_COLUMNS + ['Node'],
        ].drop_duplicates()
        keys = duplicate_table.itertuples(index=False, name=None)
        raise ValueError(
            'The gain timing table must contain exactly one row per '
            f'route/node; duplicate rows: {_format_keys(keys)}'
        )

    route_keys = set(
        route_table[ROUTE_KEY_COLUMNS].itertuples(index=False, name=None)
    )
    timing_route_keys = set(
        gain_timing_table[ROUTE_KEY_COLUMNS].itertuples(
            index=False,
            name=None,
        )
    )
    unexpected_routes = sorted(timing_route_keys - route_keys, key=str)
    if unexpected_routes:
        raise ValueError(
            'The gain timing table contains routes absent from the route '
            f'table: {_format_keys(unexpected_routes)}'
        )


def load_timing_from_dict(segment_path):
    with bz2.BZ2File(segment_path, 'rb') as input_file:
        return pickle.load(input_file)


def _validate_timing_array(array, segment_id, route, label):
    if not isinstance(array, np.ndarray):
        raise ValueError(
            f'Timing dictionary {label} for segment {segment_id}, route '
            f'{route} must be a numpy array'
        )
    if array.ndim != 1:
        raise ValueError(
            f'Timing dictionary {label} for segment {segment_id}, route '
            f'{route} must be one-dimensional'
        )
    if array.size == 0:
        raise ValueError(
            f'Timing dictionary {label} for segment {segment_id}, route '
            f'{route} must not be empty'
        )
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(
            f'Timing dictionary {label} for segment {segment_id}, route '
            f'{route} must be numeric'
        )


def _validate_segment_timing_dict(
    segment_route_table,
    segment_timing_table,
    timing_dict,
    segment_id,
):
    if not isinstance(timing_dict, Mapping):
        raise ValueError(
            f'Timing dictionary for segment {segment_id} must be a mapping'
        )

    expected_routes = set(segment_route_table['Route'])
    actual_routes = set(timing_dict)
    missing_routes = sorted(expected_routes - actual_routes, key=str)
    unexpected_routes = sorted(actual_routes - expected_routes, key=str)
    if missing_routes or unexpected_routes:
        details = []
        if missing_routes:
            details.append(f"missing: {', '.join(map(str, missing_routes))}")
        if unexpected_routes:
            details.append(
                f"unexpected: {', '.join(map(str, unexpected_routes))}"
            )
        raise ValueError(
            'Timing dictionary route coverage does not match the route '
            f"table for segment {segment_id} ({'; '.join(details)})"
        )

    for route in segment_route_table['Route']:
        route_entry = timing_dict[route]
        if not isinstance(route_entry, Mapping):
            raise ValueError(
                f'Timing dictionary entry for segment {segment_id}, route '
                f'{route} must be a mapping'
            )
        route_timing = route_entry.get('Timing')
        if not isinstance(route_timing, Mapping):
            raise ValueError(
                f'Timing dictionary entry for segment {segment_id}, route '
                f'{route} must contain a Timing mapping'
            )
        if 'WGD' not in route_timing:
            raise ValueError(
                f'Timing dictionary for segment {segment_id}, route {route} '
                'is missing the WGD array'
            )

        route_timing_metadata = segment_timing_table.loc[
            segment_timing_table['Route'] == route
        ]
        metadata_nodes = set(route_timing_metadata['Node'])
        dictionary_nodes = set(route_timing) - {'WGD'}
        missing_nodes = sorted(metadata_nodes - dictionary_nodes, key=str)
        unexpected_nodes = sorted(dictionary_nodes - metadata_nodes, key=str)
        if missing_nodes or unexpected_nodes:
            details = []
            if missing_nodes:
                details.append(
                    f"missing: {', '.join(map(str, missing_nodes))}"
                )
            if unexpected_nodes:
                details.append(
                    f"unexpected: {', '.join(map(str, unexpected_nodes))}"
                )
            raise ValueError(
                'Timing dictionary node coverage does not match the gain '
                f'timing table for segment {segment_id}, route {route} '
                f"({'; '.join(details)})"
            )

        expected_length = None
        for label in ['WGD', *route_timing_metadata['Node'].tolist()]:
            array = route_timing[label]
            _validate_timing_array(array, segment_id, route, label)
            if expected_length is None:
                expected_length = array.size
            elif array.size != expected_length:
                raise ValueError(
                    'All timing arrays for segment '
                    f'{segment_id}, route {route} must have the same length'
                )


def _normalized_route_probabilities(segment_route_table):
    probabilities = segment_route_table['Probability'].to_numpy(dtype=float)
    normalizer = probabilities.sum()
    if (
        not np.isfinite(probabilities).all()
        or not np.isfinite(normalizer)
        or normalizer <= 0
    ):
        return None
    return probabilities / normalizer


def produce_timing_segment_tables(
    segment_route_table,
    segment_timing_table,
    timing_dict,
    segment_id,
    n_samples=100,
):
    """Draw routes once per posterior sample and timings jointly per route."""
    _validate_segment_timing_dict(
        segment_route_table,
        segment_timing_table,
        timing_dict,
        segment_id,
    )
    probabilities = _normalized_route_probabilities(segment_route_table)
    if probabilities is None:
        return None, None

    sampled_routes = np.random.choice(
        segment_route_table['Route'].to_numpy(),
        size=n_samples,
        replace=True,
        p=probabilities,
    )
    segment_metadata = segment_route_table.iloc[0][
        SEGMENT_METADATA_COLUMNS
    ].to_dict()
    timing_metadata_by_route = {
        route: route_table.set_index('Node')['Node_Phasing'].to_dict()
        for route, route_table in segment_timing_table.groupby(
            'Route',
            sort=False,
        )
    }

    gain_draws = []
    route_draws = []
    for posterior_sample_index, route in enumerate(sampled_routes):
        route_timing = timing_dict[route]['Timing']
        timing_sample_index = np.random.randint(
            0,
            route_timing['WGD'].size,
        )
        wgd_timing = route_timing['WGD'][timing_sample_index]

        route_draws.append({
            **segment_metadata,
            'Posterior_Sample_Index': posterior_sample_index,
            'Route': route,
            'WGD_Timing': wgd_timing,
        })

        node_phasing = timing_metadata_by_route.get(route, {})
        nodes = list(node_phasing)
        gain_timings = np.asarray([
            route_timing[node][timing_sample_index]
            for node in nodes
        ])
        timing_order = np.argsort(gain_timings, kind='stable')
        for gain_index, timing_position in enumerate(timing_order):
            node = nodes[timing_position]
            gain_draws.append({
                **segment_metadata,
                'Posterior_Sample_Index': posterior_sample_index,
                'Route': route,
                'Node': node,
                'Node_Phasing': node_phasing[node],
                'Gain_Timing': gain_timings[timing_position],
                'WGD_Timing': wgd_timing,
                'Gain_Index': gain_index,
            })

    return (
        pd.DataFrame(gain_draws, columns=GAIN_DRAW_COLUMNS),
        pd.DataFrame(route_draws, columns=ROUTE_DRAW_COLUMNS),
    )


def _validate_n_posterior_samples(n_posterior_samples):
    if (
        isinstance(n_posterior_samples, (bool, np.bool_))
        or not isinstance(n_posterior_samples, Integral)
        or n_posterior_samples <= 0
    ):
        raise ValueError('n_posterior_samples must be a positive integer')


def get_sample_posterior_tables(
    route_table_path,
    timing_table_path,
    input_dir,
    sample_id: str,
    n_posterior_samples: int = 100,
    apply_penalty: bool = False,
):
    """Return gain draws and the one-row-per-draw route ledger."""
    if not isinstance(apply_penalty, (bool, np.bool_)):
        raise ValueError('apply_penalty must be a boolean')
    _validate_n_posterior_samples(n_posterior_samples)

    route_table = _read_table(route_table_path, 'route table')
    gain_timing_table = _read_table(
        timing_table_path,
        'gain timing table',
    )
    _validate_route_table(route_table, sample_id)
    _validate_gain_timing_table(gain_timing_table, route_table, sample_id)
    if apply_penalty:
        route_table = apply_non_parsimony_penalty(route_table)

    timing_dict_dir = Path(input_dir) / f'{sample_id}_timing_dicts'
    gain_frames = []
    route_frames = []
    for (_, segment_id), segment_route_table in route_table.groupby(
        SEGMENT_KEY_COLUMNS,
        sort=False,
    ):
        segment_timing_table = gain_timing_table.loc[
            (gain_timing_table['Sample_ID'] == str(sample_id))
            & (gain_timing_table['Segment_ID'] == segment_id)
        ]
        timing_dict_path = (
            timing_dict_dir / f'{segment_id}_timing_dict.bz2'
        )
        if not timing_dict_path.is_file():
            raise FileNotFoundError(
                f'Missing timing dictionary for segment {segment_id}: '
                f'{timing_dict_path}'
            )
        timing_dict = load_timing_from_dict(timing_dict_path)

        gain_frame, route_frame = produce_timing_segment_tables(
            segment_route_table,
            segment_timing_table,
            timing_dict,
            segment_id,
            n_samples=n_posterior_samples,
        )
        if gain_frame is None:
            warnings.warn(
                f'WARNING {segment_id} has invalid route probabilities; '
                'skipping in posterior.',
                stacklevel=2,
            )
            continue
        gain_frames.append(gain_frame)
        route_frames.append(route_frame)

    if gain_frames:
        gain_draw_table = pd.concat(gain_frames, ignore_index=True)
        gain_draw_table = gain_draw_table.sort_values(
            by=SEGMENT_KEY_COLUMNS
            + ['Posterior_Sample_Index', 'Gain_Index'],
            kind='stable',
        ).reset_index(drop=True)
    else:
        gain_draw_table = pd.DataFrame(columns=GAIN_DRAW_COLUMNS)

    if route_frames:
        route_draw_table = pd.concat(route_frames, ignore_index=True)
        route_draw_table = route_draw_table.sort_values(
            by=SEGMENT_KEY_COLUMNS + ['Posterior_Sample_Index'],
            kind='stable',
        ).reset_index(drop=True)
    else:
        route_draw_table = pd.DataFrame(columns=ROUTE_DRAW_COLUMNS)

    return gain_draw_table, route_draw_table


def _validate_draw_tables(gain_draw_table, route_draw_table):
    _require_columns(gain_draw_table, GAIN_DRAW_COLUMNS, 'gain draw table')
    _require_columns(route_draw_table, ROUTE_DRAW_COLUMNS, 'route draw table')
    if route_draw_table.empty:
        if not gain_draw_table.empty:
            raise ValueError(
                'The gain draw table cannot contain rows when the route draw '
                'table is empty'
            )
        return

    draw_key_columns = SEGMENT_KEY_COLUMNS + ['Posterior_Sample_Index']
    _validate_identifiers(
        route_draw_table,
        SEGMENT_KEY_COLUMNS + ['Route'],
        'route draw table',
    )
    duplicate_draws = route_draw_table.duplicated(
        subset=draw_key_columns,
        keep=False,
    )
    if duplicate_draws.any():
        raise ValueError(
            'The route draw table must contain exactly one row per segment '
            'and Posterior_Sample_Index'
        )

    duplicate_gains = gain_draw_table.duplicated(
        subset=draw_key_columns + ['Gain_Index'],
        keep=False,
    )
    if duplicate_gains.any():
        raise ValueError(
            'The gain draw table must contain at most one row per segment, '
            'Posterior_Sample_Index, and Gain_Index'
        )

    if gain_draw_table.empty:
        return
    _validate_identifiers(
        gain_draw_table,
        SEGMENT_KEY_COLUMNS + ['Route'],
        'gain draw table',
    )
    ledger_routes = route_draw_table.set_index(draw_key_columns)['Route']
    gain_keys = pd.MultiIndex.from_frame(gain_draw_table[draw_key_columns])
    matching_routes = ledger_routes.reindex(gain_keys)
    if matching_routes.isna().any():
        raise ValueError(
            'Every gain draw must reference a draw in the route draw table'
        )
    if not np.array_equal(
        matching_routes.to_numpy(),
        gain_draw_table['Route'].to_numpy(),
    ):
        raise ValueError(
            'Gain draw routes must match the corresponding route draw rows'
        )


def get_segment_posterior_table_summary(
    segment_gain_draw_table,
    segment_route_draw_table,
):
    """Summarize one segment using every route draw as the denominator."""
    _validate_draw_tables(
        segment_gain_draw_table,
        segment_route_draw_table,
    )
    segment_keys = segment_route_draw_table[
        SEGMENT_KEY_COLUMNS
    ].drop_duplicates()
    if len(segment_keys) > 1:
        raise ValueError('Segment summary inputs must contain one segment')

    n_samples = len(segment_route_draw_table)
    summary_rows = []
    for gain_index, gain_index_table in segment_gain_draw_table.groupby(
        'Gain_Index',
        sort=True,
    ):
        gain_timings = gain_index_table['Gain_Timing'].to_numpy(dtype=float)
        proportion = (
            gain_index_table['Posterior_Sample_Index'].nunique() / n_samples
        )

        wgd_timings = gain_index_table['WGD_Timing'].to_numpy(dtype=float)
        if np.isnan(wgd_timings).any():
            wgd_timing_median = np.nan
            wgd_timing_low_ci = np.nan
            wgd_timing_high_ci = np.nan
            pre_wgd_probability = np.nan
            post_wgd_probability = np.nan
        else:
            wgd_timing_median = np.median(wgd_timings)
            wgd_timing_low_ci = np.percentile(wgd_timings, 2.5)
            wgd_timing_high_ci = np.percentile(wgd_timings, 97.5)
            pre_wgd_probability = np.mean(wgd_timings > gain_timings)
            post_wgd_probability = 1 - pre_wgd_probability

        summary_rows.append({
            'Gain_Index': gain_index + 1,
            'Proportion': proportion,
            'Timing_Median': np.median(gain_timings),
            'Timing_Low_CI': np.percentile(gain_timings, 2.5),
            'Timing_High_CI': np.percentile(gain_timings, 97.5),
            'Pre_WGD_Probability': pre_wgd_probability,
            'Post_WGD_Probability': post_wgd_probability,
            'WGD_Timing_Median': wgd_timing_median,
            'WGD_Timing_Low_CI': wgd_timing_low_ci,
            'WGD_Timing_High_CI': wgd_timing_high_ci,
        })
    return pd.DataFrame(summary_rows, columns=SUMMARY_VALUE_COLUMNS)


def get_sample_posterior_table_summary(
    gain_draw_table,
    route_draw_table,
    min_proportion_threshold=0.8,
):
    """Summarize gain draws with route-ledger draw counts as denominators."""
    if (
        isinstance(min_proportion_threshold, (bool, np.bool_))
        or not isinstance(min_proportion_threshold, Real)
        or not 0 <= min_proportion_threshold <= 1
    ):
        raise ValueError(
            'min_proportion_threshold must be a number between 0 and 1'
        )
    _validate_draw_tables(gain_draw_table, route_draw_table)
    if route_draw_table.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    segment_metadata = route_draw_table[
        SEGMENT_METADATA_COLUMNS
    ].drop_duplicates()
    inconsistent_metadata = segment_metadata.duplicated(
        subset=SEGMENT_KEY_COLUMNS,
        keep=False,
    )
    if inconsistent_metadata.any():
        raise ValueError(
            'Route draw rows for a segment must have consistent metadata'
        )

    summary_frames = []
    for segment_key, segment_route_draw_table in route_draw_table.groupby(
        SEGMENT_KEY_COLUMNS,
        sort=False,
    ):
        segment_gain_draw_table = gain_draw_table.loc[
            (gain_draw_table['Sample_ID'] == segment_key[0])
            & (gain_draw_table['Segment_ID'] == segment_key[1])
        ]
        segment_summary = get_segment_posterior_table_summary(
            segment_gain_draw_table,
            segment_route_draw_table,
        )
        if segment_summary.empty:
            continue
        segment_summary['Sample_ID'] = segment_key[0]
        segment_summary['Segment_ID'] = segment_key[1]
        summary_frames.append(segment_summary)

    if not summary_frames:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    summary_table = pd.concat(summary_frames, ignore_index=True)
    summary_table = pd.merge(
        segment_metadata,
        summary_table,
        on=SEGMENT_KEY_COLUMNS,
        how='inner',
        validate='one_to_many',
    )
    summary_table = summary_table.loc[
        summary_table['Proportion'] >= min_proportion_threshold
    ]
    return summary_table[SUMMARY_COLUMNS].reset_index(drop=True)
