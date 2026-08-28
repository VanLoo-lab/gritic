import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from gritic import posteriortablegen
from gritic.intervaltools import IntervalSpec
from gritic.tableschemas import (
    GAIN_DRAW_COLUMNS,
    GAIN_TIMING_TABLE_COLUMNS,
    ROUTE_DRAW_COLUMNS,
    ROUTE_TABLE_COLUMNS,
    SUMMARY_COLUMNS,
)
from gritic.timingio import write_timing_archive


def make_route_row(
    *,
    sample_id='sample',
    segment_id='segment',
    route='route-a',
    probability=1.0,
    penalized_probability=1.0,
    average_events=2.0,
    chromosome='1',
    segment_start=10,
    segment_end=20,
):
    return {
        'Sample_ID': sample_id,
        'Segment_ID': segment_id,
        'Route': route,
        'Chromosome': chromosome,
        'Segment_Start': segment_start,
        'Segment_End': segment_end,
        'Major_CN': 2,
        'Minor_CN': 1,
        'Total_CN': 3,
        'N_Mutations': 12,
        'Mutation_Rate': 0.01,
        'Probability': probability,
        'Penalized_Probability': penalized_probability,
        'Average_N_Events': average_events,
        'Average_Pre_WGD_Losses': 0.0,
        'Average_Post_WGD_Losses': 0.0,
        'Time': 0.1,
        'Density': 1.0,
        'WGD_Status': True,
        'WGD_Timing': 0.5,
        'WGD_Timing_CI_Low': 0.4,
        'WGD_Timing_CI_High': 0.6,
    }


def make_gain_timing_row(
    *,
    sample_id='sample',
    segment_id='segment',
    route='route-a',
    node=10,
    node_phasing='Major',
):
    return {
        'Sample_ID': sample_id,
        'Segment_ID': segment_id,
        'Route': route,
        'Node': node,
        'Node_Phasing': node_phasing,
        'Timing': 0.5,
        'Timing_CI_Low': 0.2,
        'Timing_CI_High': 0.8,
    }


def route_table(rows):
    return pd.DataFrame(rows, columns=ROUTE_TABLE_COLUMNS)


def timing_table(rows):
    return pd.DataFrame(rows, columns=GAIN_TIMING_TABLE_COLUMNS)


def make_route_draw(
    *,
    sample_id='sample',
    segment_id='segment',
    posterior_index=0,
    route='route-a',
    wgd_timing=0.5,
    chromosome='1',
):
    row = make_route_row(
        sample_id=sample_id,
        segment_id=segment_id,
        route=route,
        chromosome=chromosome,
    )
    return {
        column: row[column]
        for column in ROUTE_DRAW_COLUMNS
        if column in row
    } | {
        'Posterior_Sample_Index': posterior_index,
        'Route': route,
        'WGD_Timing': wgd_timing,
    }


def make_gain_draw(
    *,
    sample_id='sample',
    segment_id='segment',
    posterior_index=0,
    route='route-a',
    node=10,
    node_phasing='Major',
    gain_timing=0.25,
    wgd_timing=0.5,
    gain_index=1,
    chromosome='1',
):
    route_draw = make_route_draw(
        sample_id=sample_id,
        segment_id=segment_id,
        posterior_index=posterior_index,
        route=route,
        wgd_timing=wgd_timing,
        chromosome=chromosome,
    )
    return {
        column: route_draw[column]
        for column in GAIN_DRAW_COLUMNS
        if column in route_draw
    } | {
        'Posterior_Sample_Index': posterior_index,
        'Route': route,
        'Node': node,
        'Node_Phasing': node_phasing,
        'Gain_Timing': gain_timing,
        'WGD_Timing': wgd_timing,
        'Gain_Index': gain_index,
    }


class PenalizedProbabilityTest(unittest.TestCase):
    def test_penalty_is_normalized_independently_within_each_segment(self):
        table = pd.DataFrame([
            make_route_row(
                segment_id='s1',
                route='a',
                probability=0.75,
                average_events=1,
            ),
            make_route_row(
                segment_id='s1',
                route='b',
                probability=0.25,
                average_events=2,
            ),
            make_route_row(
                segment_id='s2',
                route='c',
                probability=2,
                average_events=3,
            ),
            make_route_row(
                segment_id='s2',
                route='d',
                probability=1,
                average_events=3,
            ),
        ])
        original = table.copy(deep=True)

        result = posteriortablegen.add_penalized_probability(table)

        pd.testing.assert_frame_equal(table, original)
        for _, segment in result.groupby(['Sample_ID', 'Segment_ID']):
            self.assertAlmostEqual(segment['Penalized_Probability'].sum(), 1)
        coefficient = posteriortablegen.NON_PARSIMONY_PENALTY_COEFFICIENT
        expected_s1 = np.array([
            0.75 * np.exp(-coefficient),
            0.25 * np.exp(-2 * coefficient),
        ])
        expected_s1 /= expected_s1.sum()
        np.testing.assert_allclose(
            result.loc[result['Segment_ID'] == 's1', 'Penalized_Probability'],
            expected_s1,
        )
        np.testing.assert_allclose(
            result.loc[result['Segment_ID'] == 's2', 'Penalized_Probability'],
            [2 / 3, 1 / 3],
        )

    def test_log_space_penalty_remains_stable_for_large_event_counts(self):
        table = pd.DataFrame([
            make_route_row(route='a', probability=0.5, average_events=10_000),
            make_route_row(route='b', probability=0.5, average_events=10_001),
        ])
        result = posteriortablegen.add_penalized_probability(table)
        self.assertTrue(np.isfinite(result['Penalized_Probability']).all())
        self.assertAlmostEqual(result['Penalized_Probability'].sum(), 1)
        self.assertGreater(
            result.loc[0, 'Penalized_Probability'],
            result.loc[1, 'Penalized_Probability'],
        )

    def test_zero_or_nonfinite_segment_weights_produce_nan_probabilities(self):
        for probabilities, event_counts in (
            ([0.0, 0.0], [1.0, 2.0]),
            ([np.nan, 1.0], [1.0, 2.0]),
            ([0.5, 0.5], [np.inf, 2.0]),
        ):
            with self.subTest(probabilities=probabilities, events=event_counts):
                table = pd.DataFrame([
                    make_route_row(
                        route='a',
                        probability=probabilities[0],
                        average_events=event_counts[0],
                    ),
                    make_route_row(
                        route='b',
                        probability=probabilities[1],
                        average_events=event_counts[1],
                    ),
                ])
                result = posteriortablegen.add_penalized_probability(table)
                self.assertTrue(result['Penalized_Probability'].isna().all())

    def test_zero_probability_route_stays_zero_when_segment_is_valid(self):
        table = pd.DataFrame([
            make_route_row(route='a', probability=1, average_events=3),
            make_route_row(route='b', probability=0, average_events=0),
        ])
        result = posteriortablegen.add_penalized_probability(table)
        np.testing.assert_array_equal(
            result['Penalized_Probability'].to_numpy(),
            [1.0, 0.0],
        )

    def test_rejects_negative_probabilities_and_event_counts(self):
        cases = (
            ({'probability': -0.1}, 'probabilities'),
            ({'average_events': -1}, 'Average_N_Events'),
        )
        for overrides, message in cases:
            with self.subTest(overrides=overrides):
                row = make_route_row(**overrides)
                with self.assertRaisesRegex(ValueError, message):
                    posteriortablegen.add_penalized_probability(
                        pd.DataFrame([row])
                    )

    def test_requires_keys_and_rejects_duplicate_routes(self):
        with self.assertRaisesRegex(ValueError, 'missing columns: Route'):
            posteriortablegen.add_penalized_probability(pd.DataFrame({
                'Sample_ID': ['s'],
                'Segment_ID': ['x'],
                'Probability': [1],
                'Average_N_Events': [1],
            }))
        duplicate = pd.DataFrame([
            make_route_row(),
            make_route_row(),
        ])
        with self.assertRaisesRegex(ValueError, 'duplicate routes'):
            posteriortablegen.add_penalized_probability(duplicate)


class InputTableValidationTest(unittest.TestCase):
    def test_route_table_coerces_numeric_strings(self):
        table = route_table([make_route_row(
            probability='0.75',
            penalized_probability='0.25',
            average_events='2.5',
        )])
        posteriortablegen._validate_route_table(table, 'sample')
        self.assertEqual(table.loc[0, 'Probability'], 0.75)
        self.assertEqual(table.loc[0, 'Penalized_Probability'], 0.25)
        self.assertEqual(table.loc[0, 'Average_N_Events'], 2.5)

    def test_route_table_rejects_missing_columns_and_blank_identifiers(self):
        table = route_table([make_route_row()]).drop(columns=['Density'])
        with self.assertRaisesRegex(ValueError, 'missing columns: Density'):
            posteriortablegen._validate_route_table(table, 'sample')

        for column in ('Sample_ID', 'Segment_ID', 'Route'):
            with self.subTest(column=column):
                table = route_table([make_route_row()])
                table.loc[0, column] = '  '
                with self.assertRaisesRegex(ValueError, column):
                    posteriortablegen._validate_route_table(table, 'sample')

    def test_route_table_rejects_wrong_sample_and_negative_values(self):
        table = route_table([make_route_row(sample_id='other')])
        with self.assertRaisesRegex(ValueError, 'other Sample_ID'):
            posteriortablegen._validate_route_table(table, 'sample')

        cases = (
            ('Probability', -0.1, 'probabilities'),
            ('Penalized_Probability', -0.1, 'probabilities'),
            ('Average_N_Events', -0.1, 'Average_N_Events'),
        )
        for column, value, message in cases:
            with self.subTest(column=column):
                table = route_table([make_route_row()])
                table.loc[0, column] = value
                with self.assertRaisesRegex(ValueError, message):
                    posteriortablegen._validate_route_table(table, 'sample')

    def test_route_table_rejects_inconsistent_segment_metadata(self):
        first = make_route_row(route='a')
        second = make_route_row(route='b')
        second['Chromosome'] = '2'
        table = route_table([first, second])
        with self.assertRaisesRegex(ValueError, 'Chromosome'):
            posteriortablegen._validate_route_table(table, 'sample')

    def test_gain_timing_table_accepts_empty_schema_and_valid_rows(self):
        routes = route_table([make_route_row()])
        posteriortablegen._validate_gain_timing_table(
            timing_table([]),
            routes,
            'sample',
        )
        gains = timing_table([
            make_gain_timing_row(node=10, node_phasing='Major'),
            make_gain_timing_row(node=11, node_phasing='Minor'),
        ])
        posteriortablegen._validate_gain_timing_table(gains, routes, 'sample')

    def test_gain_timing_table_rejects_missing_columns_and_identifiers(self):
        routes = route_table([make_route_row()])
        gains = timing_table([make_gain_timing_row()]).drop(
            columns=['Timing_CI_High']
        )
        with self.assertRaisesRegex(ValueError, 'Timing_CI_High'):
            posteriortablegen._validate_gain_timing_table(
                gains,
                routes,
                'sample',
            )

        for column, value in (
            ('Sample_ID', ''),
            ('Segment_ID', None),
            ('Route', ' '),
            ('Node', np.nan),
        ):
            with self.subTest(column=column):
                gains = timing_table([make_gain_timing_row()])
                gains.loc[0, column] = value
                with self.assertRaisesRegex(ValueError, column):
                    posteriortablegen._validate_gain_timing_table(
                        gains,
                        routes,
                        'sample',
                    )

    def test_gain_timing_table_rejects_wrong_sample_duplicate_and_unknown_route(self):
        routes = route_table([make_route_row()])
        cases = []
        cases.append((
            timing_table([make_gain_timing_row(sample_id='other')]),
            'other Sample_ID',
        ))
        cases.append((
            timing_table([make_gain_timing_row(), make_gain_timing_row()]),
            'duplicate rows',
        ))
        cases.append((
            timing_table([make_gain_timing_row(route='missing')]),
            'routes absent',
        ))
        for gains, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    posteriortablegen._validate_gain_timing_table(
                        gains,
                        routes,
                        'sample',
                    )

    def test_read_table_preserves_leading_zero_identifiers(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'table.tsv'
            pd.DataFrame({
                'Sample_ID': ['001'],
                'Segment_ID': ['0002'],
                'Route': ['0003'],
                'Chromosome': ['01'],
            }).to_csv(path, sep='\t', index=False)
            result = posteriortablegen._read_table(path, 'test table')
        self.assertEqual(
            result.iloc[0].to_dict(),
            {
                'Sample_ID': '001',
                'Segment_ID': '0002',
                'Route': '0003',
                'Chromosome': '01',
            },
        )

    def test_read_table_rejects_headerless_empty_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'empty.tsv'
            path.touch()
            with self.assertRaisesRegex(ValueError, 'has no header'):
                posteriortablegen._read_table(path, 'route table')


class TimingDictionaryValidationTest(unittest.TestCase):
    def setUp(self):
        self.routes = route_table([
            make_route_row(route='a', probability=0.5),
            make_route_row(route='b', probability=0.5),
        ])
        self.gains = timing_table([
            make_gain_timing_row(route='a', node=10),
            make_gain_timing_row(route='a', node=11),
            make_gain_timing_row(route='b', node=20),
        ])
        self.valid = {
            'a': {'Timing': {
                'WGD': np.array([0.5, 0.6]),
                10: np.array([0.2, 0.3]),
                11: np.array([0.4, 0.1]),
            }},
            'b': {'Timing': {
                'WGD': np.array([0.7, 0.8]),
                20: np.array([0.6, 0.9]),
            }},
        }

    def validate(self, timing_dict):
        posteriortablegen._validate_segment_timing_dict(
            self.routes,
            self.gains,
            timing_dict,
            'segment',
        )

    def test_accepts_complete_numeric_one_dimensional_arrays(self):
        self.validate(self.valid)

    def test_rejects_non_mapping_root_and_route_entry(self):
        for value, message in (
            ([], 'must be a mapping'),
            ({'a': [], 'b': self.valid['b']}, 'route a must be a mapping'),
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(value)

    def test_rejects_missing_and_unexpected_routes(self):
        cases = (
            ({'a': self.valid['a']}, 'missing: b'),
            ({**self.valid, 'c': self.valid['b']}, 'unexpected: c'),
        )
        for value, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(value)

    def test_rejects_missing_timing_mapping_and_wgd_array(self):
        cases = (
            ({**self.valid, 'a': {}}, 'contain a Timing mapping'),
            (
                {**self.valid, 'a': {'Timing': {10: np.array([0.2])}}},
                'missing the WGD array',
            ),
        )
        for value, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(value)

    def test_rejects_missing_and_unexpected_node_arrays(self):
        missing = {
            **self.valid,
            'a': {'Timing': {
                'WGD': np.array([0.5, 0.6]),
                10: np.array([0.2, 0.3]),
            }},
        }
        unexpected = {
            **self.valid,
            'a': {'Timing': {
                **self.valid['a']['Timing'],
                99: np.array([0.1, 0.2]),
            }},
        }
        for value, message in (
            (missing, 'missing: 11'),
            (unexpected, 'unexpected: 99'),
        ):
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(value)

    def test_rejects_invalid_array_types_shapes_sizes_and_dtypes(self):
        invalid_arrays = (
            ([0.5, 0.6], 'numpy array'),
            (np.array([[0.5, 0.6]]), 'one-dimensional'),
            (np.array([], dtype=float), 'must not be empty'),
            (np.array(['a', 'b']), 'must be numeric'),
        )
        for invalid_array, message in invalid_arrays:
            with self.subTest(message=message):
                value = {
                    **self.valid,
                    'a': {'Timing': {
                        **self.valid['a']['Timing'],
                        10: invalid_array,
                    }},
                }
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(value)

    def test_rejects_different_array_lengths_within_route(self):
        value = {
            **self.valid,
            'a': {'Timing': {
                **self.valid['a']['Timing'],
                10: np.array([0.2]),
            }},
        }
        with self.assertRaisesRegex(ValueError, 'same length'):
            self.validate(value)


class PosteriorSegmentGenerationTest(unittest.TestCase):
    def setUp(self):
        self.routes = route_table([
            make_route_row(route='a', probability=0.25),
            make_route_row(route='b', probability=0.75),
        ])
        self.gains = timing_table([
            make_gain_timing_row(route='a', node=10, node_phasing='Major'),
            make_gain_timing_row(route='a', node=11, node_phasing='Minor'),
            make_gain_timing_row(route='b', node=20, node_phasing='Major'),
        ])
        self.timing_dict = {
            'a': {'Timing': {
                'WGD': np.array([0.5, 0.6]),
                10: np.array([0.2, 0.3]),
                11: np.array([0.2, 0.1]),
            }},
            'b': {'Timing': {
                'WGD': np.array([0.7, 0.8]),
                20: np.array([0.9, 0.4]),
            }},
        }

    def test_probability_normalization(self):
        table = pd.DataFrame({'Probability': [2.0, 1.0]})
        np.testing.assert_allclose(
            posteriortablegen._normalized_route_probabilities(table),
            [2 / 3, 1 / 3],
        )
        for values in ([0.0, 0.0], [np.nan, 1.0], [np.inf, 1.0]):
            with self.subTest(values=values):
                table = pd.DataFrame({'Probability': values})
                self.assertIsNone(
                    posteriortablegen._normalized_route_probabilities(table)
                )

    def test_draws_route_once_and_uses_one_joint_timing_index(self):
        with mock.patch.object(
            posteriortablegen.np.random,
            'choice',
            return_value=np.array(['b', 'a']),
        ) as choice, mock.patch.object(
            posteriortablegen.np.random,
            'randint',
            side_effect=[1, 0],
        ) as randint:
            gain_draws, route_draws = (
                posteriortablegen.produce_timing_segment_tables(
                    self.routes,
                    self.gains,
                    self.timing_dict,
                    'segment',
                    n_samples=2,
                )
            )

        np.testing.assert_allclose(choice.call_args.kwargs['p'], [0.25, 0.75])
        self.assertEqual(randint.call_count, 2)
        self.assertEqual(route_draws['Route'].tolist(), ['b', 'a'])
        self.assertEqual(route_draws['WGD_Timing'].tolist(), [0.8, 0.5])
        self.assertEqual(gain_draws['Route'].tolist(), ['b', 'a', 'a'])
        self.assertEqual(gain_draws['Gain_Timing'].tolist(), [0.4, 0.2, 0.2])
        self.assertEqual(gain_draws['Node'].tolist(), [20, 10, 11])
        self.assertEqual(gain_draws['Gain_Index'].tolist(), [1, 1, 2])
        self.assertEqual(
            gain_draws['Node_Phasing'].tolist(),
            ['Major', 'Major', 'Minor'],
        )
        self.assertEqual(
            gain_draws.loc[gain_draws['Route'] == 'a', 'WGD_Timing'].tolist(),
            [0.5, 0.5],
        )
        self.assertEqual(list(gain_draws.columns), GAIN_DRAW_COLUMNS)
        self.assertEqual(list(route_draws.columns), ROUTE_DRAW_COLUMNS)

    def test_invalid_probabilities_return_pair_of_none(self):
        self.routes['Probability'] = 0
        result = posteriortablegen.produce_timing_segment_tables(
            self.routes,
            self.gains,
            self.timing_dict,
            'segment',
        )
        self.assertEqual(result, (None, None))

    def test_n_posterior_samples_requires_positive_non_boolean_integer(self):
        for value in (0, -1, 1.5, True, np.bool_(False), '2'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'positive integer'):
                    posteriortablegen._validate_n_posterior_samples(value)
        for value in (1, np.int64(2)):
            posteriortablegen._validate_n_posterior_samples(value)


class SamplePosteriorArchiveIntegrationTest(unittest.TestCase):
    def test_reads_tables_and_archives_for_multiple_segments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [
                make_route_row(
                    sample_id='001',
                    segment_id='0002',
                    route='r2',
                    probability=1,
                    penalized_probability=1,
                    chromosome='01',
                ),
                make_route_row(
                    sample_id='001',
                    segment_id='0001',
                    route='r1',
                    probability=1,
                    penalized_probability=1,
                    chromosome='01',
                ),
            ]
            route_path = root / 'routes.tsv'
            route_table(rows).to_csv(route_path, sep='\t', index=False)
            gains = timing_table([
                make_gain_timing_row(
                    sample_id='001', segment_id='0002', route='r2', node=20
                ),
                make_gain_timing_row(
                    sample_id='001', segment_id='0001', route='r1', node=10
                ),
            ])
            timing_path = root / 'gains.tsv'
            gains.to_csv(timing_path, sep='\t', index=False)
            archive_dir = root / '001_timing_dicts'
            write_timing_archive({
                'r2': {'Timing': {
                    'WGD': np.array([0.6, 0.7]),
                    20: np.array([0.3, 0.4]),
                }},
            }, archive_dir, '0002')
            write_timing_archive({
                'r1': {'Timing': {
                    'WGD': np.array([0.8, 0.9]),
                    10: np.array([0.1, 0.2]),
                }},
            }, archive_dir, '0001')

            with mock.patch.object(
                posteriortablegen.np.random,
                'choice',
                side_effect=[np.array(['r2', 'r2']), np.array(['r1', 'r1'])],
            ), mock.patch.object(
                posteriortablegen.np.random,
                'randint',
                side_effect=[0, 1, 1, 0],
            ):
                gain_draws, route_draws = (
                    posteriortablegen.get_sample_posterior_tables(
                        route_path,
                        timing_path,
                        root,
                        '001',
                        n_posterior_samples=2,
                    )
                )

        self.assertEqual(route_draws['Sample_ID'].unique().tolist(), ['001'])
        self.assertEqual(route_draws['Chromosome'].unique().tolist(), ['01'])
        self.assertEqual(
            list(route_draws[['Segment_ID', 'Posterior_Sample_Index']]
                 .itertuples(index=False, name=None)),
            [('0001', 0), ('0001', 1), ('0002', 0), ('0002', 1)],
        )
        self.assertEqual(
            list(gain_draws[['Segment_ID', 'Gain_Timing']]
                 .itertuples(index=False, name=None)),
            [('0001', 0.2), ('0001', 0.1), ('0002', 0.3), ('0002', 0.4)],
        )

    def test_uses_penalized_probability_column_when_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [
                make_route_row(
                    route='a', probability=0.9, penalized_probability=0.2
                ),
                make_route_row(
                    route='b', probability=0.1, penalized_probability=0.8
                ),
            ]
            route_path = root / 'routes.tsv'
            route_table(rows).to_csv(route_path, sep='\t', index=False)
            timing_path = root / 'gains.tsv'
            timing_table([
                make_gain_timing_row(route='a', node=1),
                make_gain_timing_row(route='b', node=2),
            ]).to_csv(timing_path, sep='\t', index=False)
            write_timing_archive({
                'a': {'Timing': {
                    'WGD': np.array([0.5]),
                    1: np.array([0.2]),
                }},
                'b': {'Timing': {
                    'WGD': np.array([0.5]),
                    2: np.array([0.3]),
                }},
            }, root / 'sample_timing_dicts', 'segment')

            with mock.patch.object(
                posteriortablegen.np.random,
                'choice',
                return_value=np.array(['b']),
            ) as choice, mock.patch.object(
                posteriortablegen.np.random,
                'randint',
                return_value=0,
            ):
                posteriortablegen.get_sample_posterior_tables(
                    route_path,
                    timing_path,
                    root,
                    'sample',
                    n_posterior_samples=1,
                    apply_penalty=True,
                )
        np.testing.assert_allclose(choice.call_args.kwargs['p'], [0.2, 0.8])

    def test_incomplete_archive_pair_raises_file_not_found(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            route_path = root / 'routes.tsv'
            timing_path = root / 'gains.tsv'
            route_table([make_route_row()]).to_csv(
                route_path, sep='\t', index=False
            )
            timing_table([make_gain_timing_row()]).to_csv(
                timing_path, sep='\t', index=False
            )
            archive_dir = root / 'sample_timing_dicts'
            archive_dir.mkdir()
            (archive_dir / 'segment_timing_dict.npz').touch()
            with self.assertRaisesRegex(
                FileNotFoundError,
                'Incomplete timing archive pair',
            ):
                posteriortablegen.get_sample_posterior_tables(
                    route_path,
                    timing_path,
                    root,
                    'sample',
                )

    def test_invalid_segment_probabilities_warn_and_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            route_path = root / 'routes.tsv'
            timing_path = root / 'gains.tsv'
            route_table([make_route_row(
                probability=np.nan,
                penalized_probability=np.nan,
            )]).to_csv(route_path, sep='\t', index=False)
            timing_table([make_gain_timing_row()]).to_csv(
                timing_path, sep='\t', index=False
            )
            write_timing_archive({
                'route-a': {'Timing': {
                    'WGD': np.array([0.5]),
                    10: np.array([0.2]),
                }},
            }, root / 'sample_timing_dicts', 'segment')
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                gain_draws, route_draws = (
                    posteriortablegen.get_sample_posterior_tables(
                        route_path,
                        timing_path,
                        root,
                        'sample',
                    )
                )
        self.assertEqual(len(caught), 1)
        self.assertIn('invalid route probabilities', str(caught[0].message))
        self.assertTrue(gain_draws.empty)
        self.assertTrue(route_draws.empty)
        self.assertEqual(list(gain_draws.columns), GAIN_DRAW_COLUMNS)
        self.assertEqual(list(route_draws.columns), ROUTE_DRAW_COLUMNS)

    def test_public_argument_validation_precedes_file_reads(self):
        for kwargs, message in (
            ({'apply_penalty': 1}, 'apply_penalty must be a boolean'),
            ({'n_posterior_samples': 0}, 'positive integer'),
        ):
            with self.subTest(kwargs=kwargs):
                with self.assertRaisesRegex(ValueError, message):
                    posteriortablegen.get_sample_posterior_tables(
                        'missing-routes',
                        'missing-timings',
                        'missing-input',
                        'sample',
                        **kwargs,
                    )


class PosteriorDrawValidationAndSummaryTest(unittest.TestCase):
    def valid_tables(self):
        route_draws = pd.DataFrame([
            make_route_draw(posterior_index=index, route=f'r{index % 2}')
            for index in range(4)
        ], columns=ROUTE_DRAW_COLUMNS)
        gain_draws = pd.DataFrame([
            make_gain_draw(
                posterior_index=index,
                route=f'r{index % 2}',
                gain_timing=[0.1, 0.2, 0.3, 0.4][index],
                wgd_timing=[0.2, 0.1, 0.5, 0.3][index],
                gain_index=1,
            )
            for index in range(4)
        ] + [
            make_gain_draw(
                posterior_index=index,
                route=f'r{index % 2}',
                node=11,
                gain_timing=[0.6, 0.7][index],
                wgd_timing=[0.2, 0.1][index],
                gain_index=2,
            )
            for index in range(2)
        ], columns=GAIN_DRAW_COLUMNS)
        return gain_draws, route_draws

    def test_validates_empty_tables_and_required_columns(self):
        posteriortablegen._validate_draw_tables(
            pd.DataFrame(columns=GAIN_DRAW_COLUMNS),
            pd.DataFrame(columns=ROUTE_DRAW_COLUMNS),
        )
        with self.assertRaisesRegex(ValueError, 'missing columns: Route'):
            posteriortablegen._validate_draw_tables(
                pd.DataFrame(columns=GAIN_DRAW_COLUMNS),
                pd.DataFrame(columns=[
                    column for column in ROUTE_DRAW_COLUMNS if column != 'Route'
                ]),
            )
        with self.assertRaisesRegex(ValueError, 'cannot contain rows'):
            posteriortablegen._validate_draw_tables(
                pd.DataFrame([make_gain_draw()], columns=GAIN_DRAW_COLUMNS),
                pd.DataFrame(columns=ROUTE_DRAW_COLUMNS),
            )

    def test_rejects_duplicate_route_draws(self):
        gain_draws, route_draws = self.valid_tables()
        route_draws = pd.concat([route_draws, route_draws.iloc[[0]]])
        with self.assertRaisesRegex(ValueError, 'exactly one row'):
            posteriortablegen._validate_draw_tables(gain_draws, route_draws)

    def test_rejects_invalid_or_duplicate_gain_indexes(self):
        gain_draws, route_draws = self.valid_tables()
        for value in (0, -1, 1.5, True):
            with self.subTest(value=value):
                invalid = gain_draws.copy()
                invalid['Gain_Index'] = invalid['Gain_Index'].astype(object)
                invalid.loc[0, 'Gain_Index'] = value
                with self.assertRaisesRegex(ValueError, 'positive integers'):
                    posteriortablegen._validate_draw_tables(
                        invalid,
                        route_draws,
                    )
        duplicate = pd.concat([gain_draws, gain_draws.iloc[[0]]])
        with self.assertRaisesRegex(ValueError, 'at most one row'):
            posteriortablegen._validate_draw_tables(duplicate, route_draws)

    def test_rejects_gain_without_ledger_draw_or_with_different_route(self):
        gain_draws, route_draws = self.valid_tables()
        missing = gain_draws.copy()
        missing.loc[0, 'Posterior_Sample_Index'] = 99
        with self.assertRaisesRegex(ValueError, 'reference a draw'):
            posteriortablegen._validate_draw_tables(missing, route_draws)

        mismatch = gain_draws.copy()
        mismatch.loc[0, 'Route'] = 'different'
        with self.assertRaisesRegex(ValueError, 'routes must match'):
            posteriortablegen._validate_draw_tables(mismatch, route_draws)

    def test_segment_summary_uses_route_draws_as_denominator(self):
        gain_draws, route_draws = self.valid_tables()
        result = posteriortablegen.get_segment_posterior_table_summary(
            gain_draws,
            route_draws,
            interval=IntervalSpec(1),
        )
        first = result.loc[result['Gain_Index'] == 1].iloc[0]
        second = result.loc[result['Gain_Index'] == 2].iloc[0]
        self.assertEqual(first['Proportion'], 1)
        self.assertEqual(second['Proportion'], 0.5)
        self.assertEqual(first['Timing_Median'], 0.25)
        self.assertEqual(first['Timing_Low_CI'], 0.1)
        self.assertEqual(first['Timing_High_CI'], 0.4)
        self.assertEqual(first['Pre_WGD_Probability'], 0.5)
        self.assertEqual(first['Post_WGD_Probability'], 0.5)
        self.assertEqual(first['WGD_Timing_Median'], 0.25)

    def test_nan_wgd_timing_propagates_only_wgd_summary_fields(self):
        gain_draws, route_draws = self.valid_tables()
        gain_draws.loc[
            gain_draws['Gain_Index'] == 1,
            'WGD_Timing',
        ] = np.nan
        result = posteriortablegen.get_segment_posterior_table_summary(
            gain_draws,
            route_draws,
            interval=IntervalSpec(1),
        )
        first = result.loc[result['Gain_Index'] == 1].iloc[0]
        self.assertEqual(first['Timing_Median'], 0.25)
        for column in (
            'Pre_WGD_Probability',
            'Post_WGD_Probability',
            'WGD_Timing_Median',
            'WGD_Timing_Low_CI',
            'WGD_Timing_High_CI',
        ):
            self.assertTrue(np.isnan(first[column]), column)

    def test_segment_summary_rejects_multiple_segments(self):
        gain_draws, route_draws = self.valid_tables()
        extra_route = make_route_draw(segment_id='other', posterior_index=0)
        route_draws = pd.concat([
            route_draws,
            pd.DataFrame([extra_route], columns=ROUTE_DRAW_COLUMNS),
        ])
        with self.assertRaisesRegex(ValueError, 'one segment'):
            posteriortablegen.get_segment_posterior_table_summary(
                gain_draws,
                route_draws,
            )

    def test_sample_summary_filters_by_proportion_and_attaches_metadata(self):
        gain_draws, route_draws = self.valid_tables()
        result = posteriortablegen.get_sample_posterior_table_summary(
            gain_draws,
            route_draws,
            min_proportion_threshold=0.75,
            interval=IntervalSpec(1),
        )
        self.assertEqual(list(result.columns), SUMMARY_COLUMNS)
        self.assertEqual(result['Gain_Index'].tolist(), [1])
        self.assertEqual(result.loc[0, 'Sample_ID'], 'sample')
        self.assertEqual(result.loc[0, 'Segment_ID'], 'segment')
        self.assertEqual(result.loc[0, 'Chromosome'], '1')

    def test_sample_summary_handles_no_routes_and_no_gains(self):
        empty = posteriortablegen.get_sample_posterior_table_summary(
            pd.DataFrame(columns=GAIN_DRAW_COLUMNS),
            pd.DataFrame(columns=ROUTE_DRAW_COLUMNS),
        )
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), SUMMARY_COLUMNS)

        route_only = pd.DataFrame([
            make_route_draw(posterior_index=0),
        ], columns=ROUTE_DRAW_COLUMNS)
        empty = posteriortablegen.get_sample_posterior_table_summary(
            pd.DataFrame(columns=GAIN_DRAW_COLUMNS),
            route_only,
        )
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), SUMMARY_COLUMNS)

    def test_sample_summary_rejects_invalid_threshold_and_metadata(self):
        gain_draws, route_draws = self.valid_tables()
        for value in (-0.1, 1.1, np.nan, True, '0.8'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                    posteriortablegen.get_sample_posterior_table_summary(
                        gain_draws,
                        route_draws,
                        min_proportion_threshold=value,
                    )

        route_draws.loc[1, 'Chromosome'] = '2'
        with self.assertRaisesRegex(ValueError, 'consistent metadata'):
            posteriortablegen.get_sample_posterior_table_summary(
                gain_draws,
                route_draws,
            )


if __name__ == '__main__':
    unittest.main()
