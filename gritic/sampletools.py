import warnings
import unicodedata
from numbers import Integral, Real

import numpy as np
import pandas as pd
from scipy.stats import binom, poisson

from numba import njit,prange

import gritic.multiplicityoptimiser as multiplicityoptimiser
import gritic.dataloader as dataloader


DEFAULT_MIN_MUTATION_ALT_COUNT = 3
DEFAULT_MIN_MUTATION_COVERAGE = 10
DEFAULT_COVERAGE_VAF_QUANTILE = 0.9
DEFAULT_MIN_SUBCLONE_CCF = 0.01
DEFAULT_MAX_SUBCLONE_CCF = 0.9
DEFAULT_MIN_SUBCLONE_FRACTION = 0.1
DEFAULT_AUTOSOME_COUNT = 22

PHASE_GROUP_ORDER = {
    'non_phased': 0,
    'major': 1,
    'minor': 2,
}

_VALIDATED_INPUT_TABLES = object()

_WINDOWS_FORBIDDEN_FILENAME_CHARACTERS = frozenset('<>:"/\\|?*')
_WINDOWS_RESERVED_DEVICE_NAMES = frozenset({
    'CON',
    'PRN',
    'AUX',
    'NUL',
    'CLOCK$',
    'CONIN$',
    'CONOUT$',
    *(f'COM{number}' for number in range(1, 10)),
    *(f'LPT{number}' for number in range(1, 10)),
    'COM¹',
    'COM²',
    'COM³',
    'LPT¹',
    'LPT²',
    'LPT³',
})
_MAX_PATH_COMPONENT_UNITS = 255
_LONGEST_SAMPLE_ID_OUTPUT_SUFFIX = (
    '_posterior_timing_table_summary_penalty_False.tsv'
)


def _validate_non_negative_integer(value, parameter_name):
    if (
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or value < 0
    ):
        raise ValueError(f'{parameter_name} must be a non-negative integer')
    return int(value)


def _validate_unit_interval_number(value, parameter_name):
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or not 0 <= value <= 1
    ):
        raise ValueError(
            f'{parameter_name} must be a finite number between 0 and 1'
        )
    return float(value)


def validate_min_subclone_ccf(value):
    """Return a finite lower subclone-CCF bound in ``(0, 1]``."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or not 0 < value <= 1
    ):
        raise ValueError(
            'min_subclone_ccf must be a finite number greater than 0 and '
            'at most 1'
        )
    return float(value)


def validate_coverage_vaf_quantile(value):
    """Return a finite VAF quantile in the NumPy-supported ``[0, 1]``."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or not 0 <= value <= 1
    ):
        raise ValueError(
            'coverage_vaf_quantile must be a finite number between 0 and 1'
        )
    return float(value)


def validate_purity(value):
    """Return a finite purity in the biologically meaningful ``(0, 1]``."""
    return dataloader.validate_purity(value)


def frequency_weighted_quantile(values, frequencies, quantile):
    """Return NumPy's linear quantile of an integer-frequency sample.

    This is equivalent to ``np.quantile(np.repeat(values, frequencies), q)``
    without materializing the expanded sample.  GRITIC uses it for the
    high-VAF coverage estimate, where count groups carry the number of SNVs
    represented by each distinct read-count pair.
    """
    values = np.asarray(values, dtype=np.float64)
    frequencies = np.asarray(frequencies, dtype=np.int64)
    quantile = validate_coverage_vaf_quantile(quantile)
    if values.ndim != 1 or frequencies.ndim != 1:
        raise ValueError('values and frequencies must be one-dimensional')
    if values.size == 0 or values.size != frequencies.size:
        raise ValueError(
            'values and frequencies must have the same non-zero length'
        )
    if np.any(frequencies <= 0):
        raise ValueError('frequencies must contain positive integers')

    order = np.argsort(values, kind='stable')
    ordered_values = values[order]
    cumulative_frequencies = np.cumsum(
        frequencies[order],
        dtype=np.int64,
    )
    expanded_size = int(cumulative_frequencies[-1])
    position = (expanded_size - 1) * quantile
    lower_position = int(np.floor(position))
    upper_position = int(np.ceil(position))
    interpolation_weight = position - lower_position
    lower_index = np.searchsorted(
        cumulative_frequencies,
        lower_position,
        side='right',
    )
    upper_index = np.searchsorted(
        cumulative_frequencies,
        upper_position,
        side='right',
    )
    lower_value = ordered_values[lower_index]
    upper_value = ordered_values[upper_index]

    # Match NumPy's scalar linear interpolation, including its numerically
    # stable upper-end expression for weights at or above one half.
    difference = upper_value - lower_value
    if interpolation_weight >= 0.5:
        return upper_value - difference * (1.0 - interpolation_weight)
    return lower_value + difference * interpolation_weight


def get_normal_total_copy_number(chromosome, sex):
    """Return normal-cell copy number for an autosome or sex chromosome."""
    return dataloader.get_normal_total_copy_number(chromosome, sex)


def _canonical_phase_labels(phasing):
    """Return the compact phase labels used by the grouping tables."""
    return phasing.where(phasing.notna(), 'non_phased').astype(str)


def build_sample_group_tables(mutation_table):
    """Build sample-wide count/phase dictionaries and segment weights.

    Count groups are unique read-count pairs across the sample. Phase groups
    are unique ``(Count_Group_ID, Phasing)`` pairs across the sample. Only the
    sparse segment/phase-group association carries mutation frequencies.
    """
    if mutation_table.empty:
        raise ValueError('Cannot group an empty mutation table')

    mutation_table = mutation_table.drop(
        columns=['Count_Group_ID', 'Phase_Group_ID'],
        errors='ignore',
    ).copy()
    if 'Phasing' not in mutation_table.columns:
        mutation_table['Phasing'] = np.nan

    count_columns = ['Tumor_Ref_Count', 'Tumor_Alt_Count']
    count_group_table = (
        mutation_table[count_columns]
        .drop_duplicates()
        .sort_values(count_columns, kind='mergesort')
        .reset_index(drop=True)
    )
    count_group_table.insert(
        0,
        'Count_Group_ID',
        np.arange(len(count_group_table), dtype=np.int64),
    )

    count_lookup = count_group_table.set_index(count_columns)[
        'Count_Group_ID'
    ]
    mutation_count_index = pd.MultiIndex.from_frame(
        mutation_table[count_columns]
    )
    mutation_count_groups = count_lookup.reindex(
        mutation_count_index
    ).to_numpy(dtype=np.int64)

    phase_rows = pd.DataFrame({
        'Count_Group_ID': mutation_count_groups,
        'Phasing': _canonical_phase_labels(
            mutation_table['Phasing']
        ).to_numpy(),
    })
    phase_group_table = phase_rows.drop_duplicates().copy()
    phase_group_table['_Phase_Order'] = (
        phase_group_table['Phasing'].map(PHASE_GROUP_ORDER).fillna(3)
    )
    phase_group_table = (
        phase_group_table
        .sort_values(
            ['Count_Group_ID', '_Phase_Order', 'Phasing'],
            kind='mergesort',
        )
        .drop(columns='_Phase_Order')
        .reset_index(drop=True)
    )
    phase_group_table.insert(
        0,
        'Phase_Group_ID',
        np.arange(len(phase_group_table), dtype=np.int64),
    )

    phase_lookup = phase_group_table.set_index(
        ['Count_Group_ID', 'Phasing']
    )['Phase_Group_ID']
    mutation_phase_index = pd.MultiIndex.from_frame(phase_rows)
    mutation_table['Phase_Group_ID'] = phase_lookup.reindex(
        mutation_phase_index
    ).to_numpy(dtype=np.int64)

    segment_group_table = (
        mutation_table
        .groupby(
            ['Segment_ID', 'Phase_Group_ID'],
            sort=True,
            as_index=False,
        )
        .size()
        .rename(columns={'size': 'N_Mutations'})
        .sort_values(
            ['Segment_ID', 'Phase_Group_ID'],
            kind='mergesort',
        )
        .reset_index(drop=True)
    )
    segment_group_table['Phase_Group_ID'] = segment_group_table[
        'Phase_Group_ID'
    ].astype(np.int64)
    segment_group_table['N_Mutations'] = segment_group_table[
        'N_Mutations'
    ].astype(np.int64)

    if int(segment_group_table['N_Mutations'].sum()) != len(mutation_table):
        raise AssertionError('Segment-group weights do not cover all mutations')

    return (
        mutation_table,
        count_group_table,
        phase_group_table,
        segment_group_table,
    )


def build_likelihood_context_tables(mutation_table, sex):
    """Build deterministic observation contexts and a segment/context map."""
    segment_columns = [
        'Segment_ID',
        'Chromosome',
        'Major_CN',
        'Minor_CN',
    ]
    segment_rows = mutation_table[segment_columns].drop_duplicates()
    duplicate_segments = segment_rows.duplicated('Segment_ID', keep=False)
    if duplicate_segments.any():
        raise ValueError(
            'Chromosome and copy number must be unique within each segment'
        )
    segment_rows = segment_rows.copy()
    segment_rows['Normal_Total_CN'] = [
        get_normal_total_copy_number(chromosome, sex)
        for chromosome in segment_rows['Chromosome']
    ]

    context_columns = ['Major_CN', 'Minor_CN', 'Normal_Total_CN']
    likelihood_context_table = (
        segment_rows[context_columns]
        .drop_duplicates()
        .sort_values(context_columns, kind='mergesort')
        .reset_index(drop=True)
    )
    likelihood_context_table.insert(
        0,
        'Likelihood_Context_ID',
        np.arange(len(likelihood_context_table), dtype=np.int64),
    )
    segment_context_table = (
        segment_rows[['Segment_ID'] + context_columns]
        .merge(
            likelihood_context_table,
            on=context_columns,
            how='left',
            sort=False,
            validate='many_to_one',
        )[['Segment_ID', 'Likelihood_Context_ID']]
        .sort_values('Segment_ID', kind='mergesort')
        .reset_index(drop=True)
    )
    segment_context_table['Likelihood_Context_ID'] = segment_context_table[
        'Likelihood_Context_ID'
    ].astype(np.int64)
    return likelihood_context_table, segment_context_table


def _normalised_count_likelihoods(
    ref_counts,
    alt_counts,
    state_vafs,
):
    """Return the legacy-normalized read likelihood for each state/group."""
    total_counts = ref_counts + alt_counts
    log_probabilities = binom.logpmf(
        alt_counts[None, :],
        total_counts[None, :],
        state_vafs[:, None],
    )
    probabilities = np.exp(
        log_probabilities - np.max(log_probabilities, axis=0)
    )
    normalising_sums = probabilities.sum(axis=0)
    normalising_sums = np.where(
        np.isclose(normalising_sums, 0),
        1,
        normalising_sums,
    )
    return probabilities / normalising_sums


def build_count_group_likelihood_table(
    count_group_table,
    phase_group_table,
    segment_group_table,
    likelihood_context_table,
    segment_context_table,
    purity,
    subclone_table,
):
    """Calculate each observed ``(context, count-group)`` likelihood once."""
    observed_pairs = (
        segment_group_table[['Segment_ID', 'Phase_Group_ID']]
        .merge(
            phase_group_table,
            on='Phase_Group_ID',
            how='left',
            sort=False,
            validate='many_to_one',
        )
        .merge(
            segment_context_table,
            on='Segment_ID',
            how='left',
            sort=False,
            validate='many_to_one',
        )[['Likelihood_Context_ID', 'Count_Group_ID']]
        .drop_duplicates()
        .sort_values(
            ['Likelihood_Context_ID', 'Count_Group_ID'],
            kind='mergesort',
        )
        .reset_index(drop=True)
    )
    if observed_pairs.isna().any().any():
        raise AssertionError('Grouping tables contain an unresolved key')

    max_major_cn = int(likelihood_context_table['Major_CN'].max())
    n_subclones = 0 if subclone_table is None else len(subclone_table)
    probability_columns = [
        f'Prob_Mult_{multiplicity}'
        for multiplicity in range(1, max_major_cn + 1)
    ]
    probability_columns.extend(
        f'Prob_Subclone_{index}' for index in range(n_subclones)
    )
    probability_values = np.full(
        (len(observed_pairs), len(probability_columns)),
        np.nan,
        dtype=np.float64,
    )
    count_lookup = count_group_table.set_index('Count_Group_ID')

    for context in likelihood_context_table.itertuples(index=False):
        context_mask = observed_pairs['Likelihood_Context_ID'].eq(
            context.Likelihood_Context_ID
        )
        if not context_mask.any():
            continue
        context_count_ids = observed_pairs.loc[
            context_mask,
            'Count_Group_ID',
        ].to_numpy(dtype=np.int64)
        context_counts = count_lookup.loc[context_count_ids]
        ref_counts = context_counts['Tumor_Ref_Count'].to_numpy(
            dtype=np.int64
        )
        alt_counts = context_counts['Tumor_Alt_Count'].to_numpy(
            dtype=np.int64
        )
        total_cn = context.Major_CN + context.Minor_CN
        mult_one_vaf = purity / (
            total_cn * purity
            + context.Normal_Total_CN * (1 - purity)
        )
        clonal_multiplicities = np.arange(
            1,
            context.Major_CN + 1,
            dtype=np.float64,
        )
        if subclone_table is None:
            state_multiplicities = clonal_multiplicities
        else:
            state_multiplicities = np.concatenate((
                clonal_multiplicities,
                subclone_table['Subclone_CCF'].to_numpy(dtype=np.float64),
            ))
        state_vafs = np.minimum(1, mult_one_vaf * state_multiplicities)
        context_probabilities = _normalised_count_likelihoods(
            ref_counts,
            alt_counts,
            state_vafs,
        ).T
        context_row_indices = np.flatnonzero(context_mask.to_numpy())
        probability_values[
            context_row_indices,
            :context.Major_CN,
        ] = context_probabilities[:, :context.Major_CN]
        if n_subclones:
            probability_values[
                context_row_indices,
                max_major_cn:,
            ] = context_probabilities[:, context.Major_CN:]

    probability_table = pd.DataFrame(
        probability_values,
        columns=probability_columns,
    )
    return pd.concat(
        [observed_pairs, probability_table],
        axis=1,
    )


def validate_sample_id(sample_id):
    """Validate that ``sample_id`` is a portable filename component.

    GRITIC uses the ID both as a directory name and as a filename prefix. The
    accepted form is therefore one non-empty path component that is safe under
    common POSIX and Windows filename rules. The ID must not be ``.`` or ``..``;
    contain a path separator, a Windows-forbidden filename character, or a
    Unicode control, format, or surrogate character; end in a dot or space; or
    have a Windows device-name stem. Its longest derived GRITIC filename must
    fit both the usual 255-byte POSIX component limit and the
    255-UTF-16-code-unit Windows component limit. The value is validated, never
    normalized or sanitized.
    """
    if not isinstance(sample_id, str):
        raise ValueError('Sample_ID must be a string')
    if not sample_id:
        raise ValueError('Sample_ID must not be empty')
    if sample_id in {'.', '..'}:
        raise ValueError("Sample_ID must not be '.' or '..'")

    forbidden_characters = sorted(
        set(sample_id) & _WINDOWS_FORBIDDEN_FILENAME_CHARACTERS
    )
    if forbidden_characters:
        raise ValueError(
            'Sample_ID contains a path separator or Windows-forbidden '
            f'filename character: {forbidden_characters[0]!r}'
        )

    for character in sample_id:
        if unicodedata.category(character) in {'Cc', 'Cf', 'Cs'}:
            raise ValueError(
                'Sample_ID must not contain Unicode control, format, or '
                'surrogate characters'
            )

    if sample_id.endswith(('.', ' ')):
        raise ValueError('Sample_ID must not end in a dot or space')

    # Windows reserves device names even when followed by an extension and
    # ignores spaces immediately before that extension during device lookup.
    device_stem = sample_id.partition('.')[0].rstrip(' ').upper()
    if device_stem in _WINDOWS_RESERVED_DEVICE_NAMES:
        raise ValueError(
            f'Sample_ID uses a reserved Windows device name: {device_stem!r}'
        )

    derived_filename = sample_id + _LONGEST_SAMPLE_ID_OUTPUT_SUFFIX
    try:
        utf8_bytes = len(derived_filename.encode('utf-8'))
        utf16_code_units = len(derived_filename.encode('utf-16-le')) // 2
    except UnicodeEncodeError as error:
        raise ValueError('Sample_ID must contain valid Unicode text') from error
    if (
        utf8_bytes > _MAX_PATH_COMPONENT_UNITS
        or utf16_code_units > _MAX_PATH_COMPONENT_UNITS
    ):
        raise ValueError(
            'Sample_ID is too long: GRITIC-derived filenames must be at most '
            f'{_MAX_PATH_COMPONENT_UNITS} UTF-8 bytes and UTF-16 code units'
        )

    return sample_id


def get_major_cn_mode_from_cn_table(cn_table, *, _validated=False):
    if _validated:
        cn_table = cn_table.copy()
    else:
        cn_table = dataloader.validate_segment_coordinates(cn_table)
        cn_table = dataloader.validate_non_overlapping_segments(cn_table)
    cn_table['Segment_Width'] = (
        dataloader.get_segment_widths(cn_table)
    )
    major_cn_mode = cn_table.groupby('Major_CN')['Segment_Width'].sum().idxmax()
    return major_cn_mode

def get_major_cn_mode(sample):
    cn_table = sample.cn_table.copy()
    normalized_chromosomes = (
        cn_table['Chromosome']
        .astype(str)
        .str.replace(r'^chr', '', regex=True)
    )
    cn_table = cn_table.loc[
        normalized_chromosomes.isin(sample.autosomes)
    ].copy()
    if cn_table.empty:
        raise ValueError(
            'Cannot infer WGD status because the copy number table contains '
            'no segments on the configured numbered autosomes'
        )
    return get_major_cn_mode_from_cn_table(cn_table, _validated=True)

@njit(parallel=True, cache=True)
def log_likelihood_numba_parallel(mult_array, mult_states, group_weights):
    log_likelihood_store = np.zeros(mult_states.shape[0])
    for i in prange(mult_states.shape[0]):
        log_likelihood = 0.0
        for group_index in range(mult_array.shape[0]):
            mutation_likelihood = 0.0
            for multiplicity_index in range(mult_array.shape[1]):
                state_probability = mult_states[i, multiplicity_index]
                if state_probability < 0.0:
                    state_probability = 0.0
                elif state_probability > 1.0:
                    state_probability = 1.0
                mutation_likelihood += (
                    state_probability
                    * mult_array[group_index, multiplicity_index]
                )
            log_likelihood += (
                group_weights[group_index]
                * np.log(mutation_likelihood + 2.2e-300)
            )
        log_likelihood_store[i] = log_likelihood
    return log_likelihood_store

def evaluate_likelihood_array_numba(
    full_states,
    non_phased_array,
    reads_correction_array,
    group_weights,
    tolerance=1e-8,
):
    full_states = np.multiply(full_states, reads_correction_array)
    full_states = full_states / np.sum(
        full_states,
        axis=1,
        keepdims=True,
    )

    prob_sums = np.sum(non_phased_array, axis=1)
    probs_sums_good = (
        (prob_sums > 1.0 - tolerance)
        & (prob_sums < 1.0 + tolerance)
    )
    non_phased_array = np.ascontiguousarray(
        non_phased_array[probs_sums_good, :],
        dtype=np.float64,
    )
    group_weights = np.ascontiguousarray(
        group_weights[probs_sums_good],
        dtype=np.int64,
    )
    if non_phased_array.shape[0] ==0:
        return np.full(full_states.shape[0], np.nan)
    likelihood_array = log_likelihood_numba_parallel(
        non_phased_array,
        np.ascontiguousarray(full_states, dtype=np.float64),
        group_weights,
    )

    return likelihood_array

class MultProbabilityStore:
    def __init__(
        self,
        mult_array_store,
        reads_correction_store,
        weight_store,
        major_cn,
        minor_cn,
        n_subclones,
    ):
        def probability_array(name):
            value = mult_array_store[name]
            if value is None:
                return None
            return np.ascontiguousarray(value, dtype=np.float64)

        def integer_weights(name):
            value = weight_store[name]
            if value is None:
                return None
            return np.ascontiguousarray(value, dtype=np.int64)

        def correction_array(name):
            value = reads_correction_store[name]
            if value is None:
                return None
            return np.ascontiguousarray(value, dtype=np.float64)

        self.non_phased_array = probability_array('Non_Phased')
        self.major_array = probability_array('Major')
        self.minor_array = probability_array('Minor')
        self.combined_array = probability_array('All')

        self.non_phased_weights = integer_weights('Non_Phased')
        self.major_weights = integer_weights('Major')
        self.minor_weights = integer_weights('Minor')
        self.combined_weights = integer_weights('All')

        
        self.reads_correction_non_phased_array = correction_array(
            'Non_Phased'
        )
        self.reads_correction_major_array = correction_array('Major')
        self.reads_correction_minor_array = correction_array('Minor')
        self.reads_correction_combined_array = correction_array('All')

        self.use_non_phased = self.non_phased_array is not None and self.non_phased_array.shape[0]>0
        self.use_major = self.major_array is not None and self.major_array.shape[0]>0
        self.use_minor = self.minor_array is not None and self.minor_array.shape[0]>0

       

        if self.use_non_phased:

            assert self.non_phased_array.shape[1] == self.reads_correction_non_phased_array.size
            assert self.non_phased_array.shape[0] == self.non_phased_weights.size
            assert np.all(self.non_phased_weights > 0)
        if self.use_major:

            assert self.major_array.shape[1] == self.reads_correction_major_array.size
            assert self.major_array.shape[0] == self.major_weights.size
            assert np.all(self.major_weights > 0)
        if self.use_minor:

            assert self.minor_array.shape[1] == self.reads_correction_minor_array.size
            assert self.minor_array.shape[0] == self.minor_weights.size
            assert np.all(self.minor_weights > 0)

        if sum([self.use_non_phased,self.use_major,self.use_minor]) == 0:
            raise ValueError('No arrays provided')

        assert self.combined_array.shape[0] == self.combined_weights.size
        assert np.all(self.combined_weights > 0)
        self.n_mutations = int(np.sum(self.combined_weights))
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.n_subclones = n_subclones
    
    def __str__(self):
        return str(np.average(
            self.combined_array,
            axis=0,
            weights=self.combined_weights,
        ))
    
    def get_subclonal_correction_array(self,subclone_table):
        weights = []
        arrays = []
        if self.use_non_phased:
            non_phased_array = np.concatenate([[self.reads_correction_non_phased_array[0]],self.reads_correction_non_phased_array[-(len(subclone_table.index)):]])
            arrays.append(non_phased_array)
            weights.append(np.sum(self.non_phased_weights))
        if self.use_major:
            major_array = np.concatenate([[self.reads_correction_major_array[0]],self.reads_correction_major_array[-(len(subclone_table.index)):]])
            arrays.append(major_array)
            weights.append(np.sum(self.major_weights))
        if self.use_minor:
            minor_array = np.concatenate([[self.reads_correction_minor_array[0]],self.reads_correction_minor_array[-(len(subclone_table.index)):]])
            arrays.append(minor_array)
            weights.append(np.sum(self.minor_weights))
        arrays = np.stack(arrays)
        weights = np.array(weights)
        weights = weights/np.sum(weights)
        return np.average(arrays,axis=0,weights=weights)

    def array_likelihood(self, mult_state, array, group_weights):
        likelihood = np.multiply(mult_state, array)
        
        likelihood = np.sum(likelihood, axis=1)
        likelihood = np.sum(
            group_weights * np.log(likelihood + 2.2e-300)
        )
        return likelihood
    def evaluate_likelihood(self,mult_state):
        mult_state = np.clip(mult_state,0.0,1.0)
        mult_state = mult_state/np.sum(mult_state)
        
        log_likelihood = self.array_likelihood(
            mult_state,
            self.combined_array,
            self.combined_weights,
        )

        return log_likelihood
    
    def evaluate_likelihood_array(self,full_states):
        #to do - implement phasing
        #added three reads correction
        subclonal_indicies = np.arange(full_states.shape[1]-self.n_subclones,full_states.shape[1])
        non_phased_indicies = np.arange(0,self.major_cn)
        major_indicies = np.arange(self.major_cn,2*self.major_cn)
        minor_indicies = np.arange(2*self.major_cn,2*self.major_cn+self.minor_cn)
        assert np.all(np.concatenate((non_phased_indicies,major_indicies,minor_indicies,subclonal_indicies)) == np.arange(full_states.shape[1]))
        

        
        non_phased_indicies = np.concatenate((non_phased_indicies,subclonal_indicies))
        major_indicies = np.concatenate((major_indicies,subclonal_indicies))
        minor_indicies = np.concatenate((minor_indicies,subclonal_indicies))

        ll = np.zeros(full_states.shape[0])
        
        if self.use_non_phased:
            full_states_non_phased = full_states[:,non_phased_indicies]
            ll_non_phased = evaluate_likelihood_array_numba(
                full_states_non_phased,
                self.non_phased_array,
                self.reads_correction_non_phased_array,
                self.non_phased_weights,
            )
            
            ll += ll_non_phased
        if self.use_major:
            full_states_major = full_states[:,major_indicies]
            ll_major = evaluate_likelihood_array_numba(
                full_states_major,
                self.major_array,
                self.reads_correction_major_array,
                self.major_weights,
            )
           
            ll += ll_major
        if self.use_minor:
            full_states_minor = full_states[:,minor_indicies]
            
            ll_minor = evaluate_likelihood_array_numba(
                full_states_minor,
                self.minor_array,
                self.reads_correction_minor_array,
                self.minor_weights,
            )
            ll += ll_minor
        return ll
    

class Sample:

    @classmethod
    def _from_validated_input_tables(cls, *args, **kwargs):
        return cls(
            *args,
            **kwargs,
            _validation_token=_VALIDATED_INPUT_TABLES,
        )
   
    def __init__(
        self,
        mutation_table,
        cn_table,
        subclone_table,
        sample_id,
        purity,
        sex=None,
        merge_cn=True,
        apply_reads_correction=True,
        drop_unmatched_snvs=False,
        drop_unmatched_chromosomes=False,
        min_mutation_alt_count=DEFAULT_MIN_MUTATION_ALT_COUNT,
        min_mutation_coverage=DEFAULT_MIN_MUTATION_COVERAGE,
        coverage_vaf_quantile=DEFAULT_COVERAGE_VAF_QUANTILE,
        min_subclone_ccf=DEFAULT_MIN_SUBCLONE_CCF,
        max_subclone_ccf=DEFAULT_MAX_SUBCLONE_CCF,
        min_subclone_fraction=DEFAULT_MIN_SUBCLONE_FRACTION,
        autosome_count=DEFAULT_AUTOSOME_COUNT,
        *,
        clip_subclone_ccf=False,
        drop_unrecognized_phasing=False,
        _validation_token=None,
    ):

        self.sample_id = validate_sample_id(sample_id)
        self.purity = validate_purity(purity)
        sex = dataloader.validate_sex_karyotype(sex)
        self.autosome_count = dataloader.validate_autosome_count(
            autosome_count,
        )
        self.autosomes = dataloader.get_autosome_labels(self.autosome_count)
        self.merge_cn = merge_cn
        self.apply_reads_correction = apply_reads_correction
        self.drop_unmatched_snvs = drop_unmatched_snvs
        self.drop_unmatched_chromosomes = drop_unmatched_chromosomes
        self.drop_unrecognized_phasing = drop_unrecognized_phasing
        self.min_mutation_alt_count = _validate_non_negative_integer(
            min_mutation_alt_count,
            'min_mutation_alt_count',
        )
        self.min_mutation_coverage = _validate_non_negative_integer(
            min_mutation_coverage,
            'min_mutation_coverage',
        )
        self.coverage_vaf_quantile = validate_coverage_vaf_quantile(
            coverage_vaf_quantile
        )
        self.min_subclone_ccf = validate_min_subclone_ccf(
            min_subclone_ccf,
        )
        self.max_subclone_ccf = _validate_unit_interval_number(
            max_subclone_ccf,
            'max_subclone_ccf',
        )
        if self.min_subclone_ccf > self.max_subclone_ccf:
            raise ValueError(
                'min_subclone_ccf must be less than or equal to '
                'max_subclone_ccf'
            )
        self.min_subclone_fraction = _validate_unit_interval_number(
            min_subclone_fraction,
            'min_subclone_fraction',
        )
        if not isinstance(clip_subclone_ccf, (bool, np.bool_)):
            raise ValueError('clip_subclone_ccf must be a boolean')
        self.clip_subclone_ccf = bool(clip_subclone_ccf)
        if (
            _validation_token is not None
            and _validation_token is not _VALIDATED_INPUT_TABLES
        ):
            raise ValueError('_validation_token is reserved for internal use')
        input_tables_validated = (
            _validation_token is _VALIDATED_INPUT_TABLES
        )
        self.supplied_segment_id_map = None
        if input_tables_validated:
            self.use_supplied_segment_ids = (
                'Segment_ID' in cn_table.columns
                and 'Segment_ID' in mutation_table.columns
            )
            self.sex = (
                dataloader.infer_sex_from_copy_number_table(cn_table)
                if sex is None
                else sex
            )
            cn_table = cn_table.copy()
            mutation_table = mutation_table.copy()
        else:
            self.use_supplied_segment_ids = (
                dataloader.validate_input_table_headers(
                    cn_table.columns,
                    mutation_table.columns,
                )
            )
            cn_table, mutation_table, self.sex = (
                dataloader.validate_or_drop_input_chromosomes(
                    cn_table,
                    mutation_table,
                    self.autosome_count,
                    sex,
                    drop_unmatched_chromosomes=(
                        self.drop_unmatched_chromosomes
                    ),
                )
            )
            mutation_table = dataloader.validate_or_drop_phasing_labels(
                mutation_table,
                drop_unrecognized_phasing=(
                    self.drop_unrecognized_phasing
                ),
            )
            cn_table = dataloader.validate_segment_coordinates(cn_table)
            cn_table = dataloader.validate_non_overlapping_segments(cn_table)
            if 'Position' in mutation_table.columns:
                mutation_table = mutation_table.copy()
                mutation_table['Position'] = dataloader.parse_positions(
                    mutation_table['Position']
                )
            mutation_table = dataloader.validate_mutation_read_counts(
                mutation_table
            )
            cn_table = dataloader.validate_copy_number_values(cn_table)
            if self.use_supplied_segment_ids:
                dataloader.validate_supplied_segment_ids(
                    cn_table,
                    mutation_table,
                    allow_unmatched=self.drop_unmatched_snvs,
                )
                dataloader.validate_source_scoped_mutation_components(
                    mutation_table['Segment_ID'],
                    dataloader.get_gritic_mutation_id_components(
                        mutation_table
                    ),
                )
                if self.drop_unmatched_snvs:
                    mutation_table = (
                        dataloader.drop_unmatched_segment_id_mutations(
                            cn_table,
                            mutation_table,
                        )
                    )
        self.chromosome_order = list(self.autosomes)
        self.chromosome_order.extend(
            dataloader.SEX_CHROMOSOME_PAIR[self.sex]
        )
        self.copy_number_table = self.process_raw_copy_number_table(
            cn_table,
            _validated=True,
        )
        self.mutation_table = self.process_raw_mutation_table(
            mutation_table,
            self.copy_number_table,
            _validated=True,
        )

        self.cn_table = cn_table.copy()
        self.subclone_table = self.format_subclone_table(subclone_table)
        (
            self.mutation_table,
            self.count_group_table,
            self.phase_group_table,
            self.segment_group_table,
        ) = build_sample_group_tables(self.mutation_table)
        (
            self.likelihood_context_table,
            self.segment_context_table,
        ) = build_likelihood_context_tables(
            self.mutation_table,
            self.sex,
        )
        self.count_group_likelihood_table = (
            build_count_group_likelihood_table(
                self.count_group_table,
                self.phase_group_table,
                self.segment_group_table,
                self.likelihood_context_table,
                self.segment_context_table,
                self.purity,
                self.subclone_table,
            )
        )
        self.segments = self.get_segments(_validated=True)

    def process_raw_mutation_table(
        self,
        mutation_table,
        cn_table,
        *,
        _validated=False,
    ):
        mutation_table = mutation_table.copy()
        
        mutation_table['Chromosome'] = (
            mutation_table['Chromosome'].str.replace(
                r'^chr', '', regex=True
            )
        )
        gritic_id_components = dataloader.get_gritic_mutation_id_components(
            mutation_table
        )
        mutation_table['_GRITIC_Mutation_ID_Component'] = gritic_id_components
        if 'Mutation_ID' in mutation_table.columns:
            mutation_table['Mutation_ID'] = gritic_id_components
        else:
            mutation_table['Mutation_ID'] = pd.NA
        if 'Position' not in mutation_table.columns:
            mutation_table['Position'] = pd.NA

        if self.use_supplied_segment_ids:
            mutation_table['Source_Segment_ID'] = (
                mutation_table['Segment_ID'].astype(str)
            )
            mutation_table['Segment_ID'] = (
                mutation_table['Segment_ID'].astype(str).map(
                    self.supplied_segment_id_map
                )
            )
            zero_copy_mutations = mutation_table['Segment_ID'].isin(
                self.zero_copy_segment_ids
            )
            zero_copy_mutation_count = int(zero_copy_mutations.sum())
            if zero_copy_mutation_count:
                zero_copy_segment_count = int(
                    mutation_table.loc[
                        zero_copy_mutations,
                        'Segment_ID',
                    ].nunique()
                )
                source_segment_ids = sorted(
                    mutation_table.loc[
                        zero_copy_mutations,
                        'Source_Segment_ID',
                    ].astype(str).unique()
                )
                snv_label = (
                    'SNV' if zero_copy_mutation_count == 1 else 'SNVs'
                )
                segment_label = (
                    'segment'
                    if zero_copy_segment_count == 1
                    else 'segments'
                )
                warnings.warn(
                    f'Dropping {zero_copy_mutation_count} {snv_label} '
                    f'assigned to {zero_copy_segment_count} zero-copy '
                    '(Major_CN=0, Minor_CN=0) copy-number '
                    f'{segment_label} because GRITIC does not model '
                    'mutations in 0+0 segments. Source Segment_ID value(s): '
                    + ', '.join(source_segment_ids)
                    + '.',
                    UserWarning,
                    stacklevel=2,
                )
                mutation_table = mutation_table.loc[
                    ~zero_copy_mutations
                ].copy()
                if mutation_table.empty:
                    raise ValueError(
                        'No mutations remain after dropping SNVs assigned to '
                        'zero-copy (0+0) segments'
                    )

        mutation_table = dataloader.assign_cn_to_snv(
            mutation_table,
            cn_table,
            use_supplied_segment_ids=self.use_supplied_segment_ids,
            drop_unmatched_snvs=self.drop_unmatched_snvs,
            _validated=_validated,
        )
        if not self.use_supplied_segment_ids:
            mutation_table['Source_Segment_ID'] = (
                mutation_table['Segment_ID'].astype(str)
            )
        gritic_id_components = mutation_table.pop(
            '_GRITIC_Mutation_ID_Component'
        )
        if not (_validated and self.use_supplied_segment_ids):
            dataloader.validate_source_scoped_mutation_components(
                mutation_table['Source_Segment_ID'],
                gritic_id_components,
            )
        mutation_table['GRITIC_Mutation_ID'] = (
            dataloader.generate_gritic_mutation_ids(
                mutation_table['Source_Segment_ID'],
                gritic_id_components,
            )
        )
        dataloader.validate_phasing_copy_number_consistency(mutation_table)
       
        mutation_table = self.filter_mutation_table(mutation_table)
     
        
        mutation_table = self.format_mutation_table(mutation_table)

        
        return mutation_table
    def process_raw_copy_number_table(
        self,
        cn_table,
        *,
        _validated=False,
    ):
        cn_table = cn_table.copy()
        cn_table['Chromosome'] = cn_table['Chromosome'].str.replace(
            r'^chr', '', regex=True
        )
        if self.merge_cn:
            if self.use_supplied_segment_ids:
                cn_table, self.supplied_segment_id_map = (
                    dataloader.merge_segments(
                        cn_table,
                        return_segment_id_map=True,
                        _validated=_validated,
                    )
                )
            else:
                cn_table = dataloader.merge_segments(
                    cn_table,
                    _validated=_validated,
                )
        else:
            if self.use_supplied_segment_ids:
                supplied_segment_ids = cn_table['Segment_ID'].astype(str)
            cn_table = dataloader.generate_segment_ids(cn_table)
            if self.use_supplied_segment_ids:
                self.supplied_segment_id_map = dict(zip(
                    supplied_segment_ids,
                    cn_table['Segment_ID'],
                ))
        cn_table['Total_CN'] = cn_table['Major_CN']+cn_table['Minor_CN']
        self.zero_copy_segment_ids = frozenset(
            cn_table.loc[
                cn_table['Major_CN'].eq(0),
                'Segment_ID',
            ].astype(str)
        )
        cn_table = cn_table[cn_table['Major_CN']>0]
        return cn_table
    def format_subclone_table(self,subclone_table):
        if subclone_table is None:
            return None
        
        subclone_table = dataloader.validate_subclone_values(
            subclone_table,
            clip_subclone_ccf=self.clip_subclone_ccf,
        )
        subclone_table = dataloader.get_valid_subclones(
            subclone_table,
            min_ccf=self.min_subclone_ccf,
            max_ccf=self.max_subclone_ccf,
            min_fraction=self.min_subclone_fraction,
        )
        if subclone_table.empty:
            return None

        
        subclone_table = dataloader.filter_excess_subclones(subclone_table)
        n_snvs = np.round(
            len(self.mutation_table.index)
            * subclone_table['Subclone_Fraction']
        ).astype(int)
        subclone_table['N_SNVs'] = n_snvs
        return subclone_table.loc[
            :,
            dataloader.SUBCLONE_OUTPUT_COLUMNS,
        ].copy()
    
    def format_mutation_table(self,mutation_table):
        mutation_table = mutation_table.reset_index(drop=True)
        if 'Phasing' not in mutation_table.columns:
            mutation_table['Phasing'] = np.nan
        if 'Context' not in mutation_table.columns:
            mutation_table['Context'] = np.nan
        canonical_order = mutation_table.sort_values(
            ['Segment_ID', 'GRITIC_Mutation_ID'],
            kind='mergesort',
        )
        canonical_indices = canonical_order.groupby(
            'Segment_ID',
            sort=False,
        ).cumcount()
        mutation_table.loc[
            canonical_order.index,
            'Segment_Mutation_Index',
        ] = canonical_indices.to_numpy()
        mutation_table['Segment_Mutation_Index'] = mutation_table[
            'Segment_Mutation_Index'
        ].astype(np.int64)
        mutation_columns = [
            'Segment_ID',
            'Source_Segment_ID',
            'Mutation_ID',
            'GRITIC_Mutation_ID',
            'Segment_Mutation_Index',
            'Chromosome',
            'Segment_Start',
            'Segment_End',
            'Major_CN',
            'Minor_CN',
            'Total_CN',
            'Tumor_Ref_Count',
            'Tumor_Alt_Count',
            'Position',
        ]
        mutation_columns.extend(['Phasing','Context'])
        mutation_table = mutation_table[mutation_columns].copy()
        mutation_table['Chromosome'] = mutation_table['Chromosome'].astype(str)
        mutation_table["Gain_Type"] = mutation_table["Major_CN"].astype(str)+ "_"+ mutation_table["Minor_CN"].astype(str)
        return mutation_table
    
    def filter_mutation_table(self,mutation_table):
        if mutation_table['Tumor_Alt_Count'].min() > 10:
            warnings.warn(
                'There are no mutations with less than 10 alt reads. This '
                'may indicate a higher threshold for mutation calling than '
                'the minimum alternate-read count configured for GRITIC '
                f'({self.min_mutation_alt_count})'
            )
        if len(mutation_table.index) ==0:
            raise ValueError("There are no mutations in the mutation table. Please check the mutation table is formatted correctly")
        mutation_table['Total_Count'] = mutation_table['Tumor_Ref_Count']+mutation_table['Tumor_Alt_Count']
        

        mutation_table = mutation_table[
            mutation_table['Tumor_Alt_Count']
            >= self.min_mutation_alt_count
        ]
        mutation_table = mutation_table[
            (
                mutation_table['Tumor_Ref_Count']
                + mutation_table['Tumor_Alt_Count']
            ) >= self.min_mutation_coverage
        ]
        if mutation_table.empty:
            raise ValueError(
                'No mutations remain after applying the configured minimum '
                'alternate-read count and coverage thresholds'
            )
        
        return mutation_table
    
    def get_segments(self,min_mutations=0, *, _validated=False):
        segments = []
        sample_group_tables = {
            'count_group_table': self.count_group_table,
            'phase_group_table': self.phase_group_table,
            'segment_group_table': self.segment_group_table,
            'likelihood_context_table': self.likelihood_context_table,
            'segment_context_table': self.segment_context_table,
            'count_group_likelihood_table': (
                self.count_group_likelihood_table
            ),
        }
        for _,segment_table in self.mutation_table.groupby('Segment_ID'):
            if len(segment_table.index) >= min_mutations:
                segment = Segment(
                    segment_table,
                    self.subclone_table,
                    self.purity,
                    self.sex,
                    apply_reads_correction=self.apply_reads_correction,
                    min_mutation_alt_count=self.min_mutation_alt_count,
                    coverage_vaf_quantile=self.coverage_vaf_quantile,
                    _validated=_validated,
                    _sample_group_tables=sample_group_tables,
                )
                segments.append(segment)
        return segments
    
    def get_mutation_table(self):
        mutation_table = pd.concat([segment.mutation_table for segment in self.segments])
        return self.sort_table_by_chromosome(mutation_table)

    def get_count_group_table(self):
        return self.count_group_table.copy()

    def get_phase_group_table(self):
        return self.phase_group_table.copy()

    def get_likelihood_context_table(self):
        return self.likelihood_context_table.copy()

    def get_segment_context_table(self):
        return self.segment_context_table.copy()

    def get_count_group_likelihood_table(self):
        return self.count_group_likelihood_table.copy()

    def get_segment_group_table(self):
        return self.segment_group_table.copy()
    
    def get_subclone_table(self):
        return self.subclone_table
    
    
    def sort_table_by_chromosome(self,table):
        table =table.copy()
        table['Chromosome'] = pd.Categorical(table['Chromosome'],categories=self.chromosome_order)
        return table.sort_values(by=['Chromosome','Segment_Start']).reset_index(drop=True)

class Segment:
    
    def __init__(
        self,
        mutation_table,
        subclone_table,
        purity,
        sex,
        apply_reads_correction=True,
        min_mutation_alt_count=DEFAULT_MIN_MUTATION_ALT_COUNT,
        coverage_vaf_quantile=DEFAULT_COVERAGE_VAF_QUANTILE,
        *,
        _validated=False,
        _sample_group_tables=None,
    ):
        if _validated:
            self.mutation_table = mutation_table.copy().reset_index(drop=True)
        else:
            self.mutation_table = dataloader.validate_mutation_read_counts(
                mutation_table
            )
            self.mutation_table = dataloader.validate_segment_coordinates(
                self.mutation_table
            ).reset_index(drop=True)
        if 'Phasing' not in self.mutation_table.columns:
            self.mutation_table['Phasing'] = np.nan
        self.n_mutations = len(self.mutation_table.index)
        if subclone_table is None:
            self.subclone_table = None
        elif _validated:
            self.subclone_table = subclone_table.copy()
        else:
            self.subclone_table = dataloader.validate_subclone_values(
                subclone_table
            )
        self.sex = sex
        self.min_mutation_alt_count = _validate_non_negative_integer(
            min_mutation_alt_count,
            'min_mutation_alt_count',
        )
        self.coverage_vaf_quantile = validate_coverage_vaf_quantile(
            coverage_vaf_quantile
        )
        self.sample_clone_fractions = self.get_sample_clone_fractions()
        self.n_subclones = self.get_n_subclones()
        self.apply_reads_correction = apply_reads_correction
        self.purity = validate_purity(purity)

        self.segment_id = self.get_unique_attribute_from_table('Segment_ID')

        
        
        self.major_cn = self.get_unique_attribute_from_table('Major_CN')
        self.minor_cn = self.get_unique_attribute_from_table('Minor_CN')
        self.total_cn = self.major_cn + self.minor_cn

        self.all_possible_clonal_multiplicities = np.arange(self.major_cn)+1
        self.all_possible_subclonal_multiplicities = self.get_all_possible_subclonal_multiplicities()
        
        self.chromosome = self.get_unique_attribute_from_table('Chromosome')
        self.start = self.get_unique_attribute_from_table('Segment_Start')
        self.end = self.get_unique_attribute_from_table('Segment_End')
        self.width = int(self.end) - int(self.start)

        self.assign_mutation_groups(_sample_group_tables)
        self.assign_multiplicity_probabilities()
        
        self.multiplicity_probabilities = self.get_multiplicity_probabilities()

        
    
    def __str__(self):
        if self.n_mutations ==1:
            mutation_text = 'Mutation'
        else:
            mutation_text ='Mutations'
        return f"{self.segment_id}- {self.major_cn}+{self.minor_cn} - {self.n_mutations} {mutation_text}"
    
    def get_all_possible_subclonal_multiplicities(self):
        if self.n_subclones ==0:
            return np.zeros(0)
        return self.subclone_table['Subclone_CCF'].to_numpy()
    
    def get_n_subclones(self):
        if self.subclone_table is None:
            return 0
        return len(self.subclone_table.index)
    
    def get_sample_clone_fractions(self):
        if self.subclone_table is None:
            return np.array([1])
        subclone_fractions = self.subclone_table['Subclone_Fraction'].to_numpy()
        total_subclone_fraction = np.sum(subclone_fractions)
        clone_fraction =  1 - total_subclone_fraction
        sample_clone_fractions = np.insert(subclone_fractions,0,clone_fraction)
        return sample_clone_fractions
    
    def get_unique_attribute_from_table(self,attribute):
        if len(self.mutation_table[attribute].unique()) >1:
            raise ValueError(f"Attribute {attribute} is not unique for segment")
        return self.mutation_table[attribute].iloc[0]

    
    def get_multiplicity_names(self):
        sample_peak_names = [f'Mult_{mult}' for mult in self.all_possible_clonal_multiplicities]
        for i in range(self.n_subclones):
            sample_peak_names.append(f"Subclone_{i}")
        return sample_peak_names
    
    def get_mutation_rate(self):
        mult_proportions = multiplicityoptimiser.unconstrained_mult_optimisation(self.multiplicity_probabilities,self.n_subclones)
        if mult_proportions is None:
            return np.nan
        mult_states = np.concatenate([np.arange(self.major_cn)+1,np.ones(self.n_subclones)])
        
        mult_states = np.divide(mult_states,self.multiplicity_probabilities.reads_correction_combined_array)
        mult_correction_factor = np.sum(np.multiply(mult_proportions,mult_states))
        
        mutation_rate = self.n_mutations*mult_correction_factor
        return mutation_rate/self.width

    def assign_mutation_groups(self, sample_group_tables=None):
        """Attach this segment to the sample-wide grouping dictionaries."""
        represented_sex_chromosomes = (
            set(self.mutation_table['Chromosome'].unique())
            & {'X', 'Y', 'Z', 'W'}
        )
        if represented_sex_chromosomes and self.sex is None:
            raise ValueError(
                'sex must be specified when timing sex-chromosome segments'
            )

        if sample_group_tables is None:
            (
                self.mutation_table,
                count_group_dictionary,
                phase_group_dictionary,
                segment_group_table,
            ) = build_sample_group_tables(self.mutation_table)
            (
                likelihood_context_table,
                segment_context_table,
            ) = build_likelihood_context_tables(
                self.mutation_table,
                self.sex,
            )
            count_group_likelihood_table = (
                build_count_group_likelihood_table(
                    count_group_dictionary,
                    phase_group_dictionary,
                    segment_group_table,
                    likelihood_context_table,
                    segment_context_table,
                    self.purity,
                    self.subclone_table,
                )
            )
        else:
            required_tables = {
                'count_group_table',
                'phase_group_table',
                'segment_group_table',
                'likelihood_context_table',
                'segment_context_table',
                'count_group_likelihood_table',
            }
            missing_tables = required_tables - set(sample_group_tables)
            if missing_tables:
                raise ValueError(
                    'Missing sample grouping tables: '
                    + ', '.join(sorted(missing_tables))
                )
            count_group_dictionary = sample_group_tables[
                'count_group_table'
            ]
            phase_group_dictionary = sample_group_tables[
                'phase_group_table'
            ]
            segment_group_table = sample_group_tables[
                'segment_group_table'
            ]
            likelihood_context_table = sample_group_tables[
                'likelihood_context_table'
            ]
            segment_context_table = sample_group_tables[
                'segment_context_table'
            ]
            count_group_likelihood_table = sample_group_tables[
                'count_group_likelihood_table'
            ]

        segment_context = segment_context_table.loc[
            segment_context_table['Segment_ID'].eq(self.segment_id)
        ]
        if len(segment_context) != 1:
            raise ValueError(
                f'Expected one likelihood context for {self.segment_id}'
            )
        self.likelihood_context_id = int(
            segment_context.iloc[0]['Likelihood_Context_ID']
        )
        context_row = likelihood_context_table.loc[
            likelihood_context_table['Likelihood_Context_ID'].eq(
                self.likelihood_context_id
            )
        ]
        if len(context_row) != 1:
            raise ValueError(
                f'Likelihood context {self.likelihood_context_id} is missing'
            )
        self.normal_total_cn = int(context_row.iloc[0]['Normal_Total_CN'])
        if (
            int(context_row.iloc[0]['Major_CN']) != self.major_cn
            or int(context_row.iloc[0]['Minor_CN']) != self.minor_cn
        ):
            raise ValueError('Likelihood context copy number does not match')

        self.segment_group_table = (
            segment_group_table.loc[
                segment_group_table['Segment_ID'].eq(self.segment_id)
            ]
            .sort_values('Phase_Group_ID', kind='mergesort')
            .reset_index(drop=True)
            .copy()
        )
        if self.segment_group_table.empty:
            raise ValueError(f'No mutation groups exist for {self.segment_id}')

        phase_group_table = self.segment_group_table.merge(
            phase_group_dictionary,
            on='Phase_Group_ID',
            how='left',
            sort=False,
            validate='many_to_one',
        )
        if phase_group_table[['Count_Group_ID', 'Phasing']].isna().any().any():
            raise AssertionError('Segment references an unknown phase group')
        self.phase_group_table = phase_group_table[[
            'Segment_ID',
            'Phase_Group_ID',
            'Count_Group_ID',
            'Phasing',
            'N_Mutations',
        ]].copy()

        count_weights = (
            self.phase_group_table
            .groupby('Count_Group_ID', sort=True, as_index=False)[
                'N_Mutations'
            ]
            .sum()
        )
        segment_count_groups = count_group_dictionary.merge(
            count_weights,
            on='Count_Group_ID',
            how='inner',
            sort=False,
            validate='one_to_one',
        )
        context_likelihoods = count_group_likelihood_table.loc[
            count_group_likelihood_table['Likelihood_Context_ID'].eq(
                self.likelihood_context_id
            )
        ].drop(columns='Likelihood_Context_ID')
        segment_count_groups = segment_count_groups.merge(
            context_likelihoods,
            on='Count_Group_ID',
            how='left',
            sort=False,
            validate='one_to_one',
        ).sort_values('Count_Group_ID', kind='mergesort').reset_index(
            drop=True
        )
        probability_columns = [
            column
            for column in segment_count_groups.columns
            if column.startswith('Prob_')
        ]
        required_probability_columns = [
            f'Prob_Mult_{multiplicity}'
            for multiplicity in range(1, self.major_cn + 1)
        ]
        required_probability_columns.extend(
            f'Prob_Subclone_{index}'
            for index in range(self.n_subclones)
        )
        missing_probabilities = (
            set(required_probability_columns) - set(probability_columns)
        )
        if missing_probabilities:
            raise ValueError(
                'Missing count-group probabilities: '
                + ', '.join(sorted(missing_probabilities))
            )
        if segment_count_groups[required_probability_columns].isna().any().any():
            raise ValueError('Count-group likelihood lookup is incomplete')
        segment_count_groups.insert(0, 'Segment_ID', self.segment_id)
        self.count_group_table = segment_count_groups

        observed_phase_counts = self.mutation_table[
            'Phase_Group_ID'
        ].value_counts(sort=False).sort_index()
        expected_phase_counts = self.segment_group_table.set_index(
            'Phase_Group_ID'
        )['N_Mutations'].sort_index()
        if not observed_phase_counts.equals(expected_phase_counts):
            raise AssertionError('Mutation/segment-group mapping is inconsistent')
        if int(self.count_group_table['N_Mutations'].sum()) != self.n_mutations:
            raise AssertionError('Count-group weights do not cover the segment')
        if int(self.phase_group_table['N_Mutations'].sum()) != self.n_mutations:
            raise AssertionError('Phase-group weights do not cover the segment')

    def get_count_group_table(self):
        return self.count_group_table.copy()

    def get_phase_group_table(self):
        return self.phase_group_table.copy()

    def assign_multiplicity_probabilities(self):
        """Calculate segment-specific coverage and detection corrections."""
        mult_one_vaf = self.purity / (
            self.total_cn * self.purity
            + self.normal_total_cn * (1 - self.purity)
        )
        combined_all_possible_multiplicities = np.concatenate((
            self.all_possible_clonal_multiplicities,
            self.all_possible_subclonal_multiplicities,
        ))
        ref_counts = self.count_group_table[
            'Tumor_Ref_Count'
        ].to_numpy(dtype=np.int64)
        alt_counts = self.count_group_table[
            'Tumor_Alt_Count'
        ].to_numpy(dtype=np.int64)
        total_counts = ref_counts + alt_counts
        group_weights = self.count_group_table[
            'N_Mutations'
        ].to_numpy(dtype=np.int64)
        vaf = alt_counts / total_counts
        vaf_quantile = frequency_weighted_quantile(
            vaf,
            group_weights,
            self.coverage_vaf_quantile,
        )
        high_vaf_groups = vaf > vaf_quantile - 0.01
        highest_vaf_average_coverage = np.average(
            total_counts[high_vaf_groups],
            weights=group_weights[high_vaf_groups],
        )

        if not self.apply_reads_correction:
            highest_vaf_average_coverage = np.average(
                total_counts,
                weights=group_weights,
            )
        self.highest_vaf_average_coverage = highest_vaf_average_coverage

        alt_count_correction_factors = []
        for multiplicity in combined_all_possible_multiplicities:
            mult_vaf = np.minimum(1, mult_one_vaf * multiplicity)
            # If total depth is Poisson and alternate reads are a binomial
            # thinning of that depth, the alternate-read count is Poisson
            # with mean coverage * VAF. Evaluate the detection probability
            # exactly instead of drawing pseudo-depths.
            alt_count_correction_factor = poisson.sf(
                self.min_mutation_alt_count - 1,
                highest_vaf_average_coverage * mult_vaf,
            )
            alt_count_correction_factors.append(
                alt_count_correction_factor
            )
        self.alt_count_correction_factors = np.array(
            alt_count_correction_factors
        )
        
        
    def get_reads_correction_array(self, allele=None):
        group_weights = self.get_multiplicity_probability_weights(allele)
        if group_weights is None:
            return None
        if allele == 'minor':
            clonal_correction = self.alt_count_correction_factors[
                :self.minor_cn
            ]
        else:
            clonal_correction = self.alt_count_correction_factors[
                :self.major_cn
            ]
        if self.n_subclones:
            subclonal_correction = self.alt_count_correction_factors[
                -self.n_subclones:
            ]
            reads_correction = np.concatenate(
                (clonal_correction, subclonal_correction)
            )
        else:
            reads_correction = clonal_correction.copy()
        if not self.apply_reads_correction:
            reads_correction = np.ones_like(reads_correction)
        
        
        return reads_correction
    
    def get_multiplicity_probabilities_array(self, allele=None):

        if allele == 'minor':
            clonal_multiplicities = np.arange(1, self.minor_cn + 1)
        else:
            clonal_multiplicities = np.arange(1, self.major_cn + 1)
        
        mult_names = [f"Prob_Mult_{mult}" for mult in clonal_multiplicities]
        
        mult_names.extend([
            f"Prob_Subclone_{subclone}"
            for subclone in range(self.n_subclones)
        ])

        if allele == 'all':
            source_table = self.count_group_table
        else:
            phase_label = 'non_phased' if allele is None else allele
            phase_groups = self.phase_group_table.loc[
                self.phase_group_table['Phasing'].eq(phase_label)
            ].sort_values('Phase_Group_ID', kind='mergesort')
            if phase_groups.empty:
                return None
            count_probabilities = self.count_group_table.set_index(
                'Count_Group_ID'
            )
            source_table = count_probabilities.loc[
                phase_groups['Count_Group_ID'].to_numpy(dtype=np.int64)
            ]

        multiplicity_probabilities = source_table[
            mult_names
        ].to_numpy(dtype=np.float64)

        normalising_sums = np.sum(
            multiplicity_probabilities,
            axis=1,
        )[:, None]
        normalising_sums = np.where(
            np.isclose(normalising_sums, 0),
            1,
            normalising_sums,
        )
        multiplicity_probabilities = multiplicity_probabilities / normalising_sums

        return multiplicity_probabilities

    def get_multiplicity_probability_weights(self, allele=None):
        if allele == 'all':
            return self.count_group_table[
                'N_Mutations'
            ].to_numpy(dtype=np.int64)
        phase_label = 'non_phased' if allele is None else allele
        weights = self.phase_group_table.loc[
            self.phase_group_table['Phasing'].eq(phase_label),
            'N_Mutations',
        ].to_numpy(dtype=np.int64)
        if weights.size == 0:
            return None
        return weights
    

    def get_multiplicity_probabilities(self):
        reads_correction_store = {
            'Non_Phased': self.get_reads_correction_array(),
            'Major': self.get_reads_correction_array('major'),
            'Minor': self.get_reads_correction_array('minor'),
            'All': self.get_reads_correction_array('all'),
        }

        mult_array_store = {
            'Non_Phased': self.get_multiplicity_probabilities_array(),
            'Major': self.get_multiplicity_probabilities_array('major'),
            'Minor': self.get_multiplicity_probabilities_array('minor'),
            'All': self.get_multiplicity_probabilities_array('all'),
        }

        weight_store = {
            'Non_Phased': self.get_multiplicity_probability_weights(),
            'Major': self.get_multiplicity_probability_weights('major'),
            'Minor': self.get_multiplicity_probability_weights('minor'),
            'All': self.get_multiplicity_probability_weights('all'),
        }
        
        return MultProbabilityStore(
            mult_array_store,
            reads_correction_store,
            weight_store,
            self.major_cn,
            self.minor_cn,
            self.n_subclones,
        )
    
    def get_info_dict(self):
        
        info_dict = {'Segment_ID':self.segment_id,'Chromosome':self.chromosome,'Segment_Start':self.start,'Segment_End':self.end}
        info_dict['Major_CN'] = self.major_cn
        info_dict['Minor_CN'] = self.minor_cn
        info_dict['Total_CN'] = self.total_cn
        info_dict['N_Mutations']=self.n_mutations
        info_dict['Mutation_Rate'] = self.get_mutation_rate()
        return info_dict
    
  
   
