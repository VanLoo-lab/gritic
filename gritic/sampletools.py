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


def get_normal_total_copy_number(chromosome, sex):
    """Return normal-cell copy number for an autosome or sex chromosome."""
    return dataloader.get_normal_total_copy_number(chromosome, sex)


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

@njit(parallel=True)                           
def log_likelihood_numba_parallel(mult_array,mult_states):
    log_likelihood_store = np.zeros(mult_states.shape[0])                                
    for i in prange(mult_states.shape[0]):
        mult_state = mult_states[i]
        mult_state = np.clip(mult_state,0.0,1.0)
        log_likelihood = np.multiply(mult_state, mult_array)
    
        log_likelihood = np.sum(log_likelihood, axis=1)
        log_likelihood = np.sum(np.log(log_likelihood + 2.2e-300))
        log_likelihood_store[i] = log_likelihood
    return log_likelihood_store

@njit  
def evaluate_likelihood_array_numba(full_states,non_phased_array,reads_correction_array,tolerance=1e-8):
    
    full_states = np.multiply(full_states,reads_correction_array)
    full_states = full_states/np.sum(full_states,axis=1).reshape(-1,1)
    
    prob_sums = np.sum(non_phased_array,axis=1)
    probs_sums_good = (prob_sums>(1.0-tolerance)) & (prob_sums<(1.0+tolerance))

    non_phased_array = non_phased_array[probs_sums_good,:]
    if non_phased_array.shape[0] ==0:
        return np.nan*np.ones_like(full_states[:,0])
    likelihood_array = log_likelihood_numba_parallel(non_phased_array,full_states)

    return likelihood_array

class MultProbabilityStore:
    def __init__(self,mult_array_store,reads_correction_store,major_cn,minor_cn,n_subclones):
        self.non_phased_array = mult_array_store['Non_Phased']
        self.major_array = mult_array_store['Major']
        self.minor_array = mult_array_store['Minor']
        self.combined_array = mult_array_store['All']

        
        self.reads_correction_non_phased_array = reads_correction_store['Non_Phased']
        self.reads_correction_major_array = reads_correction_store['Major']
        self.reads_correction_minor_array = reads_correction_store['Minor']
        self.reads_correction_combined_array = reads_correction_store['All']

        self.use_non_phased = self.non_phased_array is not None and self.non_phased_array.shape[0]>0
        self.use_major = self.major_array is not None and self.major_array.shape[0]>0
        self.use_minor = self.minor_array is not None and self.minor_array.shape[0]>0

       

        if self.use_non_phased:

            assert self.non_phased_array.shape[1] == self.reads_correction_non_phased_array.size
        if self.use_major:

            assert self.major_array.shape[1] == self.reads_correction_major_array.size
        if self.use_minor:

            assert self.minor_array.shape[1] == self.reads_correction_minor_array.size

        if sum([self.use_non_phased,self.use_major,self.use_minor]) == 0:
            print(self.non_phased_array,self.major_array,self.minor_array)
            raise ValueError('No arrays provided')

        self.n_mutations = self.non_phased_array.shape[0] if self.use_non_phased else self.major_array.shape[0] if self.use_major else self.minor_array.shape[0]
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.n_subclones = n_subclones
    
    def __str__(self):
        return str(np.average(self.non_phased_array,axis=0))
    
    def get_subclonal_correction_array(self,subclone_table):
        weights = []
        arrays = []
        if self.use_non_phased:
            non_phased_array = np.concatenate([[self.reads_correction_non_phased_array[0]],self.reads_correction_non_phased_array[-(len(subclone_table.index)):]])
            arrays.append(non_phased_array)
            weights.append(self.non_phased_array.shape[0])
        if self.use_major:
            major_array = np.concatenate([[self.reads_correction_major_array[0]],self.reads_correction_major_array[-(len(subclone_table.index)):]])
            arrays.append(major_array)
            weights.append(self.major_array.shape[0])
        if self.use_minor:
            minor_array = np.concatenate([[self.reads_correction_minor_array[0]],self.reads_correction_minor_array[-(len(subclone_table.index)):]])
            arrays.append(minor_array)
            weights.append(self.minor_array.shape[0])
        arrays = np.stack(arrays)
        weights = np.array(weights)
        weights = weights/np.sum(weights)
        return np.average(arrays,axis=0,weights=weights)

    def array_likelihood(self,mult_state,array):
        likelihood = np.multiply(mult_state, array)
        
        likelihood = np.sum(likelihood, axis=1)
        likelihood = np.sum(np.log(likelihood + 2.2e-300))
        return likelihood
    def evaluate_likelihood(self,mult_state):
        mult_state = np.clip(mult_state,0.0,1.0)
        mult_state = mult_state/np.sum(mult_state)
        
        log_likelihood= self.array_likelihood(mult_state,self.combined_array)

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
            ll_non_phased = evaluate_likelihood_array_numba(full_states_non_phased,self.non_phased_array,self.reads_correction_non_phased_array)
            
            ll += ll_non_phased
        if self.use_major:
            full_states_major = full_states[:,major_indicies]
            ll_major = evaluate_likelihood_array_numba(full_states_major,self.major_array,self.reads_correction_major_array)
           
            ll += ll_major
        if self.use_minor:
            full_states_minor = full_states[:,minor_indicies]
            
            ll_minor = evaluate_likelihood_array_numba(full_states_minor,self.minor_array,self.reads_correction_minor_array)
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
        merge_cn=False,
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
                )
                segments.append(segment)
        return segments
    
    def get_mutation_table(self):
        mutation_table = pd.concat([segment.mutation_table for segment in self.segments])
        return self.sort_table_by_chromosome(mutation_table)
    
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
    ):
        if _validated:
            self.mutation_table = mutation_table.copy()
        else:
            self.mutation_table = dataloader.validate_mutation_read_counts(
                mutation_table
            )
            self.mutation_table = dataloader.validate_segment_coordinates(
                self.mutation_table
            )
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
    
    def assign_multiplicity_probabilities(self):
        #remove mult cols if they already exist
        mult_cols = [
            col
            for col in self.mutation_table.columns
            if 'Prob_' in col or 'Alt_Count_Correction_' in col
        ]
        self.mutation_table = self.mutation_table.drop(columns=mult_cols)
        represented_sex_chromosomes = (
            set(self.mutation_table['Chromosome'].unique())
            & {'X', 'Y', 'Z', 'W'}
        )
        if represented_sex_chromosomes and self.sex is None:
            raise ValueError(
                'sex must be specified when timing sex-chromosome segments'
            )
        normal_total_cn = self.mutation_table['Chromosome'].map(
            lambda chromosome: get_normal_total_copy_number(
                chromosome,
                self.sex,
            )
        )
        
        alt_count_correction_factors = []
        multiplicity_probabilities = []

        mult_one_vaf = self.purity / (self.total_cn* self.purity + normal_total_cn * (1 - self.purity))
        
        combined_all_possible_multiplicities = np.concatenate((self.all_possible_clonal_multiplicities,self.all_possible_subclonal_multiplicities))
        peak_names = self.get_multiplicity_names()
        
        vaf = self.mutation_table["Tumor_Alt_Count"] / (self.mutation_table["Tumor_Alt_Count"] + self.mutation_table["Tumor_Ref_Count"])
        highest_vaf_m_table = self.mutation_table[
            vaf > np.quantile(vaf, self.coverage_vaf_quantile) - 0.01
        ]
        highest_vaf_average_coverage = np.average(highest_vaf_m_table['Tumor_Alt_Count']+highest_vaf_m_table['Tumor_Ref_Count'])
        
        if not self.apply_reads_correction:
            highest_vaf_average_coverage = (self.mutation_table['Tumor_Alt_Count']+self.mutation_table['Tumor_Ref_Count']).mean()

        for multiplicity in combined_all_possible_multiplicities:
            # if multiplicity is bigger than total cn then the vaf can be > 1 which causes a crash
            #any mult > major cn is set to zero subsequently
            mult_vaf = np.minimum(1, mult_one_vaf * multiplicity)

            total_counts = self.mutation_table["Tumor_Alt_Count"] + self.mutation_table["Tumor_Ref_Count"]
            mult_probability = binom.logpmf(self.mutation_table["Tumor_Alt_Count"],total_counts,mult_vaf)
            multiplicity_probabilities.append(mult_probability)

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

        multiplicity_probabilities = np.array(multiplicity_probabilities)
        alt_count_correction_factors = np.array(
            alt_count_correction_factors
        )
        
        multiplicity_probabilities = np.exp(multiplicity_probabilities-np.max(multiplicity_probabilities,axis=0))
        
        normalising_sums = np.sum(multiplicity_probabilities, axis=0)

        
        normalising_sums = np.where(np.isclose(normalising_sums,0),1,normalising_sums)
        multiplicity_probabilities = multiplicity_probabilities / normalising_sums
        new_cols = {}
        for index, peak_name in enumerate(peak_names):
            new_cols[f"Prob_{peak_name}"] = multiplicity_probabilities[index, :]
            new_cols[f'Alt_Count_Correction_{peak_name}'] = (
                alt_count_correction_factors[index, :]
            )
        new_cols_table = pd.DataFrame(new_cols,index=self.mutation_table.index)
        self.mutation_table = pd.concat((self.mutation_table,new_cols_table),axis=1)
        
        
    def get_reads_correction_array(self, allele=None):
        if allele == 'minor':
            clonal_multiplicities = np.arange(1, self.minor_cn + 1)
        else:
            clonal_multiplicities = np.arange(1, self.major_cn + 1)
        mult_names = [
            f"Alt_Count_Correction_Mult_{mult}"
            for mult in clonal_multiplicities
        ]
        mult_names.extend([
            f"Alt_Count_Correction_Subclone_{subclone}"
            for subclone in range(self.n_subclones)
        ])
 
        reads_correction = self.mutation_table[mult_names].to_numpy()

        if allele is None:
            reads_correction = reads_correction[
                self.mutation_table['Phasing'].isna(),
                :,
            ]
        elif allele != 'all':
            reads_correction = reads_correction[
                self.mutation_table['Phasing'] == allele,
                :,
            ]
        
        if reads_correction.size == 0:
            return None

        
        
        
        reads_correction = np.average(reads_correction, axis=0)
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

        multiplicity_probabilities = self.mutation_table[
            mult_names
        ].to_numpy()

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

        if allele is None:
            multiplicity_probabilities = multiplicity_probabilities[
                self.mutation_table['Phasing'].isna(),
                :,
            ]
        elif allele != 'all':
            multiplicity_probabilities = multiplicity_probabilities[
                self.mutation_table['Phasing'] == allele,
                :,
            ]

        if multiplicity_probabilities.size == 0:
            return None
        
        return multiplicity_probabilities
    

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
        
        return MultProbabilityStore(
            mult_array_store,
            reads_correction_store,
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
    
  
   
