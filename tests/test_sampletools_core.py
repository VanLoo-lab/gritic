import types
import unittest
import warnings
from unittest import mock

import numpy as np
import pandas as pd
from scipy.stats import binom, poisson

from gritic import sampletools


def make_segment_mutation_table(
    *,
    segment_id='segment-a',
    chromosome='1',
    phasing=(np.nan, 'major', 'minor'),
):
    row_count = len(phasing)
    return pd.DataFrame({
        'Segment_ID': [segment_id] * row_count,
        'Chromosome': [chromosome] * row_count,
        'Segment_Start': [10] * row_count,
        'Segment_End': [110] * row_count,
        'Major_CN': [2] * row_count,
        'Minor_CN': [1] * row_count,
        'Tumor_Ref_Count': [15, 10, 5][:row_count],
        'Tumor_Alt_Count': [5, 10, 15][:row_count],
        'Phasing': list(phasing),
    })


class ScalarParameterValidationTest(unittest.TestCase):
    def test_non_negative_integer_validation(self):
        for value in (0, 7, np.int64(3)):
            with self.subTest(valid=value):
                validated = sampletools._validate_non_negative_integer(
                    value,
                    'count',
                )
                self.assertIs(type(validated), int)
                self.assertEqual(validated, int(value))

        for value in (-1, 1.0, True, np.bool_(False), '3', None):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(
                    ValueError,
                    'count must be a non-negative integer',
                ):
                    sampletools._validate_non_negative_integer(value, 'count')

    def test_closed_unit_interval_validation(self):
        for value in (0, 1, 0.125, np.float64(0.75)):
            with self.subTest(valid=value):
                self.assertEqual(
                    sampletools._validate_unit_interval_number(
                        value,
                        'fraction',
                    ),
                    float(value),
                )

        for value in (
            -0.01,
            1.01,
            True,
            np.nan,
            np.inf,
            -np.inf,
            '0.5',
            None,
        ):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                    sampletools._validate_unit_interval_number(
                        value,
                        'fraction',
                    )

    def test_minimum_subclone_ccf_uses_open_lower_bound(self):
        for value in (1e-12, 0.5, 1):
            with self.subTest(valid=value):
                self.assertEqual(
                    sampletools.validate_min_subclone_ccf(value),
                    float(value),
                )

        for value in (0, -0.1, 1.1, True, np.nan, np.inf, '0.1'):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'greater than 0'):
                    sampletools.validate_min_subclone_ccf(value)

    def test_coverage_quantile_accepts_both_endpoints(self):
        self.assertEqual(sampletools.validate_coverage_vaf_quantile(0), 0.0)
        self.assertEqual(sampletools.validate_coverage_vaf_quantile(1), 1.0)

        for value in (-0.1, 1.1, True, np.nan, np.inf, '0.9'):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                    sampletools.validate_coverage_vaf_quantile(value)

    def test_shared_biology_validators_delegate_to_dataloader(self):
        with mock.patch.object(
            sampletools.dataloader,
            'validate_purity',
            return_value=0.7,
        ) as purity_mock:
            self.assertEqual(sampletools.validate_purity('input'), 0.7)
        purity_mock.assert_called_once_with('input')

        with mock.patch.object(
            sampletools.dataloader,
            'get_normal_total_copy_number',
            return_value=1,
        ) as copy_number_mock:
            self.assertEqual(
                sampletools.get_normal_total_copy_number('X', 'XY'),
                1,
            )
        copy_number_mock.assert_called_once_with('X', 'XY')


class PortableSampleIdTest(unittest.TestCase):
    def test_ordinary_ascii_and_unicode_ids_are_preserved(self):
        for sample_id in ('sample-01_A', '.hidden', '患者-α'):
            with self.subTest(sample_id=sample_id):
                self.assertEqual(
                    sampletools.validate_sample_id(sample_id),
                    sample_id,
                )

    def test_non_string_empty_and_dot_components_are_rejected(self):
        for sample_id in (None, 17, b'sample', '', '.', '..'):
            with self.subTest(sample_id=sample_id):
                with self.assertRaises(ValueError):
                    sampletools.validate_sample_id(sample_id)

    def test_every_windows_forbidden_character_is_rejected(self):
        for character in '<>:"/\\|?*':
            with self.subTest(character=character):
                with self.assertRaisesRegex(
                    ValueError,
                    'Windows-forbidden',
                ):
                    sampletools.validate_sample_id(f'left{character}right')

    def test_unicode_control_format_and_surrogate_characters_are_rejected(self):
        characters = ('\x00', '\n', '\u200d', chr(0xD800))
        for character in characters:
            with self.subTest(code_point=ord(character)):
                with self.assertRaisesRegex(
                    ValueError,
                    'Unicode control, format, or surrogate',
                ):
                    sampletools.validate_sample_id(f'a{character}b')

    def test_trailing_dot_or_space_is_rejected(self):
        for sample_id in ('sample.', 'sample ', 'sample. '):
            with self.subTest(sample_id=sample_id):
                with self.assertRaisesRegex(ValueError, 'end in a dot or space'):
                    sampletools.validate_sample_id(sample_id)

    def test_windows_device_names_are_rejected_case_insensitively(self):
        reserved_ids = (
            'CON',
            'con.txt',
            'PRN.results',
            'AUX',
            'NUL.data',
            'COM1',
            'lpt9.anything',
            'COM²',
            'CON .txt',
        )
        for sample_id in reserved_ids:
            with self.subTest(sample_id=sample_id):
                with self.assertRaisesRegex(ValueError, 'reserved Windows'):
                    sampletools.validate_sample_id(sample_id)

    def test_ascii_component_limit_accounts_for_longest_output_suffix(self):
        available = (
            sampletools._MAX_PATH_COMPONENT_UNITS
            - len(sampletools._LONGEST_SAMPLE_ID_OUTPUT_SUFFIX)
        )
        boundary_id = 'a' * available

        self.assertEqual(
            sampletools.validate_sample_id(boundary_id),
            boundary_id,
        )
        with self.assertRaisesRegex(ValueError, 'too long'):
            sampletools.validate_sample_id(boundary_id + 'a')

    def test_component_limit_is_measured_in_encoded_units_not_characters(self):
        suffix_bytes = len(
            sampletools._LONGEST_SAMPLE_ID_OUTPUT_SUFFIX.encode('utf-8')
        )
        max_repeated_e_acute = (
            sampletools._MAX_PATH_COMPONENT_UNITS - suffix_bytes
        ) // len('é'.encode('utf-8'))
        boundary_id = 'é' * max_repeated_e_acute

        self.assertEqual(
            sampletools.validate_sample_id(boundary_id),
            boundary_id,
        )
        with self.assertRaisesRegex(ValueError, 'too long'):
            sampletools.validate_sample_id(boundary_id + 'é')


class MajorCopyNumberModeTest(unittest.TestCase):
    def test_mode_is_weighted_by_genomic_width_and_input_is_not_mutated(self):
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', '2', '2'],
            'Segment_Start': [0, 0, 40],
            'Segment_End': [60, 40, 75],
            'Major_CN': [2, 3, 3],
            'Minor_CN': [1, 1, 1],
        })
        original = copy_number_table.copy(deep=True)

        mode = sampletools.get_major_cn_mode_from_cn_table(copy_number_table)

        self.assertEqual(mode, 3)
        pd.testing.assert_frame_equal(copy_number_table, original)
        self.assertNotIn('Segment_Width', copy_number_table.columns)

    def test_sample_mode_ignores_sex_chromosomes_and_chr_prefix(self):
        sample = types.SimpleNamespace(
            cn_table=pd.DataFrame({
                'Chromosome': ['chr1', 'X'],
                'Segment_Start': [0, 0],
                'Segment_End': [10, 10_000],
                'Major_CN': [2, 8],
                'Minor_CN': [1, 1],
            }),
            autosomes=['1', '2'],
        )

        self.assertEqual(sampletools.get_major_cn_mode(sample), 2)

    def test_sample_mode_requires_at_least_one_configured_autosome(self):
        sample = types.SimpleNamespace(
            cn_table=pd.DataFrame({
                'Chromosome': ['X'],
                'Segment_Start': [0],
                'Segment_End': [100],
                'Major_CN': [2],
                'Minor_CN': [1],
            }),
            autosomes=['1', '2'],
        )

        with self.assertRaisesRegex(ValueError, 'no segments.*autosomes'):
            sampletools.get_major_cn_mode(sample)


class MultProbabilityStoreTest(unittest.TestCase):
    def test_constructor_tracks_phase_availability_and_mutation_count(self):
        major = np.array([[0.2, 0.8], [0.6, 0.4]])
        store = sampletools.MultProbabilityStore(
            {
                'Non_Phased': None,
                'Major': major,
                'Minor': None,
                'All': major,
            },
            {
                'Non_Phased': None,
                'Major': np.ones(2),
                'Minor': None,
                'All': np.ones(2),
            },
            major_cn=2,
            minor_cn=0,
            n_subclones=0,
        )

        self.assertFalse(store.use_non_phased)
        self.assertTrue(store.use_major)
        self.assertFalse(store.use_minor)
        self.assertEqual(store.n_mutations, 2)

    def test_constructor_checks_probability_and_correction_dimensions(self):
        with self.assertRaises(AssertionError):
            sampletools.MultProbabilityStore(
                {
                    'Non_Phased': np.ones((2, 3)),
                    'Major': None,
                    'Minor': None,
                    'All': np.ones((2, 3)),
                },
                {
                    'Non_Phased': np.ones(2),
                    'Major': None,
                    'Minor': None,
                    'All': np.ones(3),
                },
                major_cn=3,
                minor_cn=0,
                n_subclones=0,
            )

    def test_subclonal_correction_is_weighted_by_phase_mutation_counts(self):
        store = sampletools.MultProbabilityStore(
            {
                'Non_Phased': np.full((2, 4), 0.25),
                'Major': np.full((1, 4), 0.25),
                'Minor': np.full((3, 3), 1 / 3),
                'All': np.full((6, 4), 0.25),
            },
            {
                'Non_Phased': np.array([0.1, 0.2, 0.3, 0.4]),
                'Major': np.array([0.2, 0.4, 0.6, 0.8]),
                'Minor': np.array([0.25, 0.45, 0.65]),
                'All': np.ones(4),
            },
            major_cn=2,
            minor_cn=1,
            n_subclones=2,
        )
        subclones = pd.DataFrame(index=[0, 1])

        correction = store.get_subclonal_correction_array(subclones)

        non_phased = np.array([0.1, 0.3, 0.4])
        major = np.array([0.2, 0.6, 0.8])
        minor = np.array([0.25, 0.45, 0.65])
        expected = (2 * non_phased + major + 3 * minor) / 6
        np.testing.assert_allclose(correction, expected)

    def test_scalar_likelihood_clips_and_normalizes_state(self):
        combined = np.array([[0.2, 0.8], [0.9, 0.1]])
        store = sampletools.MultProbabilityStore(
            {
                'Non_Phased': combined,
                'Major': None,
                'Minor': None,
                'All': combined,
            },
            {
                'Non_Phased': np.ones(2),
                'Major': None,
                'Minor': None,
                'All': np.ones(2),
            },
            major_cn=2,
            minor_cn=0,
            n_subclones=0,
        )
        state = np.array([-2.0, 2.0])

        likelihood = store.evaluate_likelihood(state)

        self.assertAlmostEqual(likelihood, np.log(0.8) + np.log(0.1))
        np.testing.assert_array_equal(state, [-2.0, 2.0])


class SegmentDeterministicBehaviorTest(unittest.TestCase):
    def test_segment_metadata_clone_state_and_phase_arrays(self):
        segment = sampletools.Segment(
            make_segment_mutation_table(),
            None,
            purity=0.8,
            sex=None,
            apply_reads_correction=False,
            _validated=True,
        )

        self.assertEqual(str(segment), 'segment-a- 2+1 - 3 Mutations')
        self.assertEqual(segment.n_mutations, 3)
        self.assertEqual(segment.width, 100)
        self.assertEqual(segment.get_n_subclones(), 0)
        np.testing.assert_array_equal(segment.sample_clone_fractions, [1])
        np.testing.assert_array_equal(
            segment.all_possible_clonal_multiplicities,
            [1, 2],
        )
        self.assertEqual(segment.get_multiplicity_names(), ['Mult_1', 'Mult_2'])

        probability_columns = ['Prob_Mult_1', 'Prob_Mult_2']
        np.testing.assert_allclose(
            segment.mutation_table[probability_columns].sum(axis=1),
            1.0,
        )
        store = segment.multiplicity_probabilities
        self.assertEqual(store.non_phased_array.shape, (1, 2))
        self.assertEqual(store.major_array.shape, (1, 2))
        self.assertEqual(store.minor_array.shape, (1, 1))
        self.assertEqual(store.combined_array.shape, (3, 2))
        np.testing.assert_array_equal(
            store.reads_correction_combined_array,
            [1.0, 1.0],
        )

    def test_subclones_extend_names_multiplicities_and_clone_fractions(self):
        subclones = pd.DataFrame({
            'Cluster': ['c1', 'c2'],
            'Subclone_CCF': [0.6, 0.25],
            'Subclone_Fraction': [0.2, 0.3],
        })
        segment = sampletools.Segment(
            make_segment_mutation_table(),
            subclones,
            purity=0.8,
            sex=None,
            apply_reads_correction=False,
            _validated=True,
        )

        self.assertEqual(segment.n_subclones, 2)
        np.testing.assert_allclose(
            segment.sample_clone_fractions,
            [0.5, 0.2, 0.3],
        )
        np.testing.assert_allclose(
            segment.all_possible_subclonal_multiplicities,
            [0.6, 0.25],
        )
        self.assertEqual(
            segment.get_multiplicity_names(),
            ['Mult_1', 'Mult_2', 'Subclone_0', 'Subclone_1'],
        )
        probability_columns = [
            'Prob_Mult_1',
            'Prob_Mult_2',
            'Prob_Subclone_0',
            'Prob_Subclone_1',
        ]
        np.testing.assert_allclose(
            segment.mutation_table[probability_columns].sum(axis=1),
            1.0,
        )

    def test_probability_and_detection_correction_match_closed_form_values(self):
        segment = sampletools.Segment(
            make_segment_mutation_table(),
            None,
            purity=0.8,
            sex=None,
            apply_reads_correction=True,
            min_mutation_alt_count=3,
            coverage_vaf_quantile=0.9,
            _validated=True,
        )
        # total CN=3, normal CN=2, purity=.8 gives one-copy VAF 2/7.
        multiplicity_vafs = np.array([2.0 / 7.0, 4.0 / 7.0])
        expected_corrections = poisson.sf(2, 20 * multiplicity_vafs)
        for index, expected in enumerate(expected_corrections, start=1):
            np.testing.assert_allclose(
                segment.mutation_table[
                    f'Alt_Count_Correction_Mult_{index}'
                ],
                expected,
            )

        first_log_probabilities = binom.logpmf(5, 20, multiplicity_vafs)
        expected_first_probabilities = np.exp(
            first_log_probabilities - first_log_probabilities.max()
        )
        expected_first_probabilities /= expected_first_probabilities.sum()
        np.testing.assert_allclose(
            segment.mutation_table.loc[
                0,
                ['Prob_Mult_1', 'Prob_Mult_2'],
            ].to_numpy(dtype=float),
            expected_first_probabilities,
        )

    def test_unique_segment_attribute_rejects_mixed_values(self):
        segment = sampletools.Segment.__new__(sampletools.Segment)
        segment.mutation_table = pd.DataFrame({'Major_CN': [2, 3]})

        with self.assertRaisesRegex(ValueError, 'Major_CN is not unique'):
            segment.get_unique_attribute_from_table('Major_CN')

    def test_sex_is_required_for_sex_chromosome_segment(self):
        with self.assertRaisesRegex(ValueError, 'sex must be specified'):
            sampletools.Segment(
                make_segment_mutation_table(chromosome='X'),
                None,
                purity=0.8,
                sex=None,
                _validated=True,
            )

    def test_mutation_rate_applies_multiplicity_and_read_corrections(self):
        segment = sampletools.Segment.__new__(sampletools.Segment)
        segment.major_cn = 2
        segment.n_subclones = 1
        segment.n_mutations = 50
        segment.width = 100
        segment.multiplicity_probabilities = types.SimpleNamespace(
            reads_correction_combined_array=np.array([1.0, 0.5, 0.25]),
        )

        with mock.patch.object(
            sampletools.multiplicityoptimiser,
            'unconstrained_mult_optimisation',
            return_value=np.array([0.5, 0.25, 0.25]),
        ) as optimizer:
            mutation_rate = segment.get_mutation_rate()

        optimizer.assert_called_once_with(segment.multiplicity_probabilities, 1)
        self.assertEqual(mutation_rate, 1.25)

    def test_mutation_rate_is_nan_when_optimisation_fails(self):
        segment = sampletools.Segment.__new__(sampletools.Segment)
        segment.major_cn = 1
        segment.n_subclones = 0
        segment.n_mutations = 10
        segment.width = 100
        segment.multiplicity_probabilities = object()

        with mock.patch.object(
            sampletools.multiplicityoptimiser,
            'unconstrained_mult_optimisation',
            return_value=None,
        ):
            self.assertTrue(np.isnan(segment.get_mutation_rate()))


class SampleTableHelperTest(unittest.TestCase):
    def test_filter_mutation_table_applies_both_read_thresholds(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.min_mutation_alt_count = 3
        sample.min_mutation_coverage = 10
        mutations = pd.DataFrame({
            'Tumor_Alt_Count': [2, 3, 3, 8],
            'Tumor_Ref_Count': [20, 6, 7, 12],
            'Mutation_ID': ['low-alt', 'low-depth', 'keep-a', 'keep-b'],
        })

        filtered = sample.filter_mutation_table(mutations)

        self.assertEqual(filtered['Mutation_ID'].tolist(), ['keep-a', 'keep-b'])
        self.assertEqual(filtered['Total_Count'].tolist(), [10, 20])

    def test_filter_mutation_table_reports_when_no_rows_remain(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.min_mutation_alt_count = 3
        sample.min_mutation_coverage = 10
        mutations = pd.DataFrame({
            'Tumor_Alt_Count': [1, 2],
            'Tumor_Ref_Count': [20, 20],
        })

        with self.assertRaisesRegex(ValueError, 'No mutations remain'):
            sample.filter_mutation_table(mutations)

    def test_high_alt_count_input_emits_diagnostic_warning(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.min_mutation_alt_count = 3
        sample.min_mutation_coverage = 10
        mutations = pd.DataFrame({
            'Tumor_Alt_Count': [11, 12],
            'Tumor_Ref_Count': [9, 8],
        })

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sample.filter_mutation_table(mutations)

        self.assertEqual(len(caught), 1)
        self.assertIn('no mutations with less than 10 alt reads', str(caught[0].message))

    def test_chromosome_sort_uses_configured_biological_order(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.chromosome_order = ['1', '2', 'X', 'Y']
        table = pd.DataFrame({
            'Chromosome': ['X', '2', '1', '1'],
            'Segment_Start': [5, 20, 30, 10],
            'value': ['x', 'two', 'one-late', 'one-early'],
        })

        sorted_table = sample.sort_table_by_chromosome(table)

        self.assertEqual(
            sorted_table['value'].tolist(),
            ['one-early', 'one-late', 'two', 'x'],
        )

    def test_subclone_formatting_filters_and_adds_rounded_snv_counts(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.clip_subclone_ccf = False
        sample.min_subclone_ccf = 0.01
        sample.max_subclone_ccf = 0.9
        sample.min_subclone_fraction = 0.1
        sample.mutation_table = pd.DataFrame(index=range(11))
        subclones = pd.DataFrame({
            'Cluster': ['high', 'low'],
            'Subclone_CCF': [0.7, 0.3],
            'Subclone_Fraction': [0.3, 0.2],
        })

        formatted = sample.format_subclone_table(subclones)

        self.assertEqual(
            formatted.columns.tolist(),
            ['Cluster', 'Subclone_CCF', 'Subclone_Fraction', 'N_SNVs'],
        )
        self.assertEqual(formatted['Cluster'].tolist(), ['high', 'low'])
        self.assertEqual(formatted['N_SNVs'].tolist(), [3, 2])

    def test_subclone_formatting_returns_none_when_every_clone_is_filtered(self):
        sample = sampletools.Sample.__new__(sampletools.Sample)
        sample.clip_subclone_ccf = False
        sample.min_subclone_ccf = 0.01
        sample.max_subclone_ccf = 0.9
        sample.min_subclone_fraction = 0.1
        sample.mutation_table = pd.DataFrame(index=range(10))
        subclones = pd.DataFrame({
            'Cluster': ['too-high'],
            'Subclone_CCF': [0.95],
            'Subclone_Fraction': [0.2],
        })

        self.assertIsNone(sample.format_subclone_table(subclones))
        self.assertIsNone(sample.format_subclone_table(None))


if __name__ == '__main__':
    unittest.main()
