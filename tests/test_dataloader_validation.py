import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from gritic import dataloader


class DataloaderFixtures:
    @staticmethod
    def copy_number_table(**overrides):
        values = {
            'Chromosome': ['1', '2'],
            'Segment_Start': [0, 100],
            'Segment_End': [100, 300],
            'Major_CN': [2, 3],
            'Minor_CN': [1, 1],
            'Segment_ID': ['segment-a', 'segment-b'],
        }
        values.update(overrides)
        return pd.DataFrame(values)

    @staticmethod
    def mutation_table(**overrides):
        values = {
            'Chromosome': ['1', '2'],
            'Position': [50, 150],
            'Tumor_Ref_Count': [10, 12],
            'Tumor_Alt_Count': [5, 3],
            'Mutation_ID': ['mutation-a', 'mutation-b'],
            'Segment_ID': ['segment-a', 'segment-b'],
        }
        values.update(overrides)
        return pd.DataFrame(values)


class InputSchemaValidationTest(DataloaderFixtures, unittest.TestCase):
    def test_header_validation_selects_supplied_segment_id_mode(self):
        copy_columns = list(dataloader.COPY_NUMBER_REQUIRED_COLUMNS)
        mutation_columns = [
            *dataloader.MUTATION_REQUIRED_COLUMNS,
            'Position',
        ]
        self.assertFalse(dataloader.validate_input_table_headers(
            copy_columns,
            mutation_columns,
        ))

        copy_columns.append('Segment_ID')
        mutation_columns = [
            *dataloader.MUTATION_REQUIRED_COLUMNS,
            'Mutation_ID',
            'Segment_ID',
        ]
        self.assertTrue(dataloader.validate_input_table_headers(
            copy_columns,
            mutation_columns,
        ))

    def test_header_validation_reports_missing_required_columns_in_schema_order(self):
        with self.assertRaisesRegex(
            ValueError,
            r'Copy number table is missing required column\(s\): '
            'Segment_Start, Major_CN',
        ):
            dataloader.validate_input_table_headers(
                ['Chromosome', 'Segment_End', 'Minor_CN'],
                [*dataloader.MUTATION_REQUIRED_COLUMNS, 'Position'],
            )

        with self.assertRaisesRegex(
            ValueError,
            r'Mutation table is missing required column\(s\): '
            'Tumor_Ref_Count, Tumor_Alt_Count',
        ):
            dataloader.validate_input_table_headers(
                dataloader.COPY_NUMBER_REQUIRED_COLUMNS,
                ['Chromosome', 'Position'],
            )

    def test_header_validation_requires_identity_and_assignment_information(self):
        copy_columns = list(dataloader.COPY_NUMBER_REQUIRED_COLUMNS)
        mutation_columns = list(dataloader.MUTATION_REQUIRED_COLUMNS)
        with self.assertRaisesRegex(ValueError, 'either Mutation_ID or Position'):
            dataloader.validate_input_table_headers(
                copy_columns,
                mutation_columns,
            )

        mutation_columns.append('Mutation_ID')
        with self.assertRaisesRegex(
            ValueError,
            'Position for copy-number segment assignment',
        ):
            dataloader.validate_input_table_headers(
                copy_columns,
                mutation_columns,
            )

    def test_mutation_identity_prefers_mutation_id_and_rejects_missing_values(self):
        table = pd.DataFrame({
            'Mutation_ID': ['named', 'second'],
            'Position': [10, 20],
        })
        self.assertEqual(
            dataloader.get_gritic_mutation_id_components(table).tolist(),
            ['named', 'second'],
        )
        self.assertEqual(
            dataloader.get_gritic_mutation_id_components(
                table.drop(columns='Mutation_ID')
            ).tolist(),
            ['10', '20'],
        )

        for invalid_value in (None, '', '   '):
            with self.subTest(invalid_value=invalid_value):
                with self.assertRaisesRegex(
                    ValueError,
                    'Mutation_ID value that is not null, empty, or whitespace',
                ):
                    dataloader.get_gritic_mutation_id_components(
                        pd.DataFrame({'Mutation_ID': [invalid_value]})
                    )

    def test_mutation_identity_must_be_unique_only_within_source_segment(self):
        source_ids = pd.Series(['one', 'two'])
        components = pd.Series(['same', 'same'])
        self.assertIsNone(
            dataloader.validate_source_scoped_mutation_components(
                source_ids,
                components,
            )
        )

        with self.assertRaisesRegex(
            ValueError,
            r'duplicate pair\(s\): \(one, same\)',
        ):
            dataloader.validate_source_scoped_mutation_components(
                pd.Series(['one', 'one']),
                components,
            )

    def test_generated_mutation_ids_are_unambiguous_and_preserve_index(self):
        source_ids = pd.Series(['segment:a', 'space / percent%'], index=[4, 9])
        components = pd.Series(['mutation/b', 'x:y'], index=[4, 9])

        result = dataloader.generate_gritic_mutation_ids(
            source_ids,
            components,
        )

        self.assertEqual(result.index.tolist(), [4, 9])
        self.assertEqual(str(result.dtype), 'string')
        self.assertEqual(result.tolist(), [
            'segment%3Aa:mutation%2Fb',
            'space%20%2F%20percent%25:x%3Ay',
        ])


class ScalarAndChromosomeValidationTest(DataloaderFixtures, unittest.TestCase):
    def test_sex_autosome_and_purity_validators_accept_canonical_boundaries(self):
        for sex in (None, 'XX', 'XY', 'ZZ', 'ZW'):
            with self.subTest(sex=sex):
                self.assertEqual(dataloader.validate_sex_karyotype(sex), sex)
        self.assertEqual(dataloader.validate_autosome_count(np.int64(3)), 3)
        self.assertEqual(dataloader.get_autosome_labels(3), ('1', '2', '3'))
        self.assertEqual(dataloader.validate_purity('0.25'), 0.25)
        self.assertEqual(dataloader.validate_purity(1), 1.0)

    def test_sex_autosome_and_purity_validators_reject_invalid_values(self):
        for sex in ('xx', 'XO', 1):
            with self.subTest(sex=sex):
                with self.assertRaisesRegex(ValueError, 'XX, XY, ZZ, or ZW'):
                    dataloader.validate_sex_karyotype(sex)

        for autosome_count in (True, np.bool_(False), 0, -1, 1.0, '2', None):
            with self.subTest(autosome_count=autosome_count):
                with self.assertRaisesRegex(ValueError, 'positive integer'):
                    dataloader.validate_autosome_count(autosome_count)

        for purity in (
            False,
            np.bool_(True),
            0,
            -0.1,
            1.01,
            'nan',
            float('inf'),
            object(),
        ):
            with self.subTest(purity=purity):
                with self.assertRaisesRegex(ValueError, 'Purity must be finite'):
                    dataloader.validate_purity(purity)

    def test_normal_copy_number_uses_karyotype_specific_sex_chromosomes(self):
        expected = {
            'XX': {'1': 2, 'X': 2, 'Y': 0},
            'XY': {'1': 2, 'X': 1, 'Y': 1},
            'ZZ': {'1': 2, 'Z': 2, 'W': 0},
            'ZW': {'1': 2, 'Z': 1, 'W': 1},
        }
        for sex, chromosome_values in expected.items():
            for chromosome, value in chromosome_values.items():
                with self.subTest(sex=sex, chromosome=chromosome):
                    self.assertEqual(
                        dataloader.get_normal_total_copy_number(chromosome, sex),
                        value,
                    )
        self.assertEqual(dataloader.get_normal_total_copy_number('Y', None), 2)

    def test_chromosome_normalization_removes_one_lowercase_prefix_without_mutation(self):
        table = pd.DataFrame({
            'Chromosome': ['chr1', 'chrchr2', 'Chr3', 4],
            'value': [1, 2, 3, 4],
        })

        result = dataloader.normalize_chromosome_labels(table)

        self.assertEqual(result['Chromosome'].tolist(), ['1', 'chr2', 'Chr3', '4'])
        self.assertEqual(table['Chromosome'].tolist(), ['chr1', 'chrchr2', 'Chr3', 4])

    def test_allowed_chromosomes_include_only_present_sex_chromosomes(self):
        self.assertEqual(
            dataloader.get_allowed_chromosome_labels(2, 'XX'),
            {'1', '2', 'X'},
        )
        self.assertEqual(
            dataloader.get_allowed_chromosome_labels(2, 'ZW'),
            {'1', '2', 'Z', 'W'},
        )
        with self.assertRaisesRegex(ValueError, 'sex must be resolved'):
            dataloader.get_allowed_chromosome_labels(2, None)

    def test_invalid_chromosomes_raise_or_are_dropped_with_one_warning(self):
        table = pd.DataFrame({'Chromosome': ['chr1', '3', 'Y', '3']})
        with self.assertRaisesRegex(
            ValueError,
            r'numbered autosomes 1 through 2 or X for sex XX.*3, Y',
        ):
            dataloader.validate_or_drop_chromosome_rows(
                table,
                2,
                'XX',
                table_name='Mutation table',
            )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.validate_or_drop_chromosome_rows(
                table,
                2,
                'XX',
                table_name='Mutation table',
                drop_unmatched_chromosomes=True,
            )

        self.assertEqual(result['Chromosome'].tolist(), ['1'])
        self.assertEqual(len(caught), 1)
        self.assertRegex(str(caught[0].message), r'Dropping 3 mutation table rows')
        self.assertIn('3, Y', str(caught[0].message))

    def test_chromosome_drop_can_be_silent(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.validate_or_drop_chromosome_rows(
                pd.DataFrame({'Chromosome': ['1', '99']}),
                2,
                'XX',
                table_name='Table',
                drop_unmatched_chromosomes=True,
                warn=False,
            )
        self.assertEqual(result['Chromosome'].tolist(), ['1'])
        self.assertEqual(caught, [])

    def test_joint_chromosome_filter_infers_karyotype_and_rejects_empty_outputs(self):
        copy_number = pd.DataFrame({'Chromosome': ['chr1', 'W']})
        mutations = pd.DataFrame({'Chromosome': ['chr1', 'W']})
        filtered_cn, filtered_mutations, sex = (
            dataloader.validate_or_drop_input_chromosomes(
                copy_number,
                mutations,
                autosome_count=1,
                sex=None,
            )
        )
        self.assertEqual(sex, 'ZW')
        self.assertEqual(filtered_cn['Chromosome'].tolist(), ['1', 'W'])
        self.assertEqual(filtered_mutations['Chromosome'].tolist(), ['1', 'W'])

        cases = (
            (['99'], ['99'], 'No copy-number segments or mutations'),
            (['99'], ['1'], 'No copy-number segments remain'),
            (['1'], ['99'], 'No mutations remain'),
        )
        for copy_chromosomes, mutation_chromosomes, message in cases:
            with self.subTest(message=message):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    with self.assertRaisesRegex(ValueError, message):
                        dataloader.validate_or_drop_input_chromosomes(
                            pd.DataFrame({'Chromosome': copy_chromosomes}),
                            pd.DataFrame({'Chromosome': mutation_chromosomes}),
                            autosome_count=1,
                            sex='XX',
                            drop_unmatched_chromosomes=True,
                        )

    def test_sex_inference_handles_all_supported_systems_and_mixed_error(self):
        cases = (
            (['1'], 'XX'),
            (['chrX'], 'XX'),
            (['X', 'Y'], 'XY'),
            (['Z'], 'ZZ'),
            (['chrZ', 'chrW'], 'ZW'),
        )
        for chromosomes, expected in cases:
            with self.subTest(chromosomes=chromosomes):
                self.assertEqual(
                    dataloader.infer_sex_from_copy_number_table(
                        pd.DataFrame({'Chromosome': chromosomes})
                    ),
                    expected,
                )
        with self.assertRaisesRegex(ValueError, 'mixes X/Y and Z/W'):
            dataloader.infer_sex_from_copy_number_table(
                pd.DataFrame({'Chromosome': ['X', 'W']})
            )


class NumericTableValidationTest(DataloaderFixtures, unittest.TestCase):
    def test_positions_accept_exact_integer_spellings_and_int64_limit(self):
        maximum = np.iinfo(np.int64).max
        positions = pd.Series(['0', '2.000', '1e3', str(maximum)], index=[8, 6, 4, 2])

        result = dataloader.parse_positions(positions)

        self.assertEqual(result.tolist(), [0, 2, 1000, maximum])
        self.assertEqual(result.index.tolist(), [8, 6, 4, 2])
        self.assertEqual(result.dtype, np.dtype(np.int64))

    def test_positions_reject_fractional_nonfinite_negative_and_overflow_values(self):
        maximum = np.iinfo(np.int64).max
        positions = pd.Series([
            '-1',
            '1.5',
            'nan',
            'Infinity',
            str(maximum + 1),
            None,
            '-1',
        ])
        with self.assertRaisesRegex(
            ValueError,
            'invalid value\\(s\\): -1, 1.5, nan, Infinity, '
            f'{maximum + 1}$',
        ):
            dataloader.parse_positions(positions)

    def test_segment_coordinates_are_canonicalized_without_mutating_input(self):
        table = pd.DataFrame({
            'Segment_Start': ['0e0', '10.000'],
            'Segment_End': ['10', '2e1'],
            'other': ['a', 'b'],
        }, index=[3, 7])

        result = dataloader.validate_segment_coordinates(table)

        self.assertEqual(result['Segment_Start'].tolist(), [0, 10])
        self.assertEqual(result['Segment_End'].tolist(), [10, 20])
        self.assertEqual(result['Segment_Start'].dtype, np.dtype(np.int64))
        self.assertEqual(result.index.tolist(), [3, 7])
        self.assertEqual(table['Segment_Start'].tolist(), ['0e0', '10.000'])
        self.assertEqual(dataloader.get_segment_widths(result).tolist(), [10, 10])
        self.assertEqual(dataloader.get_segment_widths(result).dtype, object)

    def test_segment_coordinates_report_parse_errors_from_both_columns(self):
        table = pd.DataFrame({
            'Segment_Start': ['bad', '1.5'],
            'Segment_End': ['nan', str(np.iinfo(np.int64).max + 1)],
        })
        with self.assertRaisesRegex(ValueError, 'Segment_Start.*bad, 1.5.*Segment_End.*nan'):
            dataloader.validate_segment_coordinates(table)

    def test_segment_coordinates_reject_negative_empty_and_overwide_intervals(self):
        minimum = np.iinfo(np.int64).min
        maximum = np.iinfo(np.int64).max
        table = pd.DataFrame({
            'Segment_Start': [-1, 10, minimum],
            'Segment_End': [5, 10, maximum],
        })
        with self.assertRaisesRegex(
            ValueError,
            'Segment_Start must be non-negative.*'
            'Segment_End must be greater.*segment width must not exceed',
        ):
            dataloader.validate_segment_coordinates(table)

    def test_read_counts_accept_zero_in_one_allele_and_canonicalize(self):
        table = pd.DataFrame({
            'Tumor_Ref_Count': ['0', '2e1'],
            'Tumor_Alt_Count': ['4.0', '0'],
        }, index=[2, 5])

        result = dataloader.validate_mutation_read_counts(table)

        self.assertEqual(result['Tumor_Ref_Count'].tolist(), [0, 20])
        self.assertEqual(result['Tumor_Alt_Count'].tolist(), [4, 0])
        self.assertEqual(result['Tumor_Alt_Count'].dtype, np.dtype(np.int64))
        self.assertEqual(table['Tumor_Alt_Count'].tolist(), ['4.0', '0'])

    def test_read_counts_reject_invalid_values_zero_coverage_and_sum_overflow(self):
        invalid = pd.DataFrame({
            'Tumor_Ref_Count': ['-1', 'nan'],
            'Tumor_Alt_Count': ['1.5', 'bad'],
        })
        with self.assertRaisesRegex(ValueError, 'Tumor_Ref_Count.*-1, nan.*Tumor_Alt_Count.*1.5, bad'):
            dataloader.validate_mutation_read_counts(invalid)

        with self.assertRaisesRegex(ValueError, 'must be greater than 0'):
            dataloader.validate_mutation_read_counts(pd.DataFrame({
                'Tumor_Ref_Count': [0],
                'Tumor_Alt_Count': [0],
            }))

        maximum = np.iinfo(np.int64).max
        with self.assertRaisesRegex(ValueError, 'must not exceed the maximum'):
            dataloader.validate_mutation_read_counts(pd.DataFrame({
                'Tumor_Ref_Count': [maximum],
                'Tumor_Alt_Count': [1],
            }))

    def test_copy_numbers_are_canonicalized_and_preserve_other_columns(self):
        table = pd.DataFrame({
            'Major_CN': ['3.0', '2e0'],
            'Minor_CN': ['1', '2.000'],
            'Segment_ID': ['a', 'b'],
        })

        result = dataloader.validate_copy_number_values(table)

        self.assertEqual(result['Major_CN'].tolist(), [3, 2])
        self.assertEqual(result['Minor_CN'].tolist(), [1, 2])
        self.assertEqual(result['Major_CN'].dtype, np.dtype(np.int64))
        self.assertEqual(result['Segment_ID'].tolist(), ['a', 'b'])
        self.assertEqual(table['Major_CN'].tolist(), ['3.0', '2e0'])

    def test_copy_numbers_report_parse_and_biological_constraint_errors(self):
        invalid = pd.DataFrame({
            'Major_CN': ['nan', '1.5'],
            'Minor_CN': ['bad', str(np.iinfo(np.int64).max + 1)],
        })
        with self.assertRaisesRegex(ValueError, 'Major_CN.*nan, 1.5.*Minor_CN.*bad'):
            dataloader.validate_copy_number_values(invalid)

        semantic = pd.DataFrame({
            'Major_CN': [-1, 1],
            'Minor_CN': [-2, 2],
        })
        with self.assertRaisesRegex(
            ValueError,
            'Major_CN must be non-negative.*Minor_CN must be non-negative.*'
            'Minor_CN must not exceed Major_CN',
        ):
            dataloader.validate_copy_number_values(semantic)


class SegmentAndMutationAssignmentTest(DataloaderFixtures, unittest.TestCase):
    def test_supplied_segment_ids_validate_missing_duplicates_unknown_and_chromosome(self):
        valid_cn = self.copy_number_table()
        valid_mutations = self.mutation_table()
        self.assertIsNone(dataloader.validate_supplied_segment_ids(
            valid_cn,
            valid_mutations,
        ))

        cases = (
            (
                valid_cn.assign(Segment_ID=['', 'segment-b']),
                valid_mutations,
                'Every copy number table row',
                False,
            ),
            (
                valid_cn,
                valid_mutations.assign(Segment_ID=['segment-a', None]),
                'Every mutation table row',
                False,
            ),
            (
                valid_cn.assign(Segment_ID=['same', 'same']),
                valid_mutations,
                'globally unique',
                False,
            ),
            (
                valid_cn,
                valid_mutations.assign(Segment_ID=['segment-a', 'unknown']),
                'not present in the copy number table',
                False,
            ),
            (
                valid_cn,
                valid_mutations.assign(Chromosome=['2', '2']),
                'Chromosome does not match',
                False,
            ),
        )
        for copy_number, mutations, message, allow_unmatched in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    dataloader.validate_supplied_segment_ids(
                        copy_number,
                        mutations,
                        allow_unmatched=allow_unmatched,
                    )

        dataloader.validate_supplied_segment_ids(
            valid_cn.assign(Chromosome=['chr1', 'chr2']),
            valid_mutations,
        )
        dataloader.validate_supplied_segment_ids(
            valid_cn,
            valid_mutations.assign(Segment_ID=['segment-a', 'unknown']),
            allow_unmatched=True,
        )

    def test_unmatched_supplied_segment_ids_are_dropped_with_pluralized_warning(self):
        mutations = self.mutation_table(
            Chromosome=['1', '1', '2'],
            Segment_ID=['segment-a', 'missing-one', 'missing-two'],
            Mutation_ID=['kept', 'gone-a', 'gone-b'],
            Position=[10, 20, 150],
            Tumor_Ref_Count=[10, 10, 10],
            Tumor_Alt_Count=[2, 2, 2],
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.drop_unmatched_segment_id_mutations(
                self.copy_number_table(),
                mutations,
            )

        self.assertEqual(result['Mutation_ID'].tolist(), ['kept'])
        self.assertEqual(len(caught), 1)
        self.assertIn('Dropping 2 unmatched SNVs', str(caught[0].message))

    def test_coordinate_assignment_uses_half_open_intervals_and_overwrites_annotations(self):
        copy_number = self.copy_number_table(
            Chromosome=['1', '1'],
            Segment_Start=[0, 100],
            Segment_End=[100, 200],
            Major_CN=[2, 4],
            Minor_CN=[1, 2],
            Segment_ID=['ignored-a', 'ignored-b'],
        )
        mutations = pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Position': ['0', '99', '100'],
            'Mutation_ID': ['left-edge', 'left-end', 'right-edge'],
            'Segment_Start': [-1, -1, -1],
            'Major_CN': [99, 99, 99],
        })

        result = dataloader.assign_cn_to_snv(
            mutations,
            copy_number,
            use_supplied_segment_ids=False,
        )

        self.assertEqual(result['Segment_ID'].tolist(), [
            '1-0-100',
            '1-0-100',
            '1-100-200',
        ])
        self.assertEqual(result['Major_CN'].tolist(), [2, 2, 4])
        self.assertEqual(result['Minor_CN'].tolist(), [1, 1, 2])
        self.assertEqual(result['Total_CN'].tolist(), [3, 3, 6])
        self.assertEqual(mutations['Major_CN'].tolist(), [99, 99, 99])

    def test_coordinate_assignment_rejects_or_drops_unmatched_mutations(self):
        copy_number = self.copy_number_table().iloc[[0]].copy()
        mutations = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Position': [50, 100],
            'Mutation_ID': ['kept', 'unmatched'],
        })
        with self.assertRaisesRegex(ValueError, '1 unmatched SNV'):
            dataloader.assign_cn_to_snv(
                mutations,
                copy_number,
                use_supplied_segment_ids=False,
            )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.assign_cn_to_snv(
                mutations,
                copy_number,
                use_supplied_segment_ids=False,
                drop_unmatched_snvs=True,
            )
        self.assertEqual(result['Mutation_ID'].tolist(), ['kept'])
        self.assertEqual(len(caught), 1)
        self.assertIn('Dropping 1 unmatched SNV', str(caught[0].message))

    def test_supplied_id_assignment_joins_many_mutations_to_one_segment(self):
        copy_number = self.copy_number_table().iloc[[0]].copy()
        mutations = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_ID': [1, 1],
            'Mutation_ID': ['a', 'b'],
        })
        copy_number['Segment_ID'] = [1]

        result = dataloader.assign_cn_to_snv(
            mutations,
            copy_number,
            use_supplied_segment_ids=True,
        )

        self.assertEqual(result['Segment_ID'].tolist(), ['1', '1'])
        self.assertEqual(result['Total_CN'].tolist(), [3, 3])

    def test_segment_sorting_is_numeric_then_sex_then_other_and_coordinate_ordered(self):
        table = pd.DataFrame({
            'Chromosome': ['other', 'X', '10', '2', '1', '1', 'Y', 'W', 'Z'],
            'Segment_Start': [0, 0, 0, 0, 10, 0, 0, 0, 0],
            'Segment_End': [1, 1, 1, 1, 20, 10, 1, 1, 1],
        })

        result = dataloader.sort_copy_number_segments(table)

        self.assertEqual(result['Chromosome'].tolist(), [
            '1', '1', '2', '10', 'X', 'Y', 'Z', 'W', 'other',
        ])
        self.assertEqual(result.loc[:1, 'Segment_Start'].tolist(), [0, 10])

    def test_nonoverlap_validation_accepts_adjacency_and_rejects_nested_overlap(self):
        adjacent = pd.DataFrame({
            'Chromosome': ['2', '1', '1'],
            'Segment_Start': [0, 10, 0],
            'Segment_End': [100, 20, 10],
        })
        result = dataloader.validate_non_overlapping_segments(adjacent)
        self.assertEqual(result['Chromosome'].tolist(), ['1', '1', '2'])

        nested = pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Segment_Start': [0, 10, 20],
            'Segment_End': [100, 30, 40],
        })
        with self.assertRaisesRegex(
            ValueError,
            r'1:0-100 and 1:10-30.*1:0-100 and 1:20-40',
        ):
            dataloader.validate_non_overlapping_segments(nested)

    def test_generated_segment_ids_do_not_modify_source_table(self):
        table = self.copy_number_table()
        result = dataloader.generate_segment_ids(table)
        self.assertEqual(result['Segment_ID'].tolist(), ['1-0-100', '2-100-300'])
        self.assertEqual(table['Segment_ID'].tolist(), ['segment-a', 'segment-b'])

    def test_seeded_half_open_assignment_matches_every_generated_boundary_case(self):
        random_generator = np.random.default_rng(20260827)
        widths = random_generator.integers(2, 40, size=24)
        gaps = random_generator.integers(0, 5, size=24)
        starts = []
        ends = []
        next_start = 0
        for width, gap in zip(widths, gaps):
            starts.append(next_start)
            next_end = next_start + int(width)
            ends.append(next_end)
            next_start = next_end + int(gap)
        minor_copy_numbers = random_generator.integers(0, 3, size=24)
        major_copy_numbers = (
            minor_copy_numbers + random_generator.integers(0, 4, size=24)
        )
        copy_number = pd.DataFrame({
            'Chromosome': ['1'] * 24,
            'Segment_Start': starts,
            'Segment_End': ends,
            'Major_CN': major_copy_numbers,
            'Minor_CN': minor_copy_numbers,
        })
        mutation_positions = []
        expected_ids = []
        for start, end in zip(starts, ends):
            mutation_positions.extend((start, end - 1))
            expected_ids.extend((f'1-{start}-{end}', f'1-{start}-{end}'))
        mutations = pd.DataFrame({
            'Chromosome': ['1'] * len(mutation_positions),
            'Position': mutation_positions,
            'Mutation_ID': [
                f'mutation-{index}' for index in range(len(mutation_positions))
            ],
        })

        assigned = dataloader.assign_cn_to_snv(
            mutations,
            copy_number,
            use_supplied_segment_ids=False,
        )

        self.assertEqual(assigned['Segment_ID'].tolist(), expected_ids)
        self.assertEqual(assigned['Position'].tolist(), mutation_positions)

    def test_seeded_adjacent_merge_equals_maximal_copy_number_runs(self):
        random_generator = np.random.default_rng(8272026)
        widths = random_generator.integers(1, 30, size=40)
        copy_number_states = random_generator.integers(0, 5, size=40)
        state_to_copy_number = {
            0: (1, 0),
            1: (1, 1),
            2: (2, 0),
            3: (2, 1),
            4: (3, 1),
        }
        starts = np.concatenate(([0], np.cumsum(widths)[:-1])).astype(int)
        ends = np.cumsum(widths).astype(int)
        copy_numbers = [
            state_to_copy_number[int(state)] for state in copy_number_states
        ]
        source_ids = [f'source-{index}' for index in range(len(widths))]
        copy_number = pd.DataFrame({
            'Chromosome': ['1'] * len(widths),
            'Segment_Start': starts,
            'Segment_End': ends,
            'Major_CN': [value[0] for value in copy_numbers],
            'Minor_CN': [value[1] for value in copy_numbers],
            'Segment_ID': source_ids,
        })

        merged, source_id_map = dataloader.merge_segments(
            copy_number,
            return_segment_id_map=True,
        )

        expected_runs = []
        run_start = int(starts[0])
        run_source_ids = [source_ids[0]]
        previous_copy_number = copy_numbers[0]
        for index in range(1, len(widths)):
            if copy_numbers[index] != previous_copy_number:
                expected_runs.append((
                    run_start,
                    int(ends[index - 1]),
                    previous_copy_number,
                    run_source_ids,
                ))
                run_start = int(starts[index])
                run_source_ids = []
                previous_copy_number = copy_numbers[index]
            run_source_ids.append(source_ids[index])
        expected_runs.append((
            run_start,
            int(ends[-1]),
            previous_copy_number,
            run_source_ids,
        ))

        self.assertEqual(len(merged), len(expected_runs))
        self.assertEqual(
            dataloader.get_segment_widths(merged).sum(),
            int(widths.sum()),
        )
        for row, expected_run in zip(
            merged.itertuples(index=False),
            expected_runs,
        ):
            start, end, copy_number_value, run_ids = expected_run
            self.assertEqual((row.Segment_Start, row.Segment_End), (start, end))
            self.assertEqual(
                (row.Major_CN, row.Minor_CN),
                copy_number_value,
            )
            expected_merged_id = f'1-{start}-{end}'
            self.assertEqual(row.Segment_ID, expected_merged_id)
            for source_id in run_ids:
                self.assertEqual(source_id_map[source_id], expected_merged_id)


class SubcloneValidationTest(unittest.TestCase):
    @staticmethod
    def table(**overrides):
        values = {
            'Cluster': ['one', 'two'],
            'Subclone_CCF': [0.8, 0.4],
            'Subclone_Fraction': [0.3, 0.2],
        }
        values.update(overrides)
        return pd.DataFrame(values)

    def test_subclone_schema_numeric_values_and_fraction_tolerance(self):
        with self.assertRaisesRegex(
            ValueError,
            r'missing required column\(s\): Subclone_CCF, Subclone_Fraction',
        ):
            dataloader.validate_subclone_values(pd.DataFrame({'Cluster': ['a']}))

        table = self.table(
            Subclone_CCF=['0.8', '0.4'],
            Subclone_Fraction=[0.6, 0.4000000000005],
        )
        result = dataloader.validate_subclone_values(table)
        self.assertEqual(result['Subclone_CCF'].dtype, np.dtype(float))
        self.assertAlmostEqual(result['Subclone_Fraction'].sum(), 1.0)
        self.assertEqual(table['Subclone_CCF'].tolist(), ['0.8', '0.4'])

    def test_subclone_numeric_columns_reject_booleans_nonfinite_range_and_large_sum(self):
        cases = (
            ('Subclone_CCF', [True, 0.2], 'Subclone_CCF'),
            ('Subclone_CCF', [float('nan'), 0.2], 'Subclone_CCF'),
            ('Subclone_CCF', [-0.1, 0.2], 'Subclone_CCF'),
            ('Subclone_Fraction', [float('inf'), 0.2], 'Subclone_Fraction'),
            ('Subclone_Fraction', [1.1, 0.2], 'Subclone_Fraction'),
        )
        for column, values, message in cases:
            with self.subTest(column=column, values=values):
                with self.assertRaisesRegex(ValueError, message):
                    dataloader.validate_subclone_values(
                        self.table(**{column: values})
                    )

        with self.assertRaisesRegex(ValueError, 'sum to no more than 1'):
            dataloader.validate_subclone_values(
                self.table(Subclone_Fraction=[0.7, 0.4])
            )

    def test_optional_ccf_clipping_warns_and_never_relaxes_other_validation(self):
        table = self.table(Subclone_CCF=[-0.25, 1.5])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.validate_subclone_values(
                table,
                clip_subclone_ccf=True,
            )

        self.assertEqual(result['Subclone_CCF'].tolist(), [0.0, 1.0])
        self.assertEqual(len(caught), 1)
        self.assertIn('Clipping 2 Subclone_CCF values', str(caught[0].message))
        self.assertIn('1 below 0; 1 above 1', str(caught[0].message))

        with self.assertRaisesRegex(ValueError, 'clip_subclone_ccf must be a boolean'):
            dataloader.validate_subclone_values(table, clip_subclone_ccf='yes')
        with self.assertRaisesRegex(ValueError, 'Subclone_CCF must contain finite'):
            dataloader.validate_subclone_values(
                self.table(Subclone_CCF=[float('nan'), 2]),
                clip_subclone_ccf=True,
            )

    def test_valid_subclone_filter_has_inclusive_ccf_and_strict_fraction_thresholds(self):
        table = pd.DataFrame({
            'Cluster': ['low', 'middle', 'high', 'outside'],
            'Subclone_CCF': [0.1, 0.5, 0.9, 0.95],
            'Subclone_Fraction': [0.1, 0.2, 0.7, 1.0],
        })

        result = dataloader.get_valid_subclones(
            table,
            min_ccf=0.1,
            max_ccf=0.9,
            min_fraction=0.1,
        )

        self.assertEqual(result['Cluster'].tolist(), ['middle', 'high'])
        self.assertEqual(table['Cluster'].tolist(), ['low', 'middle', 'high', 'outside'])

    def test_excess_subclones_preserve_top_clone_and_weight_combine_the_rest(self):
        self.assertIsNone(dataloader.filter_excess_subclones(None))
        singleton = self.table().iloc[[0]]
        self.assertIs(dataloader.filter_excess_subclones(singleton), singleton)

        table = pd.DataFrame({
            'Cluster': ['low', 'top', 'middle'],
            'Subclone_CCF': [0.2, 0.9, 0.5],
            'Subclone_Fraction': [0.1, 0.4, 0.3],
        })
        result = dataloader.filter_excess_subclones(table)

        self.assertEqual(result['Cluster'].tolist(), ['top', 'middle'])
        self.assertEqual(result['Subclone_Fraction'].tolist(), [0.4, 0.4])
        self.assertAlmostEqual(result['Subclone_CCF'].iloc[1], 0.425)
        self.assertEqual(table['Cluster'].tolist(), ['low', 'top', 'middle'])


class DataloaderIOAndPloidyTest(DataloaderFixtures, unittest.TestCase):
    def write_tables(self, directory, copy_number, mutations):
        copy_path = Path(directory) / 'copy-number.tsv'
        mutation_path = Path(directory) / 'mutations.tsv'
        copy_number.to_csv(copy_path, sep='\t', index=False)
        mutations.to_csv(mutation_path, sep='\t', index=False)
        return copy_path, mutation_path

    def test_file_loading_preserves_identifier_text_and_discards_stale_annotations(self):
        copy_number = pd.DataFrame({
            'Chromosome': ['chr1'],
            'Segment_Start': ['000'],
            'Segment_End': ['100'],
            'Major_CN': ['2'],
            'Minor_CN': ['1'],
            'Segment_ID': ['001'],
        })
        mutations = pd.DataFrame({
            'Chromosome': ['chr1'],
            'Position': ['010'],
            'Tumor_Ref_Count': ['09'],
            'Tumor_Alt_Count': ['03'],
            'Mutation_ID': ['0007'],
            'Segment_ID': ['001'],
            'Segment_Start': [999],
            'Segment_End': [1000],
            'Major_CN': [99],
            'Minor_CN': [98],
        })
        with tempfile.TemporaryDirectory() as directory:
            paths = self.write_tables(directory, copy_number, mutations)
            loaded_cn, loaded_mutations = dataloader.load_input_tables(
                *paths,
                sex='XX',
            )

        self.assertEqual(loaded_cn['Chromosome'].tolist(), ['1'])
        self.assertEqual(loaded_cn['Segment_ID'].tolist(), ['001'])
        self.assertEqual(loaded_cn['Segment_Start'].tolist(), [0])
        self.assertEqual(loaded_mutations['Mutation_ID'].tolist(), ['0007'])
        self.assertEqual(loaded_mutations['Segment_ID'].tolist(), ['001'])
        self.assertEqual(loaded_mutations['Position'].tolist(), [10])
        self.assertEqual(loaded_mutations['Tumor_Ref_Count'].tolist(), [9])
        for column in dataloader.MUTATION_COPY_NUMBER_ANNOTATION_COLUMNS:
            self.assertNotIn(column, loaded_mutations.columns)

    def test_file_loading_position_mode_validates_and_sorts_segments(self):
        copy_number = self.copy_number_table(
            Chromosome=['2', '1'],
            Segment_Start=[100, 0],
            Segment_End=[300, 100],
        ).drop(columns='Segment_ID')
        mutations = self.mutation_table().drop(columns='Segment_ID')
        with tempfile.TemporaryDirectory() as directory:
            paths = self.write_tables(directory, copy_number, mutations)
            loaded_cn, loaded_mutations = dataloader.load_input_tables(
                *paths,
                sex='XX',
            )
        self.assertEqual(loaded_cn['Chromosome'].tolist(), ['1', '2'])
        self.assertEqual(loaded_mutations['Position'].tolist(), [50, 150])

    def test_file_loading_can_warn_and_drop_unknown_supplied_ids(self):
        mutations = self.mutation_table(
            Segment_ID=['segment-a', 'missing'],
        )
        with tempfile.TemporaryDirectory() as directory:
            paths = self.write_tables(directory, self.copy_number_table(), mutations)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                _, loaded_mutations = dataloader.load_input_tables(
                    *paths,
                    sex='XX',
                    drop_unmatched_snvs=True,
                )
        self.assertEqual(loaded_mutations['Mutation_ID'].tolist(), ['mutation-a'])
        self.assertEqual(len(caught), 1)
        self.assertIn('Dropping 1 unmatched SNV', str(caught[0].message))

    def test_ploidy_and_normal_ploidy_are_width_weighted(self):
        copy_number = pd.DataFrame({
            'Chromosome': ['1', 'X'],
            'Segment_Start': [0, 0],
            'Segment_End': [100, 300],
            'Major_CN': [3, 1],
            'Minor_CN': [1, 0],
        })

        self.assertAlmostEqual(
            dataloader.calculate_ploidy(copy_number, sex='XY'),
            1.75,
        )
        self.assertAlmostEqual(
            dataloader.calculate_normal_ploidy(copy_number, sex='XY'),
            1.25,
        )

    def test_nrpcc_returns_expected_coverage_and_tumor_ploidy(self):
        copy_number = pd.DataFrame({
            'Chromosome': ['1'],
            'Segment_Start': [0],
            'Segment_End': [100],
            'Major_CN': [3],
            'Minor_CN': [1],
        })
        mutations = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Tumor_Ref_Count': [10, 20],
            'Tumor_Alt_Count': [10, 20],
        })

        nrpcc, coverage, tumor_ploidy = dataloader.calculate_nrpcc(
            copy_number,
            mutations,
            purity=0.5,
            sex='XX',
        )

        self.assertEqual(coverage, 30)
        self.assertEqual(tumor_ploidy, 4)
        self.assertEqual(nrpcc, 5)


if __name__ == '__main__':
    unittest.main()
