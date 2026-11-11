import unittest
import warnings

import pandas as pd

from gritic import dataloader, sampletools


class InputPreparationTest(unittest.TestCase):
    @staticmethod
    def copy_number_table():
        return pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Segment_Start': [0, 100, 200],
            'Segment_End': [100, 200, 300],
            'Major_CN': [2, 0, 2],
            'Minor_CN': [1, 0, 1],
            'Segment_ID': ['left', 'zero', 'right'],
        })

    @staticmethod
    def mutation_table(*, supplied_ids=True):
        table = pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Position': [50, 150, 250],
            'Tumor_Ref_Count': [10, 10, 10],
            'Tumor_Alt_Count': [5, 5, 5],
            'Mutation_ID': ['left-snv', 'zero-snv', 'right-snv'],
        })
        if supplied_ids:
            table['Segment_ID'] = ['left', 'zero', 'right']
        return table

    def prepare(self, cn_table=None, mutation_table=None, **kwargs):
        kwargs.setdefault('min_mutation_count', 1)
        return dataloader.prepare_sample_inputs(
            self.copy_number_table() if cn_table is None else cn_table,
            self.mutation_table() if mutation_table is None else mutation_table,
            sex='XX',
            **kwargs,
        )

    def test_supplied_zero_copy_drop_merges_flanks_and_preserves_source_ids(self):
        with self.assertWarnsRegex(UserWarning, 'Dropping 1 SNV.*zero-copy'):
            prepared = self.prepare()

        self.assertEqual(prepared.sex, 'XX')
        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'],
        )
        self.assertEqual(
            prepared.supplied_segment_id_map,
            {'left': '1-0-300', 'right': '1-0-300'},
        )
        mutations = prepared.mutation_table.set_index('Mutation_ID')
        self.assertEqual(set(mutations.index), {'left-snv', 'right-snv'})
        self.assertEqual(mutations.loc['left-snv', 'Source_Segment_ID'], 'left')
        self.assertEqual(mutations.loc['right-snv', 'Source_Segment_ID'], 'right')
        self.assertEqual(
            mutations.loc['left-snv', 'GRITIC_Mutation_ID'], 'left:left-snv',
        )
        self.assertEqual(set(mutations['Segment_ID']), {'1-0-300'})

    def test_position_zero_copy_drop_is_independent_of_unmatched_policy(self):
        for drop_unmatched_snvs in (False, True):
            with self.subTest(drop_unmatched_snvs=drop_unmatched_snvs):
                with self.assertWarnsRegex(UserWarning, 'Dropping 1 SNV.*zero-copy'):
                    prepared = self.prepare(
                        mutation_table=self.mutation_table(supplied_ids=False),
                        drop_unmatched_snvs=drop_unmatched_snvs,
                    )

                self.assertIsNone(prepared.supplied_segment_id_map)
                mutations = prepared.mutation_table
                self.assertEqual(set(mutations['Position']), {50, 250})
                self.assertEqual(set(mutations['Major_CN']), {2})
                self.assertEqual(set(mutations['Source_Segment_ID']), {'1-0-300'})
                self.assertEqual(
                    set(mutations['GRITIC_Mutation_ID']),
                    {'1-0-300:left-snv', '1-0-300:right-snv'},
                )

    def test_supplied_ids_determine_zero_copy_ownership_even_with_positions(self):
        mutations = self.mutation_table()
        mutations['Position'] = [150, 50, 250]

        with self.assertWarnsRegex(UserWarning, 'zero-copy'):
            prepared = self.prepare(mutation_table=mutations)

        self.assertEqual(
            set(prepared.mutation_table['Mutation_ID']),
            {'left-snv', 'right-snv'},
        )
        self.assertEqual(
            prepared.mutation_table.set_index('Mutation_ID').loc[
                'left-snv', 'Position'
            ],
            150,
        )

    def test_zero_copy_coordinate_exclusion_uses_half_open_original_intervals(self):
        mutations = pd.DataFrame({
            'Chromosome': ['1'] * 4,
            'Position': [99, 100, 199, 200],
            'Tumor_Ref_Count': [10] * 4,
            'Tumor_Alt_Count': [5] * 4,
        })

        with self.assertWarnsRegex(UserWarning, 'Dropping 2 SNVs.*zero-copy'):
            prepared = self.prepare(mutation_table=mutations)

        self.assertEqual(set(prepared.mutation_table['Position']), {99, 200})
        self.assertEqual(
            set(prepared.mutation_table['GRITIC_Mutation_ID']),
            {'1-0-300:99', '1-0-300:200'},
        )
        self.assertEqual(
            sorted(prepared.mutation_table['Segment_Mutation_Index']), [0, 1],
        )

    def test_consecutive_zero_copy_segments_are_excluded_before_merging(self):
        cn_table = pd.DataFrame({
            'Chromosome': ['1'] * 4,
            'Segment_Start': [0, 100, 200, 300],
            'Segment_End': [100, 200, 300, 400],
            'Major_CN': [2, 0, 0, 2],
            'Minor_CN': [1, 0, 0, 1],
            'Segment_ID': ['left', 'zero-a', 'zero-b', 'right'],
        })
        for supplied_ids in (False, True):
            mutations = pd.DataFrame({
                'Chromosome': ['1'] * 4,
                'Position': [50, 150, 250, 350],
                'Tumor_Ref_Count': [10] * 4,
                'Tumor_Alt_Count': [5] * 4,
            })
            if supplied_ids:
                mutations['Segment_ID'] = cn_table['Segment_ID']
            with self.subTest(supplied_ids=supplied_ids):
                with self.assertWarnsRegex(UserWarning, 'Dropping 2 SNVs.*zero-copy'):
                    prepared = self.prepare(cn_table, mutations)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-400'],
                )
                self.assertEqual(set(prepared.mutation_table['Position']), {50, 350})

    def test_maximum_merge_gap_counts_the_dropped_interval_and_is_inclusive(self):
        for supplied_ids in (False, True):
            for gap, expected in (
                (None, ['1-0-300']),
                (0, ['1-0-100', '1-200-300']),
                (99, ['1-0-100', '1-200-300']),
                (100, ['1-0-300']),
            ):
                with self.subTest(supplied_ids=supplied_ids, gap=gap):
                    with self.assertWarnsRegex(UserWarning, 'zero-copy'):
                        prepared = self.prepare(
                            mutation_table=self.mutation_table(
                                supplied_ids=supplied_ids,
                            ),
                            max_merge_gap=gap,
                        )

                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(),
                        expected,
                    )
                    self.assertEqual(set(prepared.mutation_table['Position']), {50, 250})

    def test_merge_opt_out_preserves_positive_input_segments(self):
        with self.assertWarnsRegex(UserWarning, 'zero-copy'):
            prepared = self.prepare(merge_cn=False)

        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(),
            ['1-0-100', '1-200-300'],
        )
        self.assertEqual(
            prepared.supplied_segment_id_map,
            {'left': '1-0-100', 'right': '1-200-300'},
        )

    def test_uncovered_gap_inside_merged_interval_retains_current_assignment(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        prepared = self.prepare(
            cn_table,
            self.mutation_table(supplied_ids=False),
        )

        self.assertEqual(prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'])
        self.assertEqual(set(prepared.mutation_table['Position']), {50, 150, 250})

    def test_gap_only_mutations_do_not_seed_empty_source_intervals(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        mutations = self.mutation_table(supplied_ids=False).iloc[[1]].copy()
        for drop_unmatched_snvs in (False, True):
            with self.subTest(drop_unmatched_snvs=drop_unmatched_snvs):
                with self.assertRaisesRegex(
                    ValueError,
                    'No copy-number segments.*before merging',
                ):
                    self.prepare(
                        cn_table, mutations,
                        drop_unmatched_snvs=drop_unmatched_snvs,
                    )

    def test_gap_mutation_does_not_extend_a_supported_region_to_an_empty_flank(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        mutations = self.mutation_table(supplied_ids=False).iloc[[0, 1]].copy()

        with self.assertRaisesRegex(ValueError, '1 unmatched SNV'):
            self.prepare(cn_table, mutations)
        with self.assertWarnsRegex(UserWarning, '1 unmatched SNV'):
            prepared = self.prepare(cn_table, mutations, drop_unmatched_snvs=True)

        self.assertEqual(prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-100'])
        self.assertEqual(prepared.mutation_table['Position'].tolist(), [50])

    def test_all_read_filtered_mutations_raise_before_segment_merging(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        for supplied_ids in (False, True):
            mutations = self.mutation_table(
                supplied_ids=supplied_ids,
            ).iloc[[0, 2]].copy()
            mutations['Tumor_Alt_Count'] = 2
            with self.subTest(supplied_ids=supplied_ids):
                with self.assertRaisesRegex(
                    ValueError, 'No mutations remain after applying.*thresholds',
                ):
                    self.prepare(cn_table, mutations)

    def test_supplied_duplicate_mutation_labels_remain_distinct_after_merge(self):
        mutations = self.mutation_table()
        mutations['Mutation_ID'] = ['repeated', 'zero-snv', 'repeated']

        with self.assertWarnsRegex(UserWarning, 'zero-copy'):
            prepared = self.prepare(mutation_table=mutations)

        self.assertEqual(
            set(prepared.mutation_table['GRITIC_Mutation_ID']),
            {'left:repeated', 'right:repeated'},
        )
        self.assertEqual(
            sorted(prepared.mutation_table['Segment_Mutation_Index']), [0, 1],
        )

    def test_position_duplicate_labels_are_checked_in_final_segment_scope(self):
        mutations = self.mutation_table(supplied_ids=False)
        mutations['Mutation_ID'] = ['repeated', 'zero-snv', 'repeated']

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with self.assertRaisesRegex(ValueError, 'unique.*source Segment_ID'):
                self.prepare(mutation_table=mutations)

    def test_empty_positive_intervals_are_removed_before_flanks_merge(self):
        for supplied_ids in (False, True):
            for major_cn, minor_cn in ((1, 0), (9, 1), (2, 0)):
                cn_table = self.copy_number_table()
                cn_table['Major_CN'] = [3, major_cn, 3]
                cn_table['Minor_CN'] = [1, minor_cn, 1]
                mutations = self.mutation_table(
                    supplied_ids=supplied_ids,
                ).iloc[[0, 2]].copy()
                with self.subTest(
                    supplied_ids=supplied_ids, major_cn=major_cn,
                    minor_cn=minor_cn,
                ):
                    prepared = self.prepare(cn_table, mutations)

                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(),
                        ['1-0-300'],
                    )
                    self.assertEqual(len(prepared.mutation_table), 2)
                    self.assertEqual(
                        prepared.source_copy_number_table['Major_CN'].tolist(),
                        [3, major_cn, 3],
                    )
                    if supplied_ids:
                        self.assertEqual(
                            prepared.supplied_segment_id_map,
                            {'left': '1-0-300', 'right': '1-0-300'},
                        )
                        self.assertEqual(
                            set(prepared.mutation_table['GRITIC_Mutation_ID']),
                            {'left:left-snv', 'right:right-snv'},
                        )
                    else:
                        self.assertIsNone(prepared.supplied_segment_id_map)

    def test_read_thresholds_determine_whether_an_interval_is_empty_before_merging(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [1, 0]
        for supplied_ids in (False, True):
            for ref_count, alt_count, min_alt, min_coverage, retained in (
                (8, 2, 3, 10, False),
                (6, 3, 3, 10, False),
                (7, 3, 3, 10, True),
                (7, 3, 4, 10, False),
                (7, 3, 3, 11, False),
            ):
                mutations = self.mutation_table(supplied_ids=supplied_ids)
                mutations.loc[1, ['Tumor_Ref_Count', 'Tumor_Alt_Count']] = [
                    ref_count, alt_count,
                ]
                with self.subTest(
                    supplied_ids=supplied_ids, ref_count=ref_count,
                    alt_count=alt_count, min_alt=min_alt,
                    min_coverage=min_coverage,
                ):
                    prepared = self.prepare(
                        cn_table, mutations,
                        min_mutation_alt_count=min_alt,
                        min_mutation_coverage=min_coverage,
                    )

                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(),
                        ['1-0-100', '1-100-200', '1-200-300']
                        if retained else ['1-0-300'],
                    )
                    self.assertEqual(
                        set(prepared.mutation_table['Position']),
                        {50, 150, 250} if retained else {50, 250},
                    )

    def test_supported_intervening_states_remain_barriers_regardless_of_timing_eligibility(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [9, 1]
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                prepared = self.prepare(
                    cn_table, self.mutation_table(supplied_ids=supplied_ids),
                )

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-100-200', '1-200-300'],
                )
                self.assertEqual(prepared.copy_number_table['Major_CN'].tolist(), [2, 9, 2])

    def test_empty_intervals_are_dropped_when_merging_is_disabled(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [5, 1]
        for supplied_ids in (False, True):
            mutations = self.mutation_table(
                supplied_ids=supplied_ids,
            ).iloc[[0, 2]].copy()
            with self.subTest(supplied_ids=supplied_ids):
                prepared = self.prepare(cn_table, mutations, merge_cn=False)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-200-300'],
                )
                self.assertEqual(len(prepared.source_copy_number_table), 3)
                if supplied_ids:
                    self.assertEqual(
                        prepared.supplied_segment_id_map,
                        {'left': '1-0-100', 'right': '1-200-300'},
                    )

    def test_maximum_merge_gap_includes_empty_positive_intervals(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [5, 1]
        for supplied_ids in (False, True):
            mutations = self.mutation_table(
                supplied_ids=supplied_ids,
            ).iloc[[0, 2]].copy()
            for gap, expected in (
                (0, ['1-0-100', '1-200-300']),
                (99, ['1-0-100', '1-200-300']),
                (100, ['1-0-300']),
            ):
                with self.subTest(supplied_ids=supplied_ids, gap=gap):
                    prepared = self.prepare(cn_table, mutations, max_merge_gap=gap)

                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(),
                        expected,
                    )

    def test_coordinate_mutation_presence_uses_original_half_open_intervals(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [5, 1]
        mutations = self.mutation_table(supplied_ids=False).iloc[[0, 2]].copy()
        mutations['Position'] = [99, 200]

        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'],
        )
        self.assertEqual(prepared.mutation_table['Major_CN'].tolist(), [2, 2])
        self.assertEqual(set(prepared.mutation_table['Position']), {99, 200})

    def test_supplied_ids_determine_mutation_presence_even_when_positions_disagree(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [5, 1]
        mutations = self.mutation_table().iloc[[0, 2]].copy()
        mutations['Position'] = [150, 150]

        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'],
        )
        self.assertEqual(
            prepared.supplied_segment_id_map,
            {'left': '1-0-300', 'right': '1-0-300'},
        )
        self.assertEqual(set(prepared.mutation_table['Source_Segment_ID']), {'left', 'right'})

    def test_empty_same_cn_edges_do_not_extend_a_supported_interval(self):
        cn_table = self.copy_number_table()
        cn_table['Major_CN'] = 2
        cn_table['Minor_CN'] = 1
        for supplied_ids in (False, True):
            mutations = self.mutation_table(
                supplied_ids=supplied_ids,
            ).iloc[[1]].copy()
            with self.subTest(supplied_ids=supplied_ids):
                prepared = self.prepare(cn_table, mutations)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-100-200'],
                )
                self.assertEqual(len(prepared.source_copy_number_table), 3)

    def test_true_unmatched_mutations_still_require_the_drop_option(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        mutations = self.mutation_table(supplied_ids=False)
        mutations.loc[1, 'Position'] = 300

        with self.assertRaisesRegex(ValueError, '1 unmatched SNV'):
            self.prepare(cn_table, mutations)
        with self.assertWarnsRegex(UserWarning, '1 unmatched SNV'):
            prepared = self.prepare(cn_table, mutations, drop_unmatched_snvs=True)

        self.assertEqual(set(prepared.mutation_table['Position']), {50, 250})

    def test_chromosome_filtering_applies_to_both_original_tables(self):
        cn_table = self.copy_number_table().iloc[[0, 2]].copy()
        cn_table.loc[2, 'Chromosome'] = 'chr23'
        mutations = self.mutation_table().iloc[[0, 2]].copy()
        mutations.loc[2, 'Chromosome'] = 'chr23'

        with self.assertRaisesRegex(ValueError, 'Chromosome'):
            self.prepare(cn_table, mutations)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            prepared = self.prepare(
                cn_table, mutations, drop_unmatched_chromosomes=True,
            )

        self.assertTrue(any('unmatched Chromosome' in str(w.message) for w in caught))
        self.assertEqual(prepared.source_copy_number_table['Segment_ID'].tolist(), ['left'])
        self.assertEqual(prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-100'])
        self.assertEqual(prepared.mutation_table['Mutation_ID'].tolist(), ['left-snv'])

    def test_all_zero_copy_mutations_raise_after_contextual_warning(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    with self.assertRaisesRegex(ValueError, r'zero-copy|0\+0'):
                        self.prepare(
                            self.copy_number_table().iloc[[1]].copy(),
                            self.mutation_table(supplied_ids=supplied_ids).iloc[[1]].copy(),
                        )

                self.assertTrue(any('zero-copy' in str(w.message) for w in caught))

    def test_preparation_preserves_original_profile_without_mutating_inputs(self):
        cn_table = self.copy_number_table()
        cn_table['Chromosome'] = 'chr1'
        mutations = self.mutation_table()
        mutations['Chromosome'] = 'chr1'
        original_cn = cn_table.copy(deep=True)
        original_mutations = mutations.copy(deep=True)

        with self.assertWarnsRegex(UserWarning, 'zero-copy'):
            prepared = self.prepare(cn_table, mutations)

        pd.testing.assert_frame_equal(cn_table, original_cn)
        pd.testing.assert_frame_equal(mutations, original_mutations)
        self.assertEqual(prepared.source_copy_number_table['Chromosome'].tolist(), ['1'] * 3)
        self.assertEqual(prepared.source_copy_number_table['Major_CN'].tolist(), [2, 0, 2])
        self.assertEqual(prepared.source_copy_number_table['Segment_ID'].tolist(), ['left', 'zero', 'right'])
        self.assertEqual(dataloader.calculate_ploidy(prepared.source_copy_number_table, sex='XX'), 2.0)

    def test_sample_keeps_original_cn_profile_for_modal_major_cn(self):
        cn_table = self.copy_number_table()
        cn_table['Segment_Start'] = [0, 100, 1100]
        cn_table['Segment_End'] = [100, 1100, 1200]
        mutations = self.mutation_table()
        mutations['Position'] = [50, 150, 1150]

        with self.assertWarnsRegex(UserWarning, 'zero-copy'):
            sample = sampletools.Sample(
                mutations, cn_table, None, 'TEST', 0.8, sex='XX',
                min_mutation_count=1,
            )

        self.assertEqual(sample.copy_number_table['Segment_ID'].tolist(), ['1-0-1200'])
        self.assertEqual(sample.cn_table['Major_CN'].tolist(), [2, 0, 2])
        self.assertEqual(sampletools.get_major_cn_mode(sample), 0)

    def test_sample_keeps_empty_positive_intervals_for_original_major_cn_mode(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [5, 1]
        cn_table['Segment_Start'] = [0, 100, 1100]
        cn_table['Segment_End'] = [100, 1100, 1200]
        mutations = self.mutation_table().iloc[[0, 2]].copy()

        sample = sampletools.Sample(
            mutations, cn_table, None, 'TEST', 0.8, sex='XX',
            min_mutation_count=1,
        )

        self.assertEqual(sample.copy_number_table['Segment_ID'].tolist(), ['1-0-1200'])
        self.assertEqual(sample.cn_table['Major_CN'].tolist(), [2, 5, 2])
        self.assertEqual(sampletools.get_major_cn_mode(sample), 5)
        self.assertEqual(len(sample.segments), 1)
        self.assertEqual(sample.segments[0].n_mutations, 2)


class MinimumMutationCountTest(unittest.TestCase):
    @staticmethod
    def tables(counts=(4, 3, 4), *, supplied_ids=True):
        cn_table = pd.DataFrame({
            'Chromosome': ['1'] * 3,
            'Segment_Start': [0, 100, 200],
            'Segment_End': [100, 200, 300],
            'Major_CN': [3, 1, 3],
            'Minor_CN': [1, 1, 1],
            'Segment_ID': ['left', 'island', 'right'],
        })
        mutations = pd.DataFrame([
            {
                'Chromosome': '1',
                'Position': 100 * segment_index + mutation_index,
                'Tumor_Ref_Count': 10,
                'Tumor_Alt_Count': 5,
                'Mutation_ID': f'{segment_id}-{mutation_index}',
                'Segment_ID': segment_id,
            }
            for segment_index, (segment_id, count) in enumerate(
                zip(cn_table['Segment_ID'], counts)
            )
            for mutation_index in range(count)
        ])
        if not supplied_ids:
            mutations = mutations.drop(columns='Segment_ID')
        return cn_table, mutations

    @staticmethod
    def prepare(cn_table, mutations, **kwargs):
        return dataloader.prepare_sample_inputs(
            cn_table, mutations, sex='XX', **kwargs,
        )

    def test_default_drops_three_mutation_source_and_its_mutations_before_merging(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                cn_table, mutations = self.tables(supplied_ids=supplied_ids)
                prepared = self.prepare(cn_table, mutations)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 8)
                self.assertFalse(
                    prepared.mutation_table['Mutation_ID'].str.startswith('island-').any()
                )
                self.assertEqual(len(prepared.source_copy_number_table), 3)
                if supplied_ids:
                    self.assertEqual(
                        prepared.supplied_segment_id_map,
                        {'left': '1-0-300', 'right': '1-0-300'},
                    )

    def test_four_mutation_source_meets_default_and_remains_a_cn_barrier(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                prepared = self.prepare(*self.tables((4, 4, 4), supplied_ids=supplied_ids))

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-100-200', '1-200-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 12)

    def test_custom_threshold_is_inclusive(self):
        for threshold in (1, 3, 5):
            for supplied_ids in (False, True):
                with self.subTest(threshold=threshold, supplied_ids=supplied_ids):
                    counts = (threshold, threshold - 1, threshold)
                    prepared = self.prepare(
                        *self.tables(counts, supplied_ids=supplied_ids),
                        min_mutation_count=threshold,
                    )
                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'],
                    )
                    self.assertEqual(len(prepared.mutation_table), 2 * threshold)

        prepared = self.prepare(*self.tables(), min_mutation_count=3)
        self.assertEqual(len(prepared.copy_number_table), 3)
        self.assertEqual(len(prepared.mutation_table), 11)

    def test_zero_threshold_disables_source_count_filter(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                cn_table, mutations = self.tables((4, 0, 4), supplied_ids=supplied_ids)
                prepared = self.prepare(cn_table, mutations, min_mutation_count=0)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-100-200', '1-200-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 8)

    def test_read_filters_are_applied_before_counting_source_support(self):
        for supplied_ids in (False, True):
            for ref_count, alt_count in ((10, 2), (6, 3), (7, 3)):
                with self.subTest(
                    supplied_ids=supplied_ids, ref_count=ref_count, alt_count=alt_count,
                ):
                    cn_table, mutations = self.tables((4, 4, 4), supplied_ids=supplied_ids)
                    mutations.loc[
                        mutations['Mutation_ID'].eq('island-3'),
                        ['Tumor_Ref_Count', 'Tumor_Alt_Count'],
                    ] = [ref_count, alt_count]
                    prepared = self.prepare(cn_table, mutations)
                    island_retained = ref_count + alt_count >= 10 and alt_count >= 3
                    self.assertEqual(len(prepared.copy_number_table), 3 if island_retained else 1)
                    self.assertEqual(len(prepared.mutation_table), 12 if island_retained else 8)

    def test_adjacent_equal_cn_sources_cannot_pool_support_to_pass_threshold(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                cn_table, mutations = self.tables((3, 3, 4), supplied_ids=supplied_ids)
                cn_table['Major_CN'] = 3
                prepared = self.prepare(cn_table, mutations)
                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(), ['1-200-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 4)

    def test_count_filter_still_applies_when_adjacent_segments_are_kept(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                prepared = self.prepare(
                    *self.tables(supplied_ids=supplied_ids), merge_cn=False,
                )
                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-200-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 8)

    def test_removed_source_width_counts_toward_maximum_merge_gap(self):
        for supplied_ids in (False, True):
            for gap in (0, 99, 100):
                with self.subTest(supplied_ids=supplied_ids, gap=gap):
                    prepared = self.prepare(
                        *self.tables(supplied_ids=supplied_ids), max_merge_gap=gap,
                    )
                    self.assertEqual(
                        prepared.copy_number_table['Segment_ID'].tolist(),
                        ['1-0-300'] if gap == 100 else ['1-0-100', '1-200-300'],
                    )
                    self.assertEqual(len(prepared.mutation_table), 8)

    def test_truly_uncovered_gap_mutations_are_retained_after_flanks_merge(self):
        cn_table, mutations = self.tables(supplied_ids=False)
        cn_table = cn_table.iloc[[0, 2]].copy()
        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(prepared.copy_number_table['Segment_ID'].tolist(), ['1-0-300'])
        self.assertEqual(len(prepared.mutation_table), 11)

    def test_gap_mutations_do_not_count_as_source_support(self):
        cn_table, mutations = self.tables((3, 20, 4), supplied_ids=False)
        cn_table = cn_table.iloc[[0, 2]].copy()
        with self.assertWarnsRegex(UserWarning, '20 unmatched SNVs'):
            prepared = self.prepare(cn_table, mutations, drop_unmatched_snvs=True)

        self.assertEqual(prepared.copy_number_table['Segment_ID'].tolist(), ['1-200-300'])
        self.assertEqual(len(prepared.mutation_table), 4)

    def test_supplied_ids_control_counting_and_drops_when_positions_disagree(self):
        cn_table, mutations = self.tables()
        mutations['Position'] = 150
        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(len(prepared.mutation_table), 8)
        self.assertEqual(set(prepared.mutation_table['Source_Segment_ID']), {'left', 'right'})

    def test_half_open_boundary_belongs_to_following_source(self):
        cn_table, mutations = self.tables((4, 3, 4), supplied_ids=False)
        mutations.loc[mutations['Mutation_ID'].eq('left-3'), 'Position'] = 100
        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(),
            ['1-100-200', '1-200-300'],
        )
        self.assertEqual(len(prepared.mutation_table), 8)
        self.assertIn('left-3', set(prepared.mutation_table['Mutation_ID']))

    def test_no_source_meeting_threshold_fails_before_merging(self):
        for supplied_ids in (False, True):
            with self.subTest(supplied_ids=supplied_ids):
                with self.assertRaisesRegex(ValueError, 'No copy-number segments.*before merging'):
                    self.prepare(*self.tables((3, 3, 3), supplied_ids=supplied_ids))

    def test_sample_uses_default_threshold_and_passes_custom_threshold(self):
        cn_table, mutations = self.tables()
        for kwargs, expected in (({}, [8]), ({'min_mutation_count': 3}, [4, 3, 4])):
            with self.subTest(kwargs=kwargs):
                sample = sampletools.Sample(
                    mutations, cn_table, None, 'TEST', 0.8, sex='XX', **kwargs,
                )
                self.assertEqual([segment.n_mutations for segment in sample.segments], expected)

    def test_threshold_rejects_invalid_programmatic_values(self):
        cn_table, mutations = self.tables()
        for threshold in (-1, 1.0, True, '4', None):
            for prepare in (
                lambda: self.prepare(cn_table, mutations, min_mutation_count=threshold),
                lambda: sampletools.Sample(
                    mutations, cn_table, None, 'TEST', 0.8, sex='XX',
                    min_mutation_count=threshold,
                ),
            ):
                with self.subTest(threshold=threshold):
                    with self.assertRaisesRegex(ValueError, 'min_mutation_count'):
                        prepare()


if __name__ == '__main__':
    unittest.main()
