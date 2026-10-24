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

    def test_positive_intervening_states_remain_merge_barriers_without_mutations(self):
        # Timing eligibility and mutation presence must not filter the CN profile.
        for major_cn, minor_cn in ((1, 0), (9, 1), (2, 0)):
            cn_table = self.copy_number_table()
            cn_table['Major_CN'] = [3, major_cn, 3]
            cn_table['Minor_CN'] = [1, minor_cn, 1]
            mutations = self.mutation_table().iloc[[0, 2]].copy()
            with self.subTest(major_cn=major_cn, minor_cn=minor_cn):
                prepared = self.prepare(cn_table, mutations)

                self.assertEqual(
                    prepared.copy_number_table['Segment_ID'].tolist(),
                    ['1-0-100', '1-100-200', '1-200-300'],
                )
                self.assertEqual(len(prepared.mutation_table), 2)

    def test_read_filtering_does_not_remove_a_positive_cn_merge_barrier(self):
        cn_table = self.copy_number_table()
        cn_table.loc[1, ['Major_CN', 'Minor_CN']] = [1, 0]
        mutations = self.mutation_table()
        mutations.loc[1, 'Tumor_Alt_Count'] = 2

        prepared = self.prepare(cn_table, mutations)

        self.assertEqual(
            prepared.copy_number_table['Segment_ID'].tolist(),
            ['1-0-100', '1-100-200', '1-200-300'],
        )
        self.assertEqual(set(prepared.mutation_table['Position']), {50, 250})

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
            )

        self.assertEqual(sample.copy_number_table['Segment_ID'].tolist(), ['1-0-1200'])
        self.assertEqual(sample.cn_table['Major_CN'].tolist(), [2, 0, 2])
        self.assertEqual(sampletools.get_major_cn_mode(sample), 0)


if __name__ == '__main__':
    unittest.main()
