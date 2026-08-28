import unittest

import numpy as np
import pandas as pd

from gritic import sampletools


def make_repeated_segment_table():
    count_phase_rows = [
        (17, 3, np.nan),
        (17, 3, np.nan),
        (17, 3, 'major'),
        (10, 10, 'major'),
        (10, 10, 'minor'),
        (10, 10, 'minor'),
        (8, 12, np.nan),
        (8, 12, np.nan),
        (8, 12, np.nan),
    ]
    return pd.DataFrame({
        'Mutation_ID': [f'mutation-{index}' for index in range(9)],
        'Segment_ID': ['segment-a'] * 9,
        'Chromosome': ['1'] * 9,
        'Segment_Start': [100] * 9,
        'Segment_End': [300] * 9,
        'Major_CN': [2] * 9,
        'Minor_CN': [1] * 9,
        'Tumor_Ref_Count': [row[0] for row in count_phase_rows],
        'Tumor_Alt_Count': [row[1] for row in count_phase_rows],
        'Phasing': [row[2] for row in count_phase_rows],
    })


def make_two_segment_frequency_table():
    rows = []
    mutation_index = 0
    segment_specs = (
        ('segment-a', 99, 1),
        ('segment-b', 1, 99),
    )
    for segment_id, low_vaf_count, high_vaf_count in segment_specs:
        for ref_count, alt_count, n_mutations in (
            (90, 10, low_vaf_count),
            (5, 15, high_vaf_count),
        ):
            for _ in range(n_mutations):
                rows.append({
                    'Mutation_ID': f'mutation-{mutation_index}',
                    'Segment_ID': segment_id,
                    'Chromosome': '1',
                    'Segment_Start': 100 if segment_id == 'segment-a' else 300,
                    'Segment_End': 300 if segment_id == 'segment-a' else 500,
                    'Major_CN': 2,
                    'Minor_CN': 1,
                    'Tumor_Ref_Count': ref_count,
                    'Tumor_Alt_Count': alt_count,
                    'Phasing': np.nan,
                })
                mutation_index += 1
    return pd.DataFrame(rows)


def build_group_bundle(mutation_table, sex=None):
    (
        grouped_mutations,
        count_groups,
        phase_groups,
        segment_groups,
    ) = sampletools.build_sample_group_tables(mutation_table)
    likelihood_contexts, segment_contexts = (
        sampletools.build_likelihood_context_tables(
            grouped_mutations,
            sex,
        )
    )
    count_likelihoods = sampletools.build_count_group_likelihood_table(
        count_groups,
        phase_groups,
        segment_groups,
        likelihood_contexts,
        segment_contexts,
        purity=0.8,
        subclone_table=None,
    )
    return grouped_mutations, {
        'count_group_table': count_groups,
        'phase_group_table': phase_groups,
        'segment_group_table': segment_groups,
        'likelihood_context_table': likelihood_contexts,
        'segment_context_table': segment_contexts,
        'count_group_likelihood_table': count_likelihoods,
    }


class FrequencyWeightedQuantileTest(unittest.TestCase):
    def test_matches_expanded_numpy_linear_quantile_bit_for_bit(self):
        rng = np.random.default_rng(419_021)
        quantiles = [0.0, 0.1, 0.5, 0.9, 1.0]
        for size in range(1, 24):
            values = rng.random(size)
            frequencies = rng.integers(
                1,
                25,
                size=size,
                dtype=np.int64,
            )
            expanded = np.repeat(values, frequencies)
            for quantile in quantiles:
                with self.subTest(size=size, quantile=quantile):
                    expected = np.quantile(expanded, quantile)
                    actual = sampletools.frequency_weighted_quantile(
                        values,
                        frequencies,
                        quantile,
                    )
                    self.assertEqual(actual, expected)


class SegmentCountGroupingTest(unittest.TestCase):
    def make_segment(self, mutation_table):
        return sampletools.Segment(
            mutation_table,
            None,
            purity=0.8,
            sex=None,
            apply_reads_correction=True,
            _validated=True,
        )

    def test_group_tables_and_mutation_map_are_exact_and_dense(self):
        segment = self.make_segment(make_repeated_segment_table())

        self.assertEqual(len(segment.count_group_table), 3)
        self.assertEqual(len(segment.phase_group_table), 5)
        np.testing.assert_array_equal(
            segment.count_group_table['Count_Group_ID'],
            np.arange(3),
        )
        np.testing.assert_array_equal(
            segment.phase_group_table['Phase_Group_ID'],
            np.arange(5),
        )
        self.assertEqual(
            int(segment.count_group_table['N_Mutations'].sum()),
            len(segment.mutation_table),
        )
        self.assertEqual(
            int(segment.phase_group_table['N_Mutations'].sum()),
            len(segment.mutation_table),
        )

        mapped = segment.mutation_table.merge(
            segment.phase_group_table,
            on=['Segment_ID', 'Phase_Group_ID'],
            how='left',
            validate='many_to_one',
            suffixes=('_Mutation', '_Group'),
        ).merge(
            segment.count_group_table,
            on=['Segment_ID', 'Count_Group_ID'],
            how='left',
            validate='many_to_one',
            suffixes=('_Mutation', '_Group'),
        )
        mapped_phasing = mapped['Phasing_Mutation'].astype(
            'string'
        ).fillna('non_phased')
        np.testing.assert_array_equal(mapped_phasing, mapped['Phasing_Group'])
        np.testing.assert_array_equal(
            mapped['Tumor_Ref_Count_Mutation'],
            mapped['Tumor_Ref_Count_Group'],
        )
        np.testing.assert_array_equal(
            mapped['Tumor_Alt_Count_Mutation'],
            mapped['Tumor_Alt_Count_Group'],
        )

        observed_phase_counts = segment.mutation_table.groupby(
            'Phase_Group_ID'
        ).size()
        expected_phase_counts = segment.phase_group_table.set_index(
            'Phase_Group_ID'
        )['N_Mutations']
        pd.testing.assert_series_equal(
            observed_phase_counts.sort_index(),
            expected_phase_counts.sort_index(),
            check_names=False,
        )

    def test_group_ids_and_likelihoods_do_not_depend_on_input_order(self):
        mutation_table = make_repeated_segment_table()
        first = self.make_segment(mutation_table)
        second = self.make_segment(
            mutation_table.sample(frac=1.0, random_state=932).reset_index(
                drop=True
            )
        )

        pd.testing.assert_frame_equal(
            first.count_group_table,
            second.count_group_table,
        )
        pd.testing.assert_frame_equal(
            first.phase_group_table,
            second.phase_group_table,
        )
        first_map = first.mutation_table.set_index('Mutation_ID')[
            'Phase_Group_ID'
        ].sort_index()
        second_map = second.mutation_table.set_index('Mutation_ID')[
            'Phase_Group_ID'
        ].sort_index()
        pd.testing.assert_series_equal(first_map, second_map)

    def test_weighted_coverage_summary_matches_expanded_rows(self):
        mutation_table = make_repeated_segment_table()
        segment = self.make_segment(mutation_table)
        vaf = mutation_table['Tumor_Alt_Count'] / (
            mutation_table['Tumor_Alt_Count']
            + mutation_table['Tumor_Ref_Count']
        )
        high_rows = vaf > np.quantile(vaf, 0.9) - 0.01
        expected_coverage = (
            mutation_table.loc[high_rows, 'Tumor_Alt_Count']
            + mutation_table.loc[high_rows, 'Tumor_Ref_Count']
        ).mean()

        self.assertEqual(
            segment.highest_vaf_average_coverage,
            expected_coverage,
        )

    def test_group_weighted_likelihood_matches_expanded_group_rows(self):
        segment = self.make_segment(make_repeated_segment_table())
        store = segment.multiplicity_probabilities
        state = np.array([0.37, 0.63])
        expanded_probabilities = np.repeat(
            store.combined_array,
            store.combined_weights,
            axis=0,
        )
        expected = np.log(
            expanded_probabilities @ state + 2.2e-300
        ).sum()

        self.assertAlmostEqual(store.evaluate_likelihood(state), expected)


class SampleWideCountGroupingTest(unittest.TestCase):
    def test_shared_context_and_count_has_one_likelihood_row(self):
        mutation_table = make_two_segment_frequency_table()
        grouped_mutations, tables = build_group_bundle(mutation_table)

        self.assertEqual(len(tables['count_group_table']), 2)
        self.assertEqual(len(tables['likelihood_context_table']), 1)
        self.assertEqual(len(tables['count_group_likelihood_table']), 2)
        self.assertEqual(len(tables['segment_group_table']), 4)

        segments = []
        for _, segment_rows in grouped_mutations.groupby('Segment_ID'):
            segments.append(sampletools.Segment(
                segment_rows,
                None,
                purity=0.8,
                sex=None,
                apply_reads_correction=True,
                _validated=True,
                _sample_group_tables=tables,
            ))

        first, second = segments
        self.assertNotEqual(
            first.highest_vaf_average_coverage,
            second.highest_vaf_average_coverage,
        )
        self.assertFalse(np.allclose(
            first.alt_count_correction_factors,
            second.alt_count_correction_factors,
        ))
        probability_columns = ['Prob_Mult_1', 'Prob_Mult_2']
        pd.testing.assert_frame_equal(
            first.count_group_table[
                ['Count_Group_ID'] + probability_columns
            ].reset_index(drop=True),
            second.count_group_table[
                ['Count_Group_ID'] + probability_columns
            ].reset_index(drop=True),
        )

    def test_sex_chromosome_normal_copy_number_controls_context_reuse(self):
        mutation_table = pd.DataFrame({
            'Segment_ID': ['autosome', 'x-segment'],
            'Chromosome': ['1', 'X'],
            'Major_CN': [2, 2],
            'Minor_CN': [1, 1],
        })

        female_contexts, female_segments = (
            sampletools.build_likelihood_context_tables(
                mutation_table,
                sex='XX',
            )
        )
        self.assertEqual(len(female_contexts), 1)
        self.assertEqual(
            female_segments['Likelihood_Context_ID'].nunique(),
            1,
        )
        self.assertEqual(female_contexts.loc[0, 'Normal_Total_CN'], 2)

        male_contexts, male_segments = (
            sampletools.build_likelihood_context_tables(
                mutation_table,
                sex='XY',
            )
        )
        self.assertEqual(len(male_contexts), 2)
        self.assertEqual(
            male_segments['Likelihood_Context_ID'].nunique(),
            2,
        )
        self.assertEqual(
            set(male_contexts['Normal_Total_CN']),
            {1, 2},
        )

    def test_global_ids_and_context_ids_are_input_order_invariant(self):
        mutation_table = make_two_segment_frequency_table()
        first_mutations, first = build_group_bundle(mutation_table)
        second_mutations, second = build_group_bundle(
            mutation_table.sample(frac=1.0, random_state=175).reset_index(
                drop=True
            )
        )

        for table_name in (
            'count_group_table',
            'phase_group_table',
            'segment_group_table',
            'likelihood_context_table',
            'segment_context_table',
            'count_group_likelihood_table',
        ):
            with self.subTest(table=table_name):
                pd.testing.assert_frame_equal(
                    first[table_name],
                    second[table_name],
                )

        first_map = first_mutations.set_index('Mutation_ID')[
            'Phase_Group_ID'
        ].sort_index()
        second_map = second_mutations.set_index('Mutation_ID')[
            'Phase_Group_ID'
        ].sort_index()
        pd.testing.assert_series_equal(first_map, second_map)


if __name__ == '__main__':
    unittest.main()
