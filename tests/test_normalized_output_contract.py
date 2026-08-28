import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

from gritic import gritictimer, tableschemas


class NormalizedOutputContractTest(unittest.TestCase):
    @staticmethod
    def make_sample():
        tables = {
            'count_group': pd.DataFrame({
                'Count_Group_ID': [0, 1],
                'Tumor_Ref_Count': [10, 20],
                'Tumor_Alt_Count': [5, 5],
            }),
            'phase_group': pd.DataFrame({
                'Phase_Group_ID': [0, 1, 2],
                'Count_Group_ID': [0, 0, 1],
                'Phasing': ['non_phased', 'major', 'non_phased'],
            }),
            'likelihood_context': pd.DataFrame({
                'Likelihood_Context_ID': [0, 1],
                'Major_CN': [1, 2],
                'Minor_CN': [1, 0],
                'Normal_Total_CN': [2, 2],
            }),
            'segment_context': pd.DataFrame({
                'Segment_ID': ['A', 'B'],
                'Likelihood_Context_ID': [0, 1],
            }),
            'count_group_likelihood': pd.DataFrame({
                'Likelihood_Context_ID': [0, 1],
                'Count_Group_ID': [0, 1],
                'Prob_Mult_1': [1.0, 0.4],
                'Prob_Mult_2': [np.nan, 0.6],
            }),
            'segment_group': pd.DataFrame({
                'Segment_ID': ['A', 'A', 'B'],
                'Phase_Group_ID': [0, 1, 2],
                'N_Mutations': [1, 1, 1],
            }),
        }
        mutation_table = pd.DataFrame({
            'Segment_ID': ['A', 'A', 'B'],
            'Segment_Mutation_Index': [0, 1, 0],
            'Phase_Group_ID': [0, 1, 2],
            'Chromosome': ['1', '1', '2'],
            'Major_CN': [1, 1, 2],
            'Minor_CN': [1, 1, 0],
            'Total_CN': [2, 2, 2],
            'Gain_Type': ['1_1', '1_1', '2_0'],
            'Tumor_Ref_Count': [10, 10, 20],
            'Tumor_Alt_Count': [5, 5, 5],
            'Phasing': [pd.NA, 'major', pd.NA],
            'Context': ['ACA', 'ACC', 'ACG'],
        })

        def getter(name):
            return lambda: tables[name].copy()

        sample = SimpleNamespace(
            sample_id='NORMALIZED',
            sex=None,
            get_count_group_table=getter('count_group'),
            get_phase_group_table=getter('phase_group'),
            get_likelihood_context_table=getter('likelihood_context'),
            get_segment_context_table=getter('segment_context'),
            get_count_group_likelihood_table=getter(
                'count_group_likelihood'
            ),
            get_segment_group_table=getter('segment_group'),
            get_subclone_table=lambda: None,
        )
        return sample, mutation_table, tables

    def test_valid_contract_is_canonical_and_removes_mutation_redundancy(self):
        sample, mutation_table, _ = self.make_sample()
        outputs = gritictimer._prepare_group_output_tables(
            sample,
            mutation_table,
        )
        (
            output_mutations,
            count_groups,
            phase_groups,
            contexts,
            segment_contexts,
            likelihoods,
            segment_groups,
        ) = outputs

        self.assertTrue({
            'Tumor_Ref_Count',
            'Tumor_Alt_Count',
            'Phasing',
            'Major_CN',
            'Minor_CN',
            'Total_CN',
            'Gain_Type',
        }.isdisjoint(output_mutations.columns))
        self.assertIn('Context', output_mutations.columns)
        self.assertEqual(
            count_groups.columns.tolist(),
            tableschemas.COUNT_GROUP_TABLE_COLUMNS,
        )
        self.assertEqual(
            phase_groups.columns.tolist(),
            tableschemas.PHASE_GROUP_TABLE_COLUMNS,
        )
        self.assertEqual(
            contexts.columns.tolist(),
            tableschemas.LIKELIHOOD_CONTEXT_TABLE_COLUMNS,
        )
        self.assertEqual(
            segment_contexts.columns.tolist(),
            tableschemas.SEGMENT_CONTEXT_TABLE_COLUMNS,
        )
        self.assertEqual(
            likelihoods.columns.tolist(),
            tableschemas.COUNT_GROUP_LIKELIHOOD_BASE_COLUMNS
            + ['Prob_Mult_1', 'Prob_Mult_2'],
        )
        self.assertEqual(
            segment_groups.columns.tolist(),
            tableschemas.SEGMENT_GROUP_TABLE_COLUMNS,
        )
        for table in outputs[1:]:
            self.assertEqual(set(table['Sample_ID']), {'NORMALIZED'})

    def test_noncanonical_sample_wide_phase_ids_are_rejected(self):
        sample, mutation_table, tables = self.make_sample()
        tables['phase_group']['Phase_Group_ID'] = [1, 0, 2]

        with self.assertRaisesRegex(ValueError, 'canonical phasing order'):
            gritictimer._prepare_group_output_tables(sample, mutation_table)

    def test_segment_group_fanout_must_match_mutations(self):
        sample, mutation_table, tables = self.make_sample()
        tables['segment_group'].loc[0, 'N_Mutations'] = 2

        with self.assertRaisesRegex(ValueError, 'Mutation membership'):
            gritictimer._prepare_group_output_tables(sample, mutation_table)

    def test_segment_copy_number_must_match_its_context(self):
        sample, mutation_table, tables = self.make_sample()
        tables['likelihood_context'].loc[0, 'Minor_CN'] = 0

        with self.assertRaisesRegex(ValueError, 'copy number must match'):
            gritictimer._prepare_group_output_tables(sample, mutation_table)

    def test_sparse_likelihood_keys_must_exactly_cover_segment_usage(self):
        sample, mutation_table, tables = self.make_sample()
        tables['count_group_likelihood'].drop(index=1, inplace=True)

        with self.assertRaisesRegex(ValueError, 'exactly cover observed'):
            gritictimer._prepare_group_output_tables(sample, mutation_table)

    def test_multiplicities_above_context_major_cn_must_be_blank(self):
        sample, mutation_table, tables = self.make_sample()
        tables['count_group_likelihood'].loc[0, 'Prob_Mult_2'] = 0.1
        tables['count_group_likelihood'].loc[0, 'Prob_Mult_1'] = 0.9

        with self.assertRaisesRegex(ValueError, 'above a context Major_CN'):
            gritictimer._prepare_group_output_tables(sample, mutation_table)


if __name__ == '__main__':
    unittest.main()
