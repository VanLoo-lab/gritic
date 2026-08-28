import unittest

import pandas as pd

from gritic import sampletools


class SegmentMergingTest(unittest.TestCase):
    @staticmethod
    def copy_number_table():
        return pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Segment_Start': [0, 100, 250],
            'Segment_End': [100, 200, 350],
            'Major_CN': [2, 2, 2],
            'Minor_CN': [1, 1, 1],
            'Segment_ID': ['left', 'right', 'gapped'],
        })

    @staticmethod
    def mutation_table():
        return pd.DataFrame({
            'Chromosome': ['1', '1', '1'],
            'Tumor_Ref_Count': [12, 11, 10],
            'Tumor_Alt_Count': [5, 6, 7],
            'Mutation_ID': ['left-snv', 'right-snv', 'gapped-snv'],
            'Segment_ID': ['left', 'right', 'gapped'],
        })

    def make_sample(self, **kwargs):
        return sampletools.Sample(
            self.mutation_table(),
            self.copy_number_table(),
            None,
            'TEST',
            0.8,
            sex='XX',
            **kwargs,
        )

    def test_default_merges_only_exactly_adjacent_equal_cn_segments(self):
        sample = self.make_sample()

        self.assertEqual(
            sample.copy_number_table['Segment_ID'].tolist(),
            ['1-0-200', '1-250-350'],
        )
        self.assertEqual(
            sample.supplied_segment_id_map,
            {
                'left': '1-0-200',
                'right': '1-0-200',
                'gapped': '1-250-350',
            },
        )

    def test_programmatic_opt_out_preserves_input_segments(self):
        sample = self.make_sample(merge_cn=False)

        self.assertEqual(
            sample.copy_number_table['Segment_ID'].tolist(),
            ['1-0-100', '1-100-200', '1-250-350'],
        )
        self.assertEqual(
            sample.supplied_segment_id_map,
            {
                'left': '1-0-100',
                'right': '1-100-200',
                'gapped': '1-250-350',
            },
        )


if __name__ == '__main__':
    unittest.main()
