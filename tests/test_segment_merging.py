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

    def test_default_merges_consecutive_equal_cn_segments_across_gaps(self):
        sample = self.make_sample()

        self.assertEqual(
            sample.copy_number_table['Segment_ID'].tolist(),
            ['1-0-350'],
        )
        self.assertEqual(
            sample.supplied_segment_id_map,
            {
                'left': '1-0-350',
                'right': '1-0-350',
                'gapped': '1-0-350',
            },
        )

    def test_max_merge_gap_is_inclusive_and_zero_requires_touching(self):
        cases = (
            (0, ['1-0-200', '1-250-350']),
            (49, ['1-0-200', '1-250-350']),
            (50, ['1-0-350']),
        )
        for max_merge_gap, expected_segment_ids in cases:
            with self.subTest(max_merge_gap=max_merge_gap):
                sample = self.make_sample(max_merge_gap=max_merge_gap)
                self.assertEqual(
                    sample.copy_number_table['Segment_ID'].tolist(),
                    expected_segment_ids,
                )

    def test_max_merge_gap_rejects_invalid_programmatic_values(self):
        for max_merge_gap in (-1, 1.0, True, '50'):
            with self.subTest(max_merge_gap=max_merge_gap):
                with self.assertRaisesRegex(
                    ValueError,
                    'max_merge_gap must be None or a non-negative integer',
                ):
                    self.make_sample(max_merge_gap=max_merge_gap)

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
