import unittest
import warnings

import pandas as pd

from gritic import sampletools


class ZeroCopySegmentWarningTest(unittest.TestCase):
    @staticmethod
    def copy_number_table():
        return pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 10],
            'Segment_End': [10, 20],
            'Major_CN': [0, 1],
            'Minor_CN': [0, 0],
            'Segment_ID': ['zero', 'one'],
        })

    @staticmethod
    def mutation_table(include_positive_segment=True):
        rows = [
            {
                'Chromosome': '1',
                'Tumor_Ref_Count': 10,
                'Tumor_Alt_Count': 5,
                'Mutation_ID': 'z1',
                'Segment_ID': 'zero',
            },
            {
                'Chromosome': '1',
                'Tumor_Ref_Count': 11,
                'Tumor_Alt_Count': 4,
                'Mutation_ID': 'z2',
                'Segment_ID': 'zero',
            },
        ]
        if include_positive_segment:
            rows.append({
                'Chromosome': '1',
                'Tumor_Ref_Count': 10,
                'Tumor_Alt_Count': 5,
                'Mutation_ID': 'p1',
                'Segment_ID': 'one',
            })
        return pd.DataFrame(rows)

    def make_sample(self, mutation_table):
        return sampletools.Sample(
            mutation_table,
            self.copy_number_table(),
            None,
            'TEST',
            0.8,
            sex='XX',
        )

    @staticmethod
    def user_warning_messages(caught_warnings):
        return [
            str(caught.message)
            for caught in caught_warnings
            if issubclass(caught.category, UserWarning)
        ]

    def test_warns_before_dropping_mixed_zero_copy_mutations(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sample = self.make_sample(self.mutation_table())

        warning_messages = self.user_warning_messages(caught)
        self.assertEqual(len(warning_messages), 1)
        self.assertIn('Dropping 2 SNVs', warning_messages[0])
        self.assertIn('1 zero-copy', warning_messages[0])
        self.assertIn('Source Segment_ID value(s): zero', warning_messages[0])
        self.assertEqual(sample.mutation_table['Mutation_ID'].tolist(), ['p1'])
        self.assertEqual(
            sample.mutation_table['Source_Segment_ID'].tolist(),
            ['one'],
        )

    def test_all_zero_copy_mutations_warn_then_raise_contextual_error(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            with self.assertRaisesRegex(
                ValueError,
                'No mutations remain after dropping SNVs assigned to '
                r'zero-copy \(0\+0\) segments',
            ):
                self.make_sample(
                    self.mutation_table(include_positive_segment=False)
                )

        warning_messages = self.user_warning_messages(caught)
        self.assertEqual(len(warning_messages), 1)
        self.assertIn('Dropping 2 SNVs', warning_messages[0])
        self.assertIn('1 zero-copy', warning_messages[0])
        self.assertIn('Source Segment_ID value(s): zero', warning_messages[0])


if __name__ == '__main__':
    unittest.main()
