import unittest

import pandas as pd

from gritic import sampletools


class MinorPhasingCopyNumberTest(unittest.TestCase):
    @staticmethod
    def copy_number_table(minor_cn=0, include_segment_id=True):
        table = pd.DataFrame({
            'Chromosome': ['1'],
            'Segment_Start': [0],
            'Segment_End': [100],
            'Major_CN': [2],
            'Minor_CN': [minor_cn],
        })
        if include_segment_id:
            table['Segment_ID'] = 'source-segment'
        return table

    @staticmethod
    def mutation_table(phasing, *, supplied_segment_id=True, alt_count=5):
        table = pd.DataFrame({
            'Chromosome': ['1'],
            'Tumor_Ref_Count': [12],
            'Tumor_Alt_Count': [alt_count],
            'Mutation_ID': ['mutation-1'],
            'Phasing': [phasing],
        })
        if supplied_segment_id:
            table['Segment_ID'] = 'source-segment'
        else:
            table['Position'] = 50
        return table

    def make_sample(self, mutation_table, copy_number_table, **kwargs):
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            'TEST',
            0.8,
            sex='XX',
            **kwargs,
        )

    def test_rejects_minor_phasing_after_supplied_id_assignment(self):
        with self.assertRaisesRegex(
            ValueError,
            'Minor-phased mutations require Minor_CN > 0',
        ) as error:
            self.make_sample(
                self.mutation_table(' MINOR '),
                self.copy_number_table(),
            )

        message = str(error.exception)
        self.assertIn("'source-segment:mutation-1'", message)
        self.assertIn("'1-0-100'", message)

    def test_rejects_minor_phasing_after_position_assignment(self):
        with self.assertRaises(ValueError) as error:
            self.make_sample(
                self.mutation_table('minor', supplied_segment_id=False),
                self.copy_number_table(include_segment_id=False),
            )

        message = str(error.exception)
        self.assertIn('Minor_CN=0', message)
        self.assertIn("'1-0-100:mutation-1'", message)
        self.assertIn("'1-0-100'", message)

    def test_rejects_minor_phasing_after_segment_merging(self):
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 100],
            'Segment_End': [100, 200],
            'Major_CN': [2, 2],
            'Minor_CN': [0, 0],
            'Segment_ID': ['left', 'right'],
        })
        mutation_table = self.mutation_table('minor')
        mutation_table['Segment_ID'] = 'right'

        with self.assertRaises(ValueError) as error:
            self.make_sample(
                mutation_table,
                copy_number_table,
            )

        message = str(error.exception)
        self.assertIn("'right:mutation-1'", message)
        self.assertIn("'1-0-200'", message)

    def test_validation_precedes_mutation_count_filtering(self):
        with self.assertRaisesRegex(
            ValueError,
            'Minor-phased mutations require Minor_CN > 0',
        ):
            self.make_sample(
                self.mutation_table('minor', alt_count=0),
                self.copy_number_table(),
            )

    def test_absent_unphased_and_major_phasing_are_allowed(self):
        mutation_table = pd.concat(
            [
                self.mutation_table(None),
                self.mutation_table(' MAJOR '),
            ],
            ignore_index=True,
        )
        mutation_table.loc[1, 'Mutation_ID'] = 'mutation-2'

        sample = self.make_sample(
            mutation_table,
            self.copy_number_table(),
        )

        self.assertTrue(pd.isna(sample.mutation_table['Phasing'].iloc[0]))
        self.assertEqual(sample.mutation_table['Phasing'].iloc[1], 'major')

        absent_phasing_table = self.mutation_table(None).drop(
            columns='Phasing'
        )
        absent_phasing_sample = self.make_sample(
            absent_phasing_table,
            self.copy_number_table(),
        )
        self.assertTrue(
            absent_phasing_sample.mutation_table['Phasing'].isna().all()
        )

    def test_original_balanced_wgd_copy_number_remains_valid(self):
        sample = self.make_sample(
            self.mutation_table('minor'),
            self.copy_number_table(minor_cn=2),
        )

        mutation = sample.mutation_table.iloc[0]
        self.assertEqual(mutation['Phasing'], 'minor')
        self.assertEqual(mutation['Minor_CN'], 2)
        self.assertEqual(sample.segments[0].minor_cn, 2)


if __name__ == '__main__':
    unittest.main()
