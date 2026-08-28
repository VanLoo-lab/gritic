import unittest
import warnings

import pandas as pd

from gritic import dataloader


class AdditionalDataloaderEdgeTest(unittest.TestCase):
    @staticmethod
    def subclone_table(**overrides):
        values = {
            'Cluster': ['one', 'two'],
            'Subclone_CCF': [0.8, 0.4],
            'Subclone_Fraction': [0.3, 0.2],
        }
        table = pd.DataFrame(values)
        for column, column_values in overrides.items():
            table[column] = pd.Series(column_values, dtype=object)
        return table

    @staticmethod
    def xy_copy_number_table():
        return pd.DataFrame({
            'Chromosome': ['1', 'Y'],
            'Segment_Start': [0, 0],
            'Segment_End': [100, 100],
            'Major_CN': [3, 1],
            'Minor_CN': [1, 0],
        })

    def test_negative_segment_end_reports_its_specific_constraint(self):
        with self.assertRaisesRegex(
            ValueError,
            'Segment_End must be non-negative',
        ) as error:
            dataloader.validate_segment_coordinates(pd.DataFrame({
                'Segment_Start': [0],
                'Segment_End': [-1],
            }))

        self.assertIn(
            'Segment_End must be greater than Segment_Start',
            str(error.exception),
        )

    def test_subclone_float_conversion_exceptions_become_validation_errors(self):
        cases = (
            ('Subclone_CCF', ['not-a-number', 0.4]),
            ('Subclone_Fraction', [0.3, object()]),
        )
        for column, values in cases:
            with self.subTest(column=column, value_type=type(values[-1])):
                with self.assertRaisesRegex(
                    ValueError,
                    f'{column} must contain finite values between 0 and 1',
                ):
                    dataloader.validate_subclone_values(
                        self.subclone_table(**{column: values})
                    )

    def test_ploidy_wrappers_infer_xy_karyotype_when_sex_is_omitted(self):
        copy_number_table = self.xy_copy_number_table()

        tumor_ploidy = dataloader.calculate_ploidy(
            copy_number_table,
            sex=None,
        )
        normal_ploidy = dataloader.calculate_normal_ploidy(
            copy_number_table,
            sex=None,
        )

        self.assertEqual(tumor_ploidy, 2.5)
        self.assertEqual(normal_ploidy, 1.5)

    def test_ploidy_wrappers_reject_when_chromosome_filtering_removes_every_row(self):
        copy_number_table = pd.DataFrame({
            'Chromosome': ['99'],
            'Segment_Start': [0],
            'Segment_End': [100],
            'Major_CN': [2],
            'Minor_CN': [1],
        })
        for function in (
            dataloader.calculate_ploidy,
            dataloader.calculate_normal_ploidy,
        ):
            with self.subTest(function=function.__name__):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    with self.assertRaisesRegex(
                        ValueError,
                        'No copy-number segments remain after chromosome filtering',
                    ):
                        function(
                            copy_number_table,
                            sex='XX',
                            autosome_count=1,
                            drop_unmatched_chromosomes=True,
                        )
                self.assertEqual(len(caught), 1)
                self.assertIn('Dropping 1 copy number table row', str(caught[0].message))


if __name__ == '__main__':
    unittest.main()
