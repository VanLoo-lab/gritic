import unittest

import numpy as np
import pandas as pd

from gritic import posteriortablegen
from gritic.tableschemas import GAIN_TIMING_TABLE_COLUMNS


class PosteriorNodePhasingValidationTest(unittest.TestCase):
    @staticmethod
    def route_table():
        return pd.DataFrame({
            'Sample_ID': ['sample'],
            'Segment_ID': ['segment'],
            'Route': ['route'],
        })

    @staticmethod
    def gain_timing_table(labels):
        rows = []
        for node, label in enumerate(labels):
            rows.append({
                'Sample_ID': 'sample',
                'Segment_ID': 'segment',
                'Route': 'route',
                'Node': node,
                'Node_Phasing': label,
                'Timing': 0.5,
                'Timing_CI_Low': 0.25,
                'Timing_CI_High': 0.75,
            })
        return pd.DataFrame(rows, columns=GAIN_TIMING_TABLE_COLUMNS)

    def validate(self, labels):
        posteriortablegen._validate_gain_timing_table(
            self.gain_timing_table(labels),
            self.route_table(),
            'sample',
        )

    def test_accepts_exact_major_minor_labels(self):
        self.validate(['Major', 'Minor'])

    def test_rejects_every_other_node_phasing_value(self):
        invalid_cases = (
            ('A', "'A'"),
            ('B', "'B'"),
            ('major', "'major'"),
            ('minor', "'minor'"),
            ('Other', "'Other'"),
            ('NaN', "'NaN'"),
            ('', '<blank>'),
            ('   ', '<blank>'),
            (np.nan, '<missing/NaN>'),
            (pd.NA, '<missing/NaN>'),
            (None, '<missing/NaN>'),
        )

        for value, expected_description in invalid_cases:
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                    ValueError,
                    'exact, case-sensitive values Major or Minor',
                ) as error_context:
                    self.validate([value])
                self.assertIn(
                    expected_description,
                    str(error_context.exception),
                )


if __name__ == '__main__':
    unittest.main()
