import unittest
from fractions import Fraction

import numpy as np

from gritic import validation


class ScalarParameterValidationTest(unittest.TestCase):
    def test_integer_bounds_and_optional_values(self):
        cases = (
            ({'allow_none': True}, (None, 0, 7), (-1,)),
            ({'minimum': 1}, (1, np.int64(2)), (None, 0, -1)),
            (
                {'maximum': 1, 'allow_none': True},
                (None, 0, np.int64(1)),
                (-1, 2),
            ),
            (
                {'maximum': 2**32 - 1, 'allow_none': True},
                (None, 0, np.uint32(2**32 - 1)),
                (-1, 2**32),
            ),
        )
        for bounds, valid, invalid in cases:
            for value in valid:
                with self.subTest(bounds=bounds, valid=value):
                    result = validation.validate_integer(value, 'count', **bounds)
                    if value is None:
                        self.assertIsNone(result)
                    else:
                        self.assertIs(type(result), int)
                        self.assertEqual(result, value)
            for value in (*invalid, True, np.bool_(False), 1.0, '1', np.nan):
                with self.subTest(bounds=bounds, invalid=value):
                    with self.assertRaisesRegex(ValueError, 'count must be'):
                        validation.validate_integer(value, 'count', **bounds)

    def test_booleans_are_normalized_without_coercing_other_types(self):
        for value in (True, False, np.bool_(True), np.bool_(False)):
            with self.subTest(valid=value):
                self.assertIs(validation.validate_boolean(value, 'flag'), bool(value))
        for value in (0, 1, np.int64(1), 'True', None):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'flag must be a boolean'):
                    validation.validate_boolean(value, 'flag')

    def test_non_negative_integer_validation(self):
        for value in (0, 7, np.int64(3)):
            with self.subTest(valid=value):
                validated = validation.validate_integer(
                    value,
                    'count',
                )
                self.assertIs(type(validated), int)
                self.assertEqual(validated, int(value))

        for value in (-1, 1.0, True, np.bool_(False), '3', None):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(
                    ValueError,
                    'count must be a non-negative integer',
                ):
                    validation.validate_integer(value, 'count')

    def test_closed_unit_interval_validation(self):
        for value in (0, 1, 0.125, np.float64(0.75), Fraction(3, 4)):
            with self.subTest(valid=value):
                self.assertEqual(
                    validation.validate_proportion(
                        value,
                        'fraction',
                    ),
                    float(value),
                )

        for value in (
            -0.01,
            1.01,
            True,
            np.bool_(False),
            np.nan,
            np.inf,
            -np.inf,
            '0.5',
            None,
        ):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                    validation.validate_proportion(
                        value,
                        'fraction',
                    )

    def test_proportion_can_exclude_zero(self):
        for value in (1e-12, 0.5, 1, Fraction(3, 4)):
            with self.subTest(valid=value):
                self.assertEqual(
                    validation.validate_proportion(
                        value, 'min_subclone_ccf', allow_zero=False,
                    ),
                    float(value),
                )

        for value in (0, -0.1, 1.1, True, np.nan, np.inf, '0.1'):
            with self.subTest(invalid=value):
                with self.assertRaisesRegex(ValueError, 'greater than 0'):
                    validation.validate_proportion(
                        value, 'min_subclone_ccf', allow_zero=False,
                    )
