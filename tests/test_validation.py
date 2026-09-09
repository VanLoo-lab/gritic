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


class PortableSampleIdTest(unittest.TestCase):
    def test_ordinary_ascii_and_unicode_ids_are_preserved(self):
        for sample_id in ('sample-01_A', '.hidden', '患者-α'):
            with self.subTest(sample_id=sample_id):
                self.assertEqual(
                    validation.validate_sample_id(sample_id),
                    sample_id,
                )

    def test_non_string_empty_and_dot_components_are_rejected(self):
        for sample_id in (None, 17, b'sample', '', '.', '..'):
            with self.subTest(sample_id=sample_id):
                with self.assertRaises(ValueError):
                    validation.validate_sample_id(sample_id)

    def test_every_windows_forbidden_character_is_rejected(self):
        for character in '<>:"/\\|?*':
            with self.subTest(character=character):
                with self.assertRaisesRegex(
                    ValueError,
                    'Windows-forbidden',
                ):
                    validation.validate_sample_id(f'left{character}right')

    def test_unicode_control_format_and_surrogate_characters_are_rejected(self):
        characters = ('\x00', '\n', '\u200d', chr(0xD800))
        for character in characters:
            with self.subTest(code_point=ord(character)):
                with self.assertRaisesRegex(
                    ValueError,
                    'Unicode control, format, or surrogate',
                ):
                    validation.validate_sample_id(f'a{character}b')

    def test_trailing_dot_or_space_is_rejected(self):
        for sample_id in ('sample.', 'sample ', 'sample. '):
            with self.subTest(sample_id=sample_id):
                with self.assertRaisesRegex(ValueError, 'end in a dot or space'):
                    validation.validate_sample_id(sample_id)

    def test_windows_device_names_are_rejected_case_insensitively(self):
        reserved_ids = (
            'CON',
            'con.txt',
            'PRN.results',
            'AUX',
            'NUL.data',
            'COM1',
            'lpt9.anything',
            'COM²',
            'CON .txt',
        )
        for sample_id in reserved_ids:
            with self.subTest(sample_id=sample_id):
                with self.assertRaisesRegex(ValueError, 'reserved Windows'):
                    validation.validate_sample_id(sample_id)

    def test_ascii_component_limit_accounts_for_longest_output_suffix(self):
        available = (
            validation._MAX_PATH_COMPONENT_UNITS
            - len(validation._LONGEST_SAMPLE_ID_OUTPUT_SUFFIX)
        )
        boundary_id = 'a' * available

        self.assertEqual(
            validation.validate_sample_id(boundary_id),
            boundary_id,
        )
        with self.assertRaisesRegex(ValueError, 'too long'):
            validation.validate_sample_id(boundary_id + 'a')

    def test_component_limit_is_measured_in_encoded_units_not_characters(self):
        suffix_bytes = len(
            validation._LONGEST_SAMPLE_ID_OUTPUT_SUFFIX.encode('utf-8')
        )
        max_repeated_e_acute = (
            validation._MAX_PATH_COMPONENT_UNITS - suffix_bytes
        ) // len('é'.encode('utf-8'))
        boundary_id = 'é' * max_repeated_e_acute

        self.assertEqual(
            validation.validate_sample_id(boundary_id),
            boundary_id,
        )
        with self.assertRaisesRegex(ValueError, 'too long'):
            validation.validate_sample_id(boundary_id + 'é')
