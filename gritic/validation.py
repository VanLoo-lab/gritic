"""Shared scalar argument checks; table and model validation stay with callers."""

import argparse
import unicodedata
from numbers import Integral, Real

import numpy as np


_WINDOWS_FORBIDDEN_FILENAME_CHARACTERS = frozenset('<>:"/\\|?*')
_WINDOWS_RESERVED_DEVICE_NAMES = frozenset({
    'CON',
    'PRN',
    'AUX',
    'NUL',
    'CLOCK$',
    'CONIN$',
    'CONOUT$',
    *(f'COM{number}' for number in range(1, 10)),
    *(f'LPT{number}' for number in range(1, 10)),
    'COM¹',
    'COM²',
    'COM³',
    'LPT¹',
    'LPT²',
    'LPT³',
})
_MAX_PATH_COMPONENT_UNITS = 255
_LONGEST_SAMPLE_ID_OUTPUT_SUFFIX = (
    '_posterior_timing_table_summary_penalty_False.tsv'
)


def validate_sample_id(sample_id):
    """Validate that ``sample_id`` is a portable filename component.

    GRITIC, MUTIC, and SIGTIC use the ID in artifact filenames. The
    accepted form is therefore one non-empty path component that is safe under
    common POSIX and Windows filename rules. The ID must not be ``.`` or ``..``;
    contain a path separator, a Windows-forbidden filename character, or a
    Unicode control, format, or surrogate character; end in a dot or space; or
    have a Windows device-name stem. Its longest derived GRITIC filename must
    fit both the usual 255-byte POSIX component limit and the
    255-UTF-16-code-unit Windows component limit. The value is validated, never
    normalized or sanitized.
    """
    if not isinstance(sample_id, str):
        raise ValueError('Sample_ID must be a string')
    if not sample_id:
        raise ValueError('Sample_ID must not be empty')
    if sample_id in {'.', '..'}:
        raise ValueError("Sample_ID must not be '.' or '..'")

    forbidden_characters = sorted(
        set(sample_id) & _WINDOWS_FORBIDDEN_FILENAME_CHARACTERS
    )
    if forbidden_characters:
        raise ValueError(
            'Sample_ID contains a path separator or Windows-forbidden '
            f'filename character: {forbidden_characters[0]!r}'
        )

    for character in sample_id:
        if unicodedata.category(character) in {'Cc', 'Cf', 'Cs'}:
            raise ValueError(
                'Sample_ID must not contain Unicode control, format, or '
                'surrogate characters'
            )

    if sample_id.endswith(('.', ' ')):
        raise ValueError('Sample_ID must not end in a dot or space')

    # Windows reserves device names even when followed by an extension and
    # ignores spaces immediately before that extension during device lookup.
    device_stem = sample_id.partition('.')[0].rstrip(' ').upper()
    if device_stem in _WINDOWS_RESERVED_DEVICE_NAMES:
        raise ValueError(
            f'Sample_ID uses a reserved Windows device name: {device_stem!r}'
        )

    derived_filename = sample_id + _LONGEST_SAMPLE_ID_OUTPUT_SUFFIX
    try:
        utf8_bytes = len(derived_filename.encode('utf-8'))
        utf16_code_units = len(derived_filename.encode('utf-16-le')) // 2
    except UnicodeEncodeError as error:
        raise ValueError('Sample_ID must contain valid Unicode text') from error
    if (
        utf8_bytes > _MAX_PATH_COMPONENT_UNITS
        or utf16_code_units > _MAX_PATH_COMPONENT_UNITS
    ):
        raise ValueError(
            'Sample_ID is too long: GRITIC-derived filenames must be at most '
            f'{_MAX_PATH_COMPONENT_UNITS} UTF-8 bytes and UTF-16 code units'
        )

    return sample_id


def validate_integer(
    value,
    parameter_name,
    *,
    minimum=0,
    maximum=None,
    allow_none=False,
):
    """Return a builtin integer within inclusive bounds, rejecting booleans."""
    if value is None and allow_none:
        return None
    if maximum is not None:
        requirement = f'an integer between {minimum} and {maximum}'
    elif minimum == 0:
        requirement = 'a non-negative integer'
    elif minimum == 1:
        requirement = 'a positive integer'
    else:
        requirement = f'an integer greater than or equal to {minimum}'
    if allow_none:
        requirement = f'None or {requirement}'
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Integral)
        or value < minimum
        or (maximum is not None and value > maximum)
    ):
        raise ValueError(f'{parameter_name} must be {requirement}')
    return int(value)


def validate_proportion(value, parameter_name, *, allow_zero=True):
    """Return a finite builtin float in [0, 1], or (0, 1] if zero is excluded."""
    bounds = (
        'between 0 and 1'
        if allow_zero
        else 'greater than 0 and at most 1'
    )
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Real)
        or not 0 <= value <= 1
        or (not allow_zero and value == 0)
    ):
        raise ValueError(f'{parameter_name} must be a finite number {bounds}')
    return float(value)


def positive_unit_interval_number(value):
    """Parse an argparse value as a finite number in (0, 1]."""
    try:
        return validate_proportion(float(value), 'value', allow_zero=False)
    except (TypeError, ValueError, OverflowError) as error:
        raise argparse.ArgumentTypeError(
            'must be a finite number greater than 0 and at most 1'
        ) from error


def validate_boolean(value, parameter_name):
    """Return a builtin boolean, accepting Python and NumPy booleans only."""
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{parameter_name} must be a boolean')
    return bool(value)
