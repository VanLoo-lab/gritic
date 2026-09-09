"""Shared scalar argument checks; table and model validation stay with callers."""

from numbers import Integral, Real

import numpy as np


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


def validate_boolean(value, parameter_name):
    """Return a builtin boolean, accepting Python and NumPy booleans only."""
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{parameter_name} must be a boolean')
    return bool(value)
