from dataclasses import dataclass, field
from numbers import Real

import numpy as np


HPD = 'hpd'
EQUAL_TAILED = 'equal-tailed'
INTERVAL_METHODS = (HPD, EQUAL_TAILED)


def validate_interval_width(width):
    """Return a credible-interval width expressed as a proportion."""
    if (
        isinstance(width, (bool, np.bool_))
        or not isinstance(width, Real)
        or not np.isfinite(width)
        or not 0 < width <= 1
    ):
        raise ValueError(
            'interval width must be a finite proportion greater than 0 and at '
            'most 1'
        )
    return float(width)


@dataclass(frozen=True)
class IntervalSpec:
    width: float
    method: str = HPD

    def __post_init__(self):
        object.__setattr__(self, 'width', validate_interval_width(self.width))
        if self.method not in INTERVAL_METHODS:
            raise ValueError(
                'interval method must be one of: '
                + ', '.join(INTERVAL_METHODS)
            )


@dataclass(frozen=True)
class TimingIntervalConfig:
    """All user-configurable intervals emitted or used by GRITIC."""

    route_gain: IntervalSpec = field(
        default_factory=lambda: IntervalSpec(0.95)
    )
    tree_gain: IntervalSpec = field(
        default_factory=lambda: IntervalSpec(0.9)
    )
    wgd_overlap: IntervalSpec = field(
        default_factory=lambda: IntervalSpec(0.9)
    )
    sample_wgd: IntervalSpec = field(
        default_factory=lambda: IntervalSpec(0.9)
    )
    posterior_summary: IntervalSpec = field(
        default_factory=lambda: IntervalSpec(0.95)
    )

    def __post_init__(self):
        for name in (
            'route_gain',
            'tree_gain',
            'wgd_overlap',
            'sample_wgd',
            'posterior_summary',
        ):
            if not isinstance(getattr(self, name), IntervalSpec):
                raise TypeError(f'{name} must be an IntervalSpec')


DEFAULT_TIMING_INTERVALS = TimingIntervalConfig()


def get_interval_bounds(samples, interval):
    """Return an empirical equal-tailed or shortest contiguous HPD interval."""
    if not isinstance(interval, IntervalSpec):
        raise TypeError('interval must be an IntervalSpec')

    values = np.asarray(samples, dtype=float)
    if values.ndim != 1:
        raise ValueError('interval samples must be one-dimensional')
    if values.size == 0:
        raise ValueError('interval samples must not be empty')
    if not np.isfinite(values).all():
        return np.nan, np.nan

    if interval.method == EQUAL_TAILED:
        tail_width = (1.0 - interval.width) / 2.0
        low, high = np.quantile(
            values,
            [tail_width, 1.0 - tail_width],
        )
        return float(low), float(high)

    ordered_values = np.sort(values)
    included_count = int(np.ceil(
        interval.width * ordered_values.size
    ))
    if included_count >= ordered_values.size:
        return float(ordered_values[0]), float(ordered_values[-1])

    interval_widths = (
        ordered_values[included_count - 1:]
        - ordered_values[:ordered_values.size - included_count + 1]
    )
    start_index = int(np.argmin(interval_widths))
    end_index = start_index + included_count - 1
    return (
        float(ordered_values[start_index]),
        float(ordered_values[end_index]),
    )
