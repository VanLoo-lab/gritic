from dataclasses import dataclass, field

import numpy as np

from gritic import validation


HPD = 'hpd'
EQUAL_TAILED = 'equal-tailed'
INTERVAL_METHODS = (HPD, EQUAL_TAILED)


@dataclass(frozen=True)
class IntervalSpec:
    width: float
    method: str = HPD

    def __post_init__(self):
        object.__setattr__(
            self, 'width',
            validation.validate_proportion(
                self.width, 'interval width', allow_zero=False,
            ),
        )
        if self.method not in INTERVAL_METHODS:
            raise ValueError(
                'interval method must be one of: '
                + ', '.join(INTERVAL_METHODS)
            )


DEFAULT_POSTERIOR_SUMMARY_INTERVAL = IntervalSpec(0.95)


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
        default_factory=lambda: DEFAULT_POSTERIOR_SUMMARY_INTERVAL
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
    """Compute sample bounds with ArviZ, preserving missing and full-range results."""
    if not isinstance(interval, IntervalSpec):
        raise TypeError('interval must be an IntervalSpec')

    values = np.asarray(samples, dtype=float)
    if values.ndim != 1:
        raise ValueError('interval samples must be one-dimensional')
    if values.size == 0:
        raise ValueError('interval samples must not be empty')
    if not np.isfinite(values).all():
        return np.nan, np.nan

    # MUTIC imports the configuration but computes its continuous intervals
    # directly, so load ArviZ only when sample intervals are requested.
    from arviz_stats.base import array_stats

    # ArviZ's nearest HDI requires prob < 1. At full probability both
    # interval methods span the observed range, which ETI handles directly.
    if interval.method == EQUAL_TAILED or interval.width == 1.0:
        bounds = array_stats.eti(values, prob=interval.width, axis=0)
    else:
        bounds = array_stats.hdi(
            values, prob=interval.width, axis=0, method='nearest',
        )
    return float(bounds[0]), float(bounds[1])


#thanks to 
#https://www.geeksforgeeks.org/find-the-point-where-maximum-intervals-overlap/
#for the basic idea
def get_ids_with_maximum_overlap(segment_ci_store,segment_width_store):

    if not segment_ci_store:
        raise ValueError('At least one timing interval is required')

    segment_ids_sorted_by_ci_low = sorted(segment_ci_store.keys(),key=lambda segment_id:segment_ci_store[segment_id][0])
    segment_ids_sorted_by_ci_high= sorted(segment_ci_store.keys(),key=lambda segment_id:segment_ci_store[segment_id][1])
    
    current_segments_in = [segment_ids_sorted_by_ci_low[0]]
    segments_with_max_overlap = list(current_segments_in)
    max_segment_overlap_width = sum([segment_width_store[segment_id] for segment_id in current_segments_in])

    n_segments = len(segment_ci_store.keys())
    i = 1
    j = 0
    best_overlap_timing = segment_ci_store[segment_ids_sorted_by_ci_low[0]][0]
    while i<n_segments and j <n_segments:
        next_arrival_time = segment_ci_store[segment_ids_sorted_by_ci_low[i]][0]
        next_leaving_time = segment_ci_store[segment_ids_sorted_by_ci_high[j]][1]


        if next_arrival_time <= next_leaving_time:
            joining_segment_id = segment_ids_sorted_by_ci_low[i]

            current_segments_in.append(joining_segment_id)
            
            current_segment_overlap_width = sum([segment_width_store[segment_id] for segment_id in current_segments_in])
            if current_segment_overlap_width > max_segment_overlap_width:
                segments_with_max_overlap = list(current_segments_in)
                max_segment_overlap_width = current_segment_overlap_width
                best_overlap_timing = next_arrival_time
            i += 1
        else:
            leaving_segment_id = segment_ids_sorted_by_ci_high[j]
            current_segments_in.remove(leaving_segment_id)
            j+=1
    return segments_with_max_overlap,max_segment_overlap_width,best_overlap_timing
