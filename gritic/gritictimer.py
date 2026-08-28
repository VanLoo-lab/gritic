import os
import json
import logging
from dataclasses import dataclass
from numbers import Integral, Real

import warnings

import pandas as pd
import numpy as np
import networkx as nx


from scipy.special import logsumexp

from scipy.optimize import nnls
from scipy.linalg import null_space

from gritic.sampletools import (
    MultProbabilityStore,
    Segment,
    get_major_cn_mode,
    validate_sample_id,
)
from gritic import dataloader, timingio
from gritic.intervaltools import (
    DEFAULT_TIMING_INTERVALS,
    TimingIntervalConfig,
    get_interval_bounds,
)
from gritic.tableschemas import (
    GAIN_TIMING_TABLE_COLUMNS,
    ROUTE_PARTICLES_REPRESENTATION,
    ROUTE_TABLE_COLUMNS,
    TIMING_REPRESENTATION_COLUMN,
    UNIFORM_NO_GAIN_REPRESENTATION,
)
from gritic.timingio import write_timing_archive

from sklearn.neighbors import NearestNeighbors

import time

import gritic.distributiontools as distributiontools
import gritic.treetools as treetools
from gritic.route_tree import RouteTree
import gritic.hitandrun as hitandrun

import gritic.posteriortablegen as posteriortablegen

import pathlib


logger = logging.getLogger(__name__)


MIN_WGD_MUTATIONS = 10
ROUTE_CONDITIONAL_SAMPLE_COUNT = 1000
NO_GAIN_CLONE_SHARE_PROPOSAL_COUNT = 50_500
COMBINED_WGD_SAMPLE_COUNT = 500
_INITIAL_DENSITY_EVALUATION_BATCH_COUNT = 101
_MAX_DENSITY_RETRY_DELAY_SECONDS = 30.0
_DENSITY_RUNTIME_MULTIPLIER = 5.0
_MULT_NEGATIVE_TOLERANCE = 1e-12
# Bound simultaneous per-segment posterior output while still sharing costly
# route geometry across high-copy-number segment pairs.
SHARED_FIT_OUTPUT_MEMORY_BUDGET = 1_900_000_000
SUBCLONE_FRACTION_PRIOR_MODES = ('adjusted', 'supplied')
DEFAULT_SUBCLONE_FRACTION_PRIOR = 'adjusted'
DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR = False


_MIRROR_ALLELE = {
    'Major': 'Minor',
    'Minor': 'Major',
}


def is_uniform_no_gain_route_tree(route_tree):
    """Return whether every mutation-bearing copy spans the full [0, 1]."""
    if route_tree.major_cn != 1 or route_tree.minor_cn not in (0, 1):
        return False

    nodes = list(route_tree.non_phased_node_order)
    if not nodes or route_tree.timeable_nodes or route_tree.wgd_nodes:
        return False

    for node in nodes:
        attributes = route_tree.node_attributes[node]
        if (
            attributes['Predecessor'] is not None
            or attributes['Successors']
            or attributes['Multiplicity'] != 1
        ):
            return False
    return True


def _density_evaluation_due(
    batches_sampled,
    sample_count,
    max_samples,
    *,
    current_time=None,
    next_evaluation_time=None,
):
    """Return whether route-proposal density should be evaluated now."""
    if sample_count >= max_samples:
        return True
    if next_evaluation_time is None:
        return batches_sampled >= _INITIAL_DENSITY_EVALUATION_BATCH_COUNT
    if current_time is None:
        raise ValueError(
            'current_time is required after a density retry has been scheduled'
        )
    return current_time > next_evaluation_time


def _density_retry_delay(sampling_time, density_time):
    """Return the bounded delay before repeating a density diagnostic."""
    return min(
        max(
            density_time * _DENSITY_RUNTIME_MULTIPLIER,
            sampling_time,
        ),
        _MAX_DENSITY_RETRY_DELAY_SECONDS,
    )


def validate_subclone_fraction_prior(subclone_fraction_prior):
    if subclone_fraction_prior not in SUBCLONE_FRACTION_PRIOR_MODES:
        permitted = ', '.join(SUBCLONE_FRACTION_PRIOR_MODES)
        raise ValueError(
            f'subclone_fraction_prior must be one of: {permitted}'
        )
    return subclone_fraction_prior


def validate_unordered_balanced_route_prior(unordered_balanced_route_prior):
    if not isinstance(unordered_balanced_route_prior, (bool, np.bool_)):
        raise ValueError('unordered_balanced_route_prior must be a boolean')
    return bool(unordered_balanced_route_prior)


def get_sample_clone_fractions(subclone_table):
    if subclone_table is None:
        return np.array([1.0])

    subclone_fractions = np.asarray(
        subclone_table['Subclone_Fraction'],
        dtype=float,
    )
    return np.concatenate((
        [1.0 - np.sum(subclone_fractions)],
        subclone_fractions,
    ))


def get_clone_share_prior_alpha(
    sample_clone_fractions,
    subclonal_correction_array=None,
    subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
):
    subclone_fraction_prior = validate_subclone_fraction_prior(
        subclone_fraction_prior
    )
    clone_fractions = np.asarray(sample_clone_fractions, dtype=float)
    if clone_fractions.ndim != 1 or clone_fractions.size == 0:
        raise ValueError(
            'sample_clone_fractions must be a non-empty one-dimensional array'
        )
    if (
        not np.isfinite(clone_fractions).all()
        or (clone_fractions < 0).any()
        or not np.isclose(np.sum(clone_fractions), 1.0)
    ):
        raise ValueError(
            'sample_clone_fractions must contain finite non-negative values '
            'that sum to 1'
        )

    if subclone_fraction_prior == 'supplied':
        prior_fractions = clone_fractions
    else:
        correction = np.asarray(subclonal_correction_array, dtype=float)
        if correction.shape != clone_fractions.shape:
            raise ValueError(
                'subclonal_correction_array must have the same shape as '
                'sample_clone_fractions'
            )
        if not np.isfinite(correction).all() or (correction <= 0).any():
            raise ValueError(
                'subclonal_correction_array must contain finite positive '
                'values'
            )
        # Normalize in log space so a very small but positive detection
        # probability does not overflow the inverse correction.
        log_prior_weights = np.full(clone_fractions.shape, -np.inf)
        positive_fraction = clone_fractions > 0
        log_prior_weights[positive_fraction] = (
            np.log(clone_fractions[positive_fraction])
            - np.log(correction[positive_fraction])
        )
        prior_fractions = np.exp(
            log_prior_weights - logsumexp(log_prior_weights)
        )

    return 1.0 + prior_fractions


def _get_mirror_node_map(source_tree, target_tree):
    """Map source nodes onto the target after exchanging allele identities."""
    allele_attribute = treetools.ALLELE_ATTRIBUTE

    def node_match(source_attributes, target_attributes):
        source_allele = source_attributes.get(allele_attribute)
        return (
            _MIRROR_ALLELE.get(source_allele)
            == target_attributes.get(allele_attribute)
            and source_attributes.get('WGD_Symbol')
            == target_attributes.get('WGD_Symbol')
            and source_attributes.get('Terminal_Node')
            == target_attributes.get('Terminal_Node')
        )

    matcher = nx.algorithms.isomorphism.DiGraphMatcher(
        source_tree,
        target_tree,
        node_match=node_match,
    )
    try:
        node_map = next(matcher.isomorphisms_iter())
    except StopIteration as error:
        raise ValueError(
            'Balanced mirror routes must be structurally isomorphic after '
            'exchanging Major and Minor alleles'
        ) from error

    if set(node_map) != set(source_tree) or set(node_map.values()) != set(
        target_tree
    ):
        raise ValueError('Balanced mirror route node mapping must be bijective')
    return node_map


@dataclass(frozen=True)
class ProposalGeometry:
    """Likelihood-independent route proposals shared by matching segments."""

    mult_store: np.ndarray
    timing_store: np.ndarray
    wgd_timing_store: np.ndarray
    density: np.ndarray


class Route:
    def __init__(self, route_id, route_tree):
        self.route_id = route_id
        self.short_id = route_id[:9]
        self.route_tree = route_tree

        self.major_cn = route_tree.major_cn
        self.minor_cn = route_tree.minor_cn
        self.total_cn = self.major_cn + self.minor_cn
        self.wgd_status = route_tree.wgd_status

        self.mirror_route_id = route_tree.mirror_route_id
        self.reset_fit()

    def reset_fit(self):
        """Clear likelihood-dependent state before fitting mutation data."""
        self.log_evidence = None
        self.node_timing = None
        self.wgd_timing_store = None
        self.n_events_store = None
        self.mult_store = None
        self.clone_share_store = None
        self.n_subclones = None
        self.unphased_mirror_source = None
        self.density = None
        self.density_high = None
        self.run_time = np.nan

    def get_average_events(self, event_type):
        if self.n_events_store is None:
            return np.nan
        return np.mean(self.n_events_store[event_type])

    @property
    def is_uniform_no_gain(self):
        return is_uniform_no_gain_route_tree(self.route_tree)

    def get_node_timing(self, node):
        if self.node_timing is None:
            return np.nan
        node_index = self.route_tree.non_phased_node_order.index(node)
        return self.node_timing[node_index, :]

    def get_cumulative_timing(self, timing_periods):
        cumulative_timing = []
        node_order = self.route_tree.non_phased_node_order
        node_positions = {
            node: position for position, node in enumerate(node_order)
        }
        for node in node_order:
            predecessor = self.route_tree.node_attributes[node]['Predecessor']
            node_period = timing_periods[:, node_positions[node]]
            if predecessor is None:
                cumulative_timing.append(node_period)
            else:
                cumulative_timing.append(
                    cumulative_timing[node_positions[predecessor]] + node_period
                )
        return np.asarray(cumulative_timing)

    @staticmethod
    def get_weighted_sample_indices(
        weights,
        n_candidates,
        n_samples=ROUTE_CONDITIONAL_SAMPLE_COUNT,
    ):
        weights = np.asarray(weights, dtype=float)
        if np.isnan(weights).any():
            return None
        if weights.shape != (n_candidates,):
            raise ValueError(
                'Importance weights must contain one value per proposal'
            )
        assert (weights > -1e-80).all()
        weights = np.clip(weights, 0, 1)
        weights = weights / np.sum(weights)
        return np.random.choice(
            np.arange(n_candidates),
            size=n_samples,
            replace=True,
            p=weights,
        )

    def get_weighted_arrays(
        self,
        cumulative_timing,
        wgd_timing_store,
        mult_store,
        weights,
        n_samples=ROUTE_CONDITIONAL_SAMPLE_COUNT,
    ):
        weighted_sample = self.get_weighted_sample_indices(
            weights,
            cumulative_timing.shape[1],
            n_samples,
        )
        if weighted_sample is None:
            cumulative_timing = (
                np.ones_like(cumulative_timing)[:, :n_samples] * np.nan
            )
            wgd_timing_store = (
                np.ones_like(wgd_timing_store)[:n_samples] * np.nan
            )
            mult_store = np.ones_like(mult_store)[:n_samples, :] * np.nan
            return cumulative_timing, wgd_timing_store, mult_store
        return (
            np.asarray(cumulative_timing)[:, weighted_sample],
            wgd_timing_store[weighted_sample],
            mult_store[weighted_sample],
        )

    def sample_mults(self, wgd_timing, n_samples):
        constraints_matrix, constraints_sum = (
            self.route_tree.get_combined_constraints(wgd_timing)
        )
        start_sol = nnls(constraints_matrix, constraints_sum)[0]
        constraints_null = null_space(constraints_matrix)
        solutions = hitandrun.hit_and_run(
            constraints_null,
            start_sol,
            n_samples=n_samples,
        )

        timing = self.get_cumulative_timing(solutions)
        mult = np.matmul(solutions, self.route_tree.timing_matrix)

        unphased_sum = np.sum(mult[:, :self.major_cn], axis=1)[:, None]
        major_sum = np.sum(
            mult[:, self.major_cn:2 * self.major_cn],
            axis=1,
        )[:, None]
        minor_sum = np.sum(
            mult[:, 2 * self.major_cn:],
            axis=1,
        )[:, None]
        combined_mult_sum = np.concatenate(
            [
                np.repeat(unphased_sum, self.major_cn, axis=1),
                np.repeat(major_sum, self.major_cn, axis=1),
                np.repeat(minor_sum, self.minor_cn, axis=1),
            ],
            axis=1,
        )
        return mult / combined_mult_sum, timing

    def get_n_events_estimate(
        self,
        node_timing,
        wgd_timing,
        n_samples=300,
    ):
        random_indices = np.random.choice(
            node_timing.shape[1],
            size=n_samples,
        )
        event_arrays = self.route_tree.get_n_events_batch(
            node_timing[:, random_indices],
            wgd_timing[random_indices],
        )
        return {
            event_type: values.tolist()
            for event_type, values in zip(
                ('N_Events', 'Pre_WGD_Losses', 'Post_WGD_Losses'),
                event_arrays,
            )
        }

    @staticmethod
    def simulate_clone_share(alpha, n_samples):
        return np.random.dirichlet(alpha, size=n_samples)

    def get_density_estimate(
        self,
        samples,
        n_test_points=1000,
        radius=0.05,
    ):
        nn_finder = NearestNeighbors(radius=radius, p=1)
        nn_finder.fit(samples)
        n_test_points = min(n_test_points, samples.shape[0])
        random_mult_indices = np.random.choice(
            samples.shape[0],
            size=n_test_points,
            replace=False,
        )
        nearest_neighbors = nn_finder.radius_neighbors(
            samples[random_mult_indices, :],
            return_distance=False,
        )
        nn_size = np.asarray([neighbors.size - 1 for neighbors in nearest_neighbors])
        return np.mean(nn_size > 0.1), np.mean(nn_size > 2.1)

    def run_geometry_sampling(
        self,
        wgd_timing_distribution,
        samples_per_run=500,
        max_samples=500_000,
        density_cut_off=0.9,
    ):
        """Sample route geometry without segment-specific clone fractions.

        The density diagnostic covers only hit-and-run timing coordinates.
        Clone shares are independent Dirichlet draws, so they do not require a
        Markov-chain mixing diagnostic and can be generated per segment later.
        """
        timing_batches = []
        mult_batches = []
        wgd_timing_batches = []
        density = np.array([np.nan, np.nan])

        sampling_started = time.perf_counter()
        next_eval_time = None
        batches_sampled = 0
        sample_count = 0

        while sample_count < max_samples:
            if self.wgd_status:
                wgd_timing = np.random.choice(wgd_timing_distribution)
            else:
                wgd_timing = np.nan
            mults, timing = self.sample_mults(wgd_timing, samples_per_run)
            timing_batches.append(timing)
            mult_batches.append(mults)
            wgd_timing_batches.append(
                np.full(timing.shape[1], wgd_timing, dtype=float)
            )
            batches_sampled += 1
            sample_count += timing.shape[1]

            current_time = (
                None
                if next_eval_time is None
                else time.perf_counter()
            )
            should_evaluate = _density_evaluation_due(
                batches_sampled,
                sample_count,
                max_samples,
                current_time=current_time,
                next_evaluation_time=next_eval_time,
            )
            if should_evaluate:
                timing_test = np.concatenate(timing_batches, axis=1).T
                sampling_time = time.perf_counter() - sampling_started
                density_started = time.perf_counter()
                density = np.asarray(self.get_density_estimate(timing_test))
                density_finished = time.perf_counter()
                density_time = density_finished - density_started
                del timing_test
                next_eval_time = density_finished + _density_retry_delay(
                    sampling_time,
                    density_time,
                )
                sampling_started = density_finished
                if density[0] >= density_cut_off:
                    break

        return ProposalGeometry(
            mult_store=np.concatenate(mult_batches, axis=0),
            timing_store=np.concatenate(timing_batches, axis=1),
            wgd_timing_store=np.concatenate(wgd_timing_batches),
            density=density,
        )

    def materialize_mult_store(
        self,
        geometry,
        alpha,
        n_subclones,
        clone_share=None,
    ):
        expected_columns = 2 * self.major_cn + self.minor_cn
        base_mult_store = np.asarray(geometry.mult_store)
        if (
            base_mult_store.ndim != 2
            or base_mult_store.shape[1] != expected_columns
        ):
            raise ValueError(
                'Proposal geometry has an unexpected number of multiplicity '
                'columns'
            )
        if n_subclones == 0:
            if clone_share is not None:
                raise ValueError(
                    'clone_share must be omitted when there are no subclones'
                )
            return base_mult_store
        if alpha is None or len(alpha) != n_subclones + 1:
            raise ValueError(
                'Clone-share prior must have one clonal and one value per '
                'subclone'
            )
        if clone_share is None:
            clone_share = self.simulate_clone_share(
                alpha,
                base_mult_store.shape[0],
            )
        clone_share = np.asarray(clone_share)
        expected_shape = (base_mult_store.shape[0], n_subclones + 1)
        if clone_share.shape != expected_shape:
            raise ValueError(
                f'clone_share must have shape {expected_shape}'
            )
        return np.concatenate(
            [
                base_mult_store * clone_share[:, :1],
                clone_share[:, 1:],
            ],
            axis=1,
        )

    def _set_uniform_no_gain_event_store(self):
        node_timing = np.ones(
            (len(self.route_tree.non_phased_node_order), 1),
            dtype=float,
        )
        wgd_timing = np.array(
            [0.5 if self.wgd_status else np.nan],
            dtype=float,
        )
        event_arrays = self.route_tree.get_n_events_batch(
            node_timing,
            wgd_timing,
        )
        self.n_events_store = {
            event_type: values.tolist()
            for event_type, values in zip(
                ('N_Events', 'Pre_WGD_Losses', 'Post_WGD_Losses'),
                event_arrays,
            )
        }

    def run_uniform_no_gain_sampling(
        self,
        mult_probabilities,
        alpha,
        n_subclones,
    ):
        """Fit clone shares without constructing deterministic timing geometry."""
        if not self.is_uniform_no_gain:
            raise ValueError(
                'Uniform no-gain sampling requires only independent terminal '
                'multiplicity-one copies'
            )

        run_started = time.perf_counter()
        self.n_subclones = n_subclones
        if n_subclones == 0:
            clone_share = np.ones((1, 1), dtype=float)
        else:
            clone_share = self.simulate_clone_share(
                alpha,
                NO_GAIN_CLONE_SHARE_PROPOSAL_COUNT,
            )

        n_proposals = clone_share.shape[0]
        clonal_columns = 2 * self.major_cn + self.minor_cn
        full_states = np.concatenate(
            [
                np.repeat(
                    clone_share[:, :1],
                    clonal_columns,
                    axis=1,
                ),
                clone_share[:, 1:],
            ],
            axis=1,
        )
        ll_store = mult_probabilities.evaluate_likelihood_array(full_states)
        if ll_store.size == 0:
            raise ValueError('A route fit must contain at least one proposal')
        self.log_evidence = logsumexp(ll_store) - np.log(ll_store.size)

        if n_subclones > 0:
            weights = np.exp(ll_store - np.max(ll_store))
            weighted_sample = self.get_weighted_sample_indices(
                weights,
                n_proposals,
            )
            if weighted_sample is None:
                self.clone_share_store = np.full(
                    (ROUTE_CONDITIONAL_SAMPLE_COUNT, n_subclones + 1),
                    np.nan,
                )
            else:
                self.clone_share_store = clone_share[
                    weighted_sample, :
                ].copy()

        self._set_uniform_no_gain_event_store()
        self.density = 1.0
        self.density_high = 1.0
        self.run_time = time.perf_counter() - run_started
        return None if n_subclones == 0 else clone_share

    def _validate_mirror_route(self, mirror_route):
        if mirror_route is None:
            return False
        if self.mirror_route_id is None:
            raise ValueError('Only balanced two-component routes have mirrors')
        if mirror_route.route_id != self.mirror_route_id:
            raise ValueError(
                'The supplied mirror route does not match the expected ordered '
                'mirror route ID'
            )
        if mirror_route.mirror_route_id != self.route_id:
            raise ValueError(
                'Balanced route mirror relationships must be reciprocal'
            )
        return mirror_route is not self

    def _transform_mirror_timing_store(self, source_store, mirror_route):
        source_store = np.asarray(source_store)
        source_order = mirror_route.route_tree.non_phased_node_order
        target_order = self.route_tree.non_phased_node_order
        if source_store.ndim != 2 or source_store.shape[0] != len(source_order):
            raise ValueError(
                'Mirror timing store rows must match the source route node order'
            )

        node_map = _get_mirror_node_map(
            mirror_route.route_tree.main_tree,
            self.route_tree.main_tree,
        )
        target_positions = {
            node: position for position, node in enumerate(target_order)
        }
        transformed_store = np.empty(
            (len(target_order), source_store.shape[1]),
            dtype=source_store.dtype,
        )
        for source_position, source_node in enumerate(source_order):
            target_node = node_map[source_node]
            transformed_store[target_positions[target_node], :] = (
                source_store[source_position, :]
            )
        return transformed_store

    def _transform_mirror_mult_store(self, source_store, n_subclones):
        source_store = np.asarray(source_store)
        if self.major_cn != self.minor_cn or self.minor_cn <= 0:
            raise ValueError(
                'Mirror multiplicity transformation requires balanced positive '
                'copy number'
            )
        expected_columns = 3 * self.major_cn + n_subclones
        if source_store.ndim != 2 or source_store.shape[1] != expected_columns:
            raise ValueError(
                'Mirror multiplicity store has an unexpected number of columns'
            )

        transformed_store = source_store.copy()
        major_start = self.major_cn
        minor_start = 2 * self.major_cn
        subclone_start = 3 * self.major_cn
        transformed_store[:, major_start:minor_start] = source_store[
            :, minor_start:subclone_start
        ]
        transformed_store[:, minor_start:subclone_start] = source_store[
            :, major_start:minor_start
        ]
        return transformed_store

    def transform_mirror_geometry(self, mirror_route, geometry):
        if not self._validate_mirror_route(mirror_route):
            raise ValueError('A distinct mirror route is required')
        return ProposalGeometry(
            mult_store=self._transform_mirror_mult_store(
                geometry.mult_store,
                n_subclones=0,
            ),
            timing_store=self._transform_mirror_timing_store(
                geometry.timing_store,
                mirror_route,
            ),
            wgd_timing_store=np.asarray(
                geometry.wgd_timing_store
            ).copy(),
            density=np.asarray(geometry.density).copy(),
        )

    @staticmethod
    def _contains_phased_mutations(mult_probabilities):
        return bool(
            mult_probabilities.use_major
            or mult_probabilities.use_minor
        )

    def _reuse_fitted_unphased_mirror(self, mirror_route, n_subclones):
        """Derive an exact oriented copy of a fitted unphased mirror route."""
        if not self._validate_mirror_route(mirror_route):
            raise ValueError('A distinct fitted mirror route is required')
        required_attributes = (
            'log_evidence',
            'node_timing',
            'wgd_timing_store',
            'mult_store',
            'n_events_store',
            'density',
            'density_high',
        )
        missing_attributes = [
            attribute
            for attribute in required_attributes
            if getattr(mirror_route, attribute, None) is None
        ]
        if missing_attributes:
            raise ValueError(
                'The mirror route must be fully fitted before reuse; missing: '
                + ', '.join(missing_attributes)
            )

        self.log_evidence = mirror_route.log_evidence
        self.node_timing = self._transform_mirror_timing_store(
            mirror_route.node_timing,
            mirror_route,
        )
        self.wgd_timing_store = mirror_route.wgd_timing_store.copy()
        self.mult_store = self._transform_mirror_mult_store(
            mirror_route.mult_store,
            n_subclones,
        )
        self.n_events_store = {
            event_type: np.asarray(event_values).copy()
            for event_type, event_values in mirror_route.n_events_store.items()
        }
        self.density = mirror_route.density
        self.density_high = mirror_route.density_high
        self.unphased_mirror_source = mirror_route

    def run_sampling(
        self,
        mult_probabilities,
        subclone_table,
        wgd_timing_distribution,
        subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
        mirror_route=None,
        proposal_geometry=None,
        clone_share=None,
        shared_geometry_time=0.0,
    ):
        run_started = time.perf_counter()
        self.unphased_mirror_source = None

        subclone_fraction_prior = validate_subclone_fraction_prior(
            subclone_fraction_prior
        )
        n_subclones = (
            0 if subclone_table is None else len(subclone_table.index)
        )
        self.n_subclones = n_subclones
        if (
            not self._contains_phased_mutations(mult_probabilities)
            and mirror_route is not None
            and mirror_route is not self
            and mirror_route.log_evidence is not None
        ):
            self._reuse_fitted_unphased_mirror(mirror_route, n_subclones)
            self.run_time = (
                time.perf_counter() - run_started + shared_geometry_time
            )
            return clone_share

        alpha = None
        if n_subclones > 0:
            sample_clone_fractions = get_sample_clone_fractions(subclone_table)
            subclonal_correction_array = None
            if subclone_fraction_prior == 'adjusted':
                subclonal_correction_array = (
                    mult_probabilities.get_subclonal_correction_array(
                        subclone_table
                    )
                )
            alpha = get_clone_share_prior_alpha(
                sample_clone_fractions,
                subclonal_correction_array,
                subclone_fraction_prior,
            )

        if self.is_uniform_no_gain:
            return self.run_uniform_no_gain_sampling(
                mult_probabilities,
                alpha,
                n_subclones,
            )

        if proposal_geometry is None:
            proposal_geometry = self.run_geometry_sampling(
                wgd_timing_distribution
            )
        if n_subclones > 0 and clone_share is None:
            clone_share = self.simulate_clone_share(
                alpha,
                proposal_geometry.mult_store.shape[0],
            )
        mult_store = self.materialize_mult_store(
            proposal_geometry,
            alpha,
            n_subclones,
            clone_share=clone_share,
        )
        timing_store = proposal_geometry.timing_store
        wgd_timing_store = proposal_geometry.wgd_timing_store

        ll_store = mult_probabilities.evaluate_likelihood_array(mult_store)
        if ll_store.size == 0:
            raise ValueError('A route fit must contain at least one proposal')
        self.log_evidence = logsumexp(ll_store) - np.log(ll_store.size)

        weights = np.exp(ll_store - np.max(ll_store))
        node_timing, wgd_timing_store, mult_store = (
            self.get_weighted_arrays(
                timing_store,
                wgd_timing_store,
                mult_store,
                weights,
            )
        )
        self.n_events_store = self.get_n_events_estimate(
            node_timing,
            wgd_timing_store,
        )

        self.node_timing = node_timing
        self.wgd_timing_store = np.asarray(wgd_timing_store)
        self.mult_store = mult_store
        self.run_time = (
            time.perf_counter() - run_started + shared_geometry_time
        )
        self.density = proposal_geometry.density[0]
        self.density_high = proposal_geometry.density[1]
        return clone_share


class RouteClassifier:
    def __init__(
        self,
        major_cn,
        minor_cn,
        wgd_status,
        wgd_trees_status,
        *,
        subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
        unordered_balanced_route_prior=DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
        route_trees=None,
    ):
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.wgd_status = wgd_status
        self.subclone_fraction_prior = validate_subclone_fraction_prior(
            subclone_fraction_prior
        )
        self.unordered_balanced_route_prior = (
            validate_unordered_balanced_route_prior(
                unordered_balanced_route_prior
            )
        )
        if route_trees is None:
            route_trees = self.generate_route_trees(wgd_trees_status)
        for route_tree in route_trees.values():
            if (
                route_tree.major_cn != self.major_cn
                or route_tree.minor_cn != self.minor_cn
                or route_tree.wgd_status != self.wgd_status
            ):
                raise ValueError(
                    'Shared route trees must match the classifier copy-number '
                    'state and WGD status'
                )
        self.routes = {
            route_id: Route(route_id, route_tree)
            for route_id, route_tree in route_trees.items()
        }
        route_representations = {
            (
                UNIFORM_NO_GAIN_REPRESENTATION
                if route.is_uniform_no_gain
                else ROUTE_PARTICLES_REPRESENTATION
            )
            for route in self.routes.values()
        }
        if len(route_representations) != 1:
            raise ValueError(
                'Every route in a classifier must use the same timing '
                'representation'
            )
        self.timing_representation = route_representations.pop()

        self.route_probabilities = {}

    def get_best_timing(self):
        best_route_id = max(self.route_probabilities, key=self.route_probabilities.get)
        best_route = self.routes[best_route_id]
        best_timing = []
        for node in best_route.route_tree.timeable_nodes:
            best_timing.append(best_route.get_node_timing(node))
        return np.array(best_timing)
    
    def generate_route_trees(self, wgd_trees_status):
        possible_trees = treetools.get_nx_trees(
            self.major_cn,
            self.minor_cn,
            self.wgd_status,
            wgd_trees_status,
        )
        route_trees = {
            tree_id: RouteTree(
                tree,
                self.major_cn,
                self.minor_cn,
                self.wgd_status,
            )
            for tree_id, tree in possible_trees.items()
        }
        for route_tree in route_trees.values():
            route_tree.mirror_route_id = None
            if self.major_cn == self.minor_cn and self.minor_cn > 0:
                route_tree.mirror_route_id = treetools.get_mirror_tree_hash(
                    route_tree.main_tree
                )
        return route_trees

    def get_route_prior_weight(self, route):
        """Return the optional unordered-model prior weight for a route."""
        if not self.unordered_balanced_route_prior:
            return 1.0
        if self.major_cn != self.minor_cn or self.minor_cn <= 0:
            return 1.0
        if route.mirror_route_id == route.route_id:
            return 1.0
        return 0.5
 
    def _prepare_fit(self):
        self.route_probabilities = {}
        for route in self.routes.values():
            route.reset_fit()

    def _finalize_fit(self):
        route_ids = sorted(self.routes)
        log_scores = np.asarray([
            self.routes[route_id].log_evidence
            + np.log(self.get_route_prior_weight(self.routes[route_id]))
            for route_id in route_ids
        ])
        probabilities = np.exp(log_scores - logsumexp(log_scores))
        self.route_probabilities = dict(zip(route_ids, probabilities))

    def fit_routes(
        self,
        mult_probabilities,
        subclone_table,
        wgd_timing_distribution,
    ):
        fit_route_classifiers(
            [(self, mult_probabilities, subclone_table)],
            wgd_timing_distribution,
        )

    def _get_output_routes(self):
        if len(self.route_probabilities) == 0:
            raise ValueError("routes hasn't been fit yet")

        output_routes = []
        full_ids_by_short_id = {}
        for route_id, probability in self.route_probabilities.items():
            route = self.routes[route_id]
            output_routes.append((route, probability))
            full_ids_by_short_id.setdefault(route.short_id, []).append(route_id)

        collisions = {
            short_id: full_ids
            for short_id, full_ids in full_ids_by_short_id.items()
            if len(full_ids) > 1
        }
        if collisions:
            collision_text = '; '.join(
                f"{short_id}: {', '.join(map(str, full_ids))}"
                for short_id, full_ids in collisions.items()
            )
            raise ValueError(
                'Route identifiers must be unique after shortening to the '
                f'output Route value; collisions: {collision_text}'
            )
        return output_routes

    def get_route_table(self):
        route_table_data = {
            'Route': [],
            TIMING_REPRESENTATION_COLUMN: [],
            'Probability': [],
            'Average_N_Events': [],
            'Average_Pre_WGD_Losses': [],
            'Average_Post_WGD_Losses': [],
            'Time': [],
            'Density': [],
        }
        for route, probability in self._get_output_routes():
            route_table_data['Route'].append(route.short_id)
            route_table_data[TIMING_REPRESENTATION_COLUMN].append(
                self.timing_representation
            )
            route_table_data['Probability'].append(probability)
            route_table_data['Average_N_Events'].append(
                route.get_average_events('N_Events')
            )
            route_table_data['Average_Pre_WGD_Losses'].append(
                route.get_average_events('Pre_WGD_Losses')
            )
            route_table_data['Average_Post_WGD_Losses'].append(
                route.get_average_events('Post_WGD_Losses')
            )
            route_table_data['Time'].append(np.round(route.run_time, 3))
            route_table_data['Density'].append(route.density)

        return pd.DataFrame(route_table_data)

    def get_timing_table(
        self,
        interval=DEFAULT_TIMING_INTERVALS.route_gain,
        rounding_digits=3,
    ):
        timing_table_data = {
            'Route': [],
            'Node': [],
            'Node_Phasing': [],
            'Timing': [],
            'Timing_CI_Low': [],
            'Timing_CI_High': [],
        }
        output_routes = self._get_output_routes()
        for route, _ in output_routes:
            for node in route.route_tree.timeable_nodes:
                node_timing = np.asarray(route.get_node_timing(node))
                timing_ci_low, timing_ci_high = get_interval_bounds(
                    node_timing,
                    interval,
                )
                timing_table_data['Route'].append(route.short_id)
                timing_table_data['Node'].append(node)
                timing_table_data['Node_Phasing'].append(
                    route.route_tree.node_attributes[node]['Phasing']
                )
                timing_table_data['Timing'].append(
                    np.round(np.median(node_timing), rounding_digits)
                )
                timing_table_data['Timing_CI_Low'].append(
                    np.round(timing_ci_low, rounding_digits)
                )
                timing_table_data['Timing_CI_High'].append(
                    np.round(timing_ci_high, rounding_digits)
                )

        timing_table = pd.DataFrame(timing_table_data)
        duplicate_gains = timing_table.duplicated(
            subset=['Route', 'Node'],
            keep=False,
        )
        if duplicate_gains.any():
            duplicate_keys = (
                timing_table.loc[duplicate_gains, ['Route', 'Node']]
                .drop_duplicates()
                .astype(str)
                .agg(':'.join, axis=1)
                .tolist()
            )
            raise ValueError(
                'Each timeable gain must have a unique Route and Node; '
                f"duplicates: {', '.join(duplicate_keys)}"
            )

        route_ids = {route.short_id for route, _ in output_routes}
        unexpected_routes = set(timing_table['Route']) - route_ids
        if unexpected_routes:
            raise ValueError(
                'Every timing-table Route must occur in the route table; '
                f"unexpected routes: {', '.join(map(str, unexpected_routes))}"
            )
        return timing_table
    
    @staticmethod
    def _route_probability_mapping(values, routes, label):
        if not isinstance(values, dict):
            raise TypeError(f'{label} must be a dictionary keyed by Route')
        expected_routes = {route.short_id for route in routes}
        observed_routes = set(values)
        if observed_routes != expected_routes:
            missing_routes = sorted(expected_routes - observed_routes)
            unexpected_routes = sorted(observed_routes - expected_routes)
            details = []
            if missing_routes:
                details.append('missing: ' + ', '.join(missing_routes))
            if unexpected_routes:
                details.append('unexpected: ' + ', '.join(unexpected_routes))
            raise ValueError(
                f'{label} route coverage does not match fitted routes ('
                + '; '.join(details)
                + ')'
            )

        probabilities = {}
        for route_id, probability in values.items():
            if (
                not isinstance(probability, Real)
                or isinstance(probability, (bool, np.bool_))
                or not np.isfinite(probability)
                or probability < 0
            ):
                raise ValueError(
                    f'{label} values must be finite and nonnegative'
                )
            probabilities[route_id] = float(probability)
        total_probability = sum(probabilities.values())
        if not np.isclose(total_probability, 1.0):
            raise ValueError(f'{label} values must sum to one')
        return probabilities

    @staticmethod
    def _get_route_interval_topology(route, representative_all_phasings):
        route_tree = route.route_tree
        node_attributes = route_tree.node_attributes
        timing_nodes = list(route_tree.timeable_nodes)
        if any(
            isinstance(node, (bool, np.bool_)) or not isinstance(node, Integral)
            for node in timing_nodes
        ):
            raise ValueError(
                'Self-describing timing archives require integer timing-node '
                'identifiers'
            )
        timing_column_by_node = {
            node: timing_column
            for timing_column, node in enumerate(timing_nodes)
        }

        def endpoint_source(node):
            attributes = node_attributes[node]
            if len(attributes['Successors']) == 0:
                return timingio.TIMING_SOURCE_ONE
            if attributes['WGD']:
                return timingio.TIMING_SOURCE_WGD
            if node not in timing_column_by_node:
                raise ValueError(
                    f'Nonterminal non-WGD route node {node!r} has no stored '
                    'timing column'
                )
            return (
                timingio.TIMING_SOURCE_COLUMN_OFFSET
                + timing_column_by_node[node]
            )

        eligible_nodes = list(node_attributes)
        if representative_all_phasings:
            eligible_nodes = [
                node for node in eligible_nodes
                if node_attributes[node]['Phasing'] == 'Major'
            ]
        multiplicity_order = list(dict.fromkeys(
            node_attributes[node]['Multiplicity']
            for node in eligible_nodes
        ))
        ordered_nodes = [
            node
            for multiplicity in multiplicity_order
            for node in eligible_nodes
            if node_attributes[node]['Multiplicity'] == multiplicity
        ]

        start_sources = []
        end_sources = []
        multiplicities = []
        phasing_masks = []
        all_phasing_bits = np.uint8(
            timingio.PHASING_BIT_NON_PHASED
            | timingio.PHASING_BIT_MAJOR
            | timingio.PHASING_BIT_MINOR
        )
        for node in ordered_nodes:
            attributes = node_attributes[node]
            predecessor = attributes['Predecessor']
            start_sources.append(
                timingio.TIMING_SOURCE_ZERO
                if predecessor is None
                else endpoint_source(predecessor)
            )
            end_sources.append(endpoint_source(node))
            multiplicities.append(attributes['Multiplicity'])
            if representative_all_phasings:
                phasing_masks.append(all_phasing_bits)
            else:
                allele_bit = {
                    'Major': timingio.PHASING_BIT_MAJOR,
                    'Minor': timingio.PHASING_BIT_MINOR,
                }.get(attributes['Phasing'])
                if allele_bit is None:
                    raise ValueError(
                        'Every route node must have Major or Minor phasing'
                    )
                phasing_masks.append(np.uint8(
                    timingio.PHASING_BIT_NON_PHASED | allele_bit
                ))

        return {
            timingio.ROUTE_PARTICLE_TIMING_NODE_ID_KEY: np.asarray(
                timing_nodes,
                dtype=np.int64,
            ),
            timingio.ROUTE_PARTICLE_INTERVAL_START_SOURCE_KEY: np.asarray(
                start_sources,
                dtype=np.int64,
            ),
            timingio.ROUTE_PARTICLE_INTERVAL_END_SOURCE_KEY: np.asarray(
                end_sources,
                dtype=np.int64,
            ),
            timingio.ROUTE_PARTICLE_INTERVAL_MULTIPLICITY_KEY: np.asarray(
                multiplicities,
                dtype=np.int64,
            ),
            timingio.ROUTE_PARTICLE_INTERVAL_PHASING_KEY: np.asarray(
                phasing_masks,
                dtype=np.uint8,
            ),
        }

    @staticmethod
    def _get_route_state_columns(route, representative_all_phasings):
        major_cn = route.major_cn
        minor_cn = route.minor_cn
        n_subclones = route.n_subclones
        n_columns = 2 * major_cn + minor_cn + n_subclones
        subclone_columns = np.arange(
            n_columns - n_subclones,
            n_columns,
            dtype=np.int64,
        )
        state_columns = (
            np.concatenate([
                np.arange(major_cn, dtype=np.int64),
                subclone_columns,
            ]),
            np.concatenate([
                np.arange(major_cn, 2 * major_cn, dtype=np.int64),
                subclone_columns,
            ]),
            np.concatenate([
                np.arange(
                    2 * major_cn,
                    2 * major_cn + minor_cn,
                    dtype=np.int64,
                ),
                subclone_columns,
            ]),
        )
        if representative_all_phasings:
            representative_columns = state_columns[1]
            state_columns = (
                representative_columns,
                representative_columns,
                representative_columns,
            )

        offsets = np.zeros(
            len(timingio.ROUTE_PARTICLE_PHASING_ORDER) + 1,
            dtype=np.int64,
        )
        offsets[1:] = np.cumsum([columns.size for columns in state_columns])
        return {
            timingio.ROUTE_PARTICLE_STATE_COLUMN_OFFSETS_KEY: offsets,
            timingio.ROUTE_PARTICLE_STATE_COLUMNS_KEY: np.concatenate(
                state_columns
            ),
        }

    def get_timing_dict(
        self,
        penalized_probabilities=None,
        *,
        target_major_cn=None,
        target_minor_cn=None,
        target_wgd_status=None,
        archive_kind=timingio.ROUTE_PARTICLE_ARCHIVE_KIND_SEGMENT,
    ):
        """Return the fitted self-describing posterior payload, if needed."""
        if self.timing_representation == UNIFORM_NO_GAIN_REPRESENTATION:
            n_subclones = {
                route.n_subclones for route, _ in self._get_output_routes()
            }
            if len(n_subclones) != 1 or None in n_subclones:
                raise ValueError(
                    'Uniform no-gain routes must have one fitted subclone '
                    'count'
                )
            if n_subclones.pop() == 0:
                return None

            timing_dict = {}
            for route, _ in self._get_output_routes():
                clone_share_store = np.asarray(route.clone_share_store)
                expected_shape = (
                    ROUTE_CONDITIONAL_SAMPLE_COUNT,
                    route.n_subclones + 1,
                )
                if clone_share_store.shape != expected_shape:
                    raise ValueError(
                        'Uniform no-gain Clone_Share must have shape '
                        f'{expected_shape}'
                    )
                timing_dict[route.short_id] = {
                    'Clone_Share': clone_share_store.copy(),
                }
            return timing_dict

        output_routes = [route for route, _ in self._get_output_routes()]
        ordinary_probabilities = self._route_probability_mapping(
            {
                route.short_id: self.route_probabilities[route.route_id]
                for route in output_routes
            },
            output_routes,
            'Probability',
        )
        penalized_probabilities = self._route_probability_mapping(
            penalized_probabilities,
            output_routes,
            'Penalized_Probability',
        )

        model_major_cn = self.major_cn
        model_minor_cn = self.minor_cn
        model_wgd_status = bool(self.wgd_status)
        if target_major_cn is None:
            target_major_cn = model_major_cn
        if target_minor_cn is None:
            target_minor_cn = model_minor_cn
        if target_wgd_status is None:
            target_wgd_status = model_wgd_status
        if (
            isinstance(target_major_cn, (bool, np.bool_))
            or not isinstance(target_major_cn, Integral)
            or target_major_cn < 1
            or isinstance(target_minor_cn, (bool, np.bool_))
            or not isinstance(target_minor_cn, Integral)
            or target_minor_cn < 0
            or not isinstance(target_wgd_status, (bool, np.bool_))
        ):
            raise ValueError('Target copy-number/WGD metadata is invalid')
        representative_all_phasings = bool(
            target_wgd_status
            and target_major_cn == 2
            and target_minor_cn == 2
            and model_major_cn == 2
            and model_minor_cn == 0
            and not model_wgd_status
        )
        if archive_kind not in (
            timingio.ROUTE_PARTICLE_ARCHIVE_KIND_SEGMENT,
            timingio.ROUTE_PARTICLE_ARCHIVE_KIND_POOLED_WGD,
        ):
            raise ValueError('archive_kind must identify a segment or pooled WGD')
        target_model = (
            int(target_major_cn),
            int(target_minor_cn),
            bool(target_wgd_status),
        )
        fitted_model = (
            int(model_major_cn),
            int(model_minor_cn),
            bool(model_wgd_status),
        )
        if (
            archive_kind == timingio.ROUTE_PARTICLE_ARCHIVE_KIND_SEGMENT
            and target_model != fitted_model
        ):
            raise ValueError(
                'Segment archive target metadata must match its fitted model'
            )
        if archive_kind == timingio.ROUTE_PARTICLE_ARCHIVE_KIND_POOLED_WGD:
            expected_fitted_model = (
                2,
                0 if target_minor_cn == 2 else int(target_minor_cn),
                False,
            )
            if (
                target_major_cn != 2
                or target_minor_cn not in (0, 1, 2)
                or not target_wgd_status
                or fitted_model != expected_fitted_model
                or len(output_routes) != 1
                or ordinary_probabilities[output_routes[0].short_id] != 1.0
                or penalized_probabilities[output_routes[0].short_id] != 1.0
            ):
                raise ValueError(
                    'Pooled WGD archives require one unit-probability 2+m '
                    'target with the corresponding non-WGD fitted model'
                )

        timing_dict = {}
        for route in output_routes:
            mult_store = np.array(
                route.mult_store,
                dtype=np.float64,
                order='C',
                copy=True,
            )
            if (
                not np.all(np.isfinite(mult_store))
                or np.any(mult_store < -_MULT_NEGATIVE_TOLERANCE)
            ):
                raise ValueError(
                    'Self-describing route particles require finite Mult '
                    'values without material negative mass'
                )
            np.maximum(mult_store, 0.0, out=mult_store)
            wgd_timing = np.ascontiguousarray(
                route.wgd_timing_store,
                dtype=np.float64,
            )
            n_particles = mult_store.shape[0]
            timing_nodes = list(route.route_tree.timeable_nodes)
            if timing_nodes:
                timing = np.ascontiguousarray(np.column_stack([
                    route.get_node_timing(node) for node in timing_nodes
                ]), dtype=np.float64)
            else:
                timing = np.empty((n_particles, 0), dtype=np.float64)
            if (
                mult_store.ndim != 2
                or n_particles == 0
                or wgd_timing.shape != (n_particles,)
                or timing.shape != (n_particles, len(timing_nodes))
            ):
                raise ValueError(
                    'Aligned Timing, WGD_Timing, and Mult particles must '
                    'contain the same nonzero number of rows'
                )
            expected_mult_columns = (
                2 * model_major_cn
                + model_minor_cn
                + route.n_subclones
            )
            if (
                mult_store.shape[1] != expected_mult_columns
                or not np.all(np.isfinite(timing))
            ):
                raise ValueError(
                    'Self-describing route particles require finite Timing '
                    'and finite nonnegative Mult values with the modelled '
                    'state-column count'
                )
            subclone_start = 2 * model_major_cn + model_minor_cn
            subclone_share = np.sum(
                mult_store[:, subclone_start:],
                axis=1,
            )
            clonal_blocks = (
                mult_store[:, :model_major_cn],
                mult_store[:, model_major_cn:2 * model_major_cn],
            )
            if model_minor_cn:
                clonal_blocks += (
                    mult_store[
                        :,
                        2 * model_major_cn:subclone_start,
                    ],
                )
            if any(
                not np.allclose(
                    np.sum(clonal_block, axis=1) + subclone_share,
                    1.0,
                    rtol=1e-7,
                    atol=1e-8,
                )
                for clonal_block in clonal_blocks
            ):
                raise ValueError(
                    'Every Mult phasing block plus the shared subclone '
                    'columns must sum to one'
                )
            if model_wgd_status:
                if not np.all(np.isfinite(wgd_timing)):
                    raise ValueError(
                        'WGD route particles require finite WGD_Timing values'
                    )
            elif not np.all(np.isnan(wgd_timing)):
                raise ValueError(
                    'Non-WGD route particles require NaN WGD_Timing values'
                )

            route_samples = {
                timingio.ROUTE_PARTICLE_PROBABILITY_KEY: np.asarray(
                    [ordinary_probabilities[route.short_id]],
                    dtype=np.float64,
                ),
                timingio.ROUTE_PARTICLE_PENALIZED_PROBABILITY_KEY: np.asarray(
                    [penalized_probabilities[route.short_id]],
                    dtype=np.float64,
                ),
                timingio.ROUTE_PARTICLE_TIMING_KEY: timing.copy(),
                timingio.ROUTE_PARTICLE_WGD_TIMING_KEY: wgd_timing.copy(),
                timingio.ROUTE_PARTICLE_MULT_KEY: mult_store.copy(),
                timingio.ROUTE_PARTICLE_TARGET_MAJOR_CN_KEY: np.asarray(
                    [target_major_cn], dtype=np.int64
                ),
                timingio.ROUTE_PARTICLE_TARGET_MINOR_CN_KEY: np.asarray(
                    [target_minor_cn], dtype=np.int64
                ),
                timingio.ROUTE_PARTICLE_N_SUBCLONES_KEY: np.asarray(
                    [route.n_subclones], dtype=np.int64
                ),
                timingio.ROUTE_PARTICLE_TARGET_WGD_STATUS_KEY: np.asarray(
                    [target_wgd_status], dtype=np.uint8
                ),
                timingio.ROUTE_PARTICLE_MODEL_MAJOR_CN_KEY: np.asarray(
                    [model_major_cn], dtype=np.int64
                ),
                timingio.ROUTE_PARTICLE_MODEL_MINOR_CN_KEY: np.asarray(
                    [model_minor_cn], dtype=np.int64
                ),
                timingio.ROUTE_PARTICLE_MODEL_WGD_STATUS_KEY: np.asarray(
                    [model_wgd_status], dtype=np.uint8
                ),
                timingio.ROUTE_PARTICLE_ARCHIVE_KIND_KEY: np.asarray(
                    [archive_kind], dtype=np.uint8
                ),
            }
            route_samples.update(self._get_route_interval_topology(
                route,
                representative_all_phasings,
            ))
            route_samples.update(self._get_route_state_columns(
                route,
                representative_all_phasings,
            ))
            if set(route_samples) != timingio.ROUTE_PARTICLE_ARCHIVE_FIELDS:
                raise AssertionError(
                    'Self-describing route-particle fields are incomplete'
                )
            timing_dict[route.short_id] = route_samples
        return timing_dict

    def get_timing_tree_labels(
        self,
        route,
        wgd_info,
        gain_interval=DEFAULT_TIMING_INTERVALS.tree_gain,
        rounding_digits=3,
    ):
        node_labels = {}
        for node in route.route_tree.non_phased_node_order:
            if node in route.route_tree.timeable_nodes:
                timing_dist = np.asarray(route.get_node_timing(node))
                timing = np.round(np.median(timing_dist), rounding_digits)
                timing_ci_low, timing_ci_high = get_interval_bounds(
                    timing_dist,
                    gain_interval,
                )
                timing_ci_low = np.round(timing_ci_low, rounding_digits)
                timing_ci_high = np.round(timing_ci_high, rounding_digits)
                node_labels[node] = f"{timing} - [{timing_ci_low},{timing_ci_high}]"
            elif node in route.route_tree.wgd_nodes:
                node_labels[node] = (
                    f"{wgd_info['WGD_Timing']} - "
                    f"[{wgd_info['WGD_Timing_CI_Low']},"
                    f"{wgd_info['WGD_Timing_CI_High']}]"
                )
            else:
                node_labels[node] = ''
        return node_labels
    
    def plot_trees(
        self,
        plot_output_dir,
        seg_title,
        wgd_info,
        gain_interval=DEFAULT_TIMING_INTERVALS.tree_gain,
    ):
        os.makedirs(plot_output_dir,exist_ok=True)

        if len(self.route_probabilities.keys())==0:
            raise ValueError("The routes haven't been fit yet")  
        best_route_id = max(self.route_probabilities, key=self.route_probabilities.get)   
        for route_id,route in self.routes.items():
            
            probability = np.round(self.route_probabilities[route_id],4)
            
            route_output_path = f"{plot_output_dir}/route_{route.short_id}.pdf"
            tree_subtitle = f"Route {route.short_id} (Probability = {probability})"
            plot_title = f"{seg_title}\n{tree_subtitle}"
            if route_id==best_route_id:
                plot_title = f"{plot_title} - (Best Fit)"

            node_labels = self.get_timing_tree_labels(
                route,
                wgd_info,
                gain_interval=gain_interval,
            )
            plotting_tree = route.route_tree.main_tree.copy()
            nx.set_node_attributes(plotting_tree,node_labels,'Label')
            
            treetools.plot_tree(plotting_tree,plot_title,output_path=route_output_path)


def fit_route_classifiers(classifier_jobs, wgd_timing_distribution):
    """Fit equal-geometry classifiers while retaining one route store at a time.

    Each job is ``(classifier, multiplicity_probabilities, subclone_table)``.
    Route geometry is generated once per copy-number state, then segment-specific
    clone fractions and likelihoods are applied independently.
    """
    classifier_jobs = list(classifier_jobs)
    if not classifier_jobs:
        return

    reference_classifier = classifier_jobs[0][0]
    reference_state = (
        reference_classifier.major_cn,
        reference_classifier.minor_cn,
        reference_classifier.wgd_status,
        tuple(sorted(reference_classifier.routes)),
    )
    for classifier, _, _ in classifier_jobs:
        classifier_state = (
            classifier.major_cn,
            classifier.minor_cn,
            classifier.wgd_status,
            tuple(sorted(classifier.routes)),
        )
        if classifier_state != reference_state:
            raise ValueError(
                'Shared route geometry requires matching copy-number state, '
                'WGD status, and route IDs'
            )
        classifier._prepare_fit()

    route_ids = reference_state[-1]
    processed_route_ids = set()
    n_classifiers = len(classifier_jobs)

    for source_route_id in route_ids:
        if source_route_id in processed_route_ids:
            continue
        source_reference = reference_classifier.routes[source_route_id]
        mirror_route_id = source_reference.mirror_route_id
        if (
            mirror_route_id is not None
            and mirror_route_id not in reference_classifier.routes
        ):
            raise ValueError(
                'Every balanced ordered route must have its mirror in the '
                'route classifier'
            )

        source_geometry = None
        shared_geometry_time = 0.0
        if not source_reference.is_uniform_no_gain:
            geometry_started = time.perf_counter()
            source_geometry = source_reference.run_geometry_sampling(
                wgd_timing_distribution
            )
            shared_geometry_time = (
                time.perf_counter() - geometry_started
            ) / n_classifiers

        has_distinct_mirror = mirror_route_id not in (None, source_route_id)
        target_geometry = None
        target_geometry_time = 0.0
        if has_distinct_mirror:
            target_reference = reference_classifier.routes[mirror_route_id]
            if target_reference.mirror_route_id != source_route_id:
                raise ValueError(
                    'Balanced route mirror relationships must be reciprocal'
                )
            needs_mirror_geometry = any(
                Route._contains_phased_mutations(mult_probabilities)
                for _, mult_probabilities, _ in classifier_jobs
            )
            if needs_mirror_geometry:
                transform_started = time.perf_counter()
                target_geometry = target_reference.transform_mirror_geometry(
                    source_reference,
                    source_geometry,
                )
                target_geometry_time = (
                    time.perf_counter() - transform_started
                ) / n_classifiers

        for classifier, mult_probabilities, subclone_table in classifier_jobs:
            source_route = classifier.routes[source_route_id]
            clone_share = source_route.run_sampling(
                mult_probabilities,
                subclone_table,
                wgd_timing_distribution,
                classifier.subclone_fraction_prior,
                proposal_geometry=source_geometry,
                shared_geometry_time=shared_geometry_time,
            )
            if has_distinct_mirror:
                classifier.routes[mirror_route_id].run_sampling(
                    mult_probabilities,
                    subclone_table,
                    wgd_timing_distribution,
                    classifier.subclone_fraction_prior,
                    mirror_route=source_route,
                    proposal_geometry=target_geometry,
                    clone_share=clone_share,
                    shared_geometry_time=target_geometry_time,
                )
            del clone_share
        processed_route_ids.add(source_route_id)
        if has_distinct_mirror:
            processed_route_ids.add(mirror_route_id)
        del source_geometry, target_geometry

    for classifier, _, _ in classifier_jobs:
        classifier._finalize_fit()


def estimate_classifier_output_bytes(route_trees, n_subclones):
    """Estimate retained numeric posterior bytes for one fitted classifier."""
    total_bytes = 0
    for route_tree in route_trees.values():
        if is_uniform_no_gain_route_tree(route_tree):
            if n_subclones > 0:
                total_bytes += 8 * ROUTE_CONDITIONAL_SAMPLE_COUNT * (
                    n_subclones + 1
                )
            total_bytes += 8 * 3
            continue

        mult_columns = (
            2 * route_tree.major_cn
            + route_tree.minor_cn
            + n_subclones
        )
        conditional_values = ROUTE_CONDITIONAL_SAMPLE_COUNT * (
            len(route_tree.non_phased_node_order)
            + mult_columns
            + 1  # WGD timing
        )
        event_values = 3 * 300
        total_bytes += 8 * (
            conditional_values + event_values
        )
    return total_bytes


def add_wgd_info_to_route_table(route_table,wgd_info,wgd_status:bool):
    route_table = route_table.copy()
    route_table['WGD_Timing'] = wgd_info['WGD_Timing']
    route_table['WGD_Timing_CI_Low'] = wgd_info['WGD_Timing_CI_Low']
    route_table['WGD_Timing_CI_High'] = wgd_info['WGD_Timing_CI_High']
    route_table['WGD_Status'] = wgd_status
    return route_table


def _append_table(table, table_path, columns, column_dtypes=None):
    table = table[columns]
    if os.path.exists(table_path):
        string_converters = {
            column: str
            for column in ('Sample_ID', 'Segment_ID', 'Chromosome', 'Route')
            if column in columns
        }
        previous_table = pd.read_csv(
            table_path,
            sep='\t',
            converters=string_converters,
        )
        table = pd.concat([previous_table, table], ignore_index=True)
    if column_dtypes is not None:
        table = table.astype(column_dtypes)
    _write_tsv(table, table_path)


def _write_tsv(table, table_path):
    """Write a public TSV."""
    table.to_csv(table_path, sep='\t', index=False)


def write_route_table(route_table, route_table_path):
    _append_table(route_table, route_table_path, ROUTE_TABLE_COLUMNS)


def write_gain_timing_table(timing_table, timing_table_path):
    _append_table(
        timing_table,
        timing_table_path,
        GAIN_TIMING_TABLE_COLUMNS,
        column_dtypes={'Node': 'Int64'},
    )


def _initialize_table(table_path, columns):
    if not os.path.exists(table_path):
        _write_tsv(
            pd.DataFrame(columns=columns),
            table_path,
        )

def get_wgd_info(
    wgd_timing_distribution,
    interval=DEFAULT_TIMING_INTERVALS.sample_wgd,
    rounding_digits=3,
):
    if wgd_timing_distribution is None:
        wgd_timing = np.nan
        wgd_timing_ci_high = np.nan
        wgd_timing_ci_low = np.nan
    else:
        wgd_timing = np.round(
            np.median(wgd_timing_distribution),
            rounding_digits,
        )
        wgd_timing_ci_low, wgd_timing_ci_high = get_interval_bounds(
            wgd_timing_distribution,
            interval,
        )
        wgd_timing_ci_low = np.round(
            wgd_timing_ci_low,
            rounding_digits,
        )
        wgd_timing_ci_high = np.round(
            wgd_timing_ci_high,
            rounding_digits,
        )
    return {
        'WGD_Timing': wgd_timing,
        'WGD_Timing_CI_Low': wgd_timing_ci_low,
        'WGD_Timing_CI_High': wgd_timing_ci_high,
    }


def get_potential_wgd_segments(
    sample,
    min_wgd_mutations: int = MIN_WGD_MUTATIONS,
):
    valid_autosomes = set(sample.autosomes)
    potential_wgd_segments = []
    for segment in sample.segments:
        if segment.major_cn ==2 and segment.n_mutations >= min_wgd_mutations and segment.chromosome in valid_autosomes:
            potential_wgd_segments.append(segment)
    return potential_wgd_segments

def get_combined_distribution(
    distributions,
    n_samples=COMBINED_WGD_SAMPLE_COUNT,
    eps=1e-300,
):
    bins = np.linspace(0,1,201)
    bin_mid_points = (bins[1:]+bins[:(bins.size-1)])/2
    binned_distributions = []
    for distribution in distributions:
        distribution = np.clip(distribution,eps,1.0-eps)
        binned_distribution = np.histogram(distribution,bins=bins)[0]
        binned_distribution = binned_distribution/np.sum(binned_distribution)
        binned_distributions.append(binned_distribution)
    binned_distributions = np.array(binned_distributions)
    binned_distributions = np.log(binned_distributions+eps)
    
    combined_distribution = np.sum(binned_distributions,axis=0)

    combined_distribution = np.exp(combined_distribution-logsumexp(combined_distribution))
   
    combined_distribution_samples = np.random.choice(bin_mid_points,p=combined_distribution,size=n_samples,replace=True)
    return combined_distribution_samples


def _get_wgd_timing_model(segment):
    """Return the copy-number geometry and likelihood store for WGD timing.

    Balanced 2+2 segments represent two copies of each parental allele created
    by the same WGD.  Their WGD timing is therefore fit with one representative
    duplicated allele (2+0 geometry).  All mutations, including phased rows,
    contribute once through the original segment's combined likelihood arrays.
    The original segment and its multiplicity store are left unchanged.
    """
    mult_probabilities = segment.multiplicity_probabilities
    if segment.major_cn != 2 or segment.minor_cn != 2:
        return segment.minor_cn, mult_probabilities

    pseudo_minor_cn = 0

    # The balanced likelihood matrix is assembled column-wise and can be
    # Fortran-contiguous.  The cached Numba likelihood kernel requires the
    # mutation-by-state input to be C-contiguous; passing the transposed layout
    # can crash inside compiled code rather than raising a Python exception.
    combined_array = np.ascontiguousarray(
        mult_probabilities.combined_array,
        dtype=np.float64,
    )
    combined_correction = np.ascontiguousarray(
        mult_probabilities.reads_correction_combined_array,
        dtype=np.float64,
    )
    combined_array_store = {
        'Non_Phased': combined_array,
        'Major': None,
        'Minor': None,
        'All': combined_array,
    }
    combined_correction_store = {
        'Non_Phased': combined_correction,
        'Major': None,
        'Minor': None,
        'All': combined_correction,
    }
    wgd_mult_probabilities = MultProbabilityStore(
        combined_array_store,
        combined_correction_store,
        major_cn=segment.major_cn,
        minor_cn=pseudo_minor_cn,
        n_subclones=mult_probabilities.n_subclones,
    )
    return pseudo_minor_cn, wgd_mult_probabilities


def _get_wgd_segment_result(segment, classifier, interval_config):
    wgd_route_table = classifier.get_route_table()
    wgd_timing_table = classifier.get_timing_table(
        interval=interval_config.route_gain,
    )
    wgd_timing_table = pd.merge(
        wgd_route_table,
        wgd_timing_table,
        on='Route',
        how='inner',
        validate='one_to_many',
    )
    for key, val in segment.get_info_dict().items():
        wgd_timing_table[key] = val

    classifier_route = list(classifier.routes.values())[0]
    segment_timing = classifier_route.get_node_timing(0)
    return wgd_timing_table, segment_timing


def time_wgd_major_cn_2(
    sample,
    output_dir,
    timing_dict_dir,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
    unordered_balanced_route_prior=DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
):
    subclone_fraction_prior = validate_subclone_fraction_prior(
        subclone_fraction_prior
    )
    unordered_balanced_route_prior = validate_unordered_balanced_route_prior(
        unordered_balanced_route_prior
    )

    wgd_timing_table_path = output_dir/f"{sample.sample_id}_gain_timing_table_wgd_segments.tsv"
    potential_wgd_segments = get_potential_wgd_segments(sample)
    if not potential_wgd_segments:
        raise ValueError(
            "No eligible segments with enough mutations for WGD timing in "
            f"sample {sample.sample_id}: at least one autosomal major-copy-"
            f"number-2 segment with {MIN_WGD_MUTATIONS} or more retained "
            "mutations is required."
        )
 
    segment_ci_store = {}
    segment_width_store = {}

    jobs_by_pseudo_minor_cn = {}
    for segment in potential_wgd_segments:
        pseudo_minor_cn, mult_probabilities = _get_wgd_timing_model(segment)
        state_jobs = jobs_by_pseudo_minor_cn.setdefault(
            pseudo_minor_cn,
            {'route_trees': None, 'records': []},
        )
        classifier = RouteClassifier(
            segment.major_cn,
            pseudo_minor_cn,
            False,
            'No_WGD',
            subclone_fraction_prior=subclone_fraction_prior,
            unordered_balanced_route_prior=unordered_balanced_route_prior,
            route_trees=state_jobs['route_trees'],
        )
        if state_jobs['route_trees'] is None:
            state_jobs['route_trees'] = {
                route_id: route.route_tree
                for route_id, route in classifier.routes.items()
            }
        state_jobs['records'].append((
            segment,
            classifier,
            mult_probabilities,
            getattr(segment, 'subclone_table', None),
        ))

    results_by_segment_id = {}
    for state_jobs in jobs_by_pseudo_minor_cn.values():
        fit_route_classifiers(
            [
                (classifier, mult_probabilities, subclone_table)
                for (
                    _,
                    classifier,
                    mult_probabilities,
                    subclone_table,
                ) in state_jobs['records']
            ],
            None,
        )
        for segment, classifier, _, _ in state_jobs['records']:
            results_by_segment_id[id(segment)] = _get_wgd_segment_result(
                segment,
                classifier,
                interval_config,
            )

    wgd_timing_tables = []
    for segment in potential_wgd_segments:
        logger.info('Timing WGD segment: %s', segment)
        wgd_timing_table, segment_timing = results_by_segment_id[id(segment)]
        wgd_timing_tables.append(wgd_timing_table)

        if not np.isfinite(segment_timing).all():
            continue
        segment_timing_ci_low, segment_timing_ci_high = (
            get_interval_bounds(
                segment_timing,
                interval_config.wgd_overlap,
            )
        )
        segment_ci_store[segment.segment_id] = (
            segment_timing_ci_low,
            segment_timing_ci_high,
        )
        segment_width_store[segment.segment_id] = segment.width

    if not segment_ci_store:
        raise ValueError(
            "No eligible major-copy-number-2 segment produced a finite WGD "
            f"timing interval for sample {sample.sample_id}."
        )
    
    wgd_timing_table = pd.concat(wgd_timing_tables)
    wgd_timing_table = wgd_timing_table.drop(columns=['Node','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Probability'])
 
    max_segment_width = sum(segment_width_store.values())
    segments_with_max_overlap,overlap_width,best_overlap_timing = distributiontools.get_ids_with_maximum_overlap(segment_ci_store,segment_width_store)
    
    overlap_proportion = overlap_width/max_segment_width
    
    wgd_timing_table['Best_Overlap_Timing'] = best_overlap_timing
    wgd_timing_table['Overlap_Proportion'] = overlap_proportion
    wgd_timing_table['Intersecting'] = wgd_timing_table['Segment_ID'].isin(segments_with_max_overlap)
    
    _write_tsv(wgd_timing_table, wgd_timing_table_path)
    
   
    overlapping_segments = [segment for segment in potential_wgd_segments if segment.segment_id in segments_with_max_overlap]
    non_overlapping_segment_ids = [segment.segment_id for segment in potential_wgd_segments if segment.segment_id not in segments_with_max_overlap]
    
    cn_distributions = get_combined_segment_timing_cn_2(
        overlapping_segments,
        sample.subclone_table,
        sample.purity,
        timing_dict_dir,
        subclone_fraction_prior=subclone_fraction_prior,
        unordered_balanced_route_prior=unordered_balanced_route_prior,
    )
    wgd_timing_distribution = get_combined_distribution(cn_distributions)
    
    return wgd_timing_distribution,non_overlapping_segment_ids,overlap_proportion,best_overlap_timing


def _time_combined_wgd_segment(
    minor_cn,
    mutation_table,
    combined_width,
    subclone_table,
    sample_purity,
    apply_reads_correction,
    min_mutation_alt_count,
    coverage_vaf_quantile,
    timing_dict_dir,
    subclone_fraction_prior,
    unordered_balanced_route_prior,
):
    mutation_table = mutation_table.copy()
    mutation_table['Segment_ID'] = f'Minor_CN_{minor_cn}'
    mutation_table['Chromosome'] = np.nan
    mutation_table['Segment_Start'] = 0
    mutation_table['Segment_End'] = combined_width

    new_seg = Segment(
        mutation_table,
        subclone_table,
        sample_purity,
        sex=None,
        apply_reads_correction=apply_reads_correction,
        min_mutation_alt_count=min_mutation_alt_count,
        coverage_vaf_quantile=coverage_vaf_quantile,
    )
    pseudo_minor_cn, mult_probabilities = _get_wgd_timing_model(new_seg)
    classifier = RouteClassifier(
        new_seg.major_cn,
        pseudo_minor_cn,
        False,
        'No_WGD',
        subclone_fraction_prior=subclone_fraction_prior,
        unordered_balanced_route_prior=unordered_balanced_route_prior,
    )
    classifier.fit_routes(
        mult_probabilities,
        new_seg.subclone_table,
        None,
    )

    pooled_route_table = classifier.get_route_table()
    if (
        len(pooled_route_table) != 1
        or not np.isclose(pooled_route_table.iloc[0]['Probability'], 1.0)
    ):
        raise ValueError('Pooled WGD timing must have exactly one route')
    pooled_route_id = pooled_route_table.iloc[0]['Route']
    timing_dict = classifier.get_timing_dict(
        {pooled_route_id: 1.0},
        target_major_cn=2,
        target_minor_cn=minor_cn,
        target_wgd_status=True,
        archive_kind=timingio.ROUTE_PARTICLE_ARCHIVE_KIND_POOLED_WGD,
    )
    write_timing_archive(
        timing_dict,
        timing_dict_dir,
        f'WGD_minor_cn_{minor_cn}',
    )
    return classifier.get_best_timing()[0]



def get_combined_segment_timing_cn_2(
    overlapping_segments,
    subclone_table,
    sample_purity,
    timing_dict_dir,
    subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
    unordered_balanced_route_prior=DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
):
    subclone_fraction_prior = validate_subclone_fraction_prior(
        subclone_fraction_prior
    )
    unordered_balanced_route_prior = validate_unordered_balanced_route_prior(
        unordered_balanced_route_prior
    )

    mutation_tables = []
    combined_width_by_minor_cn = {}
    min_mutation_alt_counts = {
        segment.min_mutation_alt_count
        for segment in overlapping_segments
    }
    if len(min_mutation_alt_counts) != 1:
        raise ValueError(
            'Combined WGD segments must use one minimum mutation '
            'alternate-read count'
        )
    min_mutation_alt_count = min_mutation_alt_counts.pop()
    coverage_vaf_quantiles = {
        segment.coverage_vaf_quantile
        for segment in overlapping_segments
    }
    if len(coverage_vaf_quantiles) != 1:
        raise ValueError(
            'Combined WGD segments must use one coverage VAF quantile'
        )
    coverage_vaf_quantile = coverage_vaf_quantiles.pop()
    apply_reads_correction_values = {
        segment.apply_reads_correction
        for segment in overlapping_segments
    }
    if len(apply_reads_correction_values) != 1:
        raise ValueError(
            'Combined WGD segments must use one reads-correction setting'
        )
    apply_reads_correction = apply_reads_correction_values.pop()
    for segment in overlapping_segments:
        mutation_tables.append(segment.mutation_table)
        combined_width_by_minor_cn[segment.minor_cn] = (
            combined_width_by_minor_cn.get(segment.minor_cn, 0)
            + int(segment.width)
        )
    combined_mutation_table = pd.concat(mutation_tables)
    assert len(combined_mutation_table['Major_CN'].unique()) ==1
    segment_timing_store = []
    grouped_mutations = combined_mutation_table.groupby('Minor_CN')
    for minor_cn, minor_cn_mutation_table in grouped_mutations:
        segment_timing_store.append(
            _time_combined_wgd_segment(
                minor_cn,
                minor_cn_mutation_table,
                combined_width_by_minor_cn[minor_cn],
                subclone_table,
                sample_purity,
                apply_reads_correction,
                min_mutation_alt_count,
                coverage_vaf_quantile,
                timing_dict_dir,
                subclone_fraction_prior,
                unordered_balanced_route_prior,
            )
        )
    
    return segment_timing_store


def check_permitted_cn_state(major_cn,minor_cn,wgd_status):
    #enforcing no more than 500 routes per cn state
    if minor_cn> major_cn:
        return False
    if minor_cn <0:
        return False
    if major_cn <=0:
        return False
    if major_cn ==1:
        return True
    if major_cn ==2:
        #we don't produce timing distributions for WGD segments
        if wgd_status:
            return False
        return True
    if major_cn >=3 and major_cn <=5:
        return True
    if major_cn==6 and minor_cn<=4:
        return True
    if  major_cn ==7 and minor_cn<=3:
        return True
    if major_cn ==8 and minor_cn<=1:
        return True
    return False

def get_timeable_segments(
    sample,
    *,
    wgd_status,
    excluded_segment_ids=None,
    min_mutations=0,
):

    complex_segments = {}
    excluded_segment_ids = (
        set()
        if excluded_segment_ids is None
        else set(excluded_segment_ids)
    )
    for segment in sample.segments:
        if segment.segment_id in excluded_segment_ids:
            continue
        if check_permitted_cn_state(segment.major_cn,segment.minor_cn,wgd_status) and segment.n_mutations >= min_mutations:
            segment_cn_state = (segment.major_cn,segment.minor_cn)
            if segment_cn_state not in complex_segments:
                complex_segments[segment_cn_state] = []
            complex_segments[segment_cn_state].append(segment)

    return complex_segments


def _json_scalar(value):
    if isinstance(value, np.generic):
        value = value.item()
    if (
        value is None
        or pd.isna(value)
        or (isinstance(value, float) and not np.isfinite(value))
    ):
        return None
    return value


def write_wgd_calling_info(
    wgd_info,
    overlap_proportion,
    best_overlap_timing,
    major_cn_mode,
    wgd_status,
    output_dir,
    sample_id,
):
    calling_info = {
        **wgd_info,
        'Major_CN_Mode': major_cn_mode,
        'Overlap_Proportion': overlap_proportion,
        'WGD_Status': wgd_status,
        'Best_Overlap_Timing': best_overlap_timing,
    }
    calling_info = {
        key: _json_scalar(value)
        for key, value in calling_info.items()
    }
    output_path = output_dir / f'{sample_id}_wgd_calling_info.json'
    with output_path.open('w', encoding='utf-8') as output_file:
        json.dump(calling_info, output_file, indent=2, allow_nan=False)
        output_file.write('\n')


def _write_segment_results(
    segment,
    classifier,
    output_dir,
    timing_dict_dir,
    sample_id,
    wgd_status,
    plot_trees,
    wgd_info,
    interval_config,
    route_table_path,
    timing_table_path,
):
    segment_route_table = classifier.get_route_table()
    segment_timing_table = classifier.get_timing_table(
        interval=interval_config.route_gain,
    )
    for key, val in segment.get_info_dict().items():
        segment_route_table[key] = val
    segment_route_table = add_wgd_info_to_route_table(
        segment_route_table,
        wgd_info,
        wgd_status,
    )
    segment_route_table['Sample_ID'] = sample_id
    segment_route_table = posteriortablegen.add_penalized_probability(
        segment_route_table
    )

    penalized_probabilities = dict(zip(
        segment_route_table['Route'],
        segment_route_table['Penalized_Probability'],
    ))
    timing_dict = classifier.get_timing_dict(
        penalized_probabilities,
        target_major_cn=segment.major_cn,
        target_minor_cn=segment.minor_cn,
        target_wgd_status=wgd_status,
    )
    if timing_dict is not None:
        write_timing_archive(
            timing_dict,
            timing_dict_dir,
            segment.segment_id,
        )

    segment_timing_table['Segment_ID'] = segment.segment_id
    segment_timing_table['Sample_ID'] = sample_id

    write_route_table(segment_route_table, route_table_path)
    write_gain_timing_table(segment_timing_table, timing_table_path)

    if plot_trees:
        plot_output_dir = (
            f'{output_dir}/{sample_id}_tree_plots/{segment.segment_id}'
        )
        classifier.plot_trees(
            plot_output_dir,
            str(segment),
            wgd_info,
            gain_interval=interval_config.tree_gain,
        )


def process_segments(
    segments,
    wgd_timing_distribution,
    output_dir,
    timing_dict_dir,
    sample_id,
    wgd_status,
    plot_trees,
    wgd_info,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
    unordered_balanced_route_prior=DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
):
    subclone_fraction_prior = validate_subclone_fraction_prior(
        subclone_fraction_prior
    )
    unordered_balanced_route_prior = validate_unordered_balanced_route_prior(
        unordered_balanced_route_prior
    )

    route_table_path = output_dir/f"{sample_id}_route_table.tsv"
    timing_table_path = output_dir/f"{sample_id}_gain_timing_table.tsv"

    _initialize_table(route_table_path, ROUTE_TABLE_COLUMNS)
    _initialize_table(timing_table_path, GAIN_TIMING_TABLE_COLUMNS)

    for cn_state, state_segments in segments.items():
        state_segments = list(state_segments)
        if not state_segments:
            continue

        prototype_classifier = RouteClassifier(
            cn_state[0],
            cn_state[1],
            wgd_status,
            'Default',
            subclone_fraction_prior=subclone_fraction_prior,
            unordered_balanced_route_prior=unordered_balanced_route_prior,
        )
        route_trees = {
            route_id: route.route_tree
            for route_id, route in prototype_classifier.routes.items()
        }
        max_subclones = max(
            0
            if getattr(segment, 'subclone_table', None) is None
            else len(segment.subclone_table.index)
            for segment in state_segments
        )
        estimated_output_bytes = estimate_classifier_output_bytes(
            route_trees,
            max_subclones,
        )
        batch_size = max(
            1,
            SHARED_FIT_OUTPUT_MEMORY_BUDGET
            // max(estimated_output_bytes, 1),
        )
        batch_size = min(len(state_segments), batch_size)
        logger.info(
            'Sharing route geometry for copy-number state %s across batches '
            'of at most %d segments',
            cn_state,
            batch_size,
        )

        prototype_available = True
        for batch_start in range(0, len(state_segments), batch_size):
            segment_batch = state_segments[
                batch_start:batch_start + batch_size
            ]
            classifier_jobs = []
            for segment in segment_batch:
                logger.info('Preparing gained segment: %s', segment)
                if prototype_available:
                    classifier = prototype_classifier
                    prototype_classifier = None
                    prototype_available = False
                else:
                    classifier = RouteClassifier(
                        cn_state[0],
                        cn_state[1],
                        wgd_status,
                        'Default',
                        subclone_fraction_prior=subclone_fraction_prior,
                        unordered_balanced_route_prior=(
                            unordered_balanced_route_prior
                        ),
                        route_trees=route_trees,
                    )
                classifier_jobs.append((
                    classifier,
                    segment.multiplicity_probabilities,
                    getattr(segment, 'subclone_table', None),
                ))

            fit_route_classifiers(
                classifier_jobs,
                wgd_timing_distribution,
            )

            for segment, (classifier, _, _) in zip(
                segment_batch,
                classifier_jobs,
            ):
                logger.info('Writing gained segment: %s', segment)
                _write_segment_results(
                    segment,
                    classifier,
                    output_dir,
                    timing_dict_dir,
                    sample_id,
                    wgd_status,
                    plot_trees,
                    wgd_info,
                    interval_config,
                    route_table_path,
                    timing_table_path,
                )
            del classifier_jobs
            del classifier

def _validate_wgd_count(wgd_count):
    if wgd_count is None:
        return None
    if isinstance(wgd_count, bool) or not isinstance(wgd_count, Integral):
        raise ValueError('wgd_count must be None or an integer equal to 0 or 1')
    if wgd_count not in (0, 1):
        raise ValueError('wgd_count must be None or an integer equal to 0 or 1')
    return int(wgd_count)


def _validate_min_wgd_overlap(min_wgd_overlap):
    if (
        isinstance(min_wgd_overlap, (bool, np.bool_))
        or not isinstance(min_wgd_overlap, Real)
        or not np.isfinite(min_wgd_overlap)
        or not 0 <= min_wgd_overlap <= 1
    ):
        raise ValueError(
            'min_wgd_overlap must be a finite number between 0 and 1'
        )
    return float(min_wgd_overlap)


def _validate_random_seed(random_seed):
    if random_seed is None:
        return None
    if (
        isinstance(random_seed, (bool, np.bool_))
        or not isinstance(random_seed, Integral)
        or not 0 <= random_seed <= np.iinfo(np.uint32).max
    ):
        raise ValueError(
            'random_seed must be None or an integer between 0 and 2**32 - 1'
        )
    return int(random_seed)


def _run_sample(
    sample,
    output_dir,
    timing_dict_dir,
    plot_trees,
    min_wgd_overlap,
    wgd_count,
    interval_config,
    subclone_fraction_prior,
    unordered_balanced_route_prior,
    major_cn_mode,
):
    if wgd_count is not None:
        if major_cn_mode == 1 and wgd_count == 1:
            warnings.warn(
                'Sample was specified with WGD count 1 but major CN mode is '
                '1. Please check the WGD count of the sample'
            )
        if major_cn_mode == 2 and wgd_count == 0:
            warnings.warn(
                'Sample was specified with WGD count 0 but major CN mode is '
                '2. Please check the WGD count of the sample'
            )

        if wgd_count == 1:
            (
                wgd_timing_distribution,
                non_overlapping_segment_ids,
                overlap_proportion,
                best_overlap_timing,
            ) = time_wgd_major_cn_2(
                sample,
                output_dir,
                timing_dict_dir,
                interval_config=interval_config,
                subclone_fraction_prior=subclone_fraction_prior,
                unordered_balanced_route_prior=unordered_balanced_route_prior,
            )
            wgd_status = True
        else:
            wgd_status = False
            wgd_timing_distribution = None
            non_overlapping_segment_ids = []
            best_overlap_timing = np.nan
            overlap_proportion = np.nan
    elif major_cn_mode == 1:
        wgd_status = False
        wgd_timing_distribution = None
        non_overlapping_segment_ids = []
        best_overlap_timing = np.nan
        overlap_proportion = np.nan
    else:
        (
            wgd_timing_distribution,
            non_overlapping_segment_ids,
            overlap_proportion,
            best_overlap_timing,
        ) = time_wgd_major_cn_2(
            sample,
            output_dir,
            timing_dict_dir,
            interval_config=interval_config,
            subclone_fraction_prior=subclone_fraction_prior,
            unordered_balanced_route_prior=unordered_balanced_route_prior,
        )
        if overlap_proportion < min_wgd_overlap:
            warnings.warn(
                'The major CN mode is 2, but the overlap proportion is '
                f'less than {min_wgd_overlap}. There are a lot of copy '
                'number 2 segments, but they are not overlapping enough '
                'to be confident in a WGD call. The sample will be '
                'treated as a non-WGD sample. Proceed with caution in '
                'downstream analysis for this sample'
            )
            wgd_timing_distribution = None
            non_overlapping_segment_ids = []
            wgd_status = False
        else:
            wgd_status = True

    wgd_info = get_wgd_info(
        wgd_timing_distribution,
        interval=interval_config.sample_wgd,
    )
    write_wgd_calling_info(
        wgd_info,
        overlap_proportion,
        best_overlap_timing,
        major_cn_mode,
        wgd_status,
        output_dir,
        sample.sample_id,
    )
    timeable_complex_segments = get_timeable_segments(
        sample,
        wgd_status=wgd_status,
        excluded_segment_ids=non_overlapping_segment_ids,
    )

    route_table_path = output_dir/f'{sample.sample_id}_route_table.tsv'
    timing_table_path = output_dir/f'{sample.sample_id}_gain_timing_table.tsv'
    process_segments(
        timeable_complex_segments,
        wgd_timing_distribution,
        output_dir,
        timing_dict_dir,
        sample.sample_id,
        wgd_status,
        plot_trees=plot_trees,
        wgd_info=wgd_info,
        interval_config=interval_config,
        subclone_fraction_prior=subclone_fraction_prior,
        unordered_balanced_route_prior=unordered_balanced_route_prior,
    )

    if os.path.exists(route_table_path):
        for apply_penalty in (False, True):
            (
                sample_posterior_table,
                sample_posterior_route_draw_table,
            ) = posteriortablegen.get_sample_posterior_tables(
                route_table_path,
                timing_table_path,
                output_dir,
                sample.sample_id,
                apply_penalty=apply_penalty,
            )
            sample_posterior_table_summary = (
                posteriortablegen.get_sample_posterior_table_summary(
                    sample_posterior_table,
                    sample_posterior_route_draw_table,
                    interval=interval_config.posterior_summary,
                )
            )
            posterior_table_summary_path = (
                output_dir
                / f'{sample.sample_id}_posterior_timing_table_summary_'
                f'penalty_{apply_penalty}.tsv'
            )
            _write_tsv(
                sample_posterior_table_summary,
                posterior_table_summary_path,
            )


def process_sample(
    sample,
    output_dir,
    plot_trees=False,
    min_wgd_overlap=0.6,
    wgd_count=None,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_fraction_prior=DEFAULT_SUBCLONE_FRACTION_PRIOR,
    unordered_balanced_route_prior=DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
    random_seed=None,
):
    if not isinstance(interval_config, TimingIntervalConfig):
        raise TypeError('interval_config must be a TimingIntervalConfig')
    subclone_fraction_prior = validate_subclone_fraction_prior(
        subclone_fraction_prior
    )
    unordered_balanced_route_prior = validate_unordered_balanced_route_prior(
        unordered_balanced_route_prior
    )
    wgd_count = _validate_wgd_count(wgd_count)
    min_wgd_overlap = _validate_min_wgd_overlap(min_wgd_overlap)
    random_seed = _validate_random_seed(random_seed)
    validate_sample_id(sample.sample_id)
    major_cn_mode = get_major_cn_mode(sample)
    if major_cn_mode not in (1, 2):
        raise ValueError(
            f"Major CN mode is {major_cn_mode} for sample {sample.sample_id}; "
            "GRITIC currently supports only modal major copy numbers 1 and 2"
        )

    output_dir = pathlib.Path(output_dir, pathlib.Path(sample.sample_id))
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(
                f'Sample output path is not a directory: {output_dir}'
            )
        if any(output_dir.iterdir()):
            raise FileExistsError(
                'Sample output directory must be absent or empty before a '
                f'GRITIC run: {output_dir}'
            )
    else:
        os.makedirs(output_dir, exist_ok=False)

    if random_seed is not None:
        np.random.seed(random_seed)
        hitandrun.seed_random(random_seed)


    timing_dict_dir= output_dir/f"{sample.sample_id}_timing_dicts/"
    os.makedirs(timing_dict_dir,exist_ok=True)

    sample_mutation_table = sample.get_mutation_table()
    sample_mutation_table_path = output_dir/f"{sample.sample_id}_mutation_table.tsv"
    _write_tsv(sample_mutation_table, sample_mutation_table_path)

    sample_subclone_table = sample.get_subclone_table()
    sample_subclone_table_path = output_dir/f"{sample.sample_id}_subclone_table.tsv"
    if sample_subclone_table is None:
        sample_subclone_table = pd.DataFrame(
            columns=dataloader.SUBCLONE_OUTPUT_COLUMNS,
        )
    _write_tsv(sample_subclone_table, sample_subclone_table_path)

    _run_sample(
        sample,
        output_dir,
        timing_dict_dir,
        plot_trees,
        min_wgd_overlap,
        wgd_count,
        interval_config,
        subclone_fraction_prior,
        unordered_balanced_route_prior,
        major_cn_mode,
    )
