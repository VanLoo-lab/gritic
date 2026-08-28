import unittest
from unittest import mock

import numpy as np
import pandas as pd

from gritic import gritictimer
from gritic.sampletools import MultProbabilityStore
from gritic.tableschemas import (
    TIMING_REPRESENTATION_COLUMN,
    UNIFORM_NO_GAIN_REPRESENTATION,
)


def make_multiplicity_probabilities(*, minor_cn, n_subclones):
    """Return a small valid state-likelihood store for a 1+N segment."""
    probabilities = np.array([[0.25, 0.75]])
    if n_subclones == 0:
        probabilities = np.ones((1, 1))

    arrays = {
        'Non_Phased': probabilities,
        'Major': None,
        'Minor': None,
        'All': probabilities,
    }
    corrections = {
        'Non_Phased': np.ones(probabilities.shape[1]),
        'Major': None,
        'Minor': None,
        'All': np.ones(probabilities.shape[1]),
    }
    return MultProbabilityStore(
        arrays,
        corrections,
        major_cn=1,
        minor_cn=minor_cn,
        n_subclones=n_subclones,
    )


class UniformNoGainPredicateTest(unittest.TestCase):
    def test_one_copy_states_are_eligible_with_or_without_a_wgd_call(self):
        for wgd_status in (False, True):
            for minor_cn in (0, 1):
                with self.subTest(
                    wgd_status=wgd_status,
                    minor_cn=minor_cn,
                ):
                    classifier = gritictimer.RouteClassifier(
                        1,
                        minor_cn,
                        wgd_status,
                        'Default',
                    )
                    self.assertTrue(classifier.routes)
                    self.assertTrue(all(
                        gritictimer.is_uniform_no_gain_route_tree(
                            route.route_tree
                        )
                        for route in classifier.routes.values()
                    ))

    def test_wgd_node_is_not_mistaken_for_uniform_no_gain_geometry(self):
        classifier = gritictimer.RouteClassifier(
            2,
            0,
            True,
            'Default',
        )
        wgd_routes = [
            route
            for route in classifier.routes.values()
            if route.route_tree.wgd_nodes
        ]
        self.assertTrue(wgd_routes)
        self.assertTrue(all(
            not route.route_tree.timeable_nodes
            for route in wgd_routes
        ))
        self.assertTrue(all(
            not gritictimer.is_uniform_no_gain_route_tree(route.route_tree)
            for route in wgd_routes
        ))


class UniformNoGainFitTest(unittest.TestCase):
    def test_direct_fit_preserves_constant_geometry_state_order_and_weights(self):
        proposals = np.array([
            [0.8, 0.2],
            [0.5, 0.5],
            [0.3, 0.7],
        ])
        log_likelihoods = np.array([-2.0, -1.0, 0.0])
        sampled_indices = np.resize(
            np.array([2, 0, 1]),
            gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT,
        )

        class RecordingProbabilities:
            def evaluate_likelihood_array(self, full_states):
                self.full_states = full_states.copy()
                return log_likelihoods.copy()

        for minor_cn in (0, 1):
            with self.subTest(minor_cn=minor_cn):
                classifier = gritictimer.RouteClassifier(
                    1,
                    minor_cn,
                    False,
                    'Default',
                )
                route = next(iter(classifier.routes.values()))
                probabilities = RecordingProbabilities()
                expected_states = np.concatenate(
                    [
                        np.repeat(
                            proposals[:, :1],
                            2 + minor_cn,
                            axis=1,
                        ),
                        proposals[:, 1:],
                    ],
                    axis=1,
                )

                with mock.patch.object(
                    route,
                    'simulate_clone_share',
                    return_value=proposals.copy(),
                ) as simulate_clone_share, mock.patch.object(
                    route,
                    'get_weighted_sample_indices',
                    return_value=sampled_indices,
                ) as get_weighted_sample_indices:
                    route.run_uniform_no_gain_sampling(
                        probabilities,
                        alpha=np.array([1.8, 1.2]),
                        n_subclones=1,
                    )

                simulate_clone_share.assert_called_once_with(
                    mock.ANY,
                    gritictimer.NO_GAIN_CLONE_SHARE_PROPOSAL_COUNT,
                )
                np.testing.assert_array_equal(
                    probabilities.full_states,
                    expected_states,
                )
                np.testing.assert_allclose(
                    route.log_evidence,
                    np.log(np.mean(np.exp(log_likelihoods))),
                )
                np.testing.assert_allclose(
                    get_weighted_sample_indices.call_args.args[0],
                    np.exp(log_likelihoods - np.max(log_likelihoods)),
                )
                self.assertEqual(
                    get_weighted_sample_indices.call_args.args[1],
                    proposals.shape[0],
                )
                np.testing.assert_array_equal(
                    route.clone_share_store,
                    proposals[sampled_indices],
                )

    def test_clonal_fit_skips_geometry_and_has_no_particle_payload(self):
        classifier = gritictimer.RouteClassifier(
            1,
            1,
            False,
            'Default',
        )
        probabilities = make_multiplicity_probabilities(
            minor_cn=1,
            n_subclones=0,
        )

        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            side_effect=AssertionError(
                'uniform no-gain fitting must not sample route geometry'
            ),
        ) as geometry_sampling:
            classifier.fit_routes(probabilities, None, None)

        geometry_sampling.assert_not_called()
        self.assertEqual(
            classifier.timing_representation,
            UNIFORM_NO_GAIN_REPRESENTATION,
        )
        route = next(iter(classifier.routes.values()))
        self.assertIsNone(route.clone_share_store)
        self.assertEqual(route.density, 1.0)
        self.assertEqual(route.density_high, 1.0)
        self.assertIsNone(classifier.get_timing_dict())
        route_table = classifier.get_route_table()
        self.assertEqual(
            route_table[TIMING_REPRESENTATION_COLUMN].tolist(),
            [UNIFORM_NO_GAIN_REPRESENTATION],
        )

    def test_subclonal_fit_stores_only_likelihood_resampled_clone_shares(self):
        classifier = gritictimer.RouteClassifier(
            1,
            0,
            False,
            'Default',
            subclone_fraction_prior='supplied',
        )
        probabilities = make_multiplicity_probabilities(
            minor_cn=0,
            n_subclones=1,
        )
        subclone_table = pd.DataFrame({
            'Subclone_Fraction': [0.2],
        })
        route = next(iter(classifier.routes.values()))

        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            side_effect=AssertionError(
                'uniform no-gain fitting must not sample route geometry'
            ),
        ) as geometry_sampling, mock.patch.object(
            route,
            'simulate_clone_share',
            wraps=route.simulate_clone_share,
        ) as simulate_clone_share:
            classifier.fit_routes(
                probabilities,
                subclone_table,
                None,
            )

        geometry_sampling.assert_not_called()
        self.assertEqual(
            simulate_clone_share.call_args.args[1],
            gritictimer.NO_GAIN_CLONE_SHARE_PROPOSAL_COUNT,
        )
        self.assertEqual(
            gritictimer.NO_GAIN_CLONE_SHARE_PROPOSAL_COUNT,
            50_500,
        )
        clone_shares = route.clone_share_store
        self.assertEqual(
            clone_shares.shape,
            (gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT, 2),
        )
        self.assertTrue(np.isfinite(clone_shares).all())
        self.assertTrue((clone_shares >= 0).all())
        np.testing.assert_allclose(clone_shares.sum(axis=1), 1.0)

        timing_dict = classifier.get_timing_dict()
        self.assertEqual(set(timing_dict), {route.short_id})
        self.assertEqual(
            set(timing_dict[route.short_id]),
            {'Clone_Share'},
        )
        np.testing.assert_array_equal(
            timing_dict[route.short_id]['Clone_Share'],
            clone_shares,
        )
        self.assertFalse(np.shares_memory(
            timing_dict[route.short_id]['Clone_Share'],
            clone_shares,
        ))


if __name__ == '__main__':
    unittest.main()
