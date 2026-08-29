import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

from gritic import gritictimer, timingio


class BalancedMirrorReuseTest(unittest.TestCase):
    def get_classifier_and_mirror_pair(self):
        classifier = gritictimer.RouteClassifier(
            4,
            4,
            False,
            'No_WGD',
        )
        source = next(
            route
            for route in classifier.routes.values()
            if route.mirror_route_id != route.route_id
        )
        target = classifier.routes[source.mirror_route_id]
        self.assertEqual(target.mirror_route_id, source.route_id)
        return classifier, source, target

    def get_mirror_pair(self):
        _, source, target = self.get_classifier_and_mirror_pair()
        return source, target

    @staticmethod
    def proposal_geometry(route, offset=0.0):
        n_proposals = 3
        n_nodes = len(route.route_tree.non_phased_node_order)
        mult_blocks = []
        column_offset = 0
        for block_size in (
            route.major_cn,
            route.major_cn,
            route.minor_cn,
        ):
            block = (
                np.arange(
                    n_proposals * block_size,
                    dtype=float,
                ).reshape(n_proposals, block_size)
                + 1.0
                + offset
                + column_offset
            )
            block /= np.sum(block, axis=1, keepdims=True)
            mult_blocks.append(block)
            column_offset += block_size
        mult_store = np.concatenate(mult_blocks, axis=1)
        timing_store = (
            np.arange(n_nodes * n_proposals, dtype=float).reshape(
                n_nodes,
                n_proposals,
            )
            + 1000.0
            + offset
        )
        return gritictimer.ProposalGeometry(
            mult_store=mult_store,
            timing_store=timing_store,
            wgd_timing_store=np.full(n_proposals, np.nan),
            density=np.array([0.91, 0.82]) + offset,
        )

    @staticmethod
    def expected_mirror_mult(source_mult, copy_number):
        expected = source_mult.copy()
        expected[:, copy_number:2 * copy_number] = source_mult[
            :, 2 * copy_number:3 * copy_number
        ]
        expected[:, 2 * copy_number:3 * copy_number] = source_mult[
            :, copy_number:2 * copy_number
        ]
        return expected

    def populate_fitted_route(self, route, n_subclones=1):
        n_weighted = 4
        n_nodes = len(route.route_tree.non_phased_node_order)
        n_mult_columns = 3 * route.major_cn + n_subclones
        route.log_evidence = -3.5
        route.node_timing = np.arange(
            n_nodes * n_weighted,
            dtype=float,
        ).reshape(n_nodes, n_weighted)
        route.wgd_timing_store = np.full(n_weighted, np.nan)
        route.mult_store = np.arange(
            n_weighted * n_mult_columns,
            dtype=float,
        ).reshape(n_weighted, n_mult_columns)
        route.n_events_store = {
            'N_Events': np.array([2, 2, 2]),
            'Pre_WGD_Losses': np.array([np.nan, np.nan, np.nan]),
            'Post_WGD_Losses': np.array([np.nan, np.nan, np.nan]),
        }
        route.density = 0.93
        route.density_high = 0.81

    def test_mirror_geometry_is_transformed_in_memory(self):
        source, target = self.get_mirror_pair()
        source_geometry = self.proposal_geometry(source)

        observed = target.transform_mirror_geometry(
            source,
            source_geometry,
        )

        np.testing.assert_array_equal(
            observed.mult_store,
            self.expected_mirror_mult(
                source_geometry.mult_store,
                source.major_cn,
            ),
        )
        node_map = gritictimer._get_mirror_node_map(
            source.route_tree.main_tree,
            target.route_tree.main_tree,
        )
        source_positions = {
            node: position
            for position, node in enumerate(
                source.route_tree.non_phased_node_order
            )
        }
        target_positions = {
            node: position
            for position, node in enumerate(
                target.route_tree.non_phased_node_order
            )
        }
        for source_node, target_node in node_map.items():
            np.testing.assert_array_equal(
                observed.timing_store[target_positions[target_node]],
                source_geometry.timing_store[
                    source_positions[source_node]
                ],
            )
        np.testing.assert_array_equal(
            observed.wgd_timing_store,
            source_geometry.wgd_timing_store,
        )
        np.testing.assert_array_equal(
            observed.density,
            source_geometry.density,
        )
        self.assertFalse(
            np.shares_memory(
                observed.mult_store,
                source_geometry.mult_store,
            )
        )

    def test_fit_route_classifiers_shares_one_geometry_per_route_unit(self):
        first = gritictimer.RouteClassifier(4, 4, False, 'No_WGD')
        route_trees = {
            route_id: route.route_tree
            for route_id, route in first.routes.items()
        }
        second = gritictimer.RouteClassifier(
            4,
            4,
            False,
            'No_WGD',
            route_trees=route_trees,
        )
        first_probabilities = SimpleNamespace(
            use_major=True,
            use_minor=False,
        )
        second_probabilities = SimpleNamespace(
            use_major=True,
            use_minor=False,
        )
        sampling_calls = []
        fitting_calls = []

        def sample_geometry(route, _wgd_distribution):
            sampling_calls.append(route.route_id)
            return self.proposal_geometry(
                route,
                offset=float(len(sampling_calls) * 10000),
            )

        def fit_route(
            route,
            mult_probabilities,
            _subclone_table,
            _wgd_distribution,
            _fraction_prior,
            *,
            mirror_route=None,
            proposal_geometry=None,
            clone_share=None,
            shared_geometry_time=0.0,
        ):
            del shared_geometry_time
            fitting_calls.append((
                route.route_id,
                mult_probabilities,
                proposal_geometry,
                mirror_route,
            ))
            route.log_evidence = 0.0
            probability_id = id(mult_probabilities)
            if mirror_route is None:
                generated_share = object()
                source_clone_shares[probability_id] = generated_share
                return generated_share
            test_case.assertIs(
                clone_share,
                source_clone_shares[probability_id],
            )
            return clone_share

        test_case = self
        source_clone_shares = {}
        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            autospec=True,
            side_effect=sample_geometry,
        ), mock.patch.object(
            gritictimer.Route,
            'run_sampling',
            autospec=True,
            side_effect=fit_route,
        ):
            gritictimer.fit_route_classifiers(
                [
                    (first, first_probabilities, None),
                    (second, second_probabilities, None),
                ],
                None,
            )

        route_units = {
            frozenset((route.route_id, route.mirror_route_id))
            if route.mirror_route_id is not None
            else frozenset((route.route_id,))
            for route in first.routes.values()
        }
        self.assertEqual(len(sampling_calls), len(route_units))
        self.assertEqual(len(fitting_calls), 2 * len(first.routes))
        for route_id in first.routes:
            route_fits = [
                call for call in fitting_calls if call[0] == route_id
            ]
            self.assertEqual(len(route_fits), 2)
            self.assertIs(route_fits[0][2], route_fits[1][2])
            self.assertIsNotNone(route_fits[0][2])
        self.assertAlmostEqual(sum(first.route_probabilities.values()), 1.0)
        self.assertAlmostEqual(sum(second.route_probabilities.values()), 1.0)

    def test_run_sampling_keeps_scalar_evidence_not_full_likelihood(self):
        classifier = gritictimer.RouteClassifier(3, 1, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        geometry = self.proposal_geometry(route)
        likelihoods = np.array([-4.0, -1.0, -2.0])
        test_case = self

        class ProbabilityStore:
            use_major = False
            use_minor = False

            @staticmethod
            def evaluate_likelihood_array(mult_store):
                test_case.assertEqual(mult_store.shape[0], likelihoods.size)
                return likelihoods.copy()

        event_store = {
            'N_Events': np.array([2]),
            'Pre_WGD_Losses': np.array([np.nan]),
            'Post_WGD_Losses': np.array([np.nan]),
        }
        with mock.patch.object(
            route,
            'get_weighted_arrays',
            return_value=(
                geometry.timing_store,
                geometry.wgd_timing_store,
                geometry.mult_store,
            ),
        ), mock.patch.object(
            route,
            'get_n_events_estimate',
            return_value=event_store,
        ):
            route.run_sampling(
                ProbabilityStore(),
                None,
                None,
                proposal_geometry=geometry,
            )

        expected = np.log(np.mean(np.exp(likelihoods)))
        self.assertAlmostEqual(route.log_evidence, expected)
        self.assertFalse(hasattr(route, 'll_store'))
        np.testing.assert_array_equal(
            route.node_timing,
            geometry.timing_store,
        )
        np.testing.assert_array_equal(
            route.wgd_timing_store,
            geometry.wgd_timing_store,
        )
        np.testing.assert_array_equal(
            route.mult_store,
            geometry.mult_store,
        )

    def test_unphased_fitted_mirror_reuses_scalar_and_swaps_outputs(self):
        source, target = self.get_mirror_pair()
        n_subclones = 1
        self.populate_fitted_route(source, n_subclones=n_subclones)
        mult_probabilities = SimpleNamespace(
            use_major=False,
            use_minor=False,
        )
        subclone_table = SimpleNamespace(index=range(n_subclones))

        with mock.patch.object(
            target,
            'run_geometry_sampling',
            side_effect=AssertionError(
                'fitted unphased mirror should not sample geometry'
            ),
        ):
            target.run_sampling(
                mult_probabilities,
                subclone_table,
                None,
                mirror_route=source,
            )

        self.assertEqual(target.log_evidence, source.log_evidence)
        self.assertFalse(hasattr(target, 'll_store'))
        np.testing.assert_array_equal(
            target.mult_store,
            self.expected_mirror_mult(source.mult_store, source.major_cn),
        )
        np.testing.assert_array_equal(
            target.wgd_timing_store,
            source.wgd_timing_store,
        )

        node_map = gritictimer._get_mirror_node_map(
            source.route_tree.main_tree,
            target.route_tree.main_tree,
        )
        source_positions = {
            node: position
            for position, node in enumerate(
                source.route_tree.non_phased_node_order
            )
        }
        target_positions = {
            node: position
            for position, node in enumerate(
                target.route_tree.non_phased_node_order
            )
        }
        for source_node, target_node in node_map.items():
            np.testing.assert_array_equal(
                target.node_timing[target_positions[target_node]],
                source.node_timing[source_positions[source_node]],
            )

    def test_unphased_mirrors_have_equal_posteriors_and_archived_samples(self):
        classifier, source, target = self.get_classifier_and_mirror_pair()

        class UnphasedProbabilityStore:
            use_major = False
            use_minor = False

            @staticmethod
            def evaluate_likelihood_array(_mult_store):
                return np.array([-3.0, -1.0, -2.0])

        def sample_geometry(route, _wgd_distribution):
            geometry = self.proposal_geometry(route)
            geometry.timing_store[:] = np.linspace(
                0.1,
                0.9,
                geometry.timing_store.size,
            ).reshape(geometry.timing_store.shape)
            return geometry

        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            autospec=True,
            side_effect=sample_geometry,
        ), mock.patch.object(
            gritictimer.Route,
            'get_n_events_estimate',
            autospec=True,
            return_value={
                'N_Events': np.array([2]),
                'Pre_WGD_Losses': np.array([np.nan]),
                'Post_WGD_Losses': np.array([np.nan]),
            },
        ):
            np.random.seed(7)
            classifier.fit_routes(
                UnphasedProbabilityStore(),
                None,
                None,
            )

        self.assertEqual(
            classifier.route_probabilities[source.route_id],
            classifier.route_probabilities[target.route_id],
        )
        self.assertEqual(target.log_evidence, source.log_evidence)
        self.assertFalse(hasattr(source, 'll_store'))

        timing_dict = classifier.get_timing_dict({
            route.short_id: classifier.route_probabilities[route.route_id]
            for route in classifier.routes.values()
        })
        with tempfile.TemporaryDirectory() as directory:
            archive_path, manifest_path = timingio.write_timing_archive(
                timing_dict,
                directory,
                'balanced-segment',
            )
            archived = timingio.load_timing_archive(
                archive_path,
                manifest_path,
            )

        source_archive = archived[source.short_id]
        target_archive = archived[target.short_id]
        self.assertEqual(
            set(source_archive),
            timingio.ROUTE_PARTICLE_ARCHIVE_FIELDS,
        )
        self.assertEqual(
            set(target_archive),
            timingio.ROUTE_PARTICLE_ARCHIVE_FIELDS,
        )
        np.testing.assert_array_equal(
            target_archive['Mult'],
            self.expected_mirror_mult(
                source_archive['Mult'],
                source.major_cn,
            ),
        )
        np.testing.assert_array_equal(
            target_archive['WGD_Timing'],
            source_archive['WGD_Timing'],
        )
        node_map = gritictimer._get_mirror_node_map(
            source.route_tree.main_tree,
            target.route_tree.main_tree,
        )
        source_columns = dict(zip(
            source_archive['Timing_Node_ID'].tolist(),
            source_archive['Timing'].T,
        ))
        target_columns = dict(zip(
            target_archive['Timing_Node_ID'].tolist(),
            target_archive['Timing'].T,
        ))
        for source_node in source.route_tree.timeable_nodes:
            np.testing.assert_array_equal(
                target_columns[node_map[source_node]],
                source_columns[source_node],
            )

    def test_repeated_fit_replaces_scalar_evidence_and_stale_probabilities(self):
        classifier, source, target = self.get_classifier_and_mirror_pair()

        class ChangingUnphasedProbabilityStore:
            use_major = False
            use_minor = False

            def __init__(self, likelihoods):
                self.likelihoods = np.asarray(likelihoods, dtype=float)

            def evaluate_likelihood_array(self, _mult_store):
                return self.likelihoods.copy()

        def sample_geometry(route, _wgd_distribution):
            return self.proposal_geometry(route)

        event_store = {
            'N_Events': np.array([2]),
            'Pre_WGD_Losses': np.array([np.nan]),
            'Post_WGD_Losses': np.array([np.nan]),
        }
        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            autospec=True,
            side_effect=sample_geometry,
        ), mock.patch.object(
            gritictimer.Route,
            'get_n_events_estimate',
            autospec=True,
            return_value=event_store,
        ):
            classifier.fit_routes(
                ChangingUnphasedProbabilityStore([-3.0, -1.0, -2.0]),
                None,
                None,
            )
            first_evidence = source.log_evidence
            classifier.route_probabilities['stale-route'] = 1.0

            classifier.fit_routes(
                ChangingUnphasedProbabilityStore([-12.0, -8.0, -10.0]),
                None,
                None,
            )

        expected = np.log(np.mean(np.exp([-12.0, -8.0, -10.0])))
        self.assertNotIn('stale-route', classifier.route_probabilities)
        self.assertNotEqual(source.log_evidence, first_evidence)
        self.assertAlmostEqual(source.log_evidence, expected)
        self.assertEqual(target.log_evidence, source.log_evidence)

    def test_phased_mirror_reuses_geometry_but_recomputes_evidence(self):
        source, target = self.get_mirror_pair()
        source_geometry = self.proposal_geometry(source)
        target_geometry = target.transform_mirror_geometry(
            source,
            source_geometry,
        )
        self.populate_fitted_route(source, n_subclones=0)

        class PhasedProbabilityStore:
            use_major = True
            use_minor = False

            def __init__(self):
                self.observed_mult_store = None

            def evaluate_likelihood_array(self, mult_store):
                self.observed_mult_store = mult_store.copy()
                return np.array([-4.0, -1.0, -2.0])

        mult_probabilities = PhasedProbabilityStore()
        with mock.patch.object(
            target,
            'run_geometry_sampling',
            side_effect=AssertionError('mirror geometry should be supplied'),
        ), mock.patch.object(
            target,
            'get_weighted_arrays',
            wraps=target.get_weighted_arrays,
        ) as get_weighted_arrays, mock.patch.object(
            target,
            'get_n_events_estimate',
            return_value={
                'N_Events': np.array([2]),
                'Pre_WGD_Losses': np.array([np.nan]),
                'Post_WGD_Losses': np.array([np.nan]),
            },
        ):
            target.run_sampling(
                mult_probabilities,
                None,
                None,
                mirror_route=source,
                proposal_geometry=target_geometry,
            )

        expected_mult = self.expected_mirror_mult(
            source_geometry.mult_store,
            source.major_cn,
        )
        np.testing.assert_array_equal(
            mult_probabilities.observed_mult_store,
            expected_mult,
        )
        expected_evidence = np.log(np.mean(np.exp([-4.0, -1.0, -2.0])))
        self.assertAlmostEqual(target.log_evidence, expected_evidence)
        self.assertNotEqual(target.log_evidence, source.log_evidence)
        self.assertIsNone(target.unphased_mirror_source)
        self.assertFalse(hasattr(target, 'll_store'))
        get_weighted_arrays.assert_called_once()
        np.testing.assert_allclose(
            get_weighted_arrays.call_args.args[3],
            np.exp(np.array([-3.0, 0.0, -1.0])),
        )

    def test_phased_mirror_pair_reuses_one_clone_share_draw(self):
        classifier, first, second = self.get_classifier_and_mirror_pair()
        source_route_id = min(first.route_id, second.route_id)
        source = classifier.routes[source_route_id]
        target = classifier.routes[source.mirror_route_id]
        classifier.routes = {
            source.route_id: source,
            target.route_id: target,
        }
        classifier.subclone_fraction_prior = 'supplied'
        source_geometry = self.proposal_geometry(source)
        clone_share = np.array([
            [0.8, 0.2],
            [0.7, 0.3],
            [0.6, 0.4],
        ])
        likelihood_inputs = []

        class PhasedProbabilityStore:
            use_major = True
            use_minor = False

            @staticmethod
            def evaluate_likelihood_array(mult_store):
                likelihood_inputs.append(mult_store.copy())
                return np.array([-3.0, -1.0, -2.0])

        subclone_table = pd.DataFrame({
            'Subclone_Fraction': [0.2],
        })
        with mock.patch.object(
            source,
            'run_geometry_sampling',
            return_value=source_geometry,
        ), mock.patch.object(
            gritictimer.Route,
            'simulate_clone_share',
            autospec=True,
            return_value=clone_share,
        ) as simulate_clone_share:
            gritictimer.fit_route_classifiers(
                [(classifier, PhasedProbabilityStore(), subclone_table)],
                None,
            )

        simulate_clone_share.assert_called_once()
        self.assertEqual(len(likelihood_inputs), 2)
        source_mult, target_mult = likelihood_inputs
        np.testing.assert_allclose(
            source_mult[:, :-1],
            source_geometry.mult_store * clone_share[:, :1],
        )
        np.testing.assert_array_equal(
            source_mult[:, -1],
            clone_share[:, -1],
        )
        np.testing.assert_array_equal(
            target_mult[:, -1],
            source_mult[:, -1],
        )
        np.testing.assert_allclose(
            target_mult[:, :-1],
            self.expected_mirror_mult(
                source_mult[:, :-1],
                source.major_cn,
            ),
        )

    def test_wgd_decorated_mirror_node_map_preserves_wgd_structure(self):
        classifier = gritictimer.RouteClassifier(
            3,
            3,
            True,
            'Default',
        )
        source = next(
            route
            for route in classifier.routes.values()
            if (
                route.mirror_route_id != route.route_id
                and route.route_tree.wgd_nodes
            )
        )
        target = classifier.routes[source.mirror_route_id]
        node_map = gritictimer._get_mirror_node_map(
            source.route_tree.main_tree,
            target.route_tree.main_tree,
        )

        source_tree = source.route_tree.main_tree
        target_tree = target.route_tree.main_tree
        for source_node, target_node in node_map.items():
            self.assertEqual(
                source_tree.nodes[source_node]['WGD_Symbol'],
                target_tree.nodes[target_node]['WGD_Symbol'],
            )
            self.assertEqual(
                gritictimer._MIRROR_ALLELE[
                    source_tree.nodes[source_node][
                        gritictimer.treetools.ALLELE_ATTRIBUTE
                    ]
                ],
                target_tree.nodes[target_node][
                    gritictimer.treetools.ALLELE_ATTRIBUTE
                ],
            )


if __name__ == '__main__':
    unittest.main()
