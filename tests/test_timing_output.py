import unittest
from unittest import mock

import numpy as np

from gritic import gritictimer, timingio
from gritic.tableschemas import ROUTE_PARTICLES_REPRESENTATION


class TimingOutputTest(unittest.TestCase):
    @staticmethod
    def fitted_classifier(*, wgd_status, n_particles=2):
        classifier = gritictimer.RouteClassifier(
            2,
            0,
            wgd_status,
            'WGD' if wgd_status else 'No_WGD',
        )
        route = next(iter(classifier.routes.values()))
        route.node_timing = np.full(
            (len(route.route_tree.non_phased_node_order), n_particles),
            0.4,
            dtype=np.float64,
        )
        route.wgd_timing_store = np.full(
            n_particles,
            0.5 if wgd_status else np.nan,
            dtype=np.float64,
        )
        route.mult_store = np.tile(
            np.array([0.5, 0.5, 0.5, 0.5], dtype=np.float64),
            (n_particles, 1),
        )
        route.n_subclones = 0
        classifier.route_probabilities = {route.route_id: 1.0}
        classifier._get_output_routes = mock.Mock(
            return_value=[(route, 1.0)]
        )
        return classifier, route

    def test_serializes_only_aligned_conditional_particles(self):
        classifier = gritictimer.RouteClassifier(2, 0, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        node = route.route_tree.timeable_nodes[0]
        route.node_timing = np.vstack([
            np.array([0.15, 0.25, 0.35]),
            np.ones(3),
            np.ones(3),
        ])
        route.wgd_timing_store = np.full(3, np.nan)
        route.mult_store = np.array([
            [0.1, 0.9, 0.1, 0.9],
            [0.2, 0.8, 0.2, 0.8],
            [0.3, 0.7, 0.3, 0.7],
        ])
        route.n_subclones = 0
        route.log_evidence = -101.0
        classifier.route_probabilities = {route.route_id: 1.0}
        classifier._get_output_routes = mock.Mock(
            return_value=[(route, 1.0)]
        )
        with mock.patch.object(
            gritictimer.np.random,
            'choice',
            side_effect=AssertionError(
                'Timing serialization must not resample fitted particles'
            ),
        ):
            output = classifier.get_timing_dict(
                {route.short_id: 1.0}
            )[route.short_id]

        self.assertEqual(set(output), timingio.ROUTE_PARTICLE_ARCHIVE_FIELDS)
        self.assertFalse(hasattr(route, 'll_store'))
        np.testing.assert_array_equal(
            output['WGD_Timing'],
            route.wgd_timing_store,
        )
        np.testing.assert_array_equal(
            output['Timing'][:, 0],
            route.node_timing[0],
        )
        np.testing.assert_array_equal(output['Timing_Node_ID'], [node])
        np.testing.assert_array_equal(
            output['Mult'],
            route.mult_store,
        )
        self.assertEqual(
            output['Mult'].shape[0],
            output['WGD_Timing'].shape[0],
        )
        self.assertEqual(
            output['Mult'].shape[0],
            output['Timing'].shape[0],
        )
        self.assertFalse(np.shares_memory(
            output['WGD_Timing'],
            route.wgd_timing_store,
        ))
        self.assertFalse(np.shares_memory(
            output['Timing'],
            route.node_timing,
        ))
        self.assertFalse(np.shares_memory(output['Mult'], route.mult_store))
        np.testing.assert_array_equal(output['Interval_Start_Source'], [0, 3, 3])
        np.testing.assert_array_equal(output['Interval_End_Source'], [3, 1, 1])
        np.testing.assert_array_equal(output['Interval_Multiplicity'], [2, 1, 1])
        np.testing.assert_array_equal(output['Interval_Phasing'], [3, 3, 3])
        np.testing.assert_array_equal(output['State_Column_Offsets'], [0, 2, 4, 4])
        np.testing.assert_array_equal(output['State_Columns'], [0, 1, 2, 3])
        for key in (
            'Probability',
            'Penalized_Probability',
            'Timing',
            'WGD_Timing',
            'Mult',
        ):
            self.assertEqual(output[key].dtype, np.dtype(np.float64))
        for key in (
            'Timing_Node_ID',
            'Interval_Start_Source',
            'Interval_End_Source',
            'Interval_Multiplicity',
            'State_Column_Offsets',
            'State_Columns',
            'Target_Major_CN',
            'Target_Minor_CN',
            'N_Subclones',
            'Model_Major_CN',
            'Model_Minor_CN',
        ):
            self.assertEqual(output[key].dtype, np.dtype(np.int64))

    def test_pooled_two_plus_two_explicitly_aliases_every_phase_to_major(self):
        classifier = gritictimer.RouteClassifier(2, 0, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        route.node_timing = np.vstack([
            np.array([0.2, 0.4]),
            np.ones(2),
            np.ones(2),
        ])
        route.wgd_timing_store = np.full(2, np.nan)
        route.mult_store = np.array([
            [0.20, 0.70, 0.20, 0.70, 0.10],
            [0.45, 0.45, 0.45, 0.45, 0.10],
        ])
        route.n_subclones = 1
        classifier.route_probabilities = {route.route_id: 1.0}

        output = classifier.get_timing_dict(
            {route.short_id: 1.0},
            target_major_cn=2,
            target_minor_cn=2,
            target_wgd_status=True,
            archive_kind=timingio.ROUTE_PARTICLE_ARCHIVE_KIND_POOLED_WGD,
        )[route.short_id]

        np.testing.assert_array_equal(output['Interval_Phasing'], [7, 7, 7])
        np.testing.assert_array_equal(output['State_Column_Offsets'], [0, 3, 6, 9])
        np.testing.assert_array_equal(
            output['State_Columns'],
            [2, 3, 4, 2, 3, 4, 2, 3, 4],
        )
        np.testing.assert_array_equal(output['Target_Major_CN'], [2])
        np.testing.assert_array_equal(output['Target_Minor_CN'], [2])
        np.testing.assert_array_equal(output['Target_WGD_Status'], [1])
        np.testing.assert_array_equal(output['Model_Major_CN'], [2])
        np.testing.assert_array_equal(output['Model_Minor_CN'], [0])
        np.testing.assert_array_equal(output['Model_WGD_Status'], [0])
        np.testing.assert_array_equal(output['Archive_Kind'], [1])

    def test_mult_serialization_clips_only_numerical_negative_noise(self):
        classifier = gritictimer.RouteClassifier(2, 0, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        route.node_timing = np.array([
            [0.4],
            [1.0],
            [1.0],
        ])
        route.wgd_timing_store = np.array([np.nan])
        route.n_subclones = 0
        classifier.route_probabilities = {route.route_id: 1.0}

        numerical_negative = -5e-13
        route.mult_store = np.array([[
            numerical_negative,
            1.0 - numerical_negative,
            numerical_negative,
            1.0 - numerical_negative,
        ]])
        output = classifier.get_timing_dict(
            {route.short_id: 1.0}
        )[route.short_id]

        np.testing.assert_array_equal(
            output['Mult'],
            [[0.0, 1.0 - numerical_negative, 0.0, 1.0 - numerical_negative]],
        )
        self.assertEqual(route.mult_store[0, 0], numerical_negative)

        route.mult_store[0, 0] = -2e-12
        with self.assertRaisesRegex(ValueError, 'material negative mass'):
            classifier.get_timing_dict({route.short_id: 1.0})

    def test_timing_serialization_uses_shared_boundary_tolerance(self):
        classifier, route = self.fitted_classifier(wgd_status=False)
        tolerance = timingio.ROUTE_PARTICLE_TIMING_TOLERANCE
        node = route.route_tree.timeable_nodes[0]
        node_index = route.route_tree.non_phased_node_order.index(node)
        route.node_timing[node_index] = np.array([
            -0.5 * tolerance,
            np.nextafter(1.0, np.inf),
        ])
        original_timing = route.node_timing.copy()

        output = classifier.get_timing_dict(
            {route.short_id: 1.0}
        )[route.short_id]

        np.testing.assert_array_equal(output['Timing'][:, 0], [0.0, 1.0])
        np.testing.assert_array_equal(route.node_timing, original_timing)

        for invalid_value in (
            -2.0 * tolerance,
            1.0 + 2.0 * tolerance,
            np.nan,
            np.inf,
            -np.inf,
        ):
            with self.subTest(invalid_value=invalid_value):
                route.node_timing[node_index] = [invalid_value, 0.4]
                with self.assertRaisesRegex(ValueError, 'finite Timing'):
                    classifier.get_timing_dict({route.short_id: 1.0})

    def test_wgd_timing_serialization_uses_shared_boundary_tolerance(self):
        classifier, route = self.fitted_classifier(wgd_status=True)
        tolerance = timingio.ROUTE_PARTICLE_TIMING_TOLERANCE
        route.wgd_timing_store = np.array([
            -0.5 * tolerance,
            1.0 + 0.5 * tolerance,
        ])
        original_wgd_timing = route.wgd_timing_store.copy()

        output = classifier.get_timing_dict(
            {route.short_id: 1.0}
        )[route.short_id]

        np.testing.assert_array_equal(output['WGD_Timing'], [0.0, 1.0])
        np.testing.assert_array_equal(
            route.wgd_timing_store,
            original_wgd_timing,
        )

        for invalid_value in (
            -2.0 * tolerance,
            1.0 + 2.0 * tolerance,
            np.nan,
            np.inf,
            -np.inf,
        ):
            with self.subTest(invalid_value=invalid_value):
                route.wgd_timing_store = np.array([invalid_value, 0.5])
                with self.assertRaisesRegex(ValueError, 'finite WGD_Timing'):
                    classifier.get_timing_dict({route.short_id: 1.0})


if __name__ == '__main__':
    unittest.main()
