import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from gritic import gritictimer


class TimingOutputTest(unittest.TestCase):
    def test_serializes_only_aligned_conditional_particles(self):
        node = 7
        route = gritictimer.Route.__new__(gritictimer.Route)
        route.route_tree = SimpleNamespace(
            timeable_nodes=[node],
            non_phased_node_order=[node],
        )

        route.short_id = 'route-one'
        route.node_timing = np.array([[0.15, 0.25, 0.35]])
        route.wgd_timing_store = np.array([0.45, 0.55, 0.65])
        route.mult_store = np.array([
            [0.1, 0.9],
            [0.2, 0.8],
            [0.3, 0.7],
        ])
        route.log_evidence = -101.0
        route.unphased_mirror_source = SimpleNamespace(
            short_id='source-route'
        )

        classifier = gritictimer.RouteClassifier.__new__(
            gritictimer.RouteClassifier
        )
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
            output = classifier.get_timing_dict()[route.short_id]

        self.assertEqual(set(output), {'Timing', 'Mult'})
        self.assertFalse(hasattr(route, 'll_store'))
        np.testing.assert_array_equal(
            output['Timing']['WGD'],
            route.wgd_timing_store,
        )
        np.testing.assert_array_equal(
            output['Timing'][node],
            route.node_timing[0],
        )
        np.testing.assert_array_equal(
            output['Mult'],
            route.mult_store,
        )
        self.assertEqual(
            output['Mult'].shape[0],
            output['Timing']['WGD'].shape[0],
        )
        self.assertEqual(
            output['Mult'].shape[0],
            output['Timing'][node].shape[0],
        )
        self.assertFalse(np.shares_memory(
            output['Timing']['WGD'],
            route.wgd_timing_store,
        ))
        self.assertFalse(np.shares_memory(
            output['Timing'][node],
            route.node_timing,
        ))
        self.assertFalse(np.shares_memory(output['Mult'], route.mult_store))


if __name__ == '__main__':
    unittest.main()
