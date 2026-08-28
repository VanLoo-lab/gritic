import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from gritic import gritictimer


class TimingOutputTest(unittest.TestCase):
    def test_top_level_omits_ll_while_raw_samples_remain_aligned(self):
        node = 7
        route = gritictimer.Route.__new__(gritictimer.Route)
        route.route_tree = SimpleNamespace(
            timeable_nodes=[node],
            non_phased_node_order=[node],
        )

        source_mult = np.array([
            [1.0, 11.0],
            [2.0, 12.0],
            [3.0, 13.0],
        ])
        source_timing = np.array([[0.1, 0.2, 0.3]])
        source_wgd = np.array([0.4, 0.5, 0.6])
        source_ll = -10.0 * source_mult[:, 0]
        raw_indexes = np.array([2, 0, 1])
        with mock.patch.object(
            gritictimer.np.random,
            'randint',
            return_value=raw_indexes,
        ):
            route.raw_samples = route.get_raw_samples_store(
                source_mult,
                source_timing,
                source_wgd,
                source_ll,
                n_samples=raw_indexes.size,
            )

        np.testing.assert_array_equal(
            route.raw_samples['LL'],
            -10.0 * route.raw_samples['Mult'][:, 0],
        )
        np.testing.assert_array_equal(
            route.raw_samples['Timing'][node],
            source_timing[0, raw_indexes],
        )
        np.testing.assert_array_equal(
            route.raw_samples['WGD_Timing'],
            source_wgd[raw_indexes],
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
            'randint',
            side_effect=AssertionError(
                'Timing serialization must not resample fitted particles'
            ),
        ):
            output = classifier.get_timing_dict()[route.short_id]

        self.assertNotIn('LL', output)
        self.assertFalse(hasattr(route, 'll_store'))
        self.assertIs(output['Raw_Samples'], route.raw_samples)
        self.assertIn('LL', output['Raw_Samples'])
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
