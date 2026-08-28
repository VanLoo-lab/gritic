import unittest

import networkx as nx
import numpy as np

from gritic import treetools
from gritic.route_tree import RouteTree


class RouteTreeBatchEventCountingTest(unittest.TestCase):
    @staticmethod
    def make_split_wgd_route():
        tree = nx.DiGraph()
        tree.add_edges_from([(0, 1), (0, 2), (3, 4), (3, 5)])
        nx.set_node_attributes(tree, False, "WGD_Symbol")
        nx.set_node_attributes(tree, False, "Terminal_Node")
        nx.set_node_attributes(tree, "Major", "Allele")
        nx.set_node_attributes(tree, "Minor", "Allele")
        for node in (0, 1, 2):
            tree.nodes[node]["Allele"] = "Major"
        for node in (1, 2, 4, 5):
            tree.nodes[node]["Terminal_Node"] = True
        tree.nodes[0]["WGD_Symbol"] = True
        return RouteTree(
            tree,
            major_cn=2,
            minor_cn=2,
            wgd_status=True,
        )

    @staticmethod
    def scalar_counts(route_tree, node_timing, wgd_timing):
        counts = [
            route_tree.get_n_events(node_timing[:, index], wgd_timing[index])
            for index in range(node_timing.shape[1])
        ]
        return tuple(np.asarray(values) for values in zip(*counts))

    def assert_batch_matches_scalar(self, route_tree, node_timing, wgd_timing):
        expected = self.scalar_counts(route_tree, node_timing, wgd_timing)
        actual = route_tree.get_n_events_batch(node_timing, wgd_timing)
        for expected_values, actual_values in zip(expected, actual):
            np.testing.assert_equal(actual_values, expected_values)

    def test_wgd_batch_matches_scalar_across_generated_routes(self):
        rng = np.random.default_rng(807)
        trees = treetools.get_nx_trees(
            major_cn=4,
            minor_cn=2,
            wgd_status=True,
            wgd_trees_status="Only_WGD",
        )

        for route_id, tree in trees.items():
            with self.subTest(route_id=route_id):
                route_tree = RouteTree(
                    tree,
                    major_cn=4,
                    minor_cn=2,
                    wgd_status=True,
                )
                node_timing = rng.uniform(
                    0.0,
                    1.0,
                    size=(len(route_tree.non_phased_node_order), 128),
                )
                wgd_timing = rng.uniform(0.0, 1.0, size=128)
                self.assert_batch_matches_scalar(
                    route_tree,
                    node_timing,
                    wgd_timing,
                )

    def test_batch_preserves_plotting_tree_tolerance_boundaries(self):
        route_tree = self.make_split_wgd_route()
        node_indices = {
            node: index
            for index, node in enumerate(route_tree.non_phased_node_order)
        }
        tolerance = 1e-7
        wgd_time = 0.5
        lower_boundary = wgd_time - tolerance
        upper_boundary = wgd_time + tolerance

        node_timing = np.zeros((len(node_indices), 6))
        node_timing[node_indices[4], 0] = lower_boundary
        node_timing[node_indices[4], 1] = np.nextafter(lower_boundary, np.inf)
        node_timing[node_indices[3], 2] = upper_boundary
        node_timing[node_indices[4], 2] = 0.8
        node_timing[node_indices[3], 3] = np.nextafter(upper_boundary, np.inf)
        node_timing[node_indices[4], 3] = 0.8
        node_timing[node_indices[4], 4] = np.nan
        wgd_timing = np.full(6, wgd_time)
        wgd_timing[5] = np.nan

        n_events, pre_wgd_losses, post_wgd_losses = (
            route_tree.get_n_events_batch(node_timing, wgd_timing)
        )

        np.testing.assert_array_equal(post_wgd_losses, [0, 1, 2, 1, 0, 0])
        np.testing.assert_array_equal(pre_wgd_losses, np.zeros(6, dtype=int))
        np.testing.assert_array_equal(n_events, [2, 3, 4, 3, 2, 2])
        self.assert_batch_matches_scalar(route_tree, node_timing, wgd_timing)

    def test_non_wgd_batch_matches_scalar_and_marks_losses_unidentifiable(self):
        tree = next(
            iter(
                treetools.get_nx_trees(
                    major_cn=3,
                    minor_cn=0,
                    wgd_status=False,
                    wgd_trees_status="No_WGD",
                ).values()
            )
        )
        route_tree = RouteTree(
            tree,
            major_cn=3,
            minor_cn=0,
            wgd_status=False,
        )
        node_timing = np.zeros((len(route_tree.non_phased_node_order), 9))
        wgd_timing = np.linspace(0.1, 0.9, num=9)

        n_events, pre_wgd_losses, post_wgd_losses = (
            route_tree.get_n_events_batch(node_timing, wgd_timing)
        )

        np.testing.assert_array_equal(n_events, np.full(9, 3))
        self.assertTrue(np.isnan(pre_wgd_losses).all())
        self.assertTrue(np.isnan(post_wgd_losses).all())
        self.assert_batch_matches_scalar(route_tree, node_timing, wgd_timing)

    def test_empty_batch_returns_empty_arrays(self):
        route_tree = self.make_split_wgd_route()
        node_timing = np.empty((len(route_tree.non_phased_node_order), 0))

        counts = route_tree.get_n_events_batch(node_timing, np.empty(0))

        self.assertEqual([values.shape for values in counts], [(0,), (0,), (0,)])

    def test_batch_validates_input_shapes(self):
        route_tree = self.make_split_wgd_route()
        n_nodes = len(route_tree.non_phased_node_order)

        with self.assertRaisesRegex(ValueError, "2-D"):
            route_tree.get_n_events_batch(np.zeros(n_nodes), np.zeros(1))
        with self.assertRaisesRegex(ValueError, "one row per route node"):
            route_tree.get_n_events_batch(np.zeros((n_nodes - 1, 2)), np.zeros(2))
        with self.assertRaisesRegex(ValueError, "1-D"):
            route_tree.get_n_events_batch(np.zeros((n_nodes, 2)), np.zeros((1, 2)))
        with self.assertRaisesRegex(ValueError, "same number of samples"):
            route_tree.get_n_events_batch(np.zeros((n_nodes, 2)), np.zeros(3))


if __name__ == "__main__":
    unittest.main()
