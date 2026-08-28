import unittest

import networkx as nx
import numpy as np

from gritic.route_tree import RouteTree


def add_node(tree, node, allele, *, terminal=True, wgd=False):
    tree.add_node(
        node,
        Allele=allele,
        Terminal_Node=terminal,
        WGD_Symbol=wgd,
    )


def make_two_plus_one_tree(*, root_wgd=False):
    tree = nx.DiGraph()
    add_node(tree, 0, 'Major', terminal=False, wgd=root_wgd)
    add_node(tree, 1, 'Major')
    add_node(tree, 2, 'Major')
    add_node(tree, 3, 'Minor')
    tree.add_edges_from([(0, 1), (0, 2)])
    return tree


def make_three_plus_one_tree(*, wgd_node=2, include_minor=True):
    tree = nx.DiGraph()
    add_node(tree, 0, 'Major', terminal=False)
    add_node(tree, 1, 'Major')
    add_node(
        tree,
        2,
        'Major',
        terminal=False,
        wgd=wgd_node == 2,
    )
    add_node(tree, 3, 'Major')
    add_node(tree, 4, 'Major')
    tree.add_edges_from([(0, 1), (0, 2), (2, 3), (2, 4)])
    if include_minor:
        add_node(tree, 5, 'Minor')
    return tree


class RouteTreeConstructionTest(unittest.TestCase):
    def test_node_orders_and_timeable_nodes_follow_dfs(self):
        route = RouteTree(
            make_three_plus_one_tree(),
            major_cn=3,
            minor_cn=1,
            wgd_status=True,
        )
        self.assertEqual(route.non_phased_node_order, [0, 1, 2, 3, 4, 5])
        self.assertEqual(route.get_node_order('Major'), [0, 1, 2, 3, 4])
        self.assertEqual(route.get_node_order('Minor'), [5])
        self.assertEqual(route.timeable_nodes, [0])
        self.assertEqual(route.wgd_nodes, [2])

    def test_invalid_allele_for_node_order_is_rejected(self):
        route = RouteTree(
            make_two_plus_one_tree(),
            major_cn=2,
            minor_cn=1,
            wgd_status=False,
        )
        for allele in ('major', 'Other', ''):
            with self.subTest(allele=allele):
                with self.assertRaisesRegex(ValueError, 'None, Major or Minor'):
                    route.get_node_order(allele)

    def test_timing_matrix_maps_node_periods_to_unphased_and_allele_bins(self):
        route = RouteTree(
            make_two_plus_one_tree(),
            major_cn=2,
            minor_cn=1,
            wgd_status=False,
        )
        np.testing.assert_array_equal(route.timing_matrix, np.array([
            [0, 1, 0, 1, 0],
            [1, 0, 1, 0, 0],
            [1, 0, 1, 0, 0],
            [1, 0, 0, 0, 1],
        ]))

    def test_sum_constraints_have_one_row_per_root_to_leaf_path(self):
        route = RouteTree(
            make_three_plus_one_tree(),
            major_cn=3,
            minor_cn=1,
            wgd_status=True,
        )
        np.testing.assert_array_equal(route.sum_constraint_matrix, np.array([
            [1, 1, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ]))

    def test_wgd_constraints_stop_at_wgd_node(self):
        route = RouteTree(
            make_three_plus_one_tree(),
            major_cn=3,
            minor_cn=1,
            wgd_status=True,
        )
        np.testing.assert_array_equal(
            route.wgd_constraint_matrix,
            np.array([[1, 0, 1, 0, 0, 0]]),
        )
        matrix, values = route.get_combined_constraints(0.5)
        np.testing.assert_array_equal(matrix, np.array([
            [1, 1, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [1, 0, 1, 0, 0, 0],
        ]))
        np.testing.assert_array_equal(values, [1, 1, 1, 1, 0.5])

        feasible_periods = np.array([0.2, 0.8, 0.3, 0.5, 0.5, 1.0])
        np.testing.assert_allclose(matrix @ feasible_periods, values)

    def test_no_wgd_path_produces_only_sum_constraints(self):
        route = RouteTree(
            make_two_plus_one_tree(),
            major_cn=2,
            minor_cn=1,
            wgd_status=False,
        )
        self.assertIsNone(route.wgd_constraint_matrix)
        matrix, values = route.get_combined_constraints(np.nan)
        np.testing.assert_array_equal(matrix, route.sum_constraint_matrix)
        np.testing.assert_array_equal(values, np.ones(3))

    def test_post_wgd_loss_indices_keep_only_branches_outside_wgd_path(self):
        route = RouteTree(
            make_three_plus_one_tree(),
            major_cn=3,
            minor_cn=1,
            wgd_status=True,
        )
        np.testing.assert_array_equal(
            route._post_wgd_loss_node_indices,
            [1, 5],
        )
        np.testing.assert_array_equal(
            route._post_wgd_loss_predecessor_indices,
            [0, -1],
        )
        self.assertEqual(route._post_wgd_loss_node_indices.dtype, np.intp)
        self.assertEqual(
            route._post_wgd_loss_predecessor_indices.dtype,
            np.intp,
        )


class RouteTreeEventCountingTest(unittest.TestCase):
    def test_non_wgd_counts_binary_gains_and_unidentifiable_minor_loss(self):
        with_minor = RouteTree(
            make_two_plus_one_tree(),
            major_cn=2,
            minor_cn=1,
            wgd_status=False,
        )
        self.assertEqual(
            with_minor.get_n_events(np.zeros(4), np.nan)[0],
            1,
        )
        self.assertTrue(np.isnan(with_minor.get_n_events(np.zeros(4), np.nan)[1]))
        self.assertTrue(np.isnan(with_minor.get_n_events(np.zeros(4), np.nan)[2]))

        no_minor = RouteTree(
            make_two_plus_one_tree().subgraph([0, 1, 2]).copy(),
            major_cn=2,
            minor_cn=0,
            wgd_status=False,
        )
        self.assertEqual(no_minor.get_n_events(np.zeros(3), np.nan)[0], 2)

    def test_wgd_counts_post_wgd_loss_on_independent_minor_component(self):
        route = RouteTree(
            make_three_plus_one_tree(),
            major_cn=3,
            minor_cn=1,
            wgd_status=True,
        )
        no_crossing = np.array([0.2, 0.4, 0.3, 0.5, 0.5, 0.4])
        crossing = no_crossing.copy()
        crossing[5] = 0.6
        self.assertEqual(route.get_n_events(no_crossing, 0.5), (2, 0, 0))
        self.assertEqual(route.get_n_events(crossing, 0.5), (3, 0, 1))

    def test_wgd_single_component_adds_pre_wgd_minor_loss(self):
        route = RouteTree(
            make_three_plus_one_tree(include_minor=False),
            major_cn=3,
            minor_cn=0,
            wgd_status=True,
        )
        self.assertEqual(
            route.get_n_events(np.array([0.2, 0.4, 0.3, 0.5, 0.5]), 0.5),
            (3, 1, 0),
        )

if __name__ == '__main__':
    unittest.main()
