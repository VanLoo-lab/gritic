import itertools
import unittest
from unittest import mock

import networkx as nx

from gritic import gritictimer, treetools


class OrderedRouteIdentityTest(unittest.TestCase):
    @staticmethod
    def get_leaf_count(tree):
        return sum(tree.out_degree(node) == 0 for node in tree)

    def test_balanced_routes_are_full_ordered_cartesian_product(self):
        routes = treetools.get_nx_trees(
            major_cn=4,
            minor_cn=4,
            wgd_status=False,
            wgd_trees_status='No_WGD',
        )

        self.assertEqual(len(treetools.allPossibleFBT(4)), 2)
        self.assertEqual(len(routes), 4)

        component_pairs = {
            treetools.get_ordered_component_hashes(tree)
            for tree in routes.values()
        }
        allele_hashes = {pair[0] for pair in component_pairs}
        self.assertEqual(len(allele_hashes), 2)
        self.assertEqual(
            component_pairs,
            set(itertools.product(allele_hashes, repeat=2)),
        )

    def test_off_diagonal_routes_have_distinct_mirror_ids(self):
        routes = treetools.get_nx_trees(4, 4, False, 'No_WGD')
        off_diagonal_routes = []
        diagonal_routes = []
        for route_id, tree in routes.items():
            self.assertEqual(route_id, treetools.get_tree_hash(tree))
            major_hash, minor_hash = treetools.get_ordered_component_hashes(tree)
            mirror_id = treetools.get_mirror_tree_hash(tree)
            self.assertIn(mirror_id, routes)
            if major_hash == minor_hash:
                diagonal_routes.append((route_id, mirror_id))
            else:
                off_diagonal_routes.append((route_id, mirror_id))
                self.assertEqual(
                    treetools.get_ordered_component_hashes(routes[mirror_id]),
                    (minor_hash, major_hash),
                )

        self.assertEqual(len(off_diagonal_routes), 2)
        self.assertTrue(all(route_id != mirror_id for route_id, mirror_id in off_diagonal_routes))
        self.assertTrue(all(route_id == mirror_id for route_id, mirror_id in diagonal_routes))

    def test_ordered_wgd_dedup_keeps_allele_mirrors(self):
        routes = treetools.get_nx_trees(2, 2, True, 'Default')
        wgd_states = set()

        for route_id, tree in routes.items():
            major_tree, minor_tree = treetools.split_tree(tree)
            state = (
                any(tree.nodes[node]['WGD_Symbol'] for node in major_tree),
                any(tree.nodes[node]['WGD_Symbol'] for node in minor_tree),
            )
            wgd_states.add(state)
            mirror_id = treetools.get_mirror_tree_hash(tree)
            self.assertIn(mirror_id, routes)
            mirror_major, mirror_minor = treetools.split_tree(routes[mirror_id])
            mirror_state = (
                any(
                    routes[mirror_id].nodes[node]['WGD_Symbol']
                    for node in mirror_major
                ),
                any(
                    routes[mirror_id].nodes[node]['WGD_Symbol']
                    for node in mirror_minor
                ),
            )
            self.assertEqual(mirror_state, state[::-1])

        self.assertEqual(
            wgd_states,
            {(False, False), (True, False), (False, True), (True, True)},
        )

    def test_zero_minor_uses_explicit_empty_component_hash(self):
        routes = treetools.get_nx_trees(2, 0, False, 'No_WGD')
        self.assertEqual(len(routes), 1)
        tree = next(iter(routes.values()))

        _, minor_hash = treetools.get_ordered_component_hashes(tree)
        self.assertEqual(minor_hash, treetools.EMPTY_MINOR_COMPONENT_HASH)
        self.assertEqual(
            {attributes['Allele'] for _, attributes in tree.nodes(data=True)},
            {'Major'},
        )
        with self.assertRaisesRegex(ValueError, 'balanced two-allele'):
            treetools.get_mirror_tree_hash(tree)

    def test_plot_layout_uses_explicit_roles_not_component_size(self):
        major_tree = nx.DiGraph()
        major_tree.add_node(
            0,
            WGD_Symbol=False,
            Terminal_Node=True,
            Allele='Major',
        )
        minor_tree = treetools.convert_to_nx_tree(
            treetools.TreeNode(treetools.TreeNode(), treetools.TreeNode()),
            'Minor',
        )
        tree = nx.disjoint_union(major_tree, minor_tree)

        major_component, minor_component = treetools.split_tree(tree)
        positions = treetools.get_combined_hierarchy_pos(tree)
        major_root = next(iter(major_component))
        minor_root = next(
            node for node, degree in minor_component.in_degree() if degree == 0
        )

        self.assertEqual(major_component.number_of_nodes(), 1)
        self.assertEqual(minor_component.number_of_nodes(), 3)
        self.assertLess(positions[major_root][0], positions[minor_root][0])

    def test_plotting_nodes_inherit_allele_provenance(self):
        tree = next(iter(
            treetools.get_nx_trees(2, 0, False, 'No_WGD').values()
        ))
        node_order = list(nx.dfs_preorder_nodes(tree))
        plotting_tree = treetools.convert_to_plotting_tree(
            tree,
            wgd_timing=0.5,
            route_timing_summary=[0.75] * len(node_order),
            node_order=node_order,
        )

        self.assertGreater(plotting_tree.number_of_nodes(), tree.number_of_nodes())
        self.assertEqual(
            {attributes['Allele'] for _, attributes in plotting_tree.nodes(data=True)},
            {'Major'},
        )

    def test_unlabelled_route_nodes_are_rejected(self):
        tree = nx.DiGraph([(0, 1), (0, 2)])
        nx.set_node_attributes(tree, False, 'WGD_Symbol')

        with self.assertRaisesRegex(ValueError, 'Every route node must have Allele'):
            treetools.get_tree_hash(tree)

    def test_local_wgd_hash_collision_is_rejected(self):
        tree = next(iter(
            treetools.get_nx_trees(2, 2, False, 'No_WGD').values()
        ))

        with mock.patch.object(
            treetools,
            'get_tree_hash',
            return_value='forced-collision',
        ):
            with self.assertRaisesRegex(
                ValueError,
                'hash collision.*deduplicating Default WGD routes',
            ):
                treetools.get_wgd_trees(tree, 'Default')

    def test_final_route_store_hash_collision_is_rejected(self):
        with mock.patch.object(
            treetools,
            'get_tree_hash',
            return_value='forced-collision',
        ):
            with self.assertRaisesRegex(
                ValueError,
                r'hash collision.*building the 4\+0 route store',
            ):
                treetools.get_nx_trees(4, 0, False, 'No_WGD')

    def test_mirror_requires_equal_terminal_copy_counts(self):
        major_tree = nx.DiGraph([(0, 1), (1, 2)])
        minor_tree = nx.DiGraph([(0, 1), (0, 2)])
        for tree, allele in (
            (major_tree, 'Major'),
            (minor_tree, 'Minor'),
        ):
            nx.set_node_attributes(tree, False, 'WGD_Symbol')
            nx.set_node_attributes(tree, False, 'Terminal_Node')
            nx.set_node_attributes(tree, allele, 'Allele')
            for node in tree:
                if tree.out_degree(node) == 0:
                    tree.nodes[node]['Terminal_Node'] = True

        tree = nx.disjoint_union(major_tree, minor_tree)
        major_component, minor_component = treetools.split_tree(tree)
        self.assertEqual(
            major_component.number_of_nodes(),
            minor_component.number_of_nodes(),
        )
        self.assertNotEqual(
            self.get_leaf_count(major_component),
            self.get_leaf_count(minor_component),
        )

        with self.assertRaisesRegex(ValueError, 'terminal copy counts'):
            treetools.get_mirror_tree_hash(tree)

    def test_convert_to_nx_tree_initializes_supplied_current_node(self):
        tree = nx.DiGraph()
        tree.add_node(0, Existing='preserved')
        tree.add_node(2, WGD_Symbol=True)
        binary_tree = treetools.TreeNode(
            treetools.TreeNode(),
            treetools.TreeNode(),
        )

        result = treetools.convert_to_nx_tree(
            binary_tree,
            'Major',
            D=tree,
            current_node_id=2,
        )

        self.assertIs(result, tree)
        self.assertEqual(set(tree.successors(2)), {3, 4})
        self.assertEqual(tree.nodes[2]['Allele'], 'Major')
        self.assertTrue(tree.nodes[2]['WGD_Symbol'])
        self.assertFalse(tree.nodes[2]['Terminal_Node'])
        self.assertEqual(tree.nodes[0]['Existing'], 'preserved')
        for node in (3, 4):
            self.assertEqual(tree.nodes[node]['Allele'], 'Major')
            self.assertFalse(tree.nodes[node]['WGD_Symbol'])
            self.assertTrue(tree.nodes[node]['Terminal_Node'])

    def test_convert_to_nx_tree_rejects_conflicting_supplied_allele(self):
        tree = nx.DiGraph()
        tree.add_node(0, Allele='Minor')

        with self.assertRaisesRegex(ValueError, 'Cannot extend.*Minor.*Major'):
            treetools.convert_to_nx_tree(
                treetools.TreeNode(),
                'Major',
                D=tree,
                current_node_id=0,
            )

    def test_all_permitted_states_preserve_ordered_route_invariants(self):
        full_route_id_by_short_id = {}
        n_states = 0

        for wgd_status in (False, True):
            wgd_trees_status = 'Default' if wgd_status else 'No_WGD'
            for major_cn in range(1, 9):
                for minor_cn in range(major_cn + 1):
                    if not gritictimer.check_permitted_cn_state(
                        major_cn,
                        minor_cn,
                        wgd_status,
                    ):
                        continue

                    n_states += 1
                    state = f'{major_cn}+{minor_cn}, WGD={wgd_status}'
                    routes = treetools.get_nx_trees(
                        major_cn,
                        minor_cn,
                        wgd_status,
                        wgd_trees_status,
                    )
                    self.assertGreater(len(routes), 0, state)
                    self.assertLessEqual(len(routes), 500, state)
                    self.assertEqual(
                        len({route_id[:9] for route_id in routes}),
                        len(routes),
                        state,
                    )

                    for route_id, tree in routes.items():
                        self.assertEqual(route_id, treetools.get_tree_hash(tree))
                        expected_alleles = {'Major'}
                        if minor_cn > 0:
                            expected_alleles.add('Minor')
                        self.assertEqual(
                            {
                                attributes['Allele']
                                for _, attributes in tree.nodes(data=True)
                            },
                            expected_alleles,
                            state,
                        )

                        major_tree, minor_tree = treetools.split_tree(tree)
                        self.assertEqual(
                            self.get_leaf_count(major_tree),
                            major_cn,
                            state,
                        )
                        self.assertEqual(
                            self.get_leaf_count(minor_tree),
                            minor_cn,
                            state,
                        )

                        if major_cn == minor_cn and minor_cn > 0:
                            self.assertIn(
                                treetools.get_mirror_tree_hash(tree),
                                routes,
                                state,
                            )

                        short_id = route_id[:9]
                        previous_route_id = full_route_id_by_short_id.setdefault(
                            short_id,
                            route_id,
                        )
                        self.assertEqual(previous_route_id, route_id, state)

        self.assertEqual(n_states, 59)


if __name__ == '__main__':
    unittest.main()
