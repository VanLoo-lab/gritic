import tempfile
import unittest
from pathlib import Path
from unittest import mock

import networkx as nx

from gritic import treetools


def add_route_node(
    tree,
    node,
    allele,
    *,
    wgd=False,
    terminal=True,
):
    tree.add_node(
        node,
        Allele=allele,
        WGD_Symbol=wgd,
        Terminal_Node=terminal,
    )


def make_forest(*, include_minor=True, wgd_node=None):
    """Three-leaf Major component and optional two-leaf Minor component."""
    tree = nx.DiGraph()
    add_route_node(tree, 'm0', 'Major', terminal=False)
    add_route_node(tree, 'm1', 'Major')
    add_route_node(
        tree,
        'm2',
        'Major',
        wgd=wgd_node == 'm2',
        terminal=False,
    )
    add_route_node(tree, 'm3', 'Major')
    add_route_node(tree, 'm4', 'Major')
    tree.add_edges_from([
        ('m0', 'm1'),
        ('m0', 'm2'),
        ('m2', 'm3'),
        ('m2', 'm4'),
    ])
    if include_minor:
        add_route_node(tree, 'n0', 'Minor', terminal=False)
        add_route_node(tree, 'n1', 'Minor')
        add_route_node(tree, 'n2', 'Minor')
        tree.add_edges_from([('n0', 'n1'), ('n0', 'n2')])
    return tree


class FullBinaryTreeGenerationTest(unittest.TestCase):
    def test_unlabelled_full_binary_tree_counts(self):
        expected_counts = {0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 3}
        for leaf_count, expected in expected_counts.items():
            with self.subTest(leaf_count=leaf_count):
                self.assertEqual(
                    len(treetools.allPossibleFBT(leaf_count)),
                    expected,
                )

    def test_converted_trees_have_requested_leaves_and_full_branching(self):
        for leaf_count in range(1, 6):
            for binary_tree in treetools.allPossibleFBT(leaf_count):
                with self.subTest(leaf_count=leaf_count):
                    tree = treetools.convert_to_nx_tree(binary_tree, 'Major')
                    self.assertTrue(treetools.check_tree(tree))
                    terminals = [
                        node
                        for node, attributes in tree.nodes(data=True)
                        if attributes['Terminal_Node']
                    ]
                    self.assertEqual(len(terminals), leaf_count)
                    self.assertEqual(
                        tree.number_of_nodes(),
                        2 * leaf_count - 1,
                    )
                    self.assertEqual(
                        set(nx.get_node_attributes(tree, 'Allele').values()),
                        {'Major'},
                    )

    def test_convert_to_nx_trees_applies_requested_allele(self):
        trees = treetools.convert_to_nx_trees(
            treetools.allPossibleFBT(3),
            'Minor',
        )
        self.assertEqual(len(trees), 1)
        self.assertEqual(
            set(nx.get_node_attributes(trees[0], 'Allele').values()),
            {'Minor'},
        )

    def test_check_tree_rejects_unary_and_three_way_branches(self):
        unary = nx.DiGraph([(0, 1)])
        three_way = nx.DiGraph([(0, 1), (0, 2), (0, 3)])
        self.assertFalse(treetools.check_tree(unary))
        self.assertFalse(treetools.check_tree(three_way))
        self.assertTrue(treetools.check_tree(nx.DiGraph()))


class AlleleComponentValidationTest(unittest.TestCase):
    def test_split_tree_returns_exact_allele_components(self):
        tree = make_forest()
        major, minor = treetools.split_tree(tree)
        self.assertEqual(set(major), {'m0', 'm1', 'm2', 'm3', 'm4'})
        self.assertEqual(set(minor), {'n0', 'n1', 'n2'})
        self.assertEqual(set(major.edges), {
            ('m0', 'm1'),
            ('m0', 'm2'),
            ('m2', 'm3'),
            ('m2', 'm4'),
        })
        self.assertEqual(set(minor.edges), {('n0', 'n1'), ('n0', 'n2')})

    def test_split_tree_represents_absent_minor_as_empty_digraph(self):
        major, minor = treetools.split_tree(make_forest(include_minor=False))
        self.assertEqual(major.number_of_nodes(), 5)
        self.assertIsInstance(minor, nx.DiGraph)
        self.assertEqual(minor.number_of_nodes(), 0)

    def test_node_phasing_is_explicit_for_every_node(self):
        tree = make_forest()
        phasing = treetools.get_node_phasing_tree(tree)
        self.assertEqual(set(phasing), set(tree))
        self.assertEqual({phasing[node] for node in tree if node[0] == 'm'}, {
            'Major'
        })
        self.assertEqual({phasing[node] for node in tree if node[0] == 'n'}, {
            'Minor'
        })

    def test_rejects_missing_or_unknown_allele_labels(self):
        for value in (None, 'major', 'Other'):
            with self.subTest(value=value):
                tree = make_forest()
                if value is None:
                    del tree.nodes['m1']['Allele']
                else:
                    tree.nodes['m1']['Allele'] = value
                with self.assertRaisesRegex(ValueError, 'Major or Minor'):
                    treetools.split_tree(tree)

    def test_rejects_no_major_component(self):
        tree = nx.DiGraph()
        add_route_node(tree, 0, 'Minor')
        with self.assertRaisesRegex(ValueError, 'contain a Major'):
            treetools.split_tree(tree)

    def test_rejects_disconnected_same_allele_component(self):
        tree = make_forest(include_minor=False)
        add_route_node(tree, 'isolated', 'Major')
        with self.assertRaisesRegex(ValueError, 'one connected component'):
            treetools.split_tree(tree)

    def test_rejects_edges_between_alleles(self):
        tree = make_forest()
        tree.add_edge('m1', 'n0')
        with self.assertRaisesRegex(ValueError, 'cannot connect'):
            treetools.split_tree(tree)


class WgdDecorationAndPathTest(unittest.TestCase):
    def test_convert_to_wgd_tree_does_not_mutate_source(self):
        tree = make_forest()
        converted = treetools.convert_to_wgd_tree(tree, ['m0', 'n0'])
        self.assertFalse(tree.nodes['m0']['WGD_Symbol'])
        self.assertFalse(tree.nodes['n0']['WGD_Symbol'])
        self.assertTrue(converted.nodes['m0']['WGD_Symbol'])
        self.assertTrue(converted.nodes['n0']['WGD_Symbol'])

    def test_wgd_combination_status_rules(self):
        tree = make_forest()
        self.assertTrue(
            treetools.valid_wgd_node_combination(tree, [], 'Default')
        )
        self.assertFalse(
            treetools.valid_wgd_node_combination(tree, [], 'Only_WGD')
        )
        self.assertFalse(
            treetools.valid_wgd_node_combination(
                tree,
                [],
                'Only_WGD_2+2',
            )
        )
        self.assertTrue(
            treetools.valid_wgd_node_combination(tree, ['m0'], 'Only_WGD')
        )
        self.assertFalse(
            treetools.valid_wgd_node_combination(tree, ['m0'], 'No_WGD')
        )
        self.assertFalse(
            treetools.valid_wgd_node_combination(
                tree,
                ['m0'],
                'Only_WGD_2+2',
            )
        )
        self.assertFalse(
            treetools.valid_wgd_node_combination(
                tree,
                ['m0', 'm2'],
                'Default',
            )
        )
        self.assertTrue(
            treetools.valid_wgd_node_combination(
                tree,
                ['m2', 'n0'],
                'Only_WGD_2+2',
            )
        )

    def test_get_wgd_trees_enumerates_permitted_non_descendant_decorations(self):
        tree = make_forest()
        internal_nodes = {'m0', 'm2', 'n0'}
        only_wgd = treetools.get_wgd_trees(tree, 'Only_WGD')
        observed_sets = {
            frozenset(
                node
                for node, attributes in route.nodes(data=True)
                if attributes['WGD_Symbol']
            )
            for route in only_wgd
        }
        self.assertNotIn(frozenset(), observed_sets)
        self.assertTrue(
            {frozenset([node]) for node in internal_nodes}.issubset(
                observed_sets
            )
        )
        self.assertNotIn(frozenset(['m0', 'm2']), observed_sets)
        self.assertIn(frozenset(['m2', 'n0']), observed_sets)

    def test_possible_paths_cover_every_root_to_leaf_path_in_forest(self):
        paths = treetools.get_possible_paths(make_forest())
        self.assertEqual(paths, [
            ['m0', 'm1'],
            ['m0', 'm2', 'm3'],
            ['m0', 'm2', 'm4'],
            ['n0', 'n1'],
            ['n0', 'n2'],
        ])

    def test_wgd_paths_stop_at_first_wgd_and_ignore_paths_without_one(self):
        tree = make_forest(wgd_node='m2')
        self.assertEqual(treetools.get_wgd_paths(tree), [['m0', 'm2']])
        tree.nodes['n0']['WGD_Symbol'] = True
        self.assertEqual(
            treetools.get_wgd_paths(tree),
            [['m0', 'm2'], ['n0']],
        )
        self.assertEqual(treetools.get_wgd_paths(make_forest()), [])

    def test_path_helpers_accept_single_node_components(self):
        tree = nx.DiGraph()
        add_route_node(tree, 0, 'Major')
        self.assertEqual(treetools.get_possible_paths(tree), [[0]])
        self.assertEqual(treetools.get_wgd_paths(tree), [])
        tree.nodes[0]['WGD_Symbol'] = True
        self.assertEqual(treetools.get_wgd_paths(tree), [[0]])


class NodeAttributesTest(unittest.TestCase):
    def test_attributes_capture_topology_multiplicity_phasing_and_wgd_order(self):
        tree = make_forest(wgd_node='m2')
        attributes = treetools.get_node_attributes(tree, wgd_status=True)
        self.assertIsNone(attributes['m0']['Predecessor'])
        self.assertEqual(attributes['m0']['Successors'], ['m1', 'm2'])
        self.assertEqual(attributes['m0']['Multiplicity'], 3)
        self.assertEqual(attributes['m0']['WGD_Ordering'], 'Pre')
        self.assertEqual(attributes['m1']['Multiplicity'], 1)
        self.assertEqual(attributes['m1']['WGD_Ordering'], 'Calculate')
        self.assertEqual(attributes['m2']['Predecessor'], 'm0')
        self.assertEqual(attributes['m2']['Multiplicity'], 2)
        self.assertEqual(attributes['m2']['WGD_Ordering'], 'WGD')
        self.assertEqual(attributes['m3']['WGD_Ordering'], 'Post')
        self.assertEqual(set(attributes['m3']['Ancestors']), {'m2', 'm0'})
        self.assertEqual(attributes['n0']['Phasing'], 'Minor')
        self.assertEqual(attributes['n0']['WGD_Ordering'], 'Calculate')

    def test_non_wgd_status_marks_every_order_as_na(self):
        attributes = treetools.get_node_attributes(
            make_forest(wgd_node='m2'),
            wgd_status=False,
        )
        self.assertEqual(
            {value['WGD_Ordering'] for value in attributes.values()},
            {'NA'},
        )


class PlottingTreeConversionTest(unittest.TestCase):
    @staticmethod
    def crossing_tree(*, root_wgd=False):
        tree = nx.DiGraph()
        add_route_node(tree, 10, 'Major', wgd=root_wgd, terminal=False)
        add_route_node(tree, 20, 'Major')
        add_route_node(tree, 30, 'Major')
        add_route_node(tree, 40, 'Minor')
        tree.add_edges_from([(10, 20), (10, 30)])
        return tree

    def test_inserts_wgd_and_loss_nodes_at_each_unrepresented_crossing(self):
        source = self.crossing_tree()
        plotting = treetools.convert_to_plotting_tree(
            source,
            wgd_timing=0.5,
            route_timing_summary=[0.2, 0.6, 0.4, 0.7],
            node_order=[10, 20, 30, 40],
        )
        self.assertEqual(set(source), {10, 20, 30, 40})
        self.assertEqual(set(source.edges), {(10, 20), (10, 30)})
        loss_nodes = [
            node
            for node, attributes in plotting.nodes(data=True)
            if attributes['Loss_Symbol']
        ]
        inserted_wgds = [
            node
            for node, attributes in plotting.nodes(data=True)
            if attributes['WGD_Symbol']
        ]
        self.assertEqual(len(loss_nodes), 2)
        self.assertEqual(len(inserted_wgds), 2)
        self.assertEqual(
            {plotting.nodes[node]['Allele'] for node in loss_nodes},
            {'Major', 'Minor'},
        )
        for loss_node in loss_nodes:
            predecessor = next(plotting.predecessors(loss_node))
            self.assertTrue(plotting.nodes[predecessor]['WGD_Symbol'])
            self.assertEqual(plotting.out_degree(predecessor), 2)
        self.assertNotIn((10, 20), plotting.edges)
        major_wgd = next(
            node
            for node in inserted_wgds
            if plotting.nodes[node]['Allele'] == 'Major'
        )
        self.assertIn((10, major_wgd), plotting.edges)
        self.assertIn((major_wgd, 20), plotting.edges)
        minor_wgd = next(
            node
            for node in inserted_wgds
            if plotting.nodes[node]['Allele'] == 'Minor'
        )
        self.assertEqual(plotting.in_degree(minor_wgd), 0)
        self.assertIn((minor_wgd, 40), plotting.edges)

    def test_existing_wgd_anywhere_on_component_path_prevents_insertion(self):
        tree = self.crossing_tree(root_wgd=True)
        plotting = treetools.convert_to_plotting_tree(
            tree,
            wgd_timing=0.5,
            route_timing_summary=[0.2, 0.8, 0.9, 0.4],
            node_order=[10, 20, 30, 40],
        )
        self.assertEqual(plotting.number_of_nodes(), tree.number_of_nodes())
        self.assertFalse(any(
            attributes['Loss_Symbol']
            for _, attributes in plotting.nodes(data=True)
        ))

    def test_crossing_tolerance_boundaries_are_inclusive_only_above_low_bound(self):
        tolerance = 1e-7
        for current_timing, expected_loss in (
            (0.5 - tolerance, False),
            (0.5 - tolerance / 2, True),
        ):
            with self.subTest(current_timing=current_timing):
                tree = nx.DiGraph()
                add_route_node(tree, 0, 'Major')
                plotting = treetools.convert_to_plotting_tree(
                    tree,
                    wgd_timing=0.5,
                    route_timing_summary=[current_timing],
                    node_order=[0],
                    tol=tolerance,
                )
                self.assertEqual(
                    sum(
                        attributes['Loss_Symbol']
                        for _, attributes in plotting.nodes(data=True)
                    ),
                    int(expected_loss),
                )


class TreeIdentityAndLayoutTest(unittest.TestCase):
    def test_new_node_id_handles_sparse_integer_labels(self):
        tree = nx.DiGraph()
        tree.add_nodes_from([0, 1, 3, 4])
        self.assertEqual(treetools._get_new_node_id(tree), 5)
        tree = nx.DiGraph()
        tree.add_nodes_from([0, 2, 3])
        self.assertEqual(treetools._get_new_node_id(tree), 4)

    def test_combined_hierarchy_layout_offsets_minor_component(self):
        tree = nx.DiGraph()
        add_route_node(tree, 'm', 'Major')
        add_route_node(tree, 'n', 'Minor')
        positions = treetools.get_combined_hierarchy_pos(tree)
        self.assertEqual(positions['m'], (0.5, 0))
        self.assertEqual(positions['n'], (1.5, 0))

    def test_hierarchy_layout_rejects_non_tree(self):
        cycle = nx.DiGraph([(0, 1), (1, 0)])
        with self.assertRaisesRegex(TypeError, 'not a tree'):
            treetools.hierarchy_pos(cycle)

    def test_undirected_hierarchy_without_root_uses_a_valid_random_root(self):
        tree = nx.Graph([(0, 1), (1, 2)])
        with mock.patch.object(
            treetools.random,
            'choice',
            return_value=1,
        ) as choose:
            positions = treetools.hierarchy_pos(tree)
        choose.assert_called_once_with([0, 1, 2])
        self.assertEqual(set(positions), {0, 1, 2})
        self.assertEqual(positions[1], (0.5, 0))
        self.assertEqual(positions[0], (0.25, -0.2))
        self.assertEqual(positions[2], (0.75, -0.2))

    def test_write_tree_removes_terminal_nodes_from_written_copy(self):
        tree = nx.DiGraph()
        add_route_node(tree, 0, 'Major', terminal=False)
        add_route_node(tree, 1, 'Major')
        add_route_node(tree, 2, 'Major')
        tree.add_edges_from([(0, 1), (0, 2)])
        for node in tree:
            tree.nodes[node]['Full_Timing'] = 0.5
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'tree.graphml'
            treetools.write_tree(tree, path)
            written = nx.read_graphml(path)
        self.assertEqual(list(written), ['0'])
        self.assertNotIn('Full_Timing', written.nodes['0'])
        self.assertNotIn('Terminal_Node', written.nodes['0'])


if __name__ == '__main__':
    unittest.main()
