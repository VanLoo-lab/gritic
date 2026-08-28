import unittest

import networkx as nx
import numpy as np

from gritic import gritictimer, treetools
from gritic.route_tree import RouteTree


class RouteComponentPhasingTest(unittest.TestCase):
    @staticmethod
    def make_binary_tree(allele='Major'):
        tree = nx.DiGraph()
        tree.add_edges_from([(0, 1), (0, 2)])
        nx.set_node_attributes(tree, False, 'WGD_Symbol')
        nx.set_node_attributes(tree, False, 'Terminal_Node')
        nx.set_node_attributes(tree, allele, 'Allele')
        tree.nodes[1]['Terminal_Node'] = True
        tree.nodes[2]['Terminal_Node'] = True
        return tree

    def assert_component_labels_match_split(self, tree):
        node_phasing = treetools.get_node_phasing_tree(tree)
        major_tree, minor_tree = treetools.split_tree(tree)

        self.assertEqual(
            {node_phasing[node] for node in major_tree.nodes},
            {'Major'},
        )
        self.assertEqual(
            {node_phasing[node] for node in minor_tree.nodes},
            {'Minor'} if minor_tree.number_of_nodes() else set(),
        )
        return major_tree, minor_tree

    def test_single_component_is_major(self):
        tree = self.make_binary_tree()

        major_tree, minor_tree = self.assert_component_labels_match_split(tree)
        route_tree = RouteTree(
            tree,
            major_cn=2,
            minor_cn=0,
            wgd_status=False,
        )

        self.assertEqual(set(major_tree.nodes), set(tree.nodes))
        self.assertEqual(minor_tree.number_of_nodes(), 0)
        self.assertEqual(
            {
                attributes['Phasing']
                for attributes in route_tree.node_attributes.values()
            },
            {'Major'},
        )

    def test_balanced_components_use_model_major_minor_order(self):
        tree = nx.disjoint_union(
            self.make_binary_tree('Major'),
            self.make_binary_tree('Minor'),
        )

        major_tree, minor_tree = self.assert_component_labels_match_split(tree)

        self.assertEqual(major_tree.number_of_nodes(), 3)
        self.assertEqual(minor_tree.number_of_nodes(), 3)
        self.assertTrue(set(major_tree.nodes).isdisjoint(minor_tree.nodes))

    def test_balanced_gain_table_exports_major_minor_labels(self):
        classifier = gritictimer.RouteClassifier(
            major_cn=2,
            minor_cn=2,
            wgd_status=False,
            wgd_trees_status='No_WGD',
        )
        self.assertEqual(len(classifier.routes), 1)
        route_id, route = next(iter(classifier.routes.items()))
        route.node_timing = np.zeros((
            len(route.route_tree.non_phased_node_order),
            3,
        ))
        classifier.route_probabilities[route_id] = 1.0

        timing_table = classifier.get_timing_table()

        self.assertEqual(set(timing_table['Node_Phasing']), {'Major', 'Minor'})

    def test_component_roles_come_from_explicit_allele_labels(self):
        major_tree = nx.empty_graph(1, create_using=nx.DiGraph())
        nx.set_node_attributes(major_tree, False, 'WGD_Symbol')
        nx.set_node_attributes(major_tree, True, 'Terminal_Node')
        nx.set_node_attributes(major_tree, 'Major', 'Allele')
        tree = nx.disjoint_union(
            major_tree,
            self.make_binary_tree('Minor'),
        )

        major_tree, minor_tree = self.assert_component_labels_match_split(tree)

        self.assertEqual(major_tree.number_of_nodes(), 1)
        self.assertEqual(minor_tree.number_of_nodes(), 3)


if __name__ == '__main__':
    unittest.main()
