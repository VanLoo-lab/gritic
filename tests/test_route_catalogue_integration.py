"""Inference consumes catalogue identity without generating or hashing routes."""

from unittest import mock

import pytest

from gritic import gritictimer, treetools
from gritic.route_catalogue import get_routes


@pytest.mark.parametrize('major,minor,wgd,mode', [
    (1, 0, False, 'No_WGD'),
    (1, 1, True, 'Default'),
    (3, 1, True, 'Default'),
    (4, 4, True, 'Only_WGD'),
    (2, 2, True, 'Only_WGD_2+2'),
])
def test_classifier_loads_catalogue_without_enumeration_or_hashing(major, minor, wgd, mode):
    entries = get_routes(major, minor, wgd, mode)
    with mock.patch.object(treetools, 'get_nx_trees', side_effect=AssertionError('enumeration')), \
            mock.patch.object(treetools, 'get_tree_hash', side_effect=AssertionError('hashing')), \
            mock.patch.object(treetools, 'get_mirror_tree_hash', side_effect=AssertionError('mirror hashing')):
        classifier = gritictimer.RouteClassifier(major, minor, wgd, mode)
    assert list(classifier.routes) == [entry.route_id for entry in entries]
    for entry in entries:
        route = classifier.routes[entry.route_id]
        tree = route.route_tree.main_tree
        assert route.short_id == entry.short_id
        assert route.mirror_route_id == entry.mirror_route_id
        assert list(tree.nodes) == [node.node_id for node in entry.nodes]
        assert set(tree.edges) == {(node.parent_id, node.node_id) for node in entry.nodes
                                   if node.parent_id is not None}
        assert route.route_tree.timeable_nodes == [node.node_id for node in entry.nodes
                                                   if not node.terminal and not node.wgd]
        for node in entry.nodes:
            assert tree.nodes[node.node_id] == {
                'Allele': node.allele,
                'WGD_Symbol': node.wgd,
                'Terminal_Node': node.terminal,
            }


def test_separate_classifiers_cannot_mutate_catalogue_or_each_others_trees():
    first = gritictimer.RouteClassifier(3, 1, True, 'Default')
    second = gritictimer.RouteClassifier(3, 1, True, 'Default')
    route_id = next(iter(first.routes))
    first_tree = first.routes[route_id].route_tree.main_tree
    second_tree = second.routes[route_id].route_tree.main_tree
    node_id = next(iter(first_tree))
    first_tree.nodes[node_id]['Allele'] = 'changed'
    assert second_tree.nodes[node_id]['Allele'] == 'Major'
    entry = next(entry for entry in get_routes(3, 1, True) if entry.route_id == route_id)
    assert entry.nodes[0].allele == 'Major'
