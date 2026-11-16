"""The frozen catalogue preserves current route identity and node alignment."""

import copy
from dataclasses import FrozenInstanceError
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

from gritic.route_catalogue import catalogue_metadata, get_route, get_routes
from gritic.route_catalogue.catalogue import _read_catalogue


STATES = [tuple(state) for state in catalogue_metadata()["copy_number_states"]]
MODES = ("Default", "No_WGD", "Only_WGD", "Only_WGD_2+2")


@pytest.mark.parametrize("major,minor", STATES)
def test_all_routes_and_modes_match_reference_enumeration_exactly(major, minor):
    from gritic import treetools

    for wgd in (False, True):
        for mode in MODES:
            expected = treetools.get_nx_trees(major, minor, wgd, mode)
            actual = get_routes(major, minor, wgd, mode)
            assert [route.route_id for route in actual] == list(expected), (major, minor, wgd, mode)
            for route in actual:
                tree = expected[route.route_id]
                assert [node.node_id for node in route.nodes] == list(tree)
                assert [(node.parent_id, node.node_id) for node in route.nodes if node.parent_id is not None] == [
                    (next(iter(tree.predecessors(node))), node) for node in tree if tree.in_degree(node)
                ]
                for node in route.nodes:
                    assert (node.allele, node.wgd, node.terminal) == (
                        tree.nodes[node.node_id]["Allele"], tree.nodes[node.node_id]["WGD_Symbol"],
                        tree.nodes[node.node_id]["Terminal_Node"],
                    )
                mirror = treetools.get_mirror_tree_hash(tree) if major == minor else None
                assert route.mirror_route_id == mirror


def test_every_permitted_inference_state_and_pooled_wgd_model_is_covered():
    from gritic.gritictimer import check_permitted_cn_state

    for wgd in (False, True):
        for major in range(1, 10):
            for minor in range(major + 1):
                if check_permitted_cn_state(major, minor, wgd):
                    assert get_routes(major, minor, wgd)
    for minor in range(3):
        assert get_routes(2, minor, True)


def test_resolve_saved_and_full_identifiers_without_enumerating():
    for major, minor in STATES:
        for route in get_routes(major, minor, True):
            assert get_route(major, minor, True, route.route_id) is route
            assert get_route(major, minor, True, route.short_id) is route
    with pytest.raises(ValueError, match="no unique match"):
        get_route(3, 1, True, "000000000")
    with pytest.raises(ValueError, match="nine-character"):
        get_route(3, 1, True, "789")
    wgd_route = next(route for route in get_routes(3, 1, True) if any(n.wgd for n in route.nodes))
    with pytest.raises(ValueError, match="no unique match"):
        get_route(3, 1, False, wgd_route.short_id)


@pytest.mark.parametrize("arguments", [
    (9, 0, False), (8, 2, True), (0, 0, False), (3, 4, True), (True, 0, False),
    (3.0, 1, False), (3, -1, True), (3, 1, 1), (3, 1, False, "WGD"),
])
def test_unknown_or_invalid_model_is_rejected(arguments):
    with pytest.raises(ValueError):
        get_routes(*arguments)


def test_catalogue_routes_are_immutable_and_metadata_is_detached():
    route = get_routes(3, 1, True)[0]
    with pytest.raises(FrozenInstanceError):
        route.route_id = "changed"
    with pytest.raises(FrozenInstanceError):
        route.nodes[0].parent_id = 99
    metadata = catalogue_metadata()
    metadata["copy_number_states"][0][0] = 99
    metadata["route_count"] = 99
    assert catalogue_metadata()["route_count"] == 2796
    assert not {"version", "identity_version", "networkx_version"} & catalogue_metadata().keys()
    assert catalogue_metadata()["copy_number_states"][0] == [1, 0]


def test_catalogue_import_does_not_load_inference_graphs_or_arrays():
    root = Path(__file__).resolve().parents[1]
    script = (
        "import sys; sys.path.insert(0, sys.argv[1]); "
        "from gritic.route_catalogue import get_routes; assert get_routes(8,1,True); "
        "assert not any(name in {'gritic.gritictimer','gritic.treetools'} "
        "or name.split('.')[0] in {'networkx','numpy','scipy','matplotlib'} "
        "for name in sys.modules)"
    )
    result = subprocess.run([sys.executable, "-I", "-c", script, str(root)],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


@pytest.fixture
def document():
    source = Path(__file__).resolve().parents[1] / "gritic/route_catalogue/data/routes.json.gz"
    return json.loads(gzip.decompress(source.read_bytes()))


@pytest.mark.parametrize("mutation,match", [
    (lambda d: d.update(format="other"), "Invalid route catalogue format"),
    (lambda d: d.update(route_count=0), "counts"),
    (lambda d: d["states"][0]["routes"].append(copy.deepcopy(d["states"][0]["routes"][0])), "collide"),
    (lambda d: d["states"][0]["routes"][0]["nodes"][0].__setitem__(1, 99), "parent"),
    (lambda d: d["states"][0]["routes"][0]["nodes"][0].__setitem__(3, True), "WGD"),
    (lambda d: d["states"][1]["routes"][0].update(mirror_route_id="0" * 32), "mirror"),
])
def test_invalid_catalogue_data_is_rejected(document, mutation, match):
    mutation(document)
    with pytest.raises(ValueError, match=match):
        _read_catalogue(document)
