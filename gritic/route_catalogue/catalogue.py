"""Validated lookup of a frozen catalogue generated with GRITIC's identity rules."""

from dataclasses import dataclass
from functools import cache
import gzip
from importlib.resources import files
import json
from numbers import Integral
import re
from types import MappingProxyType


_STATUSES = frozenset({"Default", "No_WGD", "Only_WGD", "Only_WGD_2+2"})
_HASH = re.compile(r"[0-9a-f]{32}\Z")
_SOURCE = "data/routes.json.gz"


@dataclass(frozen=True, slots=True)
class CatalogueNode:
    node_id: int
    parent_id: int | None
    allele: str
    wgd: bool
    terminal: bool


@dataclass(frozen=True, slots=True)
class CatalogueRoute:
    route_id: str
    major_cn: int
    minor_cn: int
    nodes: tuple[CatalogueNode, ...]
    mirror_route_id: str | None = None

    @property
    def short_id(self) -> str:
        """The nine-character Route identifier written by current GRITIC."""
        return self.route_id[:9]


def _integer(value, label, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{label} must be an integer of at least {minimum}")
    return int(value)


def _validate_route(route):
    if not isinstance(route.route_id, str) or _HASH.fullmatch(route.route_id) is None:
        raise ValueError("Catalogue route IDs must be full lowercase topology hashes")
    nodes = {node.node_id: node for node in route.nodes}
    if not nodes or len(nodes) != len(route.nodes):
        raise ValueError("Catalogue node IDs must be nonempty and unique")
    children = {node_id: [] for node_id in nodes}
    roots = {"Major": [], "Minor": []}
    for node in route.nodes:
        _integer(node.node_id, "Catalogue node ID")
        if node.allele not in roots or type(node.wgd) is not bool or type(node.terminal) is not bool:
            raise ValueError("Catalogue node attributes are invalid")
        if node.parent_id is None:
            roots[node.allele].append(node.node_id)
        else:
            _integer(node.parent_id, "Catalogue parent ID")
            if node.parent_id not in nodes or nodes[node.parent_id].allele != node.allele:
                raise ValueError("Catalogue parent must exist in the same allele")
            children[node.parent_id].append(node.node_id)
    seen = set()
    for allele, copies in (("Major", route.major_cn), ("Minor", route.minor_cn)):
        if len(roots[allele]) != int(copies > 0):
            raise ValueError("Catalogue allele root count disagrees with copy number")
        terminal_count = 0
        stack = [(node_id, False) for node_id in roots[allele]]
        while stack:
            node_id, wgd_ancestor = stack.pop()
            if node_id in seen:
                raise ValueError("Catalogue contains a cycle")
            seen.add(node_id)
            node = nodes[node_id]
            if len(children[node_id]) != (0 if node.terminal else 2):
                raise ValueError("Catalogue nodes must be terminal or bifurcating")
            if node.wgd and (node.terminal or wgd_ancestor):
                raise ValueError("Catalogue WGD nodes must be internal and non-ancestral")
            terminal_count += node.terminal
            stack.extend((child, wgd_ancestor or node.wgd) for child in children[node_id])
        if terminal_count != copies:
            raise ValueError("Catalogue leaf count disagrees with copy number")
    if seen != set(nodes):
        raise ValueError("Catalogue contains disconnected nodes")
    balanced = route.major_cn == route.minor_cn and route.minor_cn > 0
    if balanced != (route.mirror_route_id is not None):
        raise ValueError("Only balanced routes must have a mirror route ID")


def _read_catalogue(document):
    if not isinstance(document, dict) or document.get("format") != "gritic_route_catalogue":
        raise ValueError("Invalid route catalogue format")
    result = {}
    total = 0
    for state in document["states"]:
        major = _integer(state["major_cn"], "Major copy number", 1)
        minor = _integer(state["minor_cn"], "Minor copy number")
        if minor > major or (major, minor) in result:
            raise ValueError("Catalogue copy-number states must be valid and unique")
        routes = tuple(CatalogueRoute(
            row["route_id"], major, minor,
            tuple(CatalogueNode(*node) for node in row["nodes"]), row["mirror_route_id"],
        ) for row in state["routes"])
        if not routes:
            raise ValueError("Catalogue states must contain routes")
        ids = {route.route_id: route for route in routes}
        if len(ids) != len(routes) or len({route.short_id for route in routes}) != len(routes):
            raise ValueError("Catalogue route IDs collide within a copy-number state")
        for route in routes:
            _validate_route(route)
            if route.mirror_route_id is not None:
                mirror = ids.get(route.mirror_route_id)
                if mirror is None or mirror.mirror_route_id != route.route_id:
                    raise ValueError("Catalogue mirror relationships must exist and be reciprocal")
        result[(major, minor)] = routes
        total += len(routes)
    if not result or document["state_count"] != len(result) or document["route_count"] != total:
        raise ValueError("Catalogue counts disagree with its contents")
    return MappingProxyType(result)


@cache
def _catalogue():
    payload = gzip.decompress(files("gritic.route_catalogue").joinpath(_SOURCE).read_bytes())
    document = json.loads(payload)
    return _read_catalogue(document), document


def catalogue_metadata():
    """Return detached provenance and coverage metadata, without route records."""
    _, document = _catalogue()
    metadata = {key: value for key, value in document.items() if key != "states"}
    metadata["copy_number_states"] = [[state["major_cn"], state["minor_cn"]]
                                      for state in document["states"]]
    return json.loads(json.dumps(metadata))


def get_routes(major_cn, minor_cn, wgd_status, wgd_trees_status="Default"):
    """Return routes in original generation order for a supported model state.

    WGD restrictions apply only when ``wgd_status`` is true, as in GRITIC's
    generator. No inference, hashing, graph construction, or dynamic fallback
    occurs during lookup. Unsupported states and modes are errors.
    """
    major = _integer(major_cn, "major_cn", 1)
    minor = _integer(minor_cn, "minor_cn")
    if type(wgd_status) is not bool:
        raise ValueError("wgd_status must be a boolean")
    if not isinstance(wgd_trees_status, str) or wgd_trees_status not in _STATUSES:
        raise ValueError(f"wgd_trees_status must be one of {sorted(_STATUSES)}")
    catalogue, _ = _catalogue()
    if (major, minor) not in catalogue:
        raise ValueError(f"Copy-number state {major}+{minor} is not in the route catalogue")
    routes = catalogue[(major, minor)]
    if not wgd_status or wgd_trees_status == "No_WGD":
        return tuple(route for route in routes if not any(node.wgd for node in route.nodes))
    if wgd_trees_status == "Only_WGD":
        return tuple(route for route in routes if any(node.wgd for node in route.nodes))
    if wgd_trees_status == "Only_WGD_2+2":
        return tuple(route for route in routes if sum(node.wgd for node in route.nodes) >= 2)
    return routes


def get_route(major_cn, minor_cn, wgd_status, route_id):
    """Resolve an exact full or nine-character saved ID within a model state."""
    if not isinstance(route_id, str) or len(route_id) not in (9, 32):
        raise ValueError("route_id must be a full hash or a nine-character saved Route ID")
    matches = [route for route in get_routes(major_cn, minor_cn, wgd_status)
               if route_id in (route.route_id, route.short_id)]
    if len(matches) != 1:
        raise ValueError(
            f"Route {route_id!r} has no unique match in the catalogue "
            f"for {major_cn}+{minor_cn}, WGD={wgd_status}"
        )
    return matches[0]
