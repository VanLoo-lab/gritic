"""Rebuild the frozen route catalogue from the GRITIC reference enumerator."""

import argparse
import gzip
import json
from pathlib import Path
import sys


def build_catalogue():
    from gritic import treetools

    # These are GRITIC's supported inference states, plus 2+m WGD models
    # used by route construction and pooled-WGD interpretation.
    states = [(major, minor) for major in range(1, 9) for minor in range(major + 1)
              if major <= 5 or major == 6 and minor <= 4
              or major == 7 and minor <= 3 or major == 8 and minor <= 1]
    records = []
    for major, minor in states:
        routes = []
        trees = treetools.get_nx_trees(major, minor, True, "Default")
        for route_id, tree in trees.items():
            nodes = []
            for node, attributes in tree.nodes(data=True):
                parents = list(tree.predecessors(node))
                if len(parents) > 1:
                    raise ValueError("A route node cannot have multiple parents")
                nodes.append([int(node), int(parents[0]) if parents else None,
                              attributes["Allele"], bool(attributes["WGD_Symbol"]),
                              bool(attributes["Terminal_Node"])])
            routes.append({
                "route_id": route_id,
                "nodes": nodes,
                "mirror_route_id": treetools.get_mirror_tree_hash(tree) if major == minor else None,
            })
        records.append({"major_cn": major, "minor_cn": minor, "routes": routes})
        print(f"{major}+{minor}: {len(routes)} routes", flush=True)
    return {
        "format": "gritic_route_catalogue",
        "hash_algorithm": "ordered Major/Minor NetworkX Weisfeiler-Lehman hashes, combined with MD5",
        "node_order": "GRITIC reference graph insertion order",
        "state_count": len(records),
        "route_count": sum(len(record["routes"]) for record in records),
        "states": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify byte-for-byte reproducibility without writing")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    from gritic.route_catalogue.catalogue import _read_catalogue

    document = build_catalogue()
    _read_catalogue(document)
    payload = json.dumps(document, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    compressed = gzip.compress(payload, compresslevel=9, mtime=0)
    destination = root / "gritic" / "route_catalogue" / "data" / "routes.json.gz"
    if args.check:
        if destination.read_bytes() != compressed:
            raise SystemExit("Generated route catalogue differs from committed data")
        print(f"Catalogue matches: {document['route_count']} routes")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(compressed)
        print(f"Saved {destination}: {document['route_count']} routes, {len(compressed):,} bytes")


if __name__ == "__main__":
    main()
