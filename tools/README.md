# Route catalogue

`gritic.route_catalogue` supplies immutable evolutionary route topologies.
The module and its compressed data are included in the GRITIC distribution.
Reading the catalogue loads only standard-library modules; it does not import
GRITIC's inference, graph, numerical, or plotting libraries.

The catalogue contains **2,796 routes across 31 copy-number states**.
It covers GRITIC's permitted timing states, plus the `2+0`, `2+1`, and `2+2`
WGD states used for route construction and pooled-WGD models. Topologies are
stored once with all valid WGD annotations. Lookup selects the requested
WGD status and route restriction.

## Read a route

```python
from gritic.route_catalogue import catalogue_metadata, get_route, get_routes

routes = get_routes(3, 1, True)
route = get_route(3, 1, True, routes[0].short_id)
for node in route.nodes:
    print(node.node_id, node.parent_id, node.allele, node.wgd, node.terminal)
```

`get_routes(major_cn, minor_cn, wgd_status, wgd_trees_status="Default")`
returns a tuple in the original GRITIC generation order. The supported modes
are `Default`, `No_WGD`, `Only_WGD`, and `Only_WGD_2+2`. The latter two require
at least one and at least two WGD nodes, respectively. As in GRITIC's reference
enumerator, restrictions apply only when `wgd_status` is true.

`get_route(major_cn, minor_cn, wgd_status, route_id)` resolves an exact full
hash or the nine-character `Route` identifier currently written by GRITIC.
Lookup validates the model state and requires a unique match. Unsupported
states, unknown modes, and unmatched identifiers raise `ValueError`; lookup
does not generate candidates or guess an alternative identity convention.

Each frozen `CatalogueRoute` contains `route_id`, `major_cn`, `minor_cn`,
`nodes`, and `mirror_route_id`, plus its derived `short_id` property. Each
frozen `CatalogueNode` contains `node_id`, `parent_id`, `allele`, `wgd`, and
`terminal`. Node numbers and order are the reference generator's values, so
they align with saved gain timing tables and `Timing_Node_ID` arrays. Mirror
IDs preserve the ordered Major/Minor relationship of balanced states.

The catalogue contains structural information only. Posterior timings,
intervals, route probabilities, sample WGD summaries, and model/target
metadata come from each sample's output files. In particular, consumers of
pooled-WGD archives must use their fitted-model metadata when resolving a
route; target copy numbers may differ from the fitted geometry.

## Identity and reproducibility

The data format is `gritic_route_catalogue`, stored in `data/routes.json.gz`.
The Major and Minor components receive directed Weisfeiler-Lehman hashes,
which are combined in allele order with MD5 and the fixed `ordered-alleles`
namespace. Consumers load these frozen identifiers directly. The catalogue
has no separate schema or identity version fields.

`catalogue_metadata()` describes coverage, hash construction, and node order.
There are no alternate hash conventions or compatibility tables in the runtime.

Generate data from the suite root using the existing environment:

```bash
.venv/bin/python gritic/tools/build_route_catalogue.py
.venv/bin/python gritic/tools/build_route_catalogue.py --check
.venv/bin/python -m pytest gritic/tests/test_route_catalogue_data.py
```

The generator writes a deterministic compressed JSON resource. It checks
structural validity, full and shortened ID uniqueness within each state,
and reciprocal mirror mappings. Tests compare every supported state and
restriction against the reference enumerator, including node IDs, edges,
attributes, hash values, mirror IDs, and ordering. Additional tests verify
runtime isolation from inference and graph libraries.
