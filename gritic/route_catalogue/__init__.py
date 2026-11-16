"""Immutable, precomputed evolutionary route topology without inference imports."""

from .catalogue import (
    CatalogueNode,
    CatalogueRoute,
    catalogue_metadata,
    get_route,
    get_routes,
)

__all__ = [
    "CatalogueNode", "CatalogueRoute",
    "catalogue_metadata", "get_route", "get_routes",
]
