"""Immutable, precomputed evolutionary route topology without inference imports."""

from .catalogue import (
    CATALOGUE_VERSION,
    IDENTITY_VERSION,
    CatalogueNode,
    CatalogueRoute,
    catalogue_metadata,
    get_route,
    get_routes,
)

__all__ = [
    "CATALOGUE_VERSION", "IDENTITY_VERSION", "CatalogueNode", "CatalogueRoute",
    "catalogue_metadata", "get_route", "get_routes",
]
