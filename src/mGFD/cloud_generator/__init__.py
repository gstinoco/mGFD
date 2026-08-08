"""
CloudGenerator — Point Cloud Generation System for mGFD

Overview:
    This package contains the core modular components for the Cloud Generation system.
    It is organized into specialized modules to separate concerns and improve maintainability.

Public API:
    generate_cloud_natural
    generate_cloud_regular
    reduce_points_by_region
"""

from .core.generator import generate_cloud_natural, generate_cloud_regular
from .core.reduction import reduce_points_by_region

__all__ = [
    "generate_cloud_natural",
    "generate_cloud_regular",
    "reduce_points_by_region",
]
