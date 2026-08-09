"""
CloudGenerator.core.point_generation.jit_geometry — Numba JIT Geometric Helpers

Overview:
    This module provides highly optimized, Numba JIT-compiled geometric algorithms.
    It includes a ray-casting point-in-polygon checker, bypassing Python overhead
    and the inability to compile external libraries like Shapely or matplotlib.path
    directly inside Numba routines.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx
    
    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI.
        Coordination of Scientific Research, CIC-UMSNH.
        Aula CIMNE-Morelia.
        SIIIA-MATH: Soluciones de Ingeniería.

Date:
    August, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.
import numba as nb                                                                                                                      # JIT compilation.

@nb.njit(cache=True, fastmath=True)                                                                                                     # Decorator for JIT compilation.
def _is_point_in_polygon_jit(px: float, py: float, poly: np.ndarray) -> bool:
    """
    _is_point_in_polygon_jit
    Numba JIT-compiled ray-casting algorithm to determine if a point is inside a polygon.
    
    Input:
        px          float           X coordinate of the point.
        py          float           Y coordinate of the point.
        poly        n x 2 ndarray   Coordinates defining the polygon boundary.
        
    Output:
        inside      bool            True if the point is strictly inside, False otherwise.
    """
    inside = False                                                                                                                      # Initialize inclusion flag.
    n      = poly.shape[0]                                                                                                              # Number of vertices.
    j      = n - 1                                                                                                                      # Start with the last vertex.

    for i in range(n):                                                                                                                  # Loop over all polygon edges.
        xi, yi = poly[i, 0], poly[i, 1]                                                                                                 # Current vertex.
        xj, yj = poly[j, 0], poly[j, 1]                                                                                                 # Previous vertex.

        intersect = ((yi > py) != (yj > py)) and (px < (xj - xi) * (py - yi) / (yj - yi) + xi)                                          # Ray-casting condition.
        if intersect:                                                                                                                   # If the ray crosses the edge.
            inside = not inside                                                                                                         # Toggle the inside flag.

        j = i                                                                                                                           # Move to the next edge.

    return inside                                                                                                                       # Return final classification.
