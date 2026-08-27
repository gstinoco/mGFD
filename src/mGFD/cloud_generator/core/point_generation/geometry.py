"""
CloudGenerator.core.point_generation.geometry — Geometry Utilities

Overview:
    This module provides fast geometric operations for the point generation,
    specifically a fast point-in-polygon checker.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx
    
    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México.
        Aula CIMNE-Morelia. México.
        SIIIA-MATH: Soluciones de Ingeniería. México.

    Based on the theoretical concepts presented in:
        "mGFD: A meshless generalized finite difference method",
        Gerardo Tinoco-Guerrero, Francisco Javier Domínguez-Mota, José Alberto Guzmán-Torres, 
        Gabriela Pedraza-Jiménez, José Gerardo Tinoco-Ruiz,
        Computers & Mathematics with Applications, Volume 195 (2025) 396-418.
        https://doi.org/10.1016/j.camwa.2025.07.034

Date:
    March, 2026.
Last Modification:
    August, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Callable                                                                                                             # Type hinting.
from matplotlib.path import Path                                                                                                        # Point-in-polygon checking fallback.
from shapely.geometry import Polygon                                                                                                    # Geometric objects and operations.

def create_fast_polygon_checker(polygon: Polygon) -> Callable:
    """
    create_fast_polygon_checker
    Returns a fast point-in-polygon checker function using matplotlib Path.
    
    Input:
        polygon         Polygon         Shapely Polygon.
    
    Output:
        contains        Callable        Checker function.
    """
    ext_path   = Path(np.array(polygon.exterior.coords))                                                                                # Create a Path from the polygon exterior.
    hole_paths = [Path(np.array(interior.coords)) for interior in polygon.interiors]                                                    # Create Paths for all the polygon holes.
    
    def contains(x: float, y: float) -> bool:                                                                                           # Define the inner checker function.
        """
        contains
        
        Checks if a given coordinate (x, y) is inside the polygon using JIT ray casting.
        
        Input:
            x                           float           X coordinate of the test point.
            y                           float           Y coordinate of the test point.
            
        Output:
            is_inside                   bool            True if the point is strictly inside the domain, False otherwise.
        """
        pt = (x, y)                                                                                                                     # Create a coordinate tuple.
        if not ext_path.contains_point(pt): return False                                                                                # Return False if outside the exterior.
        for hp in hole_paths:                                                                                                           # Iterate through all the hole paths.
            if hp.contains_point(pt): return False                                                                                      # Return False if inside any hole.
        return True                                                                                                                     # Return True if inside exterior and outside holes.
    
    return contains                                                                                                                     # Return the checker function.
