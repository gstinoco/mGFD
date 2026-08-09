"""
CloudGenerator.core.point_generation.boundary — Boundary Generation

Overview:
    This module implements boundary generation algorithms to create nodes
    along the perimeter of geometric regions.

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

from typing import List, Tuple                                                                                                          # Type hinting.

def generate_boundary_points(contour: List[Tuple[float, float]], cloud_size: float) -> np.ndarray:
    """
    generate_boundary_points
    Generate points along the boundary of a contour.
    
    Input:
        contour             List                List of (x, y) coordinate tuples defining the boundary.
        cloud_size          float               Desired spacing between points.
    
    Output:
        boundary_points     m x 2 ndarray       Array of (x, y) coordinates of boundary points.
    """
    if len(contour) < 2:                                                                                                                # Check if contour has at least two points.
        return np.empty((0, 2))                                                                                                         # Return empty array if not.
    
    boundary_points = []                                                                                                                # Initialize list for boundary points.
    
    for i in range(len(contour) - 1):                                                                                                   # Iterate over contour segments.
        x1, y1 = contour[i]                                                                                                             # Start point of the segment.
        x2, y2 = contour[i + 1]                                                                                                         # End point of the segment.
        
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                                                                                 # Calculate segment length.
        
        if distance > 0:                                                                                                                # If segment has positive length.
            # Use ceil to ensure we don't under-sample the boundary
            # This prevents large gaps when distance is slightly less than a multiple of cloud_size
            num_points = max(1, int(np.ceil(distance / cloud_size)))                                                                    # Calculate number of points for segment.
            
            for j in range(num_points):                                                                                                 # Interpolate points along the segment.
                t = j / num_points                                                                                                      # Interpolation parameter.
                x = x1 + t * (x2 - x1)                                                                                                  # Interpolated x coordinate.
                y = y1 + t * (y2 - y1)                                                                                                  # Interpolated y coordinate.
                boundary_points.append([x, y])                                                                                          # Add point to the list.
    
    if boundary_points:                                                                                                                 # Check if any points were generated.
        boundary_points   = np.array(boundary_points)                                                                                   # Convert list to NumPy array.
        rounded_points    = np.round(boundary_points / (cloud_size * 0.1)) * (cloud_size * 0.1)                                         # Round points to avoid duplicates due to float precision.
        _, unique_indices = np.unique(rounded_points, axis=0, return_index=True)                                                        # Find unique indices.
        boundary_points   = boundary_points[unique_indices]                                                                             # Keep only unique points.
        
        return boundary_points                                                                                                          # Return the filtered boundary points.
    
    return np.empty((0, 2))                                                                                                             # Return empty array if no points.
