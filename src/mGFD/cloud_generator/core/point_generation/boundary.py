"""
CloudGenerator.core.point_generation.boundary — Boundary Generation

Overview:
    This module implements boundary generation algorithms to create nodes
    along the perimeter of geometric regions.

Public API:
    generate_boundary_points    Generate nodes along domain boundary.

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
import numba as nb                                                                                                                      # JIT compilation. # Core numerical operations.

from typing import List, Tuple                                                                                                          # Type hinting.

@nb.njit(cache=True, fastmath=True)                                                                                                     # Decorator for JIT compilation.
def _generate_boundary_points_jit(contour: np.ndarray, cloud_size: float) -> np.ndarray:
    """
    _generate_boundary_points_jit
    Numba JIT-compiled helper to interpolate points along a boundary contour.
    
    Input:
        contour         n x 2 ndarray   Coordinates defining the boundary.
        cloud_size      float           Desired spacing between points.
        
    Output:
        points          m x 2 ndarray   Array of (x, y) interpolated coordinates.
    """
    n_pts    = contour.shape[0]                                                                                                         # Number of contour points.
    capacity = n_pts * 10                                                                                                               # Initial capacity for dynamic array.
    out      = np.zeros((capacity, 2), dtype=np.float64)                                                                                # Allocate output array.
    count    = 0                                                                                                                        # Counter for output points.

    for i in range(n_pts - 1):                                                                                                          # Iterate over contour segments.
        x1, y1   = contour[i, 0], contour[i, 1]                                                                                         # Start point of the segment.
        x2, y2   = contour[i + 1, 0], contour[i + 1, 1]                                                                                 # End point of the segment.
        
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                                                                                 # Calculate segment length.
        
        if distance > 0:                                                                                                                # If segment has positive length.
            num_points = int(np.ceil(distance / cloud_size))                                                                            # Calculate number of points for segment.
            
            if num_points < 1:                                                                                                          # Ensure at least one point.
                num_points = 1                                                                                                          # Fallback to 1 point.
            
            for j in range(num_points):                                                                                                 # Interpolate points along the segment.
                if count >= capacity:                                                                                                   # If we exceed capacity.
                    capacity *= 2                                                                                                       # Double capacity.
                    new_out   = np.zeros((capacity, 2), dtype=np.float64)                                                               # Allocate new array.
                    
                    for k in range(count):                                                                                              # Copy existing data.
                        new_out[k, 0] = out[k, 0]                                                                                       # Copy x.
                        new_out[k, 1] = out[k, 1]                                                                                       # Copy y.
                    
                    out = new_out                                                                                                       # Replace old array.
                
                t             = j / num_points                                                                                          # Interpolation parameter.
                out[count, 0] = x1 + t * (x2 - x1)                                                                                      # Interpolated x coordinate.
                out[count, 1] = y1 + t * (y2 - y1)                                                                                      # Interpolated y coordinate.
                count        += 1                                                                                                       # Increment counter.
                
    return out[:count]                                                                                                                  # Return valid slice of output array.

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
    
    contour_arr     = np.array(contour, dtype=np.float64)                                                                               # Convert contour to NumPy array.
    boundary_points = _generate_boundary_points_jit(contour_arr, cloud_size)                                                            # Call JIT compiled interpolator.
    
    if len(boundary_points) > 0:                                                                                                        # Check if any points were generated.
        rounded_points    = np.round(boundary_points / (cloud_size * 0.1)) * (cloud_size * 0.1)                                         # Round points to avoid duplicates.
        _, unique_indices = np.unique(rounded_points, axis=0, return_index=True)                                                        # Find unique indices.
        boundary_points   = boundary_points[unique_indices]                                                                             # Keep only unique points.
        
        return boundary_points                                                                                                          # Return the filtered boundary points.
    
    return np.empty((0, 2))                                                                                                             # Return empty array if no points.
