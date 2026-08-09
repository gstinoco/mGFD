"""
CloudGenerator.core.point_generation.relaxation — Node Relaxation

Overview:
    This module implements Lloyd's relaxation algorithm using Voronoi diagrams
    to regularize point clouds, transforming random distributions into more
    uniform, honeycomb-like structures.

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
import numba as nb                                                                                                                      # JIT compilation.                                                                                                                      # Core numerical operations.

from scipy.spatial import Voronoi                                                                                                       # Voronoi diagram computation.
from shapely.geometry import Polygon                                                                                                    # Geometric objects and operations.

from mGFD.cloud_generator.core.point_generation.geometry import create_fast_polygon_checker                                             # Fast polygon checker.
from mGFD.cloud_generator.core.point_generation.jit_geometry import _is_point_in_polygon_jit                                            # Fast JIT geometric operations.                                             # Fast polygon checker.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Decorator for JIT compilation.
def _compute_relaxation_step_jit(int_pts: np.ndarray, point_region: np.ndarray, flat_regions: np.ndarray, offsets: np.ndarray, vertices: np.ndarray, poly_coords: np.ndarray) -> tuple:
    """
    _compute_relaxation_step_jit
    Numba JIT-compiled helper to compute centroids and update node positions efficiently.
    
    Input:
        int_pts         m x 2 ndarray   Current interior points.
        point_region    m ndarray       Region index for each point.
        flat_regions    n ndarray       Flattened CSR-like array of region vertices indices.
        offsets         k ndarray       Offsets for each region in the flat array.
        vertices        v x 2 ndarray   Voronoi vertices coordinates.
        poly_coords     p x 2 ndarray   Domain polygon coordinates for inclusion check.
        
    Output:
        new_int_pts     m x 2 ndarray   Updated interior points.
        max_movement    float           Maximum displacement distance.
    """
    n            = len(int_pts)                                                                                                         # Number of interior points.
    new_int_pts  = np.zeros_like(int_pts)                                                                                               # Allocate updated points array.
    movements    = np.zeros(n, dtype=np.float64)                                                                                        # Allocate array for individual movements.

    for i in nb.prange(n):                                                                                                              # type: ignore
        reg_idx = point_region[i]                                                                                                       # Get Voronoi region index.
        start   = offsets[reg_idx]                                                                                                      # Start index in flattened array.
        end     = offsets[reg_idx + 1]                                                                                                  # End index in flattened array.
        
        if end - start == 0:                                                                                                            # If region is empty.
            new_int_pts[i] = int_pts[i]                                                                                                 # Keep original position.
        else:                                                                                                                           # If region has vertices.
            is_unbounded = False                                                                                                        # Flag for unbounded region.
            sum_x        = 0.0                                                                                                          # Accumulator for X.
            sum_y        = 0.0                                                                                                          # Accumulator for Y.
            count        = 0                                                                                                            # Number of vertices in region.
            
            for j in range(start, end):                                                                                                 # Iterate over region vertices.
                v_idx = flat_regions[j]                                                                                                 # Vertex index.
                
                if v_idx == -1:                                                                                                         # Check for unbounded vertex (-1).
                    is_unbounded = True                                                                                                 # Set flag.
                    break                                                                                                               # Stop processing this region.
                    
                sum_x += vertices[v_idx, 0]                                                                                             # Accumulate X.
                sum_y += vertices[v_idx, 1]                                                                                             # Accumulate Y.
                count += 1                                                                                                              # Increment counter.
                
            if is_unbounded or count == 0:                                                                                              # If region is invalid or unbounded.
                new_int_pts[i] = int_pts[i]                                                                                             # Keep original position.
            else:                                                                                                                       # If region is valid.
                cx = sum_x / count                                                                                                      # Compute centroid X.
                cy = sum_y / count                                                                                                      # Compute centroid Y.
                
                if _is_point_in_polygon_jit(cx, cy, poly_coords):                                                                       # Check if centroid is inside polygon.
                    new_int_pts[i, 0] = cx                                                                                              # Update point X.
                    new_int_pts[i, 1] = cy                                                                                              # Update point Y.
                    dist              = np.sqrt((cx - int_pts[i, 0])**2 + (cy - int_pts[i, 1])**2)                                      # Calculate movement.
                    movements[i]      = dist                                                                                            # Store movement for later max reduction.
                else:                                                                                                                   # If centroid is outside polygon.
                    new_int_pts[i] = int_pts[i]                                                                                         # Keep original position.
                    
    max_movement = np.max(movements)                                                                                                    # Safely extract maximum movement outside prange.
    return new_int_pts, max_movement                                                                                                    # Return updated points and tracker.

def lloyd_relaxation(interior_points: np.ndarray, boundary_points: np.ndarray, polygon: Polygon, iterations: int = 5, tolerance: float = 1e-4) -> np.ndarray:
    """
    lloyd_relaxation
    Apply Lloyd's relaxation using Voronoi diagrams to regularize the point cloud into a honeycomb-like pattern.
    Only interior points are moved; boundary points are fixed.
    
    Input:
        interior_points     m x 2 ndarray       Array of interior (x, y) points.
        boundary_points     m x 2 ndarray       Array of boundary (x, y) points.
        polygon             Polygon             Shapely polygon to restrict the points.
        iterations          int                 Number of relaxation steps. Default is 5.
        tolerance           float               Minimum movement threshold for early stopping. Default is 1e-4.
    
    Output:
        int_pts             m x 2 ndarray       Relaxed interior points.
    """
    if len(interior_points) == 0:                                                                                                       # Ensure there are points to relax.
        return interior_points                                                                                                          # Return unchanged if empty.

    # 1. Add dummy points                                                                                                               #
    minx, miny, maxx, maxy = polygon.bounds                                                                                             # Extract domain bounds.
    dx, dy                 = maxx - minx, maxy - miny                                                                                   # Calculate width and height offsets.
    dummy_points           = np.array([                                                                                                 # Define the four far-away points.
        [minx - dx, miny - dy], [maxx + dx, miny - dy],                                                                                 # Bottom-left and Bottom-right.
        [maxx + dx, maxy + dy], [minx - dx, maxy + dy]                                                                                  # Top-right and Top-left.
    ])                                                                                                                                  #
    
    int_pts     = np.array(interior_points, dtype=np.float64)                                                                           # Copy interior points array.
    bnd_pts     = np.array(boundary_points, dtype=np.float64) if len(boundary_points) > 0 else np.empty((0, 2))                         # Ensure boundary is a valid array.
    poly_coords = np.array(polygon.exterior.coords, dtype=np.float64)                                                                   # Extract polygon exterior coordinates.
    
    for step in range(iterations):                                                                                                      # Relaxation loop over iterations.
        # 2. Combine all points                                                                                                         #
        pts = np.vstack([int_pts, bnd_pts, dummy_points]) if len(bnd_pts) > 0 else np.vstack([int_pts, dummy_points])                   # Merge all coordinates for Voronoi diagram.
        
        try:                                                                                                                            # Attempt Voronoi calculation.
            vor = Voronoi(pts)                                                                                                          # Compute the Voronoi diagram.
        except Exception as e:                                                                                                          # If computation fails.
            break                                                                                                                       # Stop relaxation process.
            
        # 3. Flatten regions for JIT processing                                                                                         #
        region_lens = np.array([len(r) for r in vor.regions], dtype=np.int32)                                                           # Compute length of each region list.
        offsets     = np.zeros(len(vor.regions) + 1, dtype=np.int32)                                                                    # Allocate offsets array.
        offsets[1:] = np.cumsum(region_lens)                                                                                            # Calculate offsets via cumulative sum.
        
        if len(vor.regions) > 0:                                                                                                        # If regions exist.
            flat_regions = np.concatenate([np.array(r, dtype=np.int32) for r in vor.regions])                                           # Flatten regions list.
        else:                                                                                                                           # If no regions.
            flat_regions = np.empty(0, dtype=np.int32)                                                                                  # Create empty array.
            
        point_region = np.array(vor.point_region[:len(int_pts)], dtype=np.int32)                                                        # Slice point regions for interior points only.
        
        # 4. JIT Accelerated Centroid Calculation                                                                                       #
        int_pts, max_movement = _compute_relaxation_step_jit(int_pts, point_region, flat_regions, offsets, vor.vertices, poly_coords)   # Delegate to Numba helper.
        
        if max_movement < tolerance:                                                                                                    # Check early termination condition.
            break                                                                                                                       # Exit relaxation loop early.

    return int_pts                                                                                                                      # Return the relaxed points.
