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

from scipy.spatial import Voronoi                                                                                                       # Voronoi diagram computation.
from shapely.geometry import Polygon                                                                                                    # Geometric objects and operations.

from mGFD.cloud_generator.core.point_generation.geometry import create_fast_polygon_checker                                             # Fast polygon checker.

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

    # 1. Add dummy points
    minx, miny, maxx, maxy = polygon.bounds                                                                                             # Extract domain bounds.
    dx, dy                 = maxx - minx, maxy - miny                                                                                   # Calculate width and height offsets.
    dummy_points           = np.array([                                                                                                 # Define the four far-away points.
        [minx - dx, miny - dy], [maxx + dx, miny - dy],                                                                                 # Bottom-left and Bottom-right.
        [maxx + dx, maxy + dy], [minx - dx, maxy + dy]                                                                                  # Top-right and Top-left.
    ])
    
    int_pts       = np.array(interior_points)                                                                                           # Copy interior points array.
    bnd_pts       = np.array(boundary_points) if len(boundary_points) > 0 else np.empty((0, 2))                                         # Ensure boundary is a valid array.
    
    fast_contains = create_fast_polygon_checker(polygon)                                                                                # Instantiate polygon inclusion checker.

    for step in range(iterations):                                                                                                      # Relaxation loop over iterations.
        # 2. Combine all points
        pts = np.vstack([int_pts, bnd_pts, dummy_points]) if len(bnd_pts) > 0 else np.vstack([int_pts, dummy_points])                   # Merge all coordinates for Voronoi diagram.
        
        try:                                                                                                                            # Attempt Voronoi calculation.
            vor = Voronoi(pts)                                                                                                          # Compute the Voronoi diagram.
        except Exception as e:                                                                                                          # If computation fails.
            break                                                                                                                       # Stop relaxation process.
            
        new_int_pts  = []                                                                                                               # List for updated points.
        max_movement = 0.0                                                                                                              # Tracker for maximum point displacement.
        
        # 3. Update interior points
        for i in range(len(int_pts)):                                                                                                   # Iterate over the movable points.
            region_index = vor.point_region[i]                                                                                          # Get the cell index for this point.
            region       = vor.regions[region_index]                                                                                    # Get the vertices list for this cell.
            
            if -1 in region or len(region) == 0:                                                                                        # If the cell is unbounded or empty.
                new_int_pts.append(int_pts[i])                                                                                          # Do not move this point.
                continue                                                                                                                # Move to next point.
                
            # 4. Get Voronoi vertices
            cell_vertices = vor.vertices[region]                                                                                        # Extract cell geometric vertices.
            
            # 5. Calculate centroid
            centroid = np.mean(cell_vertices, axis=0)                                                                                   # Calculate the center of mass of the cell.
            
            # 6. Ensure centroid is inside polygon
            if fast_contains(centroid[0], centroid[1]):                                                                                 # Verify the centroid is inside the domain.
                new_int_pts.append(centroid)                                                                                            # Update point position to centroid.
                dist = np.sqrt((centroid[0] - int_pts[i][0])**2 + (centroid[1] - int_pts[i][1])**2)                                     # Compute movement distance.
                if dist > max_movement:                                                                                                 # Update maximum movement tracker.
                    max_movement = dist                                                                                                 # Store new max movement.
            else:                                                                                                                       # If centroid wandered out of the polygon.
                new_int_pts.append(int_pts[i])                                                                                          # Keep the point at its previous position.
                
        int_pts = np.array(new_int_pts)                                                                                                 # Update working array for next iteration.
        
        if max_movement < tolerance:                                                                                                    # Check early termination condition.
            break                                                                                                                       # Exit relaxation loop early.

    return int_pts                                                                                                                      # Return the relaxed points.
