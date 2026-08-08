"""
CloudGenerator.point_generation — Point Generation

Overview:
    This module implements the core algorithms for generating point clouds using
    both Regular (Grid-based) and Natural (Poisson Disk Sampling) distributions.
    It handles the complex logic of boundary and interior point generation for
    multi-region geometries.

Data conventions:
    None.

Public API:
    create_fast_polygon_checker                     Returns a fast point-in-polygon checker function.
    generate_boundary_points                        Generate points along the boundary of a contour.
    generate_interior_points                        Generate points inside a polygon using a vectorized grid-based approach.
    poisson_disk_sampling                           Generate points using Poisson Disk Sampling algorithm.
    generate_interior_points_poisson                Generate interior points using Poisson Disk Sampling.
    generate_region_cloud                           Generate a complete cloud for a single region.
    generate_region_cloud_with_uniform_density      Generate a complete cloud for a single region with uniform density.
    lloyd_relaxation                                Apply Lloyd's relaxation using Voronoi diagrams.
    generate_region_cloud_poisson                   Generate a complete cloud for a single region using Poisson Disk Sampling.
    generate_region_cloud_with_holes                Generate a cloud for the main region considering interior regions as holes.
    generate_region_cloud_with_holes_poisson        Generate a cloud with holes using Poisson Disk Sampling.
    generate_interior_regions_clouds                Generate clouds for interior regions.
    generate_interior_regions_clouds_poisson        Generate clouds for interior regions using Poisson Disk Sampling.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
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

import numpy as np                                                                  # Numerical arrays and mathematical operations.
from shapely.geometry import Point, Polygon                                         # Geometric objects and operations.
import cv2                                                                          # OpenCV for fast rasterization and masking.
import logging                                                                      # Event logging.
import random                                                                       # Random number generation.
from scipy.spatial import Voronoi                                                   # Voronoi diagram computation.
import os                                                                           # Operating system interfaces.
from concurrent.futures import ProcessPoolExecutor                                  # Parallel execution.
from matplotlib.path import Path                                                    # Point-in-polygon checking fallback.
from .utils import calculate_cloud_size, create_closed_contour                      # Utility functions for cloud generation.

def create_fast_polygon_checker(polygon: Polygon):
    """
    Returns a fast point-in-polygon checker function using matplotlib Path.
    
    Input:
        polygon         Polygon         Shapely Polygon.
    
    Output:
        contains        function        Checker function.
    """
    ext_path = Path(np.array(polygon.exterior.coords))                              # Create a Path from the polygon exterior.
    hole_paths = [Path(np.array(interior.coords)) for interior in polygon.interiors]# Create Paths for all the polygon holes.
    
    def contains(x, y):                                                             # Define the inner checker function.
        pt = (x, y)                                                                 # Create a coordinate tuple.
        if not ext_path.contains_point(pt): return False                            # Return False if outside the exterior.
        for hp in hole_paths:                                                       # Iterate through all the hole paths.
            if hp.contains_point(pt): return False                                  # Return False if inside any hole.
        return True                                                                 # Return True if inside exterior and outside holes.
    return contains                                                                 # Return the checker function.

def generate_boundary_points(contour: list[tuple[float, float]], cloud_size: float) -> np.ndarray:
    """
    Generate points along the boundary of a contour.
    
    Input:
        contour             list        List of (x, y) coordinate tuples defining the boundary.
        cloud_size          float       Desired spacing between points.
    
    Output:
        boundary_points     ndarray     Array of (x, y) coordinates of boundary points.
    """
    if len(contour) < 2:                                                            # Check if contour has at least two points.
        return np.empty((0, 2))                                                     # Return empty array if not.
    
    boundary_points = []                                                            # Initialize list for boundary points.
    
    for i in range(len(contour) - 1):                                               # Iterate over contour segments.
        x1, y1 = contour[i]                                                         # Start point of the segment.
        x2, y2 = contour[i + 1]                                                     # End point of the segment.
        
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                             # Calculate segment length.
        
        if distance > 0:                                                            # If segment has positive length.
            # Use ceil to ensure we don't under-sample the boundary
            # This prevents large gaps when distance is slightly less than a multiple of cloud_size
            num_points = max(1, int(np.ceil(distance / cloud_size)))                # Calculate number of points for segment.
            
            for j in range(num_points):                                             # Interpolate points along the segment.
                t = j / num_points                                                  # Interpolation parameter.
                x = x1 + t * (x2 - x1)                                              # Interpolated x coordinate.
                y = y1 + t * (y2 - y1)                                              # Interpolated y coordinate.
                boundary_points.append([x, y])                                      # Add point to the list.
    
    if boundary_points:                                                             # Check if any points were generated.
        boundary_points = np.array(boundary_points)                                 # Convert list to NumPy array.
        rounded_points = np.round(boundary_points / (cloud_size * 0.1)) * (cloud_size * 0.1) # Round points to avoid duplicates due to float precision.
        _, unique_indices = np.unique(rounded_points, axis=0, return_index=True)    # Find unique indices.
        boundary_points = boundary_points[unique_indices]                           # Keep only unique points.
        return boundary_points                                                      # Return the filtered boundary points.
    
    return np.empty((0, 2))                                                         # Return empty array if no points.

def generate_interior_points(polygon: Polygon, cloud_size: float) -> np.ndarray:
    """
    Generate points inside a polygon using a vectorized grid-based approach.
    
    Optimized using OpenCV mask for fast point-in-polygon checks (replacing Matplotlib Path).
    
    Input:
        polygon             Polygon     The polygon to fill with points.
        cloud_size          float       Desired spacing between points.
    
    Output:
        points              ndarray     Array of (x, y) coordinates of generated points.
    """
    try:                                                                            # Start of the main generation block.
        bounds = polygon.bounds                                                     # Get bounding box of the polygon.
        x_min, y_min, x_max, y_max = bounds                                         # Extract bounding box limits.
        
        # Create grid of points
        x_range = np.arange(x_min, x_max + cloud_size, cloud_size)                  # Create range for x coordinates.
        y_range = np.arange(y_min, y_max + cloud_size, cloud_size)                  # Create range for y coordinates.
        
        if len(x_range) == 0 or len(y_range) == 0:                                  # Check if grids are valid.
            return np.empty((0, 2))                                                 # Return empty array if not.
            
        xx, yy = np.meshgrid(x_range, y_range)                                      # Create 2D mesh grid.
        grid_points = np.vstack((xx.ravel(), yy.ravel())).T                         # Flatten and stack into point list.
        
        # Create a mask using OpenCV
        width = len(x_range)                                                        # Mask width in pixels.
        height = len(y_range)                                                       # Mask height in pixels.
        mask = np.zeros((height, width), dtype=np.uint8)                            # Initialize the mask array.
        
        # Function to convert world coords to pixel coords
        def to_pixel_coords(coords):                                                # Helper for coordinate transformation.
            pixel_coords = []                                                       # List for pixel coordinates.
            for x, y in coords:                                                     # Iterate over world coordinates.
                px = int(round((x - x_min) / cloud_size))                           # Map x to pixel index.
                py = int(round((y - y_min) / cloud_size))                           # Map y to pixel index.
                pixel_coords.append([px, py])                                       # Add to pixel list.
            return np.array(pixel_coords, dtype=np.int32)                           # Return as integer NumPy array.

        # Draw exterior
        ext_coords = list(polygon.exterior.coords)                                  # Extract exterior coordinates.
        ext_pixels = to_pixel_coords(ext_coords)                                    # Convert to pixel coordinates.
        cv2.fillPoly(mask, [ext_pixels], 1)                                         # Fill the exterior polygon in the mask.
        
        # Exclude exterior boundary (set to 0) to avoid overlapping with explicit boundary nodes
        cv2.polylines(mask, [ext_pixels], isClosed=True, color=0, thickness=1)      # Clear boundary pixels in the mask.
        
        # Draw interiors (holes)
        for interior in polygon.interiors:                                          # Iterate over hole geometries.
            int_coords = list(interior.coords)                                      # Extract hole coordinates.
            int_pixels = to_pixel_coords(int_coords)                                # Convert to pixel coordinates.
            cv2.fillPoly(mask, [int_pixels], 0)                                     # Clear hole pixels in the mask.
            
        # Select points
        selected_mask = mask.astype(bool).ravel()                                   # Flatten mask and convert to boolean.
        points = grid_points[selected_mask]                                         # Filter grid points using the mask.
        
        return points if len(points) > 0 else np.empty((0, 2))                      # Return points or empty array.
        
    except Exception as e:                                                          # Fallback on OpenCV error.
        logging.error(f"Error in vectorized point generation: {e}")                 # Log the error.
        # Fallback to Path based iterative method
        logging.info("Falling back to Path iterative method")                       # Log fallback action.
        points = []                                                                 # Initialize point list.
        # Re-generate grid points if needed (though they are local in try block)
        x_range = np.arange(x_min, x_max + cloud_size, cloud_size)                  # Re-create x range.
        y_range = np.arange(y_min, y_max + cloud_size, cloud_size)                  # Re-create y range.
        fast_contains = create_fast_polygon_checker(polygon)                        # Create path-based checker.
        
        for x in x_range:                                                           # Iterate over x grid.
            for y in y_range:                                                       # Iterate over y grid.
                if fast_contains(x, y):                                             # Check if point is inside.
                    points.append([x, y])                                           # Add to list if inside.
        return np.array(points) if points else np.empty((0, 2))                     # Return array of points.

from scipy.stats.qmc import PoissonDisk                                             # Import SciPy Poisson Disk generator.
from scipy.spatial import cKDTree                                                   # Import KDTree for distance checks.

def poisson_disk_sampling(polygon: Polygon, radius: float, k: int = 30, boundary_points: list[tuple[float, float]] = None) -> list[tuple[float, float]]:
    """
    Generate points using Poisson Disk Sampling algorithm via SciPy QMC engine for extreme speed.
    
    Input:
        polygon             Polygon     Target polygon to fill with points.
        radius              float       Minimum distance between any two generated points.
        k                   int         Ignored in SciPy version, kept for API compatibility. Default is 30.
        boundary_points     list        Pre-existing boundary points to check distances against.
    
    Output:
        final_points        list        List of (x, y) coordinate tuples representing generated points.
    """
    minx, miny, maxx, maxy = polygon.bounds                                         # Get bounding box limits.
    
    # SciPy PoissonDisk operates in [0, 1] space.
    # We must scale our dimensions and the radius down to [0, 1].
    width = maxx - minx                                                             # Calculate domain width.
    height = maxy - miny                                                            # Calculate domain height.
    
    # If the polygon is just a point or a line, we can't fill it
    if width <= 0 or height <= 0:                                                   # Ensure domain is 2D.
        return []                                                                   # Return empty list if not.
        
    scale_factor = max(width, height)                                               # Determine scale factor for normalization.
    scaled_radius = radius / scale_factor                                           # Normalize the sampling radius.
    
    try:                                                                            # Try sampling generation.
        # Initialize SciPy engine
        engine = PoissonDisk(d=2, radius=scaled_radius)                             # Create Poisson generator engine.
        # Generate samples in [0, 1]^2
        samples = engine.fill_space()                                               # Generate the normalized points.
    except Exception as e:                                                          # Fallback if sampling fails.
        logging.error(f"SciPy PoissonDisk failed: {e}. Falling back to random uniform points.") # Log failure.
        # Very rough fallback if radius is too large for the bounding box
        return []                                                                   # Return empty list.

    # Scale back to world coordinates
    samples[:, 0] = samples[:, 0] * scale_factor + minx                             # Scale and shift x coordinates.
    samples[:, 1] = samples[:, 1] * scale_factor + miny                             # Scale and shift y coordinates.
    
    # Filter points inside polygon using vectorized OpenCV mask for extreme speed
    try:                                                                            # Attempt fast masking.
        # Create a mask with a resolution high enough to capture the polygon accurately
        resolution = max(radius * 0.1, (maxx - minx) / 2000.0) # at most 2000px wide  # Calculate optimal pixel resolution.
        if resolution <= 0: resolution = 1.0                                        # Prevent zero or negative resolution.
        
        width_px = int(np.ceil((maxx - minx) / resolution)) + 1                     # Calculate mask width in pixels.
        height_px = int(np.ceil((maxy - miny) / resolution)) + 1                    # Calculate mask height in pixels.
        
        mask = np.zeros((height_px, width_px), dtype=np.uint8)                      # Initialize the OpenCV mask array.
        
        def to_pixel_coords(coords):                                                # Coordinate transformation helper.
            pixel_coords = []                                                       # List to store pixel coordinates.
            for x, y in coords:                                                     # Iterate over world coordinates.
                px = int(round((x - minx) / resolution))                            # Convert x to pixel index.
                py = int(round((y - miny) / resolution))                            # Convert y to pixel index.
                pixel_coords.append([px, py])                                       # Add to list.
            return np.array(pixel_coords, dtype=np.int32)                           # Return integer array.
            
        # Draw exterior
        ext_pixels = to_pixel_coords(list(polygon.exterior.coords))                 # Convert exterior to pixels.
        cv2.fillPoly(mask, [ext_pixels], 1)                                         # Fill the mask with exterior polygon.
        # Exclude exterior boundary to avoid overlapping
        cv2.polylines(mask, [ext_pixels], isClosed=True, color=0, thickness=1)      # Unset pixels exactly on the boundary.
        
        # Draw holes
        for interior in polygon.interiors:                                          # Iterate over internal holes.
            int_pixels = to_pixel_coords(list(interior.coords))                     # Convert hole to pixels.
            cv2.fillPoly(mask, [int_pixels], 0)                                     # Clear hole region in the mask.
            
        # Map samples to pixel coordinates and check mask
        sample_px = np.round((samples[:, 0] - minx) / resolution).astype(int)       # Get pixel x for samples.
        sample_py = np.round((samples[:, 1] - miny) / resolution).astype(int)       # Get pixel y for samples.
        
        # Filter bounds to avoid IndexError
        valid_idx = (sample_px >= 0) & (sample_px < width_px) & (sample_py >= 0) & (sample_py < height_px) # Mask bounds check.
        
        # Of the valid indices, check the mask
        in_polygon = np.zeros(len(samples), dtype=bool)                             # Initialize boolean inclusion mask.
        in_polygon[valid_idx] = mask[sample_py[valid_idx], sample_px[valid_idx]] > 0# Check pixel values.
        
        interior_points = samples[in_polygon].tolist()                              # Select samples inside the mask.
        interior_points = [(p[0], p[1]) for p in interior_points]                   # Convert to list of tuples.
        
    except Exception as e:                                                          # If fast masking fails.
        logging.error(f"Vectorized masking failed in Poisson: {e}. Falling back to Path checker.") # Log fallback.
        fast_contains = create_fast_polygon_checker(polygon)                        # Create path checker.
        interior_points = [                                                         # List comprehension for fallback filtering.
            (x, y) for x, y in samples                                              # Iterate over samples.
            if fast_contains(x, y)                                                  # Check inclusion.
        ]
    
    if not interior_points:                                                         # Check if any points survived filtering.
        return []                                                                   # Return empty list if not.
        
    if boundary_points and len(boundary_points) > 0:                                # If boundary points exist for distance checking.
        # Filter points that are too close to the boundary points using cKDTree
        boundary_array = np.array(boundary_points)                                  # Convert boundary to NumPy array.
        interior_array = np.array(interior_points)                                  # Convert interior to NumPy array.
        
        # Build KDTree for boundary points
        tree = cKDTree(boundary_array)                                              # Initialize the KDTree.
        
        # Query distances. k=1 gets the closest boundary point to each interior point
        distances, _ = tree.query(interior_array, k=1)                              # Get distances to nearest boundary point.
        
        # Keep only interior points that are at least `radius` away from any boundary
        valid_indices = distances >= (radius * 0.8) # Relaxing the radius slightly to avoid extreme culling near edges # Filter condition.
        
        final_points = interior_array[valid_indices].tolist()                       # Filter array using mask.
        return [(p[0], p[1]) for p in final_points]                                 # Convert and return valid points.
        
    return interior_points                                                          # Return points directly if no boundary checks needed.

def generate_interior_points_poisson(polygon: Polygon, cloud_size: float) -> list[tuple[float, float]]:
    """
    Generate interior points using Poisson Disk Sampling for more natural distribution.
    
    Input:
        polygon             Polygon     Shapely polygon representing the region.
        cloud_size          float       Target spacing between points (used as radius).
    
    Output:
        interior_points     list        List of (x, y) tuples representing interior points.
    """
    # Use cloud_size as the minimum distance between points
    radius = cloud_size * 0.8  # Slightly smaller to get more points                # Set the sampling radius.
    
    try:                                                                            # Start generation.
        interior_points = poisson_disk_sampling(polygon, radius)                    # Perform Poisson sampling.
        logging.info(f"Generated {len(interior_points)} interior points using Poisson Disk Sampling") # Log generation count.
        return interior_points                                                      # Return the interior points.
    except Exception as e:                                                          # On failure.
        logging.error(f"Error in Poisson Disk Sampling, falling back to grid method: {e}") # Log fallback.
        # Fallback to original method
        return generate_interior_points(polygon, cloud_size)                        # Return grid-based points.

def generate_region_cloud(region_points: list[tuple[float, float]]) -> tuple[np.ndarray, np.ndarray, float] | tuple[None, None, None]:
    """
    Generate a complete cloud for a single region.
    This function generates the SAME result every time for the same input.
    
    Input:
        region_points       list        List of region points.
        
    Output:
        boundary_points     ndarray     Boundary points.
        interior_points     ndarray     Interior points.
        cloud_size          float       Cloud size.
    """
    try:                                                                            # Start region generation.
        # Calculate cloud size for this region
        cloud_size = calculate_cloud_size(region_points)                            # Determine appropriate point spacing.
        
        # Create closed contour
        contour = create_closed_contour(region_points)                              # Obtain closed path.
        
        # Generate boundary points
        boundary_points = generate_boundary_points(contour, cloud_size)             # Generate points on the perimeter.
        
        # Create polygon for interior point generation
        polygon = Polygon(contour)                                                  # Instantiate Shapely polygon.
        if not polygon.is_valid:                                                    # Check if polygon topology is valid.
            polygon = polygon.buffer(0)                                             # Attempt to fix self-intersections.
        
        # Generate interior points
        interior_points = generate_interior_points(polygon, cloud_size)             # Generate points inside the region.
        
        logging.info(f"Generated {len(boundary_points)} boundary points and {len(interior_points)} interior points") # Log generated counts.
        
        return boundary_points, interior_points, cloud_size                         # Return all data.
        
    except Exception as e:                                                          # On any failure.
        logging.error(f"Error generating region cloud: {e}")                        # Log the error.
        return None, None, None                                                     # Return None values.

def generate_region_cloud_with_uniform_density(region_points: list[tuple[float, float]], cloud_size: float) -> tuple[np.ndarray, np.ndarray, float] | tuple[None, None, None]:
    """
    Generate a complete cloud for a single region using a specified cloud size.
    This ensures uniform density across different regions.
    
    Input:
        region_points       list        List of region points.
        cloud_size          float       Cloud size to enforce.
        
    Output:
        boundary_points     ndarray     Boundary points.
        interior_points     ndarray     Interior points.
        cloud_size          float       Cloud size.
    """
    try:                                                                            # Start fixed density generation.
        # Create closed contour
        contour = create_closed_contour(region_points)                              # Obtain closed path.
        
        # Generate boundary points using the specified cloud size
        boundary_points = generate_boundary_points(contour, cloud_size)             # Generate points on the perimeter.
        
        # Create polygon for interior point generation
        polygon = Polygon(contour)                                                  # Instantiate Shapely polygon.
        if not polygon.is_valid:                                                    # Check if polygon topology is valid.
            polygon = polygon.buffer(0)                                             # Attempt to fix self-intersections.
        
        # Generate interior points using the specified cloud size
        interior_points = generate_interior_points(polygon, cloud_size)             # Generate points inside the region.
        
        logging.info(f"Generated {len(boundary_points)} boundary points and {len(interior_points)} interior points with uniform density") # Log counts.
        
        return boundary_points, interior_points, cloud_size                         # Return generated data.
        
    except Exception as e:                                                          # On failure.
        logging.error(f"Error generating region cloud with uniform density: {e}")   # Log the error.
        return None, None, None                                                     # Return None values.

def lloyd_relaxation(interior_points: np.ndarray, boundary_points: np.ndarray, polygon: Polygon, iterations: int = 5, tolerance: float = 1e-4) -> np.ndarray:
    """
    Apply Lloyd's relaxation using Voronoi diagrams to regularize the point cloud into a honeycomb-like pattern.
    Only interior points are moved; boundary points are fixed.
    
    Input:
        interior_points     ndarray     Array of interior (x, y) points.
        boundary_points     ndarray     Array of boundary (x, y) points.
        polygon             Polygon     Shapely polygon to restrict the points.
        iterations          int         Number of relaxation steps. Default is 5.
        tolerance           float       Minimum movement threshold for early stopping. Default is 1e-4.
    
    Output:
        int_pts             ndarray     Relaxed interior points.
    """
    if len(interior_points) == 0:                                                   # Ensure there are points to relax.
        return interior_points                                                      # Return unchanged if empty.

    # To avoid boundary issues in Voronoi, add dummy points far away (bounding box corners)
    minx, miny, maxx, maxy = polygon.bounds                                         # Extract domain bounds.
    dx, dy = maxx - minx, maxy - miny                                               # Calculate width and height offsets.
    dummy_points = np.array([                                                       # Define the four far-away points.
        [minx - dx, miny - dy], [maxx + dx, miny - dy],                             # Bottom-left and Bottom-right.
        [maxx + dx, maxy + dy], [minx - dx, maxy + dy]                              # Top-right and Top-left.
    ])
    
    int_pts = np.array(interior_points)                                             # Copy interior points array.
    bnd_pts = np.array(boundary_points) if len(boundary_points) > 0 else np.empty((0, 2)) # Ensure boundary is a valid array.
    
    fast_contains = create_fast_polygon_checker(polygon)                            # Instantiate polygon inclusion checker.

    for step in range(iterations):                                                  # Relaxation loop over iterations.
        # Combine all points: [interior, boundary, dummy]
        pts = np.vstack([int_pts, bnd_pts, dummy_points]) if len(bnd_pts) > 0 else np.vstack([int_pts, dummy_points]) # Merge all coordinates for Voronoi diagram.
        
        try:                                                                        # Attempt Voronoi calculation.
            vor = Voronoi(pts)                                                      # Compute the Voronoi diagram.
        except Exception as e:                                                      # If computation fails.
            logging.warning(f"Voronoi computation failed during Lloyd relaxation: {e}") # Log warning.
            break                                                                   # Stop relaxation process.
            
        new_int_pts = []                                                            # List for updated points.
        max_movement = 0.0                                                          # Tracker for maximum point displacement.
        
        # Update only interior points (indices 0 to len(int_pts)-1)
        for i in range(len(int_pts)):                                               # Iterate over the movable points.
            region_index = vor.point_region[i]                                      # Get the cell index for this point.
            region = vor.regions[region_index]                                      # Get the vertices list for this cell.
            
            if -1 in region or len(region) == 0:                                    # If the cell is unbounded or empty.
                new_int_pts.append(int_pts[i])                                      # Do not move this point.
                continue                                                            # Move to next point.
                
            # Get vertices of this Voronoi region
            cell_vertices = vor.vertices[region]                                    # Extract cell geometric vertices.
            
            # Simple centroid of the polygon formed by vertices
            centroid = np.mean(cell_vertices, axis=0)                               # Calculate the center of mass of the cell.
            
            # Ensure centroid is strictly inside polygon using fast checker
            if fast_contains(centroid[0], centroid[1]):                             # Verify the centroid is inside the domain.
                new_int_pts.append(centroid)                                        # Update point position to centroid.
                dist = np.sqrt((centroid[0] - int_pts[i][0])**2 + (centroid[1] - int_pts[i][1])**2) # Compute movement distance.
                if dist > max_movement:                                             # Update maximum movement tracker.
                    max_movement = dist                                             # Store new max movement.
            else:                                                                   # If centroid wandered out of the polygon.
                new_int_pts.append(int_pts[i]) # keep original                      # Keep the point at its previous position.
                
        int_pts = np.array(new_int_pts)                                             # Update working array for next iteration.
        
        if max_movement < tolerance:                                                # Check early termination condition.
            logging.info(f"Lloyd relaxation converged early after {step+1} iterations (max movement: {max_movement:.6f})") # Log convergence.
            break                                                                   # Exit relaxation loop early.

    return int_pts                                                                  # Return the relaxed points.

def generate_region_cloud_poisson(region_points: list[tuple[float, float]], cloud_size: float) -> tuple[np.ndarray, np.ndarray, float] | tuple[None, None, None]:
    """
    Generate a complete cloud for a single region using Poisson Disk Sampling for more natural distribution.
    
    Input:
        region_points       list        List of region points.
        cloud_size          float       Cloud size.
        
    Output:
        boundary_points     ndarray     Boundary points.
        interior_points     ndarray     Interior points.
        cloud_size          float       Cloud size.
    """
    try:                                                                            # Start generation.
        # Create closed contour
        contour = create_closed_contour(region_points)                              # Obtain closed path.
        
        # Generate boundary points using the specified cloud size
        boundary_points = generate_boundary_points(contour, cloud_size)             # Generate perimeter nodes.
        
        # Create polygon for interior point generation
        polygon = Polygon(contour)                                                  # Instantiate Shapely polygon.
        if not polygon.is_valid:                                                    # Ensure topology is valid.
            polygon = polygon.buffer(0)                                             # Fix self-intersections.
        
        # Generate interior points using Poisson Disk Sampling
        interior_points = generate_interior_points_poisson(polygon, cloud_size)     # Populate interior with Poisson disk.
        
        # Apply Lloyd's relaxation to regularize the mesh (Honeycomb effect)
        if len(interior_points) > 0:                                                # Check if any interior points exist.
            interior_points = lloyd_relaxation(np.array(interior_points), np.array(boundary_points), polygon, iterations=5).tolist() # Apply smoothing.
        
        logging.info(f"Generated {len(boundary_points)} boundary points and {len(interior_points)} interior points using Poisson Disk Sampling") # Log counts.
        
        return boundary_points, interior_points, cloud_size                         # Return the generated data.
        
    except Exception as e:                                                          # On failure.
        logging.error(f"Error generating region cloud with Poisson Disk Sampling: {e}") # Log error details.
        # Fallback to uniform density method
        return generate_region_cloud_with_uniform_density(region_points, cloud_size)# Use simpler grid-based generator.

def generate_region_cloud_with_holes(main_region_points: list[tuple[float, float]], hole_regions_list: list[list[tuple[float, float]]], cloud_size: float = None, inside_regions: bool = False) -> tuple[np.ndarray, np.ndarray, float] | tuple[None, None, None]:
    """
    Generate a cloud for the main region (region 1) considering interior regions as holes.
    
    Input:
        main_region_points  list        Points for main region.
        hole_regions_list   list        List of hole regions.
        cloud_size          float       Cloud size.
        inside_regions      bool        Flag to process interior regions as holes.
        
    Output:
        boundary_points     ndarray     Boundary points.
        interior_points     ndarray     Interior points.
        cloud_size          float       Cloud size.
    """
    try:                                                                            # Start grid region with holes logic.
        # Calculate cloud size for the main region if not provided
        if cloud_size is None:                                                      # Check if explicit spacing is needed.
            cloud_size = calculate_cloud_size(main_region_points)                   # Compute base cloud size dynamically.
        
        # Create closed contour for main region
        main_contour = create_closed_contour(main_region_points)                    # Obtain exterior path.
        
        # Generate boundary points for main region
        boundary_points = generate_boundary_points(main_contour, cloud_size)        # Populate exterior boundary.
        
        # Create main polygon
        main_polygon = Polygon(main_contour)                                        # Create exterior geometric polygon.
        if not main_polygon.is_valid:                                               # Ensure it is valid.
            main_polygon = main_polygon.buffer(0)                                   # Clean geometry.
        
        # Create hole polygons and generate boundary points for holes
        hole_polygons = []                                                          # List to store inner void geometries.
        hole_boundary_count = 0                                                     # Counter for hole boundary nodes.
        
        for hole_region in hole_regions_list:                                       # Iterate over each hole definition.
            hole_contour = create_closed_contour(hole_region)                       # Create closed path for the hole.
            
            # Generate boundary points for this hole (it is an interior boundary)
            # Only if inside_regions is False (meaning we are NOT filling the holes with other regions)
            # If inside_regions is True, these boundaries belong to the interior regions (regions 2, 3, etc.)
            if not inside_regions:                                                  # If holes are empty spaces.
                hole_points = generate_boundary_points(hole_contour, cloud_size)    # Generate nodes along the hole border.
                if len(hole_points) > 0:                                            # Ensure valid points were generated.
                    boundary_points = np.vstack((boundary_points, hole_points))     # Add them to the main boundary array.
                    hole_boundary_count += len(hole_points)                         # Increment the tracking counter.
                
            hole_polygon = Polygon(hole_contour)                                    # Create geometric polygon for the hole.
            if not hole_polygon.is_valid:                                           # Ensure validity.
                hole_polygon = hole_polygon.buffer(0)                               # Clean geometry.
            hole_polygons.append(hole_polygon)                                      # Store the hole polygon.
        
        if hole_boundary_count > 0:                                                 # If hole boundaries were generated.
            logging.info(f"Generated {hole_boundary_count} boundary points from {len(hole_regions_list)} hole regions (inside_regions={inside_regions})") # Log count.
        
        # Create polygon with holes
        if hole_polygons:                                                           # If there are any holes.
            # Create a polygon with holes
            polygon_with_holes = main_polygon                                       # Start with the main body.
            for hole in hole_polygons:                                              # Iterate over the holes.
                # Subtract each hole from the main polygon
                polygon_with_holes = polygon_with_holes.difference(hole)            # Cut the hole out of the main shape.
        else:                                                                       # If no holes exist.
            polygon_with_holes = main_polygon                                       # Proceed with the solid main polygon.
        
        # Generate interior points avoiding holes
        interior_points = generate_interior_points(polygon_with_holes, cloud_size)  # Fill the remaining space using grid points.
        
        logging.info(f"Generated main region with holes: {len(boundary_points)} boundary points and {len(interior_points)} interior points") # Log generated point stats.
        logging.info(f"Excluded {len(hole_polygons)} hole regions from main region")# Log hole counts.
        
        return boundary_points, interior_points, cloud_size                         # Return generated data.
        
    except Exception as e:                                                          # If execution fails.
        logging.error(f"Error generating region cloud with holes: {e}")             # Log error details.
        return None, None, None                                                     # Return failure values.

def generate_region_cloud_with_holes_poisson(main_region_points: list[tuple[float, float]], hole_regions_list: list[list[tuple[float, float]]], cloud_size: float = None, inside_regions: bool = False) -> tuple[np.ndarray, np.ndarray, float] | tuple[None, None, None]:
    """
    Generate a cloud for the main region (region 1) considering interior regions as holes,
    using Poisson Disk Sampling for more natural distribution.
    
    Input:
        main_region_points  list        Points for main region.
        hole_regions_list   list        List of hole regions.
        cloud_size          float       Cloud size.
        inside_regions      bool        Flag to process interior regions as holes.
        
    Output:
        boundary_points     ndarray     Boundary points.
        interior_points     ndarray     Interior points.
        cloud_size          float       Cloud size.
    """
    try:                                                                            # Start poisson region with holes logic.
        # Calculate cloud size for the main region if not provided
        if cloud_size is None:                                                      # Check if explicit spacing is needed.
            cloud_size = calculate_cloud_size(main_region_points)                   # Compute base cloud size dynamically.
        
        # Create closed contour for main region
        main_contour = create_closed_contour(main_region_points)                    # Obtain exterior path.
        
        # Generate boundary points for main region
        boundary_points = generate_boundary_points(main_contour, cloud_size)        # Populate exterior boundary.
        
        # Create main polygon
        main_polygon = Polygon(main_contour)                                        # Create exterior geometric polygon.
        if not main_polygon.is_valid:                                               # Ensure it is valid.
            main_polygon = main_polygon.buffer(0)                                   # Clean geometry.
        
        # Create hole polygons and generate boundary points for holes
        hole_polygons = []                                                          # List to store inner void geometries.
        hole_boundary_count = 0                                                     # Counter for hole boundary nodes.
        
        for hole_region in hole_regions_list:                                       # Iterate over each hole definition.
            hole_contour = create_closed_contour(hole_region)                       # Create closed path for the hole.
            
            # Generate boundary points for this hole (it is an interior boundary)
            # Only if inside_regions is False (meaning we are NOT filling the holes with other regions)
            # If inside_regions is True, these boundaries belong to the interior regions (regions 2, 3, etc.)
            if not inside_regions:                                                  # If holes are empty spaces.
                hole_points = generate_boundary_points(hole_contour, cloud_size)    # Generate nodes along the hole border.
                if len(hole_points) > 0:                                            # Ensure valid points were generated.
                    boundary_points = np.vstack((boundary_points, hole_points))     # Add them to the main boundary array.
                    hole_boundary_count += len(hole_points)                         # Increment the tracking counter.
                
            hole_polygon = Polygon(hole_contour)                                    # Create geometric polygon for the hole.
            if not hole_polygon.is_valid:                                           # Ensure validity.
                hole_polygon = hole_polygon.buffer(0)                               # Clean geometry.
            hole_polygons.append(hole_polygon)                                      # Store the hole polygon.
        
        if hole_boundary_count > 0:                                                 # If hole boundaries were generated.
            logging.info(f"Generated {hole_boundary_count} boundary points from {len(hole_regions_list)} hole regions (inside_regions={inside_regions})") # Log count.
        
        # Create polygon with holes
        if hole_polygons:                                                           # If there are any holes.
            # Create a polygon with holes
            polygon_with_holes = main_polygon                                       # Start with the main body.
            for hole in hole_polygons:                                              # Iterate over the holes.
                # Subtract each hole from the main polygon
                polygon_with_holes = polygon_with_holes.difference(hole)            # Cut the hole out of the main shape.
        else:                                                                       # If no holes exist.
            polygon_with_holes = main_polygon                                       # Proceed with the solid main polygon.
        
        # Generate interior points avoiding holes using Poisson Disk Sampling
        interior_points = generate_interior_points_poisson(polygon_with_holes, cloud_size) # Fill space with Poisson sampling.
        
        # Apply Lloyd's relaxation
        if len(interior_points) > 0:                                                # Check if any interior points exist.
            interior_points = lloyd_relaxation(np.array(interior_points), np.array(boundary_points), polygon_with_holes, iterations=5).tolist() # Apply smoothing.
        
        logging.info(f"Generated main region with holes using Poisson: {len(boundary_points)} boundary points and {len(interior_points)} interior points") # Log generated point stats.
        logging.info(f"Excluded {len(hole_polygons)} hole regions from main region")# Log hole counts.
        
        return boundary_points, interior_points, cloud_size                         # Return generated data.
        
    except Exception as e:                                                          # If Poisson execution fails.
        logging.error(f"Error generating region cloud with holes using Poisson: {e}") # Log error details.
        # Fallback to grid-based method with holes
        return generate_region_cloud_with_holes(main_region_points, hole_regions_list) # Use simpler grid-based generator.

def generate_region_task(region_points: list[tuple[float, float]], region_id: int, main_cloud_size: float) -> tuple[np.ndarray, np.ndarray, int] | None:
    """
    Helper function for parallel region generation.
    Must be at top level for ProcessPoolExecutor (pickling).
    
    Input:
        region_points       list        Points for this region.
        region_id           int         ID of region.
        main_cloud_size     float       Base cloud size.
        
    Output:
        tuple containing (boundary_points, interior_points, region_id) or None.
    """
    try:                                                                            # Start isolated grid task execution.
        poly = Polygon(region_points)                                               # Create a geometry for checking dimensions.
        minx, miny, maxx, maxy = poly.bounds                                        # Get bounding box.
        min_dim = min(maxx - minx, maxy - miny)                                     # Compute the shortest dimension span.
        
        # Shrink cloud size for very small regions
        actual_cloud_size = main_cloud_size                                         # Start with the main spacing.
        if min_dim < main_cloud_size * 2:                                           # If region is extremely thin or small.
            actual_cloud_size = min_dim / 3.0                                       # Scale down spacing to ensure points fit.
            logging.info(f"Region {region_id} is small (min_dim={min_dim:.6f}). Adapted cloud size to {actual_cloud_size:.6f}") # Log down-scaling.

        boundary_points, interior_points, _ = generate_region_cloud_with_uniform_density(region_points, actual_cloud_size) # Delegate to grid generator.
        
        # Fallback to representative point if interior is empty
        if interior_points is not None and len(interior_points) == 0:               # If it failed to place interior points.
            logging.warning(f"Region {region_id} has 0 interior points. Using representative point fallback.") # Log fallback action.
            rep = poly.representative_point()                                       # Get a guaranteed internal coordinate.
            interior_points = np.array([[rep.x, rep.y]])                            # Use it as the single interior point.

        if boundary_points is not None and interior_points is not None:             # Validate generated data.
            logging.info(f"Generated cloud for interior region {region_id} with uniform density (cloud_size: {actual_cloud_size:.6f})") # Log completion.
            return (boundary_points, interior_points, region_id)                    # Return packaged results.
        else:                                                                       # If generation returned None.
            logging.warning(f"Failed to generate cloud for interior region {region_id}") # Log failure.
            return None                                                             # Return failure.
    except Exception as e:                                                          # If execution crashes.
        logging.error(f"Error generating cloud for interior region {region_id}: {e}") # Log error.
        return None                                                                 # Return failure.

def generate_interior_regions_clouds(regions: list[list[tuple[float, float]]], main_cloud_size: float) -> list[tuple[np.ndarray, np.ndarray, int]]:
    """
    Generate clouds for interior regions (regions 2 onwards) using the same cloud size as the main region.
    This ensures uniform density across all regions.
    
    Input:
        regions             list        List of region points.
        main_cloud_size     float       Base cloud size.
        
    Output:
        interior_clouds     list        List of tuples with results.
    """
    interior_clouds = []                                                            # Initialize the results container.
    
    # Process regions 2 onwards (skip region 1 which is index 0)
    # Using ProcessPoolExecutor for parallelization
    
    # Determine max workers (leave one core free or use all if few)
    max_workers = max(1, os.cpu_count() - 1)                                        # Determine optimal thread count.
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:                  # Start multiprocessing pool.
        futures = []                                                                # List to track running tasks.
        for i in range(1, len(regions)):                                            # Iterate over inner regions only.
            region_points = regions[i]                                              # Extract contour.
            region_id = i + 1                                                       # Compute region ID.
            futures.append(executor.submit(generate_region_task, region_points, region_id, main_cloud_size)) # Submit job.
        
        for future in futures:                                                      # Wait and collect all task results.
            try:                                                                    # Process individual task result.
                result = future.result()                                            # Fetch the job output.
                if result:                                                          # Ensure valid result.
                    interior_clouds.append(result)                                  # Add it to the main collection.
            except Exception as e:                                                  # Handle task crashes.
                logging.error(f"Error retrieving parallel task result: {e}")        # Log error.
    
    return interior_clouds                                                          # Return compiled results.

def generate_region_task_poisson(region_points: list[tuple[float, float]], region_id: int, main_cloud_size: float) -> tuple[np.ndarray, np.ndarray, int] | None:
    """
    Helper for parallel Poisson interior regions.
    
    Input:
        region_points       list        Points for this region.
        region_id           int         ID of region.
        main_cloud_size     float       Base cloud size.
        
    Output:
        tuple containing (boundary_points, interior_points, region_id) or None.
    """
    try:                                                                            # Start isolated Poisson task execution.
        poly = Polygon(region_points)                                               # Create a geometry for checking dimensions.
        minx, miny, maxx, maxy = poly.bounds                                        # Get bounding box.
        min_dim = min(maxx - minx, maxy - miny)                                     # Compute the shortest dimension span.
        
        actual_cloud_size = main_cloud_size                                         # Start with the main spacing.
        if min_dim < main_cloud_size * 2:                                           # If region is extremely thin or small.
            actual_cloud_size = min_dim / 3.0                                       # Scale down spacing to ensure points fit.
            logging.info(f"Region {region_id} is small (min_dim={min_dim:.6f}). Adapted Poisson cloud size to {actual_cloud_size:.6f}") # Log down-scaling.
            
        boundary_points, interior_points, _ = generate_region_cloud_poisson(region_points, actual_cloud_size) # Delegate to Poisson generator.
        
        if interior_points is not None and len(interior_points) == 0:               # If it failed to place interior points.
            logging.warning(f"Region {region_id} has 0 Poisson interior points. Using representative point fallback.") # Log fallback action.
            rep = poly.representative_point()                                       # Get a guaranteed internal coordinate.
            interior_points = np.array([[rep.x, rep.y]])                            # Use it as the single interior point.
            
        if boundary_points is not None and interior_points is not None:             # Validate generated data.
            return (boundary_points, interior_points, region_id)                    # Return packaged results.
        return None                                                                 # Return failure if data is missing.
    except Exception as e:                                                          # If execution crashes.
        logging.error(f"Error in poisson interior task {region_id}: {e}")           # Log error.
        return None                                                                 # Return failure.

def generate_interior_regions_clouds_poisson(regions: list[list[tuple[float, float]]], main_cloud_size: float) -> list[tuple[np.ndarray, np.ndarray, int]]:
    """
    Parallel generation for interior regions using Poisson.
    
    Input:
        regions             list        List of region points.
        main_cloud_size     float       Base cloud size.
        
    Output:
        interior_clouds     list        List of tuples with results.
    """
    interior_clouds = []                                                            # Initialize the results container.
    max_workers = max(1, os.cpu_count() - 1)                                        # Determine optimal thread count.
    with ProcessPoolExecutor(max_workers=max_workers) as executor:                  # Start multiprocessing pool.
        futures = []                                                                # List to track running tasks.
        for i in range(1, len(regions)):                                            # Iterate over inner regions only.
            futures.append(executor.submit(generate_region_task_poisson, regions[i], i + 1, main_cloud_size)) # Submit job.
        
        for future in futures:                                                      # Wait and collect all task results.
            try:                                                                    # Process individual task result.
                result = future.result()                                            # Fetch the job output.
                if result:                                                          # Ensure valid result.
                    interior_clouds.append(result)                                  # Add it to the main collection.
            except Exception as e:                                                  # Handle task crashes.
                logging.error(f"Error in parallel poisson task: {e}")               # Log error.
    return interior_clouds                                                          # Return compiled results.
