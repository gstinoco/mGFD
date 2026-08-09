"""
CloudGenerator.core.point_generation.poisson_generation — Poisson-Based Point Generation

Overview:
    This module implements the regular (Poisson Disk Sampling) point generation algorithms.
    It relies on SciPy's QMC engine for fast generation and Lloyd's relaxation for uniformity.

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
import os                                                                                                                               # Operating system interfaces.
import cv2                                                                                                                              # OpenCV for fast rasterization and masking.
import numpy as np                                                                                                                      # Core numerical operations.

from scipy.spatial import cKDTree                                                                                                       # Import KDTree for distance checks.
from shapely.geometry import Polygon                                                                                                    # Geometric objects and operations.
from scipy.stats.qmc import PoissonDisk                                                                                                 # Import SciPy Poisson Disk generator.
from typing import List, Tuple, Optional, Any                                                                                           # Type hinting.
from concurrent.futures import ProcessPoolExecutor                                                                                      # Parallel execution.

from mGFD.cloud_generator.core.point_generation.relaxation import lloyd_relaxation                                                      # Node relaxation.
from mGFD.cloud_generator.utils.utils import calculate_cloud_size, create_closed_contour                                                # Utility functions.
from mGFD.cloud_generator.core.point_generation.boundary import generate_boundary_points                                                # Boundary generation.
from mGFD.cloud_generator.core.point_generation.geometry import create_fast_polygon_checker                                             # Geometric operations.
from mGFD.cloud_generator.core.point_generation.regular_generation import (                                                             # Fallback regular generation methods.
    generate_interior_points,
    generate_region_cloud_with_uniform_density,
    generate_region_cloud_with_holes
)

def poisson_disk_sampling(polygon: Polygon, radius: float, k: int = 30, boundary_points: Optional[List[Tuple[float, float]]] = None) -> List[Tuple[float, float]]:
    """
    poisson_disk_sampling
    Generate points using Poisson Disk Sampling algorithm via SciPy QMC engine for extreme speed.
    
    Input:
        polygon             Polygon     Target polygon to fill with points.
        radius              float       Minimum distance between any two generated points.
        k                   int         Ignored in SciPy version, kept for API compatibility. Default is 30.
        boundary_points     List        Pre-existing boundary points to check distances against.
    
    Output:
        final_points        List        List of (x, y) coordinate tuples representing generated points.
    """
    minx, miny, maxx, maxy = polygon.bounds                                                                                             # Get bounding box limits.
    
    # 1. Scale dimensions and radius to [0, 1] space
    # We must scale our dimensions and the radius down to [0, 1].
    width  = maxx - minx                                                                                                                # Calculate domain width.
    height = maxy - miny                                                                                                                # Calculate domain height.
    
    # If the polygon is just a point or a line, we can't fill it
    if width <= 0 or height <= 0:                                                                                                       # Ensure domain is 2D.
        return []                                                                                                                       # Return empty list if not.
        
    scale_factor  = max(width, height)                                                                                                  # Determine scale factor for normalization.
    scaled_radius = radius / scale_factor                                                                                               # Normalize the sampling radius.
    
    try:                                                                                                                                # Try sampling generation.
        # 2. Initialize SciPy engine
        engine  = PoissonDisk(d=2, radius=scaled_radius)                                                                                # Create Poisson generator engine.
        # Generate samples in [0, 1]^2
        samples = engine.fill_space()                                                                                                   # Generate the normalized points.
    except Exception as e:                                                                                                              # Fallback if sampling fails.
        # Very rough fallback if radius is too large for the bounding box
        return []                                                                                                                       # Return empty list.

    # 3. Scale back to world coordinates
    samples[:, 0] = samples[:, 0] * scale_factor + minx                                                                                 # Scale and shift x coordinates.
    samples[:, 1] = samples[:, 1] * scale_factor + miny                                                                                 # Scale and shift y coordinates.
    
    # 4. Filter points inside polygon using vectorized OpenCV mask for extreme speed
    try:                                                                                                                                # Attempt fast masking.
        # Create a mask with a resolution high enough to capture the polygon accurately
        resolution = max(radius * 0.1, (maxx - minx) / 2000.0)                                                                          # Calculate optimal pixel resolution.
        if resolution <= 0: resolution = 1.0                                                                                            # Prevent zero or negative resolution.
        
        width_px  = int(np.ceil((maxx - minx) / resolution)) + 1                                                                        # Calculate mask width in pixels.
        height_px = int(np.ceil((maxy - miny) / resolution)) + 1                                                                        # Calculate mask height in pixels.
        
        mask = np.zeros((height_px, width_px), dtype=np.uint8)                                                                          # Initialize the OpenCV mask array.
        
        def to_pixel_coords(coords: List[Tuple[float, float]]) -> np.ndarray:                                                           # Coordinate transformation helper.
            pixel_coords = []                                                                                                           # List to store pixel coordinates.
            for x, y in coords:                                                                                                         # Iterate over world coordinates.
                px = int(round((x - minx) / resolution))                                                                                # Convert x to pixel index.
                py = int(round((y - miny) / resolution))                                                                                # Convert y to pixel index.
                pixel_coords.append([px, py])                                                                                           # Add to list.
            return np.array(pixel_coords, dtype=np.int32)                                                                               # Return integer array.
            
        # Draw exterior
        ext_pixels = to_pixel_coords(list(polygon.exterior.coords))                                                                     # Convert exterior to pixels.
        cv2.fillPoly(mask, [ext_pixels], 1)                                                                                             # type: ignore
        # Exclude exterior boundary to avoid overlapping
        cv2.polylines(mask, [ext_pixels], isClosed=True, color=0, thickness=1)                                                          # type: ignore
        
        # Draw holes
        for interior in polygon.interiors:                                                                                              # Iterate over internal holes.
            int_pixels = to_pixel_coords(list(interior.coords))                                                                         # Convert hole to pixels.
            cv2.fillPoly(mask, [int_pixels], 0)                                                                                         # type: ignore
            
        # Map samples to pixel coordinates and check mask
        sample_px = np.round((samples[:, 0] - minx) / resolution).astype(int)                                                           # Get pixel x for samples.
        sample_py = np.round((samples[:, 1] - miny) / resolution).astype(int)                                                           # Get pixel y for samples.
        
        # Filter bounds to avoid IndexError
        valid_idx = (sample_px >= 0) & (sample_px < width_px) & (sample_py >= 0) & (sample_py < height_px)                              # Mask bounds check.
        
        # Of the valid indices, check the mask
        in_polygon            = np.zeros(len(samples), dtype=bool)                                                                      # Initialize boolean inclusion mask.
        in_polygon[valid_idx] = mask[sample_py[valid_idx], sample_px[valid_idx]] > 0                                                    # Check pixel values.
        
        interior_points = samples[in_polygon].tolist()                                                                                  # Select samples inside the mask.
        interior_points = [(p[0], p[1]) for p in interior_points]                                                                       # Convert to list of tuples.
        
    except Exception as e:                                                                                                              # If fast masking fails.
        # 5. Fallback if fast masking fails
        fast_contains   = create_fast_polygon_checker(polygon)                                                                          # Create path checker.
        interior_points = [                                                                                                             # List comprehension for fallback filtering.
            (x, y) for x, y in samples                                                                                                  # Iterate over samples.
            if fast_contains(x, y)                                                                                                      # Check inclusion.
        ]
    
    if not interior_points:                                                                                                             # Check if any points survived filtering.
        return []                                                                                                                       # Return empty list if not.
        
    if boundary_points and len(boundary_points) > 0:                                                                                    # If boundary points exist for distance checking.
        # 6. Filter points that are too close to the boundary points using cKDTree
        boundary_array = np.array(boundary_points)                                                                                      # Convert boundary to NumPy array.
        interior_array = np.array(interior_points)                                                                                      # Convert interior to NumPy array.
        
        # Build KDTree for boundary points
        tree = cKDTree(boundary_array)                                                                                                  # Initialize the KDTree.
        
        # Query distances. k=1 gets the closest boundary point to each interior point
        distances, _ = tree.query(interior_array, k=1)                                                                                  # Get distances to nearest boundary point.
        
        # Keep only interior points that are at least `radius` away from any boundary
        valid_indices = distances >= (radius * 0.8)                                                                                     # Filter condition.
        
        final_points = interior_array[valid_indices].tolist()                                                                           # Filter array using mask.
        return [(p[0], p[1]) for p in final_points]                                                                                     # Convert and return valid points.
        
    return interior_points                                                                                                              # Return points directly if no boundary checks needed.

def generate_interior_points_poisson(polygon: Polygon, cloud_size: float) -> List[Tuple[float, float]]:
    """
    generate_interior_points_poisson
    Generate interior points using Poisson Disk Sampling for more natural distribution.
    
    Input:
        polygon             Polygon     Shapely polygon representing the region.
        cloud_size          float       Target spacing between points (used as radius).
    
    Output:
        interior_points     List        List of (x, y) tuples representing interior points.
    """
    # Use cloud_size as the minimum distance between points
    radius = cloud_size * 0.8                                                                                                           # Set the sampling radius.
    
    try:                                                                                                                                # Start generation.
        interior_points = poisson_disk_sampling(polygon, radius)                                                                        # Perform Poisson sampling.
        
        return interior_points                                                                                                          # Return the interior points.
    except Exception as e:                                                                                                              # On failure.
        # Fallback to original method
        
        return list(map(tuple, generate_interior_points(polygon, cloud_size)))                                                          # Return grid-based points.

def generate_region_cloud_poisson(region_points: List[Tuple[float, float]], cloud_size: float) -> Tuple[Any, Any, Any]:
    """
    generate_region_cloud_poisson
    Generate a complete cloud for a single region using Poisson Disk Sampling for more natural distribution.
    
    Input:
        region_points       List                List of region points.
        cloud_size          float               Cloud size.
        
    Output:
        boundary_points     m x 2 ndarray       Boundary points.
        interior_points     m x 2 ndarray       Interior points.
        cloud_size          float               Cloud size.
    """
    try:                                                                                                                                # Start generation.
        # 1. Create closed contour
        contour = create_closed_contour(region_points)                                                                                  # Obtain closed path.
        
        # 2. Generate boundary points
        boundary_points = generate_boundary_points(contour, cloud_size)                                                                 # Generate perimeter nodes.
        
        # 3. Create polygon for interior point generation
        polygon = Polygon(contour)                                                                                                      # Instantiate Shapely polygon.
        if not polygon.is_valid:                                                                                                        # Ensure topology is valid.
            polygon = polygon.buffer(0)                                                                                                 # Fix self-intersections.
        
        # 4. Generate interior points
        interior_points = generate_interior_points_poisson(polygon, cloud_size)                                                         # Populate interior with Poisson disk.
        
        # 5. Apply Lloyd's relaxation
        if len(interior_points) > 0:                                                                                                    # Check if any interior points exist.
            interior_points = lloyd_relaxation(np.array(interior_points), np.array(boundary_points), polygon, iterations=5).tolist()    # Apply smoothing.

        return boundary_points, interior_points, cloud_size                                                                             # Return the generated data.
        
    except Exception as e:                                                                                                              # On failure.
        
        # Fallback to uniform density method
        return generate_region_cloud_with_uniform_density(region_points, cloud_size)                                                    # Use simpler grid-based generator.

def generate_region_cloud_with_holes_poisson(main_region_points: List[Tuple[float, float]], hole_regions_list: List[List[Tuple[float, float]]], cloud_size: Optional[float] = None, inside_regions: bool = False) -> Tuple[Any, Any, Any]:
    """
    generate_region_cloud_with_holes_poisson
    Generate a cloud for the main region (region 1) considering interior regions as holes,
    using Poisson Disk Sampling for more natural distribution.
    
    Input:
        main_region_points  List                Points for main region.
        hole_regions_list   List                List of hole regions.
        cloud_size          float               Cloud size.
        inside_regions      bool                Flag to process interior regions as holes.
        
    Output:
        boundary_points     m x 2 ndarray       Boundary points.
        interior_points     m x 2 ndarray       Interior points.
        cloud_size          float               Cloud size.
    """
    try:                                                                                                                                # Start poisson region with holes logic.
        # 1. Calculate cloud size for the main region if not provided
        if cloud_size is None:                                                                                                          # Check if explicit spacing is needed.
            cloud_size = calculate_cloud_size(main_region_points)                                                                       # Compute base cloud size dynamically.
        
        # 2. Create closed contour for main region
        main_contour = create_closed_contour(main_region_points)                                                                        # Obtain exterior path.
        
        # 3. Generate boundary points for main region
        boundary_points = generate_boundary_points(main_contour, cloud_size)                                                            # Populate exterior boundary.
        
        # 4. Create main polygon
        main_polygon = Polygon(main_contour)                                                                                            # Create exterior geometric polygon.
        if not main_polygon.is_valid:                                                                                                   # Ensure it is valid.
            main_polygon = main_polygon.buffer(0)                                                                                       # Clean geometry.
        
        # 5. Create hole polygons and generate boundary points for holes
        hole_polygons       = []                                                                                                        # List to store inner void geometries.
        hole_boundary_count = 0                                                                                                         # Counter for hole boundary nodes.
        
        for hole_region in hole_regions_list:                                                                                           # Iterate over each hole definition.
            hole_contour = create_closed_contour(hole_region)                                                                           # Create closed path for the hole.
            
            # Generate boundary points for this hole (it is an interior boundary)
            # Only if inside_regions is False (meaning we are NOT filling the holes with other regions)
            # If inside_regions is True, these boundaries belong to the interior regions (regions 2, 3, etc.)
            if not inside_regions:                                                                                                      # If holes are empty spaces.
                hole_points = generate_boundary_points(hole_contour, cloud_size)                                                        # Generate nodes along the hole border.
                if len(hole_points) > 0:                                                                                                # Ensure valid points were generated.
                    boundary_points      = np.vstack((boundary_points, hole_points))                                                    # Add them to the main boundary array.
                    hole_boundary_count += len(hole_points)                                                                             # Increment the tracking counter.
                
            hole_polygon = Polygon(hole_contour)                                                                                        # Create geometric polygon for the hole.
            if not hole_polygon.is_valid:                                                                                               # Ensure validity.
                hole_polygon = hole_polygon.buffer(0)                                                                                   # Clean geometry.
            hole_polygons.append(hole_polygon)                                                                                          # Store the hole polygon.
        
        if hole_boundary_count > 0:                                                                                                     # If hole boundaries were generated.
            pass                                                                                                                        # No action needed.
                                                                                                                                        # Log count.
        
        # 6. Create polygon with holes
        if hole_polygons:                                                                                                               # If there are any holes.
            # Create a polygon with holes
            polygon_with_holes = main_polygon                                                                                           # Start with the main body.
            for hole in hole_polygons:                                                                                                  # Iterate over the holes.
                # Subtract each hole from the main polygon
                polygon_with_holes = polygon_with_holes.difference(hole)                                                                # Cut the hole out of the main shape.
        else:                                                                                                                           # If no holes exist.
            polygon_with_holes = main_polygon                                                                                           # Proceed with the solid main polygon.
        
        # 7. Generate interior points avoiding holes using Poisson Disk Sampling
        interior_points = generate_interior_points_poisson(polygon_with_holes, cloud_size)                                              # Fill space with Poisson sampling.
        
        # Apply Lloyd's relaxation
        if len(interior_points) > 0:                                                                                                    # Check if any interior points exist.
            interior_points = lloyd_relaxation(np.array(interior_points), np.array(boundary_points), polygon_with_holes, iterations=5).tolist()
                                                                                                                                        # Apply smoothing.

        return boundary_points, interior_points, cloud_size                                                                             # Return generated data.
        
    except Exception as e:                                                                                                              # If Poisson execution fails.
        # Fallback to grid-based method with holes
        return generate_region_cloud_with_holes(main_region_points, hole_regions_list)                                                  # Use simpler grid-based generator.

def generate_region_task_poisson(region_points: List[Tuple[float, float]], region_id: int, main_cloud_size: float) -> Optional[Tuple[np.ndarray, np.ndarray, int]]:
    """
    generate_region_task_poisson
    Helper for parallel Poisson interior regions.
    
    Input:
        region_points       List        Points for this region.
        region_id           int         ID of region.
        main_cloud_size     float       Base cloud size.
        
    Output:
        tuple containing (boundary_points, interior_points, region_id) or None.
    """
    try:                                                                                                                                # Start isolated Poisson task execution.
        poly                   = Polygon(region_points)                                                                                 # Create a geometry for checking dimensions.
        minx, miny, maxx, maxy = poly.bounds                                                                                            # Get bounding box.
        min_dim                = min(maxx - minx, maxy - miny)                                                                          # Compute the shortest dimension span.
        
        actual_cloud_size = main_cloud_size                                                                                             # Start with the main spacing.
        if min_dim < main_cloud_size * 2:                                                                                               # If region is extremely thin or small.
            actual_cloud_size = min_dim / 3.0                                                                                           # Scale down spacing to ensure points fit.
            
        boundary_points, interior_points, _ = generate_region_cloud_poisson(region_points, actual_cloud_size)                           # Delegate to Poisson generator.
        
        if interior_points is not None and len(interior_points) == 0:                                                                   # If it failed to place interior points.
            rep             = poly.representative_point()                                                                               # Get a guaranteed internal coordinate.
            interior_points = np.array([[rep.x, rep.y]])                                                                                # Use it as the single interior point.
            
        if boundary_points is not None and interior_points is not None:                                                                 # Validate generated data.
            return (boundary_points, interior_points, region_id)                                                                        # Return packaged results.
        return None                                                                                                                     # Return failure if data is missing.
    except Exception as e:                                                                                                              # If execution crashes.
        return None                                                                                                                     # Return failure.

def generate_interior_regions_clouds_poisson(regions: List[List[Tuple[float, float]]], main_cloud_size: float) -> List[Tuple[np.ndarray, np.ndarray, int]]:
    """
    generate_interior_regions_clouds_poisson
    Parallel generation for interior regions using Poisson.
    
    Input:
        regions             List        List of region points.
        main_cloud_size     float       Base cloud size.
        
    Output:
        interior_clouds     List        List of tuples with results.
    """
    interior_clouds = []                                                                                                                # Initialize the results container.
    max_workers     = max(1, (os.cpu_count() or 1) - 1)                                                                                 # Determine optimal thread count.
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:                                                                      # Start multiprocessing pool.
        futures = []                                                                                                                    # List to track running tasks.
        for i in range(1, len(regions)):                                                                                                # Iterate over inner regions only.
            futures.append(executor.submit(generate_region_task_poisson, regions[i], i + 1, main_cloud_size))                           # Submit job.
        
        for future in futures:                                                                                                          # Wait and collect all task results.
            try:                                                                                                                        # Process individual task result.
                result = future.result()                                                                                                # Fetch the job output.
                if result:                                                                                                              # Ensure valid result.
                    interior_clouds.append(result)                                                                                      # Add it to the main collection.
            except Exception as e:                                                                                                      # Handle task crashes.
                pass                                                                                                                    # Ignored error.
    
    return interior_clouds                                                                                                              # Return compiled results.
