"""
CloudGenerator.core.point_generation.regular_generation — Regular Point Generation

Overview:
    This module implements the regular (grid-based) point generation algorithms.
    It uses a vectorized OpenCV masking approach for fast interior point filling.

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
import numba as nb                                                                                                                      # JIT compilation.                                                                                                                      # Core numerical operations.

from shapely.geometry import Polygon                                                                                                    # Geometric objects and operations.
from typing import List, Tuple, Optional, Any                                                                                           # Type hinting.
from concurrent.futures import ProcessPoolExecutor                                                                                      # Parallel execution.

from mGFD.cloud_generator.utils.utils import calculate_cloud_size, create_closed_contour                                                # Utility functions.
from mGFD.cloud_generator.core.point_generation.boundary import generate_boundary_points                                                # Boundary generation.
from mGFD.cloud_generator.core.point_generation.geometry import create_fast_polygon_checker                                             # Fast polygon checker.
from mGFD.cloud_generator.core.point_generation.jit_geometry import _is_point_in_polygon_jit                                            # Fast JIT geometric operations.                                             # Geometric operations.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Decorator for JIT compilation.
def _generate_interior_fallback_jit(x_range: np.ndarray, y_range: np.ndarray, poly_coords: np.ndarray) -> np.ndarray:
    """
    _generate_interior_fallback_jit
    Numba JIT-compiled grid point generator using Ray-Casting.
    
    Input:
        x_range         nx ndarray      Range of X coordinates.
        y_range         ny ndarray      Range of Y coordinates.
        poly_coords     p x 2 ndarray   Coordinates defining the polygon boundary.
        
    Output:
        points          m x 2 ndarray   Array of (x, y) generated points.
    """
    nx       = len(x_range)                                                                                                             # Number of X coordinates.
    ny       = len(y_range)                                                                                                             # Number of Y coordinates.
    capacity = nx * ny                                                                                                                  # Initial capacity.
    out      = np.zeros((capacity, 2), dtype=np.float64)                                                                                # Allocate output array.
    valid    = np.zeros(capacity, dtype=np.bool_)                                                                                       # Track valid grid points.

    for i in nb.prange(nx):                                                                                                             # type: ignore
        x = x_range[i]                                                                                                                  # Extract x coordinate.
        
        for j in range(ny):                                                                                                             # Iterate over y grid.
            y   = y_range[j]                                                                                                            # Extract y coordinate.
            idx = i * ny + j                                                                                                            # Global index.
            
            if _is_point_in_polygon_jit(x, y, poly_coords):                                                                             # Check if point is inside.
                out[idx, 0] = x                                                                                                         # Store x coordinate.
                out[idx, 1] = y                                                                                                         # Store y coordinate.
                valid[idx]  = True                                                                                                      # Mark index as valid.
                
    count = 0                                                                                                                           # Counter for output points.
    for idx in range(capacity):                                                                                                         # Sequential loop for valid points count.
        if valid[idx]:                                                                                                                  # If point is inside.
            count += 1                                                                                                                  # Increment valid count.
            
    final_out = np.zeros((count, 2), dtype=np.float64)                                                                                  # Allocate condensed array.
    
    idx_out = 0                                                                                                                         # Iterator for final output array.
    for idx in range(capacity):                                                                                                         # Sequential loop for transfer.
        if valid[idx]:                                                                                                                  # If point is inside.
            final_out[idx_out, 0] = out[idx, 0]                                                                                         # Copy x.
            final_out[idx_out, 1] = out[idx, 1]                                                                                         # Copy y.
            idx_out              += 1                                                                                                   # Increment index.

    return final_out                                                                                                                    # Return valid slice of output array.

def generate_interior_points(polygon: Any, cloud_size: float) -> np.ndarray:
    """
    generate_interior_points
    Generate points inside a polygon using a vectorized grid-based approach.
    
    Optimized using OpenCV mask for fast point-in-polygon checks (replacing Matplotlib Path).
    
    Input:
        polygon             Polygon             The polygon to fill with points.
        cloud_size          float               Desired spacing between points.
    
    Output:
        points              m x 2 ndarray       Array of (x, y) coordinates of generated points.
    """
    try:                                                                                                                                # Start of the main generation block.
        bounds                     = polygon.bounds                                                                                     # Get bounding box of the polygon.
        x_min, y_min, x_max, y_max = bounds                                                                                             # Extract bounding box limits.
        
        # 1. Create grid of points
        x_range = np.arange(x_min, x_max + cloud_size, cloud_size)                                                                      # Create range for x coordinates.
        y_range = np.arange(y_min, y_max + cloud_size, cloud_size)                                                                      # Create range for y coordinates.
        
        if len(x_range) == 0 or len(y_range) == 0:                                                                                      # Check if grids are valid.
            return np.empty((0, 2))                                                                                                     # Return empty array if not.
            
        xx, yy      = np.meshgrid(x_range, y_range)                                                                                     # Create 2D mesh grid.
        grid_points = np.vstack((xx.ravel(), yy.ravel())).T                                                                             # Flatten and stack into point list.
        
        # 2. Create a mask using OpenCV
        width  = len(x_range)                                                                                                           # Mask width in pixels.
        height = len(y_range)                                                                                                           # Mask height in pixels.
        mask   = np.zeros((height, width), dtype=np.uint8)                                                                              # Initialize the mask array.
        
        # 3. Function to convert world coords to pixel coords
        def to_pixel_coords(coords: List[Tuple[float, float]]) -> np.ndarray:                                                           # Helper for coordinate transformation.
            pixel_coords = []                                                                                                           # List for pixel coordinates.
            
            for x, y in coords:                                                                                                         # Iterate over world coordinates.
                px = int(round((x - x_min) / cloud_size))                                                                               # Map x to pixel index.
                py = int(round((y - y_min) / cloud_size))                                                                               # Map y to pixel index.
                pixel_coords.append([px, py])                                                                                           # Add to pixel list.
            
            return np.array(pixel_coords, dtype=np.int32)                                                                               # Return as integer NumPy array.

        # 4. Draw exterior
        ext_coords = list(polygon.exterior.coords)                                                                                      # Extract exterior coordinates.
        ext_pixels = to_pixel_coords(ext_coords)                                                                                        # Convert to pixel coordinates.
        cv2.fillPoly(mask, [ext_pixels], 1)                                                                                             # type: ignore
        
        # Exclude exterior boundary (set to 0) to avoid overlapping with explicit boundary nodes
        cv2.polylines(mask, [ext_pixels], isClosed=True, color=0, thickness=1)                                                          # type: ignore
        
        # 5. Draw interiors (holes)
        for interior in polygon.interiors:                                                                                              # Iterate over hole geometries.
            int_coords = list(interior.coords)                                                                                          # Extract hole coordinates.
            int_pixels = to_pixel_coords(int_coords)                                                                                    # Convert to pixel coordinates.
            cv2.fillPoly(mask, [int_pixels], 0)                                                                                         # type: ignore
            
        # 6. Select points
        selected_mask = mask.astype(bool).ravel()                                                                                       # Flatten mask and convert to boolean.
        points        = grid_points[selected_mask]                                                                                      # Filter grid points using the mask.
        
        return points if len(points) > 0 else np.empty((0, 2))                                                                          # Return points or empty array.
        
    except Exception as e:                                                                                                              # Fallback on OpenCV error.
        # 7. Fallback to JIT Ray-Casting method                                                                                         #
        x_range     = np.arange(x_min, x_max + cloud_size, cloud_size, dtype=np.float64)                                                # Re-create x range.
        y_range     = np.arange(y_min, y_max + cloud_size, cloud_size, dtype=np.float64)                                                # Re-create y range.
        poly_coords = np.array(polygon.exterior.coords, dtype=np.float64)                                                               # Extract polygon exterior coordinates.
        
        points = _generate_interior_fallback_jit(x_range, y_range, poly_coords)                                                         # Use JIT compiled generation.
        
        return points if len(points) > 0 else np.empty((0, 2))                                                                          # Return array of points.

def generate_region_cloud_with_uniform_density(region_points: List[Tuple[float, float]], cloud_size: float) -> Tuple[Any, Any, Any]:
    """
    generate_region_cloud_with_uniform_density
    Generate a complete cloud for a single region using a specified cloud size.
    This ensures uniform density across different regions.
    
    Input:
        region_points       List                List of region points.
        cloud_size          float               Cloud size to enforce.
        
    Output:
        boundary_points     m x 2 ndarray       Boundary points.
        interior_points     m x 2 ndarray       Interior points.
        cloud_size          float               Cloud size.
    """
    try:                                                                                                                                # Start fixed density generation.
        # 1. Create closed contour
        contour = create_closed_contour(region_points)                                                                                  # Obtain closed path.
        
        # 2. Generate boundary points using the specified cloud size
        boundary_points = generate_boundary_points(contour, cloud_size)                                                                 # Generate points on the perimeter.
        
        # 3. Create polygon for interior point generation
        polygon = Polygon(contour)                                                                                                      # Instantiate Shapely polygon.
        
        if not polygon.is_valid:                                                                                                        # Check if polygon topology is valid.
            polygon = polygon.buffer(0)                                                                                                 # Attempt to fix self-intersections.
        
        # 4. Generate interior points using the specified cloud size
        interior_points = generate_interior_points(polygon, cloud_size)                                                                 # Generate points inside the region.

        return boundary_points, interior_points, cloud_size                                                                             # Return generated data.
        
    except Exception as e:                                                                                                              # On failure.
        
        return None, None, None                                                                                                         # Return None values.

def generate_region_cloud_with_holes(main_region_points: List[Tuple[float, float]], hole_regions_list: List[List[Tuple[float, float]]], cloud_size: Optional[float] = None, inside_regions: bool = False) -> Tuple[Any, Any, Any]:
    """
    generate_region_cloud_with_holes
    Generate a cloud for the main region (region 1) considering interior regions as holes.
    
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
    try:                                                                                                                                # Start grid region with holes logic.
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
        
        # 6. Create polygon with holes
        if hole_polygons:                                                                                                               # If there are any holes.
            # Create a polygon with holes
            polygon_with_holes = main_polygon                                                                                           # Start with the main body.
            
            for hole in hole_polygons:                                                                                                  # Iterate over the holes.
                # Subtract each hole from the main polygon
                polygon_with_holes = polygon_with_holes.difference(hole)                                                                # Cut the hole out of the main shape.
        else:                                                                                                                           # If no holes exist.
            polygon_with_holes = main_polygon                                                                                           # Proceed with the solid main polygon.
        
        # 7. Generate interior points avoiding holes
        interior_points = generate_interior_points(polygon_with_holes, cloud_size)                                                      # Fill the remaining space using grid points.

        return boundary_points, interior_points, cloud_size                                                                             # Return generated data.
        
    except Exception as e:                                                                                                              # If execution fails.
        
        return None, None, None                                                                                                         # Return failure values.

def generate_region_task(region_points: List[Tuple[float, float]], region_id: int, main_cloud_size: float) -> Optional[Tuple[np.ndarray, np.ndarray, int]]:
    """
    generate_region_task
    Helper function for parallel region generation.
    Must be at top level for ProcessPoolExecutor (pickling).
    
    Input:
        region_points       List        Points for this region.
        region_id           int         ID of region.
        main_cloud_size     float       Base cloud size.
        
    Output:
        tuple containing (boundary_points, interior_points, region_id) or None.
    """
    try:                                                                                                                                # Start isolated grid task execution.
        poly                   = Polygon(region_points)                                                                                 # Create a geometry for checking dimensions.
        minx, miny, maxx, maxy = poly.bounds                                                                                            # Get bounding box.
        min_dim                = min(maxx - minx, maxy - miny)                                                                          # Compute the shortest dimension span.
        
        # 1. Shrink cloud size for very small regions
        actual_cloud_size = main_cloud_size                                                                                             # Start with the main spacing.
        
        if min_dim < main_cloud_size * 2:                                                                                               # If region is extremely thin or small.
            actual_cloud_size = min_dim / 3.0                                                                                           # Scale down spacing to ensure points fit.

        boundary_points, interior_points, _ = generate_region_cloud_with_uniform_density(region_points, actual_cloud_size)              # Delegate to grid generator.
        
        # 2. Fallback to representative point if interior is empty
        if interior_points is not None and len(interior_points) == 0:                                                                   # If it failed to place interior points.
            rep             = poly.representative_point()                                                                               # Get a guaranteed internal coordinate.
            interior_points = np.array([[rep.x, rep.y]])                                                                                # Use it as the single interior point.

        if boundary_points is not None and interior_points is not None:                                                                 # Validate generated data.
            return (boundary_points, interior_points, region_id)                                                                        # Return packaged results.
        else:                                                                                                                           # If generation returned None.
            return None                                                                                                                 # Return failure.
    
    except Exception as e:                                                                                                              # If execution crashes.
        return None                                                                                                                     # Return failure.

def generate_interior_regions_clouds(regions: List[List[Tuple[float, float]]], main_cloud_size: float) -> List[Tuple[np.ndarray, np.ndarray, int]]:
    """
    generate_interior_regions_clouds
    Generate clouds for interior regions (regions 2 onwards) using the same cloud size as the main region.
    This ensures uniform density across all regions.
    
    Input:
        regions             List        List of region points.
        main_cloud_size     float       Base cloud size.
        
    Output:
        interior_clouds     List        List of tuples with results.
    """
    interior_clouds = []                                                                                                                # Initialize the results container.
    
    # 1. Process regions 2 onwards (skip region 1 which is index 0)
    # Using ProcessPoolExecutor for parallelization
    
    # 2. Determine max workers (leave one core free or use all if few)
    max_workers = max(1, (os.cpu_count() or 1) - 1)                                                                                     # Determine optimal thread count.
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:                                                                      # Start multiprocessing pool.
        futures = []                                                                                                                    # List to track running tasks.
        
        for i in range(1, len(regions)):                                                                                                # Iterate over inner regions only.
            region_points = regions[i]                                                                                                  # Extract contour.
            region_id     = i + 1                                                                                                       # Compute region ID.
            futures.append(executor.submit(generate_region_task, region_points, region_id, main_cloud_size))                            # Submit job.
        
        for future in futures:                                                                                                          # Wait and collect all task results.
            try:                                                                                                                        # Process individual task result.
                result = future.result()                                                                                                # Fetch the job output.
                if result:                                                                                                              # Ensure valid result.
                    interior_clouds.append(result)                                                                                      # Add it to the main collection.
            except Exception as e:                                                                                                      # Handle task crashes.
                pass                                                                                                                    # Ignored error.
    
    return interior_clouds                                                                                                              # Return compiled results.
