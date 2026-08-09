"""
CloudGenerator.classification — Classification

Overview:
    This module is responsible for classifying generated nodes as either 'boundary'
    or 'interior'. It uses advanced geometric operations and spatial indexing to
    ensure accurate and efficient classification, even for complex geometries.

Data conventions:
    None.

Public API:
    classify_nodes      Classify nodes as boundary or interior using Shapely contours.

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
import shapely                                                                                                                          # Vectorized geometry operations.

from shapely.strtree import STRtree                                                                                                     # Spatial index for fast nearest-neighbor search.
from typing import List, Tuple, Optional, Any                                                                                           # Type hinting.
from shapely.geometry import Point, Polygon, LineString                                                                                 # Geometric objects for spatial operations.

from mGFD.cloud_generator.utils.utils import calculate_dynamic_boundary_refinement                                                      # Import tolerance calculator.

def classify_nodes(points: np.ndarray, regions_list: List[int], original_regions_contours: Optional[List[List[Tuple[float, float]]]] = None, cloud_size: Optional[float] = None, inside_regions: bool = False) -> List[str]:
    """
    classify_nodes
    Classify nodes as boundary or interior using Shapely contours for precise geometric operations.
    
    This approach creates LineString/Polygon objects from contours for accurate distance calculations.
    Optimized with STRtree (R-tree spatial index) for O(log N) nearest-neighbor queries,
    significantly reducing classification time for large point clouds.
    
    Input:
        points          m x 2           ndarray         Array of point coordinates.
        regions_list                    List            List of region IDs for each point.
        original_regions_contours       List            List of original contour points for each region.
        cloud_size                      float           Cloud size parameter used in generation for dynamic boundary refinement.
        inside_regions                  bool            If True, applies region-specific boundary logic (holes are not boundaries for region 1).
    
    Output:
        classifications                 List            Node classifications ("boundary" or "interior").
    """
    classifications = []                                                                                                                # Initialize the output classifications list.
    
    # 1. Convert to numpy arrays for better performance
    points_array  = np.array(points, dtype=float)                                                                                       # Convert the list of points to a NumPy float array.
    regions_array = np.array(regions_list)                                                                                              # Convert the list of regions to a NumPy array.
    
    # 2. Calculate dynamic boundary refinement based on point density
    boundary_tolerance = calculate_dynamic_boundary_refinement(points_array, cloud_size)                                                # Compute distance threshold for boundary.
    
    # 3. Pre-calculate global domain boundaries (only once)
    global_min_x, global_max_x = np.min(points_array[:, 0]), np.max(points_array[:, 0])                                                 # Extract minimum and maximum X coordinates.
    global_min_y, global_max_y = np.min(points_array[:, 1]), np.max(points_array[:, 1])                                                 # Extract minimum and maximum Y coordinates.
    
    # 4. Create Shapely geometries from contours
    shapely_contours: List[Any] = []                                                                                                    # List for storing all contour geometries.
    valid_geometries = []                                                                                                               # List for valid geometries specifically for STRtree.
    geom_id_to_index = {}                                                                                                               # Dictionary to map object ID to its contour index.
    
    if original_regions_contours:                                                                                                       # Check if any original contours were provided.
        for idx, contour_points in enumerate(original_regions_contours):                                                                # Iterate over each set of contour points.
            if contour_points and len(contour_points) >= 3:                                                                             # Proceed if there are at least 3 points.
                try:                                                                                                                    # Enter the try block for parsing the coordinates.
                    valid_coords = []                                                                                                   # List to store parsed (x, y) tuples.
                    
                    for cp in contour_points:                                                                                           # Iterate through the points of the contour.
                        if len(cp) >= 2:                                                                                                # Verify that the point has at least 2 dimensions.
                            valid_coords.append((cp[0], cp[1]))                                                                         # Convert contour point coordinates to float tuples.
                    
                    if len(valid_coords) >= 3:                                                                                          # Proceed if we have enough valid coordinates.
                        if valid_coords[0] != valid_coords[-1]:                                                                         # Check if the polygon is open.
                            valid_coords.append(valid_coords[0])                                                                        # Close the polygon by appending the first point.
                        
                        geom = None                                                                                                     # Initialize geometry variable as None.
                        
                        if len(valid_coords) >= 4:                                                                                      # Four points are needed for a valid closed polygon.
                            try:                                                                                                        # Try creating a polygon.
                                polygon = Polygon(valid_coords)                                                                         # Instantiate a Polygon object.
                                if polygon.is_valid:                                                                                    # Check if the polygon doesn't self-intersect.
                                    geom = polygon.boundary                                                                             # Use its boundary LineString for distance calculation.
                                else:                                                                                                   # If polygon creation fails.
                                    geom = LineString(valid_coords)                                                                     # Fallback to creating a LineString.
                            except:                                                                                                     # If an exception occurs.
                                geom = LineString(valid_coords)                                                                         # Fallback to LineString.
                        else:                                                                                                           # If not enough points for a closed polygon.
                            geom = LineString(valid_coords)                                                                             # Fallback to LineString.
                        
                        shapely_contours.append(geom)                                                                                   # Store the resulting geometry.
                        
                        if geom:                                                                                                        # If a geometry was successfully created.
                            valid_geometries.append(geom)                                                                               # Store it in the valid list for the STRtree.
                            geom_id_to_index[id(geom)] = idx                                                                            # Map its unique object ID to its index.
                    else:                                                                                                               # If less than 3 valid coordinates were extracted.
                        shapely_contours.append(None)                                                                                   # Append a None to maintain indexing order.
                except Exception as e:                                                                                                  # Catch any unexpected errors during creation.
                    shapely_contours.append(None)                                                                                       # Append a None to maintain indexing order.
            else:                                                                                                                       # If the contour has fewer than 3 points.
                shapely_contours.append(None)                                                                                           # Append a None to maintain indexing order.
    
    # 5. Create global domain boundary as a rectangle                                                                                   #
    domain_boundary = None                                                                                                              # Initialize the global boundary variable.
    
    try:                                                                                                                                # Try to build the boundary box.
        domain_coords = [                                                                                                               # Define the rectangular box coordinates.
            (global_min_x, global_min_y),                                                                                               # Bottom-left corner.
            (global_max_x, global_min_y),                                                                                               # Bottom-right corner.
            (global_max_x, global_max_y),                                                                                               # Top-right corner.
            (global_min_x, global_max_y),                                                                                               # Top-left corner.
            (global_min_x, global_min_y)                                                                                                # Close back at bottom-left.
        ]                                                                                                                               #
        domain_boundary = LineString(domain_coords)                                                                                     # Create a LineString for the boundary box.
    except:                                                                                                                             # Catch any exceptions during boundary creation.
        domain_boundary = None                                                                                                          # Set it to None if failed.
        
    # 6. Classify each point using Shapely vectorized operations                                                                        #
    points_geoms = shapely.points(points_array[:, 0], points_array[:, 1])                                                               # Create vector of Point geometries.
    
    if domain_boundary is not None:                                                                                                     # Check distance to the domain bounding box.
        domain_distances = shapely.distance(points_geoms, domain_boundary)                                                              # Compute vectorized distances.
        is_boundary_mask = domain_distances <= boundary_tolerance                                                                       # Initialize classification mask.
    else:                                                                                                                               # If no domain boundary.
        is_boundary_mask = np.zeros(len(points_array), dtype=bool)                                                                      # Initialize classification mask as False.
        
    if shapely_contours:                                                                                                                # Check distances to internal Shapely contours.
        for idx, contour in enumerate(shapely_contours):                                                                                # Iterate over all valid contours.
            if contour is None:                                                                                                         # Skip invalid contours.
                continue                                                                                                                # Continue to next contour.
                
            dists = shapely.distance(points_geoms, contour)                                                                             # Vectorized distance from all points to contour.
            
            if inside_regions:                                                                                                          # If region-specific logic is enabled.
                if idx == 0:                                                                                                            # Contour 0 belongs to region 1.
                    mask = (regions_array == 1) & (dists <= boundary_tolerance)                                                         # Check region 1 points.
                else:                                                                                                                   # Other contours belong to respective sub-regions.
                    mask = (regions_array == (idx + 1)) & (dists <= boundary_tolerance)                                                 # Check sub-region points.
                is_boundary_mask |= mask                                                                                                # Update global mask.
            else:                                                                                                                       # If no region logic.
                is_boundary_mask |= (dists <= boundary_tolerance)                                                                       # Apply distance check globally.
                
    classifications = np.where(is_boundary_mask, "boundary", "interior").tolist()                                                       # Construct list of classification strings.
    
    # 8. Log classification statistics
    boundary_count = sum(1 for c in classifications if c == "boundary")                                                                 # Count the number of boundary nodes.
    interior_count = len(classifications) - boundary_count                                                                              # Compute the remaining interior nodes.
    
    return classifications                                                                                                              # Return the completed classification array.
