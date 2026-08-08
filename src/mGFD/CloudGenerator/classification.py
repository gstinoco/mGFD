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
import logging                                                                      # Logging module for informative and warning messages.
from shapely.geometry import Point, Polygon, LineString                             # Geometric objects for spatial operations.
from shapely.strtree import STRtree                                                 # Spatial index for fast nearest-neighbor search.
from .utils import calculate_dynamic_boundary_refinement                            # Import tolerance calculator.

def classify_nodes(points: np.ndarray, regions_list: list[int], original_regions_contours: list[list[tuple[float, float]]] = None, cloud_size: float = None, inside_regions: bool = False) -> list[str]:
    """
    Classify nodes as boundary or interior using Shapely contours for precise geometric operations.
    
    This approach creates LineString/Polygon objects from contours for accurate distance calculations.
    Optimized with STRtree (R-tree spatial index) for O(log N) nearest-neighbor queries,
    significantly reducing classification time for large point clouds.
    
    Input:
        points                      ndarray         Array of point coordinates.
        regions_list                list            List of region IDs for each point.
        original_regions_contours   list            List of original contour points for each region.
        cloud_size                  float           Cloud size parameter used in generation for dynamic boundary refinement.
        inside_regions              bool            If True, applies region-specific boundary logic (holes are not boundaries for region 1).
    
    Output:
        classifications             list            Node classifications ("boundary" or "interior").
    """
    classifications = []                                                            # Initialize the output classifications list.
    
    # 1. Convert to numpy arrays for better performance
    points_array = np.array(points, dtype=float)                                    # Convert the list of points to a NumPy float array.
    regions_array = np.array(regions_list)                                          # Convert the list of regions to a NumPy array.
    
    # 2. Calculate dynamic boundary refinement based on point density
    boundary_tolerance = calculate_dynamic_boundary_refinement(points_array, cloud_size) # Compute distance threshold for boundary.
    
    # 3. Pre-calculate global domain boundaries (only once)
    global_min_x, global_max_x = np.min(points_array[:, 0]), np.max(points_array[:, 0]) # Extract minimum and maximum X coordinates.
    global_min_y, global_max_y = np.min(points_array[:, 1]), np.max(points_array[:, 1]) # Extract minimum and maximum Y coordinates.
    
    # 4. Create Shapely geometries from contours
    shapely_contours = []                                                           # List for storing all contour geometries.
    valid_geometries = []                                                           # List for valid geometries specifically for STRtree.
    geom_id_to_index = {}                                                           # Dictionary to map object ID to its contour index.
    
    if original_regions_contours:                                                   # Check if any original contours were provided.
        for idx, contour_points in enumerate(original_regions_contours):            # Iterate over each set of contour points.
            if contour_points and len(contour_points) >= 3:                         # Proceed if there are at least 3 points.
                try:                                                                # Enter the try block for parsing the coordinates.
                    valid_coords = []                                               # List to store parsed (x, y) tuples.
                    for cp in contour_points:                                       # Iterate through the points of the contour.
                        if len(cp) >= 2:                                            # Verify that the point has at least 2 dimensions.
                            valid_coords.append((float(cp[0]), float(cp[1])))       # Convert contour point coordinates to float tuples.
                    
                    if len(valid_coords) >= 3:                                      # Proceed if we have enough valid coordinates.
                        if valid_coords[0] != valid_coords[-1]:                     # Check if the polygon is open.
                            valid_coords.append(valid_coords[0])                    # Close the polygon by appending the first point.
                        
                        geom = None                                                 # Initialize geometry variable as None.
                        if len(valid_coords) >= 4:                                  # Four points are needed for a valid closed polygon.
                            try:                                                    # Try creating a polygon.
                                polygon = Polygon(valid_coords)                     # Instantiate a Polygon object.
                                if polygon.is_valid:                                # Check if the polygon doesn't self-intersect.
                                    geom = polygon.boundary                         # Use its boundary LineString for distance calculation.
                                else:                                               # If polygon creation fails.
                                    geom = LineString(valid_coords)                 # Fallback to creating a LineString.
                            except:                                                 # If an exception occurs.
                                geom = LineString(valid_coords)                     # Fallback to LineString.
                        else:                                                       # If not enough points for a closed polygon.
                            geom = LineString(valid_coords)                         # Fallback to LineString.
                        
                        shapely_contours.append(geom)                               # Store the resulting geometry.
                        if geom:                                                    # If a geometry was successfully created.
                            valid_geometries.append(geom)                           # Store it in the valid list for the STRtree.
                            geom_id_to_index[id(geom)] = idx                        # Map its unique object ID to its index.
                    else:                                                           # If less than 3 valid coordinates were extracted.
                        shapely_contours.append(None)                               # Append a None to maintain indexing order.
                except Exception as e:                                              # Catch any unexpected errors during creation.
                    logging.warning(f"Failed to create Shapely geometry for contour: {e}") # Log a warning with the error details.
                    shapely_contours.append(None)                                   # Append a None to maintain indexing order.
            else:                                                                   # If the contour has fewer than 3 points.
                shapely_contours.append(None)                                       # Append a None to maintain indexing order.
    
    # 5. Build STRtree for efficient spatial queries
    spatial_index = None                                                            # Initialize the spatial index variable.
    if valid_geometries:                                                            # Check if there are valid geometries to index.
        try:                                                                        # Try block to build the STRtree.
            spatial_index = STRtree(valid_geometries)                               # Construct the R-tree index.
        except Exception as e:                                                      # If tree construction fails.
            logging.warning(f"Failed to build STRtree: {e}")                        # Log a warning message.
    
    # 6. Create global domain boundary as a rectangle
    domain_boundary = None                                                          # Initialize the global boundary variable.
    try:                                                                            # Try to build the boundary box.
        domain_coords = [                                                           # Define the rectangular box coordinates.
            (global_min_x, global_min_y),                                           # Bottom-left corner.
            (global_max_x, global_min_y),                                           # Bottom-right corner.
            (global_max_x, global_max_y),                                           # Top-right corner.
            (global_min_x, global_max_y),                                           # Top-left corner.
            (global_min_x, global_min_y)                                            # Close back at bottom-left.
        ]
        domain_boundary = LineString(domain_coords)                                 # Create a LineString for the boundary box.
    except:                                                                         # Catch any exceptions during boundary creation.
        domain_boundary = None                                                      # Set it to None if failed.
    
    # 7. Classify each point using Shapely geometric operations
    for i in range(len(points_array)):                                              # Iterate over the entire point cloud.
        x, y = points_array[i]                                                      # Extract the X and Y coordinates.
        region_id = regions_array[i]                                                # Extract the region ID of the current point.
        is_boundary = False                                                         # Default boundary classification state.
        
        try:                                                                        # Attempt to classify the point geometrically.
            point_geom = Point(x, y)                                                # Create a Shapely Point geometry.
            min_distance_to_contour = float('inf')                                  # Set the minimum distance to infinity.
            
            # Primary method: Check distance to Shapely contours
            if shapely_contours:                                                    # Check if there are contours loaded.
                if spatial_index:                                                   # Use the R-tree if available.
                    # query returns candidates (geometries) that *might* be near
                    candidates = spatial_index.query(point_geom)                    # Perform a spatial query to get near geometries.
                    
                    final_candidates = []                                           # List for candidates handling STRtree API versions.
                    
                    if hasattr(candidates, '__iter__') and len(candidates) > 0:     # Check if the query returned an iterable result.
                        first_item = list(candidates)[0]                            # Get the first candidate item.
                        if isinstance(first_item, (int, np.integer)):               # Newer Shapely API returns indices.
                            final_candidates = [valid_geometries[i] for i in candidates] # Resolve geometry from the valid geometries.
                        elif isinstance(first_item, (Polygon, LineString, Point)):  # Older Shapely API returns geometries.
                             final_candidates = list(candidates)                    # Convert iterator to list directly.
                        else:                                                       # If it is neither.
                             final_candidates = valid_geometries                    # Default to evaluating all geometries.
                    else:                                                           # If no specific candidates were found or API mismatch.
                        final_candidates = valid_geometries                         # Evaluate all valid geometries.
                    
                    if not final_candidates:                                        # Ensure the candidate list isn't empty.
                        final_candidates = valid_geometries                         # Fallback to checking all valid geometries.
                        
                    for contour in final_candidates:                                # Iterate over the localized geometries.
                        # Check if this contour is a valid boundary for the current point's region
                        if inside_regions:                                          # If region-specific logic is enabled.
                            contour_idx = geom_id_to_index.get(id(contour))         # Fetch the original index of the contour.
                            if contour_idx is not None:                             # Check if the mapping was found.
                                if region_id == 1 and contour_idx != 0:             # Rule: Region 1 checks contour 0 only.
                                    continue                                        # Skip other contours for Region 1.
                                if region_id > 1 and contour_idx != (region_id - 1):# Rule: Sub-regions check their own contour.
                                    continue                                        # Skip non-matching contours.
                        
                        dist = point_geom.distance(contour)                         # Calculate exact shortest distance to contour.
                        if isinstance(dist, (np.ndarray, list)):                    # In some environments, distance might be an array.
                            dist = np.min(dist)                                     # Reduce an array to its scalar minimum.
                        min_distance_to_contour = min(min_distance_to_contour, float(dist)) # Update the shortest distance.
                else:                                                               # If STRtree is not available.
                    # Fallback to checking all contours
                    for idx, contour in enumerate(shapely_contours):                # Perform an exhaustive iteration over contours.
                        if contour is not None:                                     # Skip invalid contours.
                            if inside_regions:                                      # If region-specific logic is enabled.
                                if region_id == 1 and idx != 0:                     # Rule: Region 1 checks contour 0 only.
                                    continue                                        # Skip other contours.
                                if region_id > 1 and idx != (region_id - 1):        # Rule: Sub-regions check their own contour.
                                    continue                                        # Skip non-matching contours.
                                    
                            min_distance_to_contour = min(min_distance_to_contour, point_geom.distance(contour)) # Calculate distance.
                
                if min_distance_to_contour <= boundary_tolerance:                   # Classify as boundary if close to any contour.
                    is_boundary = True                                              # Mark point as boundary.
            
            # Fallback method: Check distance to global domain boundary
            if not is_boundary and domain_boundary is not None:                     # If not boundary yet, test against the bounding box.
                distance_to_domain = point_geom.distance(domain_boundary)           # Measure distance to the domain bounding box.
                if distance_to_domain <= boundary_tolerance:                        # Check against the same distance tolerance.
                    is_boundary = True                                              # Mark point as boundary.
            
        except Exception as e:                                                      # Handle any calculation error for a given point.
            # Fallback to simple coordinate-based classification if Shapely fails
            logging.warning(f"Shapely operation failed for point ({x}, {y}): {e}")  # Log the fallback condition.
            if (abs(x - global_min_x) < boundary_tolerance or                       # Check left edge.
                abs(x - global_max_x) < boundary_tolerance or                       # Check right edge.
                abs(y - global_min_y) < boundary_tolerance or                       # Check bottom edge.
                abs(y - global_max_y) < boundary_tolerance):                        # Check top edge.
                is_boundary = True                                                  # Set boundary flag to True.
        
        classifications.append("boundary" if is_boundary else "interior")           # Add the final classification label for this point.
    
    # 8. Log classification statistics
    boundary_count = sum(1 for c in classifications if c == "boundary")             # Count the number of boundary nodes.
    interior_count = len(classifications) - boundary_count                          # Compute the remaining interior nodes.
    logging.info(f"Node classification using Shapely: {boundary_count} boundary, {interior_count} interior") # Output classification counts.
    logging.info(f"Boundary tolerance used: {boundary_tolerance}")                  # Output the boundary distance threshold.
    
    return classifications                                                          # Return the completed classification array.
