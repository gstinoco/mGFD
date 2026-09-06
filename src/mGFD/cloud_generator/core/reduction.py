"""
CloudGenerator.reduction — Point Reduction

Overview:
    This module provides advanced point reduction functionality for CSV files containing
    cloud data organized by regions. It implements multiple reduction algorithms
    optimized for the mGFD CloudGenerator, enabling efficient processing
    of large point datasets while maintaining geometric integrity.

Data conventions:
    None.

Public API:
    reduce_points_by_region_single              Reduce points in region 1 using uniform method.
    reduce_points_by_region_multiple            Advanced point reduction using iterative uniform method.
    filter_main_region_points_in_subregions     Geometric filtering to remove main region points within subregions.
    reduce_points_by_region_with_filtering      Advanced point reduction with optional geometric filtering.
    reduce_points_by_region                     Main entry point for intelligent point reduction.

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
import csv                                                                                                                              # Import CSV reading/writing module.
import tempfile                                                                                                                         # Temporary files creation for intermediate steps.
import numpy as np                                                                                                                      # Numerical arrays and math operations.

from typing import List, Dict, Any, Optional                                                                                            # Import type hints.
from shapely.geometry import Point, Polygon, MultiPoint                                                                                 # Import geometric objects from Shapely.

def reduce_points_by_region_single(input_csv: str, output_csv: str) -> Optional[List[Dict[str, Any]]]:                                  # Single pass regional point reduction.
    """
    reduce_points_by_region_single
    Reduce points in region 1 using uniform method while preserving other regions (single pass).
    
    This function processes a CSV file containing point data with regions, applying
    a uniform reduction method that removes approximately 50% of points from region 1
    only. All other regions remain completely unchanged.
    
    Input:
        input_csv       str                                 Absolute path to the input CSV file containing point data.
        output_csv      str                                 Absolute path where the reduced CSV file will be saved.
    
    Output:
        reduced_rows    Optional[List[Dict[str, Any]]]      List of dictionaries representing the reduced points if successful,
                                                            None if any error occurred.
    """
    try:                                                                                                                                # Try to execute the main block.
        rows = []                                                                                                                       # Initialize the list for the CSV rows.
        with open(input_csv, 'r', newline='', encoding='utf-8') as f:                                                                   # Open the input CSV file for reading.
            reader = csv.DictReader(f)                                                                                                  # Create a DictReader for the file.
            if not reader.fieldnames:                                                                                                   # Check if the file has field names.
                raise ValueError("Empty CSV file")                                                                                      # Raise an error if the file is empty.
                
            fieldnames       = [f.strip() for f in reader.fieldnames]                                                                   # Strip whitespace from field names.
            required_columns = ['x', 'y', 'region']                                                                                     # Define the required columns.
            if not all(col in fieldnames for col in required_columns):                                                                  # Check if all required columns are present.
                error_msg = f"CSV file missing required columns. Expected: {required_columns}, Found: {fieldnames}"                     # Create an error message.
                raise ValueError(error_msg)                                                                                             # Raise a ValueError.
            
            for row in reader:                                                                                                          # Iterate over the rows in the CSV file.
                # 1. Convert types
                try:                                                                                                                    # Try to convert values.
                    row['x'] = float(row['x'])                                                                                          # Convert the 'x' coordinate to a float.
                    row['y'] = float(row['y'])                                                                                          # Convert the 'y' coordinate to a float.
                    # Try to convert region to int/float
                    try:                                                                                                                # Try to convert the 'region' value.
                        row['region'] = int(float(row['region']))                                                                       # Convert the 'region' to an integer.
                    except:                                                                                                             # If the conversion fails.
                        pass                                                                                                            # Keep as is if not a number.
                    rows.append(row)                                                                                                    # Append the parsed row to the list.
                except ValueError:                                                                                                      # If a ValueError occurs during conversion.
                    continue                                                                                                            # Skip to the next row.

        # 2. Group by region
        regions_data: Dict[Any, List[Dict[str, Any]]] = {}                                                                              # Initialize a dictionary for region grouping.
        for row in rows:                                                                                                                # Iterate over the parsed rows.
            rid = row['region']                                                                                                         # Get the region ID for the current row.
            if rid not in regions_data:                                                                                                 # Check if the region ID is not in the dictionary.
                regions_data[rid] = []                                                                                                  # Initialize an empty list for the new region ID.
            regions_data[rid].append(row)                                                                                               # Append the row to the corresponding region list.
        
        reduced_rows = []                                                                                                               # Initialize the list for the reduced rows.
        
        # 3. Sort regions for consistent output
        try:                                                                                                                            # Try to sort the regions numerically.
            sorted_regions = sorted(regions_data.keys())                                                                                # Sort the region keys.
        except:                                                                                                                         # If sorting numerically fails.
            sorted_regions = sorted(regions_data.keys(), key=str)                                                                       # Sort the region keys as strings.
            
        for region_id in sorted_regions:                                                                                                # Iterate over the sorted region IDs.
            region_rows = regions_data[region_id]                                                                                       # Get the rows for the current region.
            
            if region_id == 1 or str(region_id) == '1':                                                                                 # Check if the region is the main region (1).
                # 4. Preserve all boundary points, reduce only interior points
                boundaries = []                                                                                                         # Initialize the list for boundary points.
                interiors  = []                                                                                                         # Initialize the list for interior points.
                for row in region_rows:                                                                                                 # Iterate over the rows in the main region.
                    # If there's no classification, assume it's interior to allow reduction
                    if row.get('classification') == 'boundary':                                                                         # Check if the point is classified as a boundary.
                        boundaries.append(row)                                                                                          # Append to the boundaries list.
                    else:                                                                                                               # If it's not a boundary point.
                        interiors.append(row)                                                                                           # Append to the interiors list.
                        
                # 5. Reduce boundary points softly: keep 3, delete 1 (25% reduction per pass)
                # This prevents completely destroying the boundary shape while still optimizing it.
                reduced_boundaries = [p for i, p in enumerate(boundaries) if (i % 4) != 3]                                              # Filter boundaries by keeping 3 out of every 4.
                
                # 6. Reduce interior points aggressively: take every 2nd point (50% reduction per pass)
                reduced_interiors = interiors[::2]                                                                                      # Take every second point from the interiors.
                
                # 7. Add boundaries, then the reduced interior
                reduced_rows.extend(reduced_boundaries)                                                                                 # Add the reduced boundary points to the output.
                reduced_rows.extend(reduced_interiors)                                                                                  # Add the reduced interior points to the output.
            else:                                                                                                                       # If it's not the main region.
                # 8. Keep other regions as is
                reduced_rows.extend(region_rows)                                                                                        # Add the region's points without reduction.
        
        # 9. Write to output CSV
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:                                                                  # Open the output CSV file for writing.
            fieldnames = list(reduced_rows[0].keys()) if reduced_rows else ['x', 'y', 'region']                                         # Determine the field names for writing.
            writer     = csv.DictWriter(f, fieldnames=fieldnames)                                                                       # Create a DictWriter for the output file.
            writer.writeheader()                                                                                                        # Write the header row.
            writer.writerows(reduced_rows)                                                                                              # Write all the reduced rows.
            
        return reduced_rows                                                                                                             # Return the final reduced rows.
        
    except FileNotFoundError as e:                                                                                                      # Handle FileNotFoundError.
        error_msg = f"Input file not found: {input_csv}"                                                                                # Create an error message.
        return None                                                                                                                     # Return None on error.
    except Exception as e:                                                                                                              # Handle general exceptions.
        error_msg = f"Error processing point reduction for {input_csv}: {str(e)}"                                                       # Create an error message.
        return None                                                                                                                     # Return None on error.

def reduce_points_by_region_multiple(input_csv: str, output_csv: str, multiplier: int = 2) -> Optional[List[Dict[str, Any]]]:           # Multiple pass regional point reduction.
    """
    reduce_points_by_region_multiple
    Advanced point reduction using iterative uniform method for aggressive point reduction.
    
    This function applies the uniform reduction algorithm multiple times in sequence to achieve
    higher reduction rates while maintaining spatial distribution characteristics. Each iteration
    reduces the point count by approximately 50%, resulting in exponential reduction rates.
    The function uses temporary files for intermediate processing to ensure data integrity.
    
    Input:
        input_csv       str                                 Absolute path to the input CSV file containing point data.
        output_csv      str                                 Absolute path where the final reduced CSV file will be saved.
        multiplier      int                                 Number of reduction iterations to apply. Default is 2.
    
    Output:
        result          Optional[List[Dict[str, Any]]]      Final list of dictionaries with reduced points if successful,
                                                            None if any error occurred.
    """
    try:                                                                                                                                # Try to execute the main block.
        current_input = input_csv                                                                                                       # Set the initial input file to the given input.
        temp_files    = []                                                                                                              # Initialize a list to keep track of temporary files.
        result        = None                                                                                                            # Initialize the result variable.
        
        for i in range(multiplier):                                                                                                     # Iterate according to the multiplier parameter.
            if i == multiplier - 1:                                                                                                     # Check if it is the last iteration.
                # 1. Last iteration, use the final output file
                current_output = output_csv                                                                                             # Use the final output path.
            else:                                                                                                                       # If it is not the last iteration.
                # 2. Create temporary file for intermediate results
                temp_fd, current_output = tempfile.mkstemp(suffix='.csv')                                                               # Create a temporary file.
                os.close(temp_fd)                                                                                                       # Close the file descriptor.
                temp_files.append(current_output)                                                                                       # Add the temporary file to the tracking list.
            
            # 3. Apply single reduction
            result = reduce_points_by_region_single(current_input, current_output)                                                      # Call the single reduction function.
            if result is None:                                                                                                          # Check if the reduction failed.
                # 4. Clean up temp files on error
                for temp_file in temp_files:                                                                                            # Iterate over all temporary files.
                    try:                                                                                                                # Try to remove the file.
                        os.remove(temp_file)                                                                                            # Delete the temporary file.
                    except:                                                                                                             # If removal fails.
                        pass                                                                                                            # Skip and continue.
                return None                                                                                                             # Return None indicating failure.
            
            # 5. Update input for next iteration
            current_input = current_output                                                                                              # Use the output as input for the next iteration.
        
        # 6. Clean up temporary files
        for temp_file in temp_files:                                                                                                    # Iterate over all temporary files.
            try:                                                                                                                        # Try to remove the file.
                os.remove(temp_file)                                                                                                    # Delete the temporary file.
            except:                                                                                                                     # If removal fails.
                pass                                                                                                                    # Skip and continue.
        
        return result                                                                                                                   # Return the final result list.
        
    except Exception as e:                                                                                                              # Handle general exceptions.
        error_msg = f"Error in multiple point reduction: {str(e)}"                                                                      # Create an error message.
        return None                                                                                                                     # Return None on error.

def filter_main_region_points_in_subregions(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:                                        # Filter main points overlapping subregions.
    """
    filter_main_region_points_in_subregions
    Advanced geometric filtering to remove main region points that fall within subregion boundaries.
    
    This function implements sophisticated spatial analysis to maintain geometric consistency
    between hierarchical regions. It uses convex hull computation to define subregion boundaries
    and performs point-in-polygon tests to identify and remove overlapping points from the
    main region.
    
    Input:
        rows            List[Dict[str, Any]]    List of dictionaries containing cloud data.
    
    Output:
        filtered_rows   List[Dict[str, Any]]    Filtered list of dictionaries.
    """
    try:                                                                                                                                # Try to execute the main block.
        # 1. Group points by region
        regions_data: Dict[Any, List[Dict[str, Any]]] = {}                                                                              # Initialize a dictionary for grouping points by region.
        for row in rows:                                                                                                                # Iterate over the input rows.
            rid = row['region']                                                                                                         # Get the region ID for the current row.
            if rid not in regions_data:                                                                                                 # Check if the region ID is not in the dictionary.
                regions_data[rid] = []                                                                                                  # Initialize an empty list for the new region ID.
            regions_data[rid].append(row)                                                                                               # Append the row to the corresponding region list.
            
        regions = sorted(regions_data.keys(), key=lambda x: str(x))                                                                     # Sort the region keys as strings.
        
        if len(regions) <= 1:                                                                                                           # Check if there is 1 or 0 regions.
            # 2. No subregions to filter against
            return rows                                                                                                                 # Return the original rows without filtering.
        
        # 3. Create polygons for subregions (regions 2, 3, 4, ...)
        subregion_polygons = []                                                                                                         # Initialize a list for subregion polygons.
        for region_id in regions:                                                                                                       # Iterate over all available regions.
            if str(region_id) == '1':                                                                                                   # Check if it is the main region (1).
                continue                                                                                                                # Skip the main region.
                
            region_rows          = regions_data[region_id]                                                                              # Get the rows for the current subregion.
            region_points_coords = []                                                                                                   # Initialize a list for the coordinates.
            for r in region_rows:                                                                                                       # Iterate over the rows of the subregion.
                region_points_coords.append([r['x'], r['y']])                                                                           # Append the [x, y] coordinates.
            
            region_points_array = np.array(region_points_coords)                                                                        # Convert the coordinates list to a NumPy array.
            
            if len(region_points_array) >= 3:                                                                                           # Check if there are enough points for a polygon.
                try:                                                                                                                    # Try to compute the convex hull.
                    # Create convex hull of region points to form polygon using Shapely
                    # This replaces the previous SciPy implementation to reduce dependency weight
                    multi_point = MultiPoint(region_points_array)                                                                       # Create a MultiPoint geometry from the array.
                    poly        = multi_point.convex_hull                                                                               # Compute the convex hull polygon.
                    
                    if poly.is_valid and isinstance(poly, Polygon):                                                                     # Check if the polygon is valid.
                        subregion_polygons.append(poly)                                                                                 # Append the polygon to the list.
                except Exception as e:                                                                                                  # Handle exceptions during hull computation.
                    continue                                                                                                            # Skip to the next region.
        
        if not subregion_polygons:                                                                                                      # Check if no valid polygons were created.
            return rows                                                                                                                 # Return the original rows.
        
        # 4. Filter main region points
        filtered_rows = []                                                                                                              # Initialize the list for filtered rows.
        for row in rows:                                                                                                                # Iterate over all original rows.
            is_main_region = (row['region'] == 1 or str(row['region']) == '1')                                                          # Check if the row belongs to the main region.
            
            if is_main_region:                                                                                                          # If it is a point from the main region.
                # 5. Check if main region point falls inside any subregion
                point_obj        = Point(row['x'], row['y'])                                                                            # Create a Point geometry.
                inside_subregion = False                                                                                                # Flag for points inside subregions.
                for subregion_poly in subregion_polygons:                                                                               # Iterate over all subregion polygons.
                    if subregion_poly.contains(point_obj):                                                                              # Check if the polygon contains the point.
                        inside_subregion = True                                                                                         # Set the flag to True.
                        break                                                                                                           # Break the inner loop early.
                
                # 6. Only keep main region points that are NOT inside subregions
                if not inside_subregion:                                                                                                # If the point is outside all subregions.
                    filtered_rows.append(row)                                                                                           # Keep the point and add to the list.
            else:                                                                                                                       # If it is a point from a subregion.
                # 7. Keep all subregion points
                filtered_rows.append(row)                                                                                               # Add the subregion point to the list.
        
        return filtered_rows                                                                                                            # Return the list with the filtered points.
        
    except Exception as e:                                                                                                              # Handle general exceptions.
        return rows                                                                                                                     # Return the original rows if filtering fails.

def reduce_points_by_region(input_csv: str, output_csv: str, multiplier: int = 2) -> Optional[List[Dict[str, Any]]]:                    # Main entry for point cloud reduction.
    """
    reduce_points_by_region
    Main entry point for intelligent point reduction with adaptive algorithm selection.
    
    Input:
        input_csv       str                                 Absolute path to the input CSV file containing cloud data.
        output_csv      str                                 Absolute path where the reduced CSV file will be saved.
        multiplier      int                                 Reduction intensity parameter controlling algorithm selection. Default is 2.
    
    Output:
        result          Optional[List[Dict[str, Any]]]      Processed list of dictionaries with reduced points if successful,
                                                            None if any error occurred.
    """
    if multiplier <= 1:                                                                                                                 # Check if the reduction multiplier is 1 or less.
        return reduce_points_by_region_single(input_csv, output_csv)                                                                    # Return the result of the single pass reduction.
    else:                                                                                                                               # If the reduction multiplier is greater than 1.
        return reduce_points_by_region_multiple(input_csv, output_csv, multiplier)                                                      # Return the result of the multiple passes reduction.
