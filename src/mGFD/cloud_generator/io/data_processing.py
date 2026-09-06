"""
CloudGenerator.data_processing — Data Processing

Overview:
    This module handles the loading, parsing, and validation of input data for the
    Cloud Generation system. It focuses on reading CSV boundary files and structuring
    the data for downstream processing.

Data conventions:
    None.

Public API:
    load_regions        Load boundary region data from a CSV file.

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
import csv                                                                                                                              # Module to read and write CSV files.

from typing import List, Tuple, Dict, Any                                                                                               # Import type hints.

def load_regions(csv_file: str) -> List[List[Tuple[float, float]]]:                                                                     # Load boundary regions from CSV.
    """
    Load boundary region data from a CSV file.
    
    Parses a CSV file containing boundary contours, validates the presence
    of 'x' and 'y' columns, and groups the coordinates by region ID if
    a 'region' column is provided.
    
    Input:
        csv_file    str                                 The absolute or relative path to the input CSV file.
        
    Output:
        regions     List[List[Tuple[float, float]]]     A list of regions, where each region is a list of (x, y) coordinate
                                                        tuples representing the boundary contour of that region.
    """
    try:                                                                                                                                # Start of error handling block.
        regions_dict: Dict[Any, List[Tuple[float, float]]] = {}                                                                         # Dictionary to group points by region ID.
        single_region_points = []                                                                                                       # Fallback list if no regions are specified.
        has_region_column    = False                                                                                                    # Boolean flag for the 'region' column.
        
        with open(csv_file, 'r', newline='', encoding='utf-8-sig') as f:                                                                # Open the CSV file safely.
            reader = csv.DictReader(f)                                                                                                  # Initialize the CSV dictionary reader.
            
            # 1. Check required columns
            if not reader.fieldnames:                                                                                                   # Check if the CSV has a header row.
                return []                                                                                                               # Return an empty list.

            fieldnames = [f.strip() for f in reader.fieldnames]                                                                         # Extract and clean field names.
            if 'x' not in fieldnames or 'y' not in fieldnames:                                                                          # Validate that 'x' and 'y' columns exist.
                return []                                                                                                               # Return an empty list.
                
            has_region_column    = 'region' in fieldnames                                                                               # Check and set flag if 'region' column exists.
            
            # 2. Parse rows
            for row in reader:                                                                                                          # Loop through every row in the CSV file.
                try:                                                                                                                    # Inner error handling block for parsing.
                    # Clean and convert coordinates
                    x_str = row['x'].strip() if row.get('x') else ''                                                                    # Extract and clean x coordinate string.
                    y_str = row['y'].strip() if row.get('y') else ''                                                                    # Extract and clean y coordinate string.
                    
                    if not x_str or not y_str:                                                                                          # Ensure both coordinates are present.
                        continue                                                                                                        # Skip the current row if missing.
                        
                    x = float(x_str)                                                                                                    # Convert x string to a float value.
                    y = float(y_str)                                                                                                    # Convert y string to a float value.
                    
                    if has_region_column and row.get('region'):                                                                         # Check if region data is available.
                        # Handle region column
                        region_val_str = str(row['region']).strip()                                                                     # Extract and clean the region identifier.
                        try:                                                                                                            # Try parsing region as a number.
                            region_val = float(region_val_str)                                                                          # Parse the region string as a float.
                            region_id: Any = int(region_val)                                                                            # Cast the float region to an integer.
                        except (ValueError, TypeError):                                                                                 # Handle parsing errors.
                            region_id = region_val_str                                                                                  # Fallback to using the raw string as region ID.
                            
                        if region_id not in regions_dict:                                                                               # Check if region ID is new.
                            regions_dict[region_id] = []                                                                                # Initialize an empty list for the new region.
                        regions_dict[region_id].append((x, y))                                                                          # Add the point tuple to its corresponding region.
                    else:                                                                                                               # If no region column exists.
                        single_region_points.append((x, y))                                                                             # Add the point tuple to the default region.
                        
                except ValueError:                                                                                                      # Handle conversion errors inside the row.
                    continue                                                                                                            # Skip the invalid row and continue parsing.

        # 3. Consolidate regions
        regions = []                                                                                                                    # Initialize the final regions list.
        if has_region_column and regions_dict:                                                                                          # If multiple regions were processed.
            try:                                                                                                                        # Try numerical sorting.
                sorted_keys = sorted(regions_dict.keys())                                                                               # Sort region keys numerically.
            except TypeError:                                                                                                           # If numerical sorting fails (mixed types).
                sorted_keys = sorted(regions_dict.keys(), key=str)                                                                      # Sort region keys lexicographically.
                
            for key in sorted_keys:                                                                                                     # Loop through the sorted region keys.
                regions.append(regions_dict[key])                                                                                       # Append the points of each sorted region.
        else:                                                                                                                           # If only one region exists.
            if single_region_points:                                                                                                    # Check if the default region contains points.
                regions.append(single_region_points)                                                                                    # Append the single region to the list.
        
        return regions                                                                                                                  # Return the list of regions.
        
    except Exception as e:                                                                                                              # Handle general exceptions.
        raise                                                                                                                           # Re-raise the exception.
