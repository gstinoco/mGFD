"""
CloudGenerator.generator — Point Cloud Generation System

Overview:
    This module provides comprehensive functionality for generating optimized point clouds
    from CSV boundary data for use with the meshless Generalized Finite Differences (mGFD) method.
    The module implements two advanced distribution algorithms with intelligent point reduction
    and multi-region support for complex geometries.

Data conventions:
    None.

Public API:
    generate_cloud_natural      Generate point cloud using Natural Distribution algorithm with Poisson Disk Sampling.
    generate_cloud_regular      Generate point cloud using Regular Distribution algorithm with grid-based approach.

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

import numpy as np                                                                                                                      # Core numerical operations.

from typing import Dict, Optional, Callable                                                                                             # Type hinting.

from mGFD.cloud_generator.io.export import export_to_csv                                                                                # Import the CSV export function.
from mGFD.cloud_generator.io.data_processing import load_regions                                                                        # Import load_regions function.
from mGFD.cloud_generator.utils.utils import calculate_cloud_size                                                                       # Import cloud size calculator.
from mGFD.cloud_generator.core.classification import classify_nodes                                                                     # Import classification function.
from mGFD.cloud_generator.viz.visualization import create_visualization                                                                 # Import visualization function.
from mGFD.cloud_generator.core.point_generation import (                                                                                # Import specific generation functions.
    generate_region_cloud_with_holes,                                                                                                   # Regular grid region generation.
    generate_region_cloud_with_holes_poisson,                                                                                           # Poisson region generation.
    generate_interior_regions_clouds,                                                                                                   # Interior regions regular generation.
    generate_interior_regions_clouds_poisson,                                                                                           # Interior regions Poisson generation.
)                                                                                                                                       # Execute statement.

def _generate_cloud_core(csv_file: str, output_file: str, method_name: str,                                                             # Cloud generation core pipeline.
                         main_region_strategy: Callable, interior_regions_strategy: Callable,                                           # Region generation strategies.
                         inside_regions: bool = False, cloud_size: Optional[float] = None,                                              # Internal region & spacing options.
                         density_multiplier: float = 1.0, save: bool = True, show: bool = False) -> Dict:                               # Execute cloud generation core.
    """
    _generate_cloud_core
    Core generation logic shared between different distribution methods.
    
    Input:
        csv_file                    str         Path to boundary CSV.
        output_file                 str         Output filename base.
        method_name                 str         Name of generation method ('Natural' or 'Regular').
        main_region_strategy        Callable    Function to generate main region cloud.
        interior_regions_strategy   Callable    Function to generate interior regions cloud.
        inside_regions              bool        Whether to process interior holes.
        cloud_size                  Optional    Override distance parameter.
        density_multiplier          float       Density scaling factor.
        save                        bool        Whether to save PNG/SVG visualization files to disk.
        show                        bool        Whether to show interactive plot window.
        
    Output:
        results                     Dict        Comprehensive results dictionary.
    """
    temp_reduced_file = None                                                                                                            # Initialize temporary file variable.
    
    try:                                                                                                                                # Start main try-except block.
        # 1. Initialization and logging
        working_csv_file = csv_file                                                                                                     # Set the working CSV file path.
        
        # 2. Data loading
        regions          = load_regions(working_csv_file)                                                                               # Parse and load regions from CSV.
        
        if not regions:                                                                                                                 # Check if regions were loaded properly.
            return {"success": False, "error": "Failed to load regions from CSV"}                                                       # Return failure if no regions.
        
        main_region  = regions[0]                                                                                                       # Define the main (outer) region.
        hole_regions = regions[1:] if len(regions) > 1 else []                                                                          # Define hole (inner) regions.
        
        # Calculate explicit cloud size with density multiplier
        
        # 3. Cloud size calculation
        if cloud_size is None:                                                                                                          # If cloud_size is not provided.
            base_cloud_size = calculate_cloud_size(main_region)                                                                         # Calculate the base cloud size.
            # Higher density multiplier = smaller cloud size = more points
            cloud_size      = base_cloud_size / density_multiplier                                                                      # Adjust cloud size with multiplier.

        # 4. Main region generation
        # Pass inside_regions to control hole boundary generation
        main_boundary, main_interior, actual_cloud_size = main_region_strategy(main_region, hole_regions, cloud_size, inside_regions)   # Generate main region points.
        
        if main_boundary is None or main_interior is None:                                                                              # Ensure points were generated successfully.
            return {"success": False, "error": "Failed to generate cloud for main region"}                                              # Return failure if not.
        
        # 5. Output compilation
        all_points  = []                                                                                                                # List to hold all points across regions.
        all_regions = []                                                                                                                # List to hold region IDs.
        
        for point in main_boundary:                                                                                                     # Iterate over boundary points.
            all_points.append(point)                                                                                                    # Append boundary point.
            all_regions.append(1)                                                                                                       # Assign region 1 ID.
        
        for point in main_interior:                                                                                                     # Iterate over interior points.
            all_points.append(point)                                                                                                    # Append interior point.
            all_regions.append(1)                                                                                                       # Assign region 1 ID.
        
        # 6. Interior regions generation
        if inside_regions and len(regions) > 1:                                                                                         # Check if interior regions should be processed.

            interior_clouds = interior_regions_strategy(regions, actual_cloud_size)                                                     # Generate inner clouds.

            for boundary_points, interior_points, region_id in interior_clouds:                                                         # Iterate over generated sub-regions.
                for point in boundary_points:                                                                                           # Iterate over sub-region boundaries.
                    all_points.append(point)                                                                                            # Append boundary point.
                    all_regions.append(region_id)                                                                                       # Assign current region ID.
                for point in interior_points:                                                                                           # Iterate over sub-region interiors.
                    all_points.append(point)                                                                                            # Append interior point.
                    all_regions.append(region_id)                                                                                       # Assign current region ID.

        # 7. Classification and export
        all_points_arr  = np.array(all_points)                                                                                          # Convert total list to NumPy array.
        classifications = classify_nodes(all_points_arr, all_regions, regions, actual_cloud_size, inside_regions=inside_regions)        # Classify generated points.
        success         = export_to_csv(all_points_arr, classifications, all_regions, output_file)                                      # Export the generated points.

        if not success:                                                                                                                 # Check if export failed.
            return {"success": False, "error": "Failed to export CSV"}                                                                  # Return failure status.
        
        output_base = os.path.splitext(output_file)[0]                                                                                  # Get base filename without extension.
        create_visualization(all_points_arr, all_regions, output_base, classifications, save=save, show=show)                           # Generate visual representation.
        
        # 8. Result packaging
        results = {                                                                                                                     # Build the success dictionary.
            "success": True,                                                                                                            # Set success flag.
            "output_file": output_file,                                                                                                 # Include output filename.
            "visualization_file": f"{output_base}.png",                                                                                 # Include PNG file path.
            "visualization_svg_file": f"{output_base}.svg",                                                                             # Include SVG file path.
            "total_nodes": len(all_points),                                                                                             # Include total number of nodes.
            "regions_generated": len(set(all_regions)),                                                                                 # Include total unique regions.
            "main_region_nodes": sum(1 for r in all_regions if r == 1),                                                                 # Include number of nodes in region 1.
            "interior_regions_generated": inside_regions and len(regions) > 1,                                                          # Include flag for interior generation.
            "cloud_size": actual_cloud_size                                                                                             # Include final cloud size used.
        }                                                                                                                               # Execute statement.
        
        if temp_reduced_file and os.path.exists(temp_reduced_file):                                                                     # Clean up temporary files if any.
            os.remove(temp_reduced_file)                                                                                                # Remove the temporary file.
            
        return results                                                                                                                  # Return the results dictionary.
        
    # 9. Exception handling
    except Exception as e:                                                                                                              # Catch any unexpected errors.
        error_msg = f"Error in cloud generation with {method_name} Distribution: {str(e)}"                                              # Format error message.
        
        return {"success": False, "error": error_msg}                                                                                   # Return failure dictionary.

def generate_cloud_natural(csv_file: str, output_file: str, inside_regions: bool = False,                                               # Natural point cloud generator.
                           cloud_size: Optional[float] = None, density_multiplier: float = 1.0,                                         # Spacing and density parameters.
                           save: bool = True, show: bool = False) -> Dict:                                                              # Generate cloud via Poisson disk.
    """
    generate_cloud_natural
    Generate point cloud using Natural Distribution algorithm with Poisson Disk Sampling.
    
    This function creates a natural, organic point cloud using advanced Poisson Disk
    Sampling algorithm. It processes CSV boundary data to generate optimized point cloud
    suitable for meshless Generalized Finite Differences (mGFD) method computations.
    
    Input:
        csv_file            str         Path to input CSV file containing boundary coordinates.
        output_file         str         Base path for output files (without extension).
        inside_regions      bool        Enable multi-region processing for interior holes.
        cloud_size          float       Override automatic cloud size calculation.
        density_multiplier  float       Multiplier for the point cloud density. Default is 1.0.
        save                bool        If True, save PNG/SVG files to disk.
        show                bool        If True, display interactive plot window.
    
    Output:
        results             Dict        Comprehensive results dictionary.
    """
    try:                                                                                                                                # Start generation block.
        # 1. Strategy definition
        main_strategy     = generate_region_cloud_with_holes_poisson                                                                    # Assign natural main strategy.
        interior_strategy = generate_interior_regions_clouds_poisson                                                                    # Assign natural interior strategy.
        
        return _generate_cloud_core(                                                                                                    # Execute core logic.
            csv_file, output_file, "Natural",                                                                                           # Set file paths and method name.
            main_strategy, interior_strategy,                                                                                           # Set appropriate generation strategies.
            inside_regions, cloud_size, density_multiplier, save, show                                                                  # Forward other parameters.
        )                                                                                                                               # Execute statement.
        
    # 2. Exception handling
    except Exception as e:                                                                                                              # Catch exceptions for natural generation.
        error_msg = f"Error in cloud generation with Natural Distribution: {str(e)}"                                                    # Format error message.
        
        return {"success": False, "error": error_msg}                                                                                   # Return failure dictionary.

def generate_cloud_regular(csv_file: str, output_file: str, inside_regions: bool = False,                                               # Regular point cloud generator.
                           cloud_size: Optional[float] = None, density_multiplier: float = 1.0,                                         # Spacing and density parameters.
                           save: bool = True, show: bool = False) -> Dict:                                                              # Generate cloud via regular grid.
    """
    generate_cloud_regular
    Generate point cloud using Regular Distribution algorithm with grid-based approach.
    
    This function creates a uniform, structured point cloud using grid-based algorithms.
    It processes CSV boundary data to generate optimized point cloud with consistent spacing
    suitable for meshless Generalized Finite Differences (mGFD) method computations.
    
    Input:
        csv_file            str         Path to input CSV file containing boundary coordinates.
        output_file         str         Base path for output files (without extension).
        inside_regions      bool        Enable multi-region processing for interior holes.
        cloud_size          float       Override automatic cloud size calculation.
        density_multiplier  float       Multiplier for the point cloud density. Default is 1.0.
        save                bool        If True, save PNG/SVG files to disk.
        show                bool        If True, display interactive plot window.
    
    Output:
        results             Dict        Comprehensive results dictionary.
    """
    try:                                                                                                                                # Start generation block.
        # 1. Strategy definition
        main_strategy     = generate_region_cloud_with_holes                                                                            # Assign regular main strategy.
        interior_strategy = generate_interior_regions_clouds                                                                            # Assign regular interior strategy.
        
        return _generate_cloud_core(                                                                                                    # Execute core logic.
            csv_file, output_file, "Regular",                                                                                           # Set file paths and method name.
            main_strategy, interior_strategy,                                                                                           # Set appropriate generation strategies.
            inside_regions, cloud_size, density_multiplier, save, show                                                                  # Forward other parameters.
        )                                                                                                                               # Execute statement.
        
    # 2. Exception handling
    except Exception as e:                                                                                                              # Catch exceptions for regular generation.
        error_msg = f"Error in cloud generation with Regular Distribution: {str(e)}"                                                    # Format error message.

        return {"success": False, "error": error_msg}                                                                                   # Return failure dictionary.
