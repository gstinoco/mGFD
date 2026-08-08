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

import os                                                                           # Operating system interfaces.
import logging                                                                      # Module for event logging.
import tempfile                                                                     # Temporary files creation.
import numpy as np                                                                  # Import numerical processing library.
import time                                                                         # Import timing functions.
from mGFD.cloud_generator.core import reduction as reduce_points                     # Import the reduction module.
from mGFD.cloud_generator.io.data_processing import load_regions                     # Import load_regions function.
from mGFD.cloud_generator.core.point_generation import (                             # Import specific generation functions.
    generate_region_cloud_with_holes,                                               # Regular grid region generation.
    generate_region_cloud_with_holes_poisson,                                       # Poisson region generation.
    generate_interior_regions_clouds,                                               # Interior regions regular generation.
    generate_interior_regions_clouds_poisson,                                       # Interior regions Poisson generation.
)
from mGFD.cloud_generator.core.classification import classify_nodes                  # Import classification function.
from mGFD.cloud_generator.viz.visualization import create_visualization              # Import visualization function.
from mGFD.cloud_generator.io.export import export_to_csv                             # Import the CSV export function.

def _generate_cloud_core(csv_file: str, output_file: str, method_name: str, main_region_strategy: callable, interior_regions_strategy: callable, 
                         inside_regions: bool = False, cloud_size: float = None, density_multiplier: float = 1.0) -> dict:
    """
    Core function for cloud generation to strictly follow DRY principle.
    Encapsulates common logic for Natural and Regular distribution methods.
    """
    temp_reduced_file = None                                                        # Initialize temporary file variable.
    try:                                                                            # Start main try-except block.
        logging.info(f"=== Starting Cloud Generation Process with {method_name} Distribution ===") # Log start of process.
        logging.info(f"Input file: {csv_file}")                                     # Log input filename.
        logging.info(f"Output file: {output_file}")                                 # Log output filename.
        logging.info(f"Interior regions: {inside_regions}")                         # Log if interior regions are used.
        
        working_csv_file = csv_file                                                 # Set the working CSV file path.
        
        regions = load_regions(working_csv_file)                                    # Parse and load regions from CSV.
        if not regions:                                                             # Check if regions were loaded properly.
            return {"success": False, "error": "Failed to load regions from CSV"}   # Return failure if no regions.
        
        logging.info(f"Loaded {len(regions)} region(s)")                            # Log the number of loaded regions.
        
        main_region = regions[0]                                                    # Define the main (outer) region.
        hole_regions = regions[1:] if len(regions) > 1 else []                      # Define hole (inner) regions.
        
        # Calculate explicit cloud size with density multiplier
        from mGFD.cloud_generator.utils.utils import calculate_cloud_size                                     # Import cloud size calculator.
        if cloud_size is None:                                                      # If cloud_size is not provided.
            base_cloud_size = calculate_cloud_size(main_region)                     # Calculate the base cloud size.
            # Higher density multiplier = smaller cloud size = more points
            cloud_size = base_cloud_size / density_multiplier                       # Adjust cloud size with multiplier.
            logging.info(f"Calculated base cloud size: {base_cloud_size}, final with multiplier: {cloud_size}") # Log the sizes.
            
        # Execute Main Region Strategy
        # Pass inside_regions to control hole boundary generation
        main_boundary, main_interior, actual_cloud_size = main_region_strategy(main_region, hole_regions, cloud_size, inside_regions) # Generate main region points.
        
        if main_boundary is None or main_interior is None:                          # Ensure points were generated successfully.
            return {"success": False, "error": "Failed to generate cloud for main region"} # Return failure if not.
        
        logging.info(f"Generated main region cloud with {method_name} Distribution: {len(main_boundary)} boundary + {len(main_interior)} interior points") # Log generated points.
        
        all_points = []                                                             # List to hold all points across regions.
        all_regions = []                                                            # List to hold region IDs.
        
        for point in main_boundary:                                                 # Iterate over boundary points.
            all_points.append(point)                                                # Append boundary point.
            all_regions.append(1)                                                   # Assign region 1 ID.
        
        for point in main_interior:                                                 # Iterate over interior points.
            all_points.append(point)                                                # Append interior point.
            all_regions.append(1)                                                   # Assign region 1 ID.
        
        # Execute Interior Regions Strategy
        if inside_regions and len(regions) > 1:                                     # Check if interior regions should be processed.
            logging.info(f"Generating interior regions with {method_name} Distribution...") # Log generation start.
            interior_clouds = interior_regions_strategy(regions, actual_cloud_size) # Generate inner clouds.
            
            for boundary_points, interior_points, region_id in interior_clouds:     # Iterate over generated sub-regions.
                for point in boundary_points:                                       # Iterate over sub-region boundaries.
                    all_points.append(point)                                        # Append boundary point.
                    all_regions.append(region_id)                                   # Assign current region ID.
                for point in interior_points:                                       # Iterate over sub-region interiors.
                    all_points.append(point)                                        # Append interior point.
                    all_regions.append(region_id)                                   # Assign current region ID.
                logging.info(f"Added region {region_id}: {len(boundary_points)} boundary + {len(interior_points)} interior points") # Log region additions.
        
        all_points = np.array(all_points)                                           # Convert total list to NumPy array.
        
        # Note: classify_nodes now imports from classification
        # We need to pass original_regions_contours (regions) to it for accurate classification
        classifications = classify_nodes(all_points, all_regions, regions, actual_cloud_size, inside_regions=inside_regions) # Classify generated points.
        
        success = export_to_csv(all_points, classifications, all_regions, output_file) # Export the generated points.
        if not success:                                                             # Check if export failed.
            return {"success": False, "error": "Failed to export CSV"}              # Return failure status.
        
        output_base = os.path.splitext(output_file)[0]                              # Get base filename without extension.
        create_visualization(all_points, all_regions, output_base, classifications) # Generate visual representation.
        
        results = {                                                                 # Build the success dictionary.
            "success": True,                                                        # Set success flag.
            "output_file": output_file,                                             # Include output filename.
            "visualization_file": f"{output_base}.png",                             # Include PNG file path.
            "visualization_svg_file": f"{output_base}.svg",                         # Include SVG file path.
            "total_nodes": len(all_points),                                         # Include total number of nodes.
            "regions_generated": len(set(all_regions)),                             # Include total unique regions.
            "main_region_nodes": sum(1 for r in all_regions if r == 1),             # Include number of nodes in region 1.
            "interior_regions_generated": inside_regions and len(regions) > 1,      # Include flag for interior generation.
            "cloud_size": actual_cloud_size                                         # Include final cloud size used.
        }
        
        logging.info(f"=== Cloud Generation with {method_name} Distribution Completed Successfully ===") # Log completion.
        
        if temp_reduced_file and os.path.exists(temp_reduced_file):                 # Clean up temporary files if any.
            os.remove(temp_reduced_file)                                            # Remove the temporary file.
            
        return results                                                              # Return the results dictionary.
        
    except Exception as e:                                                          # Catch any unexpected errors.
        error_msg = f"Error in cloud generation with {method_name} Distribution: {str(e)}" # Format error message.
        logging.error(error_msg)                                                    # Log the error.
        return {"success": False, "error": error_msg}                               # Return failure dictionary.

def generate_cloud_natural(csv_file: str, output_file: str, inside_regions: bool = False,
                           cloud_size: float = None, density_multiplier: float = 1.0) -> dict:
    """
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
    
    Output:
        results             dict        Comprehensive results dictionary.
    """
    try:                                                                            # Start generation block.
        # Define strategies
        main_strategy = generate_region_cloud_with_holes_poisson                    # Assign natural main strategy.
        interior_strategy = generate_interior_regions_clouds_poisson                # Assign natural interior strategy.
        
        return _generate_cloud_core(                                                # Execute core logic.
            csv_file, output_file, "Natural",                                       # Set file paths and method name.
            main_strategy, interior_strategy,                                       # Set appropriate generation strategies.
            inside_regions, cloud_size, density_multiplier                          # Forward other parameters.
        )
        
    except Exception as e:                                                          # Catch exceptions for natural generation.
        error_msg = f"Error in cloud generation with Natural Distribution: {str(e)}"# Format error message.
        logging.error(error_msg)                                                    # Log the error.
        return {"success": False, "error": error_msg}                               # Return failure dictionary.

def generate_cloud_regular(csv_file: str, output_file: str, inside_regions: bool = False,
                           cloud_size: float = None, density_multiplier: float = 1.0) -> dict:
    """
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
    
    Output:
        results             dict        Comprehensive results dictionary.
    """
    try:                                                                            # Start generation block.
        # Define strategies
        main_strategy = generate_region_cloud_with_holes                            # Assign regular main strategy.
        interior_strategy = generate_interior_regions_clouds                        # Assign regular interior strategy.
        
        return _generate_cloud_core(                                                # Execute core logic.
            csv_file, output_file, "Regular",                                       # Set file paths and method name.
            main_strategy, interior_strategy,                                       # Set appropriate generation strategies.
            inside_regions, cloud_size, density_multiplier                          # Forward other parameters.
        )
        
    except Exception as e:                                                          # Catch exceptions for regular generation.
        error_msg = f"Error in cloud generation with Regular Distribution: {str(e)}"# Format error message.
        logging.error(error_msg)                                                    # Log the error.
        return {"success": False, "error": error_msg}                               # Return failure dictionary.
