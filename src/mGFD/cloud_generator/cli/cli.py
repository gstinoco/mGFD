"""
CloudGenerator.cli — Command Line Interface

Overview:
    This module provides the command line interface (CLI) for the CloudGenerator suite.
    It exposes subcommands for generating new point clouds from contours and reducing
    the density of existing point clouds.

Data conventions:
    None.

Public API:
    handle_generate     Handler for the 'generate' subcommand.
    handle_reduce       Handler for the 'reduce' subcommand.
    main                Main entry point for the CLI.

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
import sys                                                                                                                              # System-specific parameters and functions.
import logging                                                                                                                          # Standard logging module.
import argparse                                                                                                                         # Module for parsing command-line arguments.
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Any                                                                                                                  # Type hinting.

from mGFD.cloud_generator.core.reduction import reduce_points_by_region                                                                 # Import reduction function.
from mGFD.cloud_generator.viz.visualization import create_visualization                                                                 # Import visualization creator.
from mGFD.cloud_generator.core.generator import generate_cloud_natural, generate_cloud_regular                                          # Import generation functions.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def handle_generate(args: argparse.Namespace) -> None:
    """
    handle_generate
    Handler for the 'generate' subcommand.
    
    Parses the arguments and delegates execution to the appropriate point cloud
    generation algorithm (natural or regular).
    
    Input:
        args        argparse.Namespace  Parsed command line arguments containing input, output, method,
                                        inside_regions, and density.
    
    Output:
        None
    """
    # 1. Validation
    if not os.path.exists(args.input):                                                                                                  # Check if the input file exists.
        logger.error(f"Error: Input file '{args.input}' not found.")                                                                    # Log an error message.
        sys.exit(1)                                                                                                                     # Exit with an error code.
        
    # 2. Header and info output
    if not args.quiet:
        logger.info("mGFD CloudGenerator - Generate")                                                                                   # Log header text for generation.
        logger.info("==============================")                                                                                   # Log separator line.
        logger.info(f"Input:    {args.input}")                                                                                          # Log the input file path.
        logger.info(f"Output:   {args.output}")                                                                                         # Log the output file path.
        logger.info(f"Method:   {args.method}")                                                                                         # Log the selected method.
        logger.info(f"Density:  {args.density}x")                                                                                       # Log the density multiplier.
        logger.info("Generating point cloud, please wait...")                                                                           # Log wait message.
        
    # 3. Cloud generation
    try:                                                                                                                                # Start generation process inside a try-except block.
        if args.method == "natural":                                                                                                    # Check if the selected method is 'natural'.
            result = generate_cloud_natural(                                                                                            # Call the natural generation function.
                csv_file=args.input,                                                                                                    # Pass the input CSV file.
                output_file=args.output,                                                                                                # Pass the output CSV file.
                inside_regions=args.inside_regions,                                                                                     # Pass the inside regions flag.
                density_multiplier=args.density                                                                                         # Pass the density multiplier.
            )
        else:                                                                                                                           # If the method is 'regular' (or anything else).
            result = generate_cloud_regular(                                                                                            # Call the regular generation function.
                csv_file=args.input,                                                                                                    # Pass the input CSV file.
                output_file=args.output,                                                                                                # Pass the output CSV file.
                inside_regions=args.inside_regions,                                                                                     # Pass the inside regions flag.
                density_multiplier=args.density                                                                                         # Pass the density multiplier.
            )
            
        # 4. Success output
        if not args.quiet:
            logger.info("\nSuccess!")                                                                                                   # Log success message.
            logger.info(f"Total points generated: {result.get('total_nodes', 'Unknown')}")                                              # Log total generated nodes.
            logger.info(f"Regions generated:      {result.get('regions_generated', 'Unknown')}")                                        # Log number of generated regions.
            logger.info(f"Main region points:     {result.get('main_region_nodes', 'Unknown')}")                                        # Log points in the main region.
            logger.info(f"Visualizations saved to {result.get('visualization_file', '')} and {result.get('visualization_svg_file', '')}") # Log visualization paths.
        
    # 5. Exception handling
    except Exception as e:                                                                                                              # Handle general exceptions.
        logger.error(f"\nError generating cloud: {e}")                                                                                  # Log error details.
        sys.exit(1)                                                                                                                     # Exit with an error code.

def handle_reduce(args: argparse.Namespace) -> None:
    """
    handle_reduce
    Handler for the 'reduce' subcommand.
    
    Parses the arguments and delegates execution to the point cloud reduction
    algorithm to decimate an existing point cloud while preserving its boundary.
    
    Input:
        args        argparse.Namespace  Parsed command line arguments containing input, output, and multiplier.
    
    Output:
        None
    """
    # 1. Validation
    if not os.path.exists(args.input):                                                                                                  # Check if the input file exists.
        logger.error(f"Error: Input file '{args.input}' not found.")                                                                    # Log an error message.
        sys.exit(1)                                                                                                                     # Exit with an error code.
        
    # 2. Header and info output
    if not args.quiet:
        logger.info("mGFD CloudGenerator - Reduce")                                                                                     # Log header text for reduction.
        logger.info("============================")                                                                                     # Log separator line.
        logger.info(f"Input:      {args.input}")                                                                                        # Log the input file path.
        logger.info(f"Output:     {args.output}")                                                                                       # Log the output file path.
        logger.info(f"Multiplier: {args.multiplier}")                                                                                   # Log the reduction multiplier.
        logger.info("Reducing point cloud, please wait...")                                                                             # Log wait message.
        
    # 3. Cloud reduction
    try:                                                                                                                                # Start reduction process inside a try-except block.
        result = reduce_points_by_region(                                                                                               # Call the point reduction function.
            input_csv=args.input,                                                                                                       # Pass the input CSV file.
            output_csv=args.output,                                                                                                     # Pass the output CSV file.
            multiplier=args.multiplier                                                                                                  # Pass the reduction multiplier.
        )
        
        # 4. Visualization and output
        if result is not None:                                                                                                          # Check if the reduction was successful.
            if not args.quiet:
                logger.info("\nSuccess!")                                                                                               # Log success message.
                logger.info(f"Total points after reduction: {len(result)}")                                                             # Log the new total number of points.
            
            points          = np.array([(row['x'], row['y']) for row in result])                                                        # Extract point coordinates.
            regions_list    = [row.get('region', 1) for row in result]                                                                  # Extract region list.
            classifications = [row.get('classification', 'interior') for row in result]                                                 # Extract classifications.
            output_base     = os.path.splitext(args.output)[0]                                                                          # Get output path without extension.
            
            create_visualization(points, regions_list, output_base, classifications)                                                    # Generate visual plots.
            
            if not args.quiet:
                logger.info(f"Visualizations saved to {output_base}.png and {output_base}.svg")                                         # Log saved visual plots.
        else:                                                                                                                           # If reduction returned None (failed).
            logger.error("\nError: Failed to reduce cloud. Check logs for details.")                                                    # Log error message.
            sys.exit(1)                                                                                                                 # Exit with an error code.
            
    # 5. Exception handling
    except Exception as e:                                                                                                              # Handle general exceptions.
        logger.error(f"\nError reducing cloud: {e}")                                                                                    # Log error details.
        sys.exit(1)                                                                                                                     # Exit with an error code.

def main() -> None:
    """
    main
    Main entry point for the CLI.
    
    Configures the argument parser, defines the 'generate' and 'reduce'
    subcommands, and dispatches the execution to the corresponding handler.
    
    Input:
        None
        
    Output:
        None
    """
    # 1. Global Logging configuration
    logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                       # Configure logging format for CLI.

    # 2. Parser initialization
    parser = argparse.ArgumentParser(                                                                                                   # Instantiate an ArgumentParser.
        description="mGFD CloudGenerator - Point cloud generation and optimization suite."                                              # Define program description.
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress non-error console output.")                                # Add quiet flag globally.
    
    subparsers          = parser.add_subparsers(dest="command", help="Available commands")                                              # Add sub-command parsers.
    subparsers.required = True                                                                                                          # Ensure a sub-command is always provided.
    
    # 2. Generate command configuration
    parser_gen = subparsers.add_parser("generate", help="Generate a point cloud from boundary contours.")                               # Create parser for 'generate'.
    parser_gen.add_argument(                                                                                                            # Add an argument to 'generate'.
        "-i", "--input",                                                                                                                # Flag names for the input file.
        required=True,                                                                                                                  # Make this argument mandatory.
        help="Input CSV file containing the boundary contours."                                                                         # Help text for the input argument.
    )
    parser_gen.add_argument(                                                                                                            # Add an argument to 'generate'.
        "-o", "--output",                                                                                                               # Flag names for the output file.
        required=True,                                                                                                                  # Make this argument mandatory.
        help="Output CSV file where the point cloud will be saved."                                                                     # Help text for the output argument.
    )
    parser_gen.add_argument(                                                                                                            # Add an argument to 'generate'.
        "-m", "--method",                                                                                                               # Flag names for the method.
        choices=["natural", "regular"],                                                                                                 # Allow only 'natural' or 'regular'.
        default="natural",                                                                                                              # Set default method to 'natural'.
        help="Point generation method: 'natural' (Poisson-Disk) or 'regular' (Grid-based)."                                             # Help text for the method argument.
    )
    parser_gen.add_argument(                                                                                                            # Add an argument to 'generate'.
        "--inside-regions",                                                                                                             # Flag name for generating inside regions.
        action="store_true",                                                                                                            # Store a boolean True if provided.
        help="If specified, generate points inside holes (islands) too."                                                                # Help text for inside regions flag.
    )
    parser_gen.add_argument(                                                                                                            # Add an argument to 'generate'.
        "-d", "--density",                                                                                                              # Flag names for the density.
        type=float,                                                                                                                     # Require a float value.
        default=1.0,                                                                                                                    # Set default density to 1.0.
        help="Density multiplier. Higher value = denser cloud."                                                                         # Help text for the density argument.
    )
    parser_gen.set_defaults(func=handle_generate)                                                                                       # Bind the handler function to 'generate'.
    
    # 3. Reduce command configuration
    parser_red = subparsers.add_parser("reduce", help="Reduce an existing point cloud density.")                                        # Create parser for 'reduce'.
    parser_red.add_argument(                                                                                                            # Add an argument to 'reduce'.
        "-i", "--input",                                                                                                                # Flag names for the input file.
        required=True,                                                                                                                  # Make this argument mandatory.
        help="Input CSV file containing an existing point cloud."                                                                       # Help text for the input argument.
    )
    parser_red.add_argument(                                                                                                            # Add an argument to 'reduce'.
        "-o", "--output",                                                                                                               # Flag names for the output file.
        required=True,                                                                                                                  # Make this argument mandatory.
        help="Output CSV file where the reduced point cloud will be saved."                                                             # Help text for the output argument.
    )
    parser_red.add_argument(                                                                                                            # Add an argument to 'reduce'.
        "--multiplier",                                                                                                                 # Flag name for the reduction multiplier.
        type=int,                                                                                                                       # Require an integer value.
        default=2,                                                                                                                      # Set default multiplier to 2.
        help="Reduction multiplier. 1 = ~50% reduction. 2 = ~75% reduction."                                                            # Help text for multiplier.
    )
    parser_red.set_defaults(func=handle_reduce)                                                                                         # Bind the handler function to 'reduce'.
    
    # 4. Argument parsing and execution
    args = parser.parse_args()                                                                                                          # Parse the command-line arguments.
    args.func(args)                                                                                                                     # Call the assigned function based on parsed arguments.
    
if __name__ == "__main__":                                                                                                              # Check if the script is run directly.
    main()                                                                                                                              # Execute the main function.
