"""
CloudGenerator.visualization — Visualization

Overview:
    This module handles the generation of high-quality visual representations of the
    generated point clouds. It uses Matplotlib to create static images (PNG) and
    vector graphics (SVG) that clearly distinguish between regions and node types.

Data conventions:
    None.

Public API:
    create_visualization    Create a visualization of the generated point cloud.
    render_neighbors_graph  Render the connectivity graph of neighbors.

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
import numpy as np                                                                                                                      # Numerical arrays and mathematical operations.
import matplotlib.pyplot as plt                                                                                                         # Interface for pyplot.

from typing import List, Optional                                                                                                       # Import type hints.
from matplotlib.collections import LineCollection                                                                                       # For fast line rendering.

def create_visualization(points: np.ndarray, regions_list: List[int], output_base: str, classifications: Optional[List[str]] = None, save: bool = True) -> bool:
    """
    create_visualization
    Create a visualization of the generated point cloud.
    
    Generates high-quality PNG and SVG images with differentiated colors
    for boundary and interior nodes across multiple regions.
    
    Input:
        points              m x 2 ndarray           Array of (x, y) coordinates for all points.
        regions_list        List[int]               Region assignments for each node.
        output_base         str                     Base path for output files (without extension).
        classifications     Optional[List[str]]     Node classifications ('boundary' or 'interior'). Default is None.
        save                bool                    If True, save to disk. If False, show interactive plot.
    
    Output:
        success             bool                    True if visualization was successful, False otherwise.
    """
    try:                                                                                                                                # Start plotting process.
        # 1. Define colors (RGB tuples scaled to 0-1 for matplotlib)
        interior_colors = [                                                                                                             # Colors for internal nodes.
            (51/255, 153/255, 255/255),                                                                                                 # Light Blue
            (77/255, 230/255, 77/255),                                                                                                  # Bright Green
            (255/255, 179/255, 51/255),                                                                                                 # Bright Orange
            (204/255, 77/255, 255/255),                                                                                                 # Bright Purple
            (230/255, 153/255, 51/255),                                                                                                 # Golden Brown
            (255/255, 128/255, 204/255),                                                                                                # Bright Pink
            (179/255, 179/255, 179/255)                                                                                                 # Light Gray
        ]                                                                                                                               # Execute statement.
        
        boundary_colors = [                                                                                                             # Colors for boundary nodes.
            (204/255, 0/255, 0/255),                                                                                                    # Dark Red
            (153/255, 0/255, 153/255),                                                                                                  # Dark Purple
            (0/255, 77/255, 204/255),                                                                                                   # Dark Blue
            (204/255, 153/255, 0/255),                                                                                                  # Dark Yellow
            (0/255, 102/255, 0/255),                                                                                                    # Dark Green
            (0/255, 153/255, 102/255),                                                                                                  # Dark Teal
            (51/255, 51/255, 51/255)                                                                                                    # Dark Gray
        ]                                                                                                                               # Execute statement.
        
        if len(points) == 0:                                                                                                            # Verify point data exists.
            return False                                                                                                                # Return False if empty.
            
        # 2. Convert lists to numpy arrays for efficient indexing
        regions_arr = np.array(regions_list)                                                                                            # Convert region mapping.
        if classifications:                                                                                                             # If classifications are provided.
            classifications_arr = np.array(classifications)                                                                             # Convert to array.
        else:                                                                                                                           # If not provided.
            classifications_arr = np.array(['interior'] * len(points))                                                                  # Default all to interior.
            
        # 3. Setup plot
        plt.figure(figsize=(10, 8))                                                                                                     # Create figure with standard screen DPI.
        ax = plt.gca()                                                                                                                  # Get current axes object.
        ax.set_aspect('equal')                                                                                                          # Force 1:1 aspect ratio.
        
        # Configure axes
        ax.grid(True, linestyle=':', alpha=0.6)                                                                                         # Add soft dotted grid.
        ax.set_xlabel('X Coordinate')                                                                                                   # Set X axis label.
        ax.set_ylabel('Y Coordinate')                                                                                                   # Set Y axis label.
        
        # 4. Get unique regions
        unique_regions = sorted(list(set(regions_list)))                                                                                # Find all unique region IDs.
        
        # 5. Plot points by region
        for region_id in unique_regions:                                                                                                # Iterate through each region.
            region_mask = (regions_arr == region_id)                                                                                    # Region mask.
            
            boundary_mask = region_mask & (classifications_arr == 'boundary')                                                           # Boundary points.
            if np.any(boundary_mask):                                                                                                   # If any boundary points exist here.
                boundary_pts = points[boundary_mask]                                                                                    # Extract boundary coordinates.
                color_idx    = (region_id - 1) % len(boundary_colors)                                                                   # Select boundary color for region.
                ax.scatter(boundary_pts[:, 0], boundary_pts[:, 1],                                                                      # Plot boundary scatter.
                          c=[boundary_colors[color_idx]], s=8, marker='.',                                                              # Configure color and size.
                          edgecolors='none', label=f'{region_id} Boundary', zorder=10)                                                  # Set legend label and stacking order.
            
            interior_mask = region_mask & (classifications_arr != 'boundary')                                                           # Interior points.
            if np.any(interior_mask):                                                                                                   # If any interior points exist here.
                interior_pts = points[interior_mask]                                                                                    # Extract interior coordinates.
                color_idx    = (region_id - 1) % len(interior_colors)                                                                   # Select interior color for region.
                ax.scatter(interior_pts[:, 0], interior_pts[:, 1],                                                                      # Plot interior scatter.
                          c=[interior_colors[color_idx]], s=5, marker='.',                                                              # Configure color and size.
                          alpha=0.8, edgecolors='none', label=f'{region_id} Interior', zorder=5)                                        # Add transparency and label.

        # 6. Add title and adjust layout
        plt.title(f'Generated Point Cloud\n{len(points)} Total Nodes', fontsize=14)                                                     # Set main title with node count.
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small', markerscale=1.5)                     # Place legend outside plot.
        
        # 7. Save output files
        if save:                                                                                                                        # If saving to disk is requested.
            png_file = f"{output_base}.png"                                                                                             # Save PNG.
            plt.savefig(png_file, format='png', bbox_inches='tight', dpi=300)                                                           # Write PNG to disk.
            
            svg_file = f"{output_base}.svg"                                                                                             # Save SVG.
            plt.savefig(svg_file, format='svg', bbox_inches='tight')                                                                    # Write SVG to disk.
        else:                                                                                                                           # If interactive plotting is requested.
            plt.show()                                                                                                                  # Display the interactive plot window.
            
        plt.close()                                                                                                                     # Close figure to prevent memory leak.
        
        return True                                                                                                                     # Return success state.
        
    except Exception as e:                                                                                                              # If error occurs during plotting.
        try:                                                                                                                            # Attempt cleanup.
            plt.close()                                                                                                                 # Close plot in case of error to free memory.
        except:                                                                                                                         # If close fails.
            pass                                                                                                                        # Ignore.
        return False                                                                                                                    # Return failure state.

def render_neighbors_graph(points: np.ndarray, neighbors_indices: np.ndarray, regions: List[int], output_base: str) -> bool:
    """
    render_neighbors_graph
    Render the connectivity graph of neighbors using LineCollection for performance.
    
    Input:
        points              m x 2 ndarray       Array of (x, y) coordinates for all nodes.
        neighbors_indices   m x nvec ndarray    NxK array of neighbor indices. -1 indicates no neighbor.
        regions             List[int]           Array of region identifiers for coloring.
        output_base         str                 Base path for output files (without extension).
        
    Output:
        success             bool                True if successful, False otherwise.
    """
    try:                                                                                                                                # Start graph plotting process.
        points_arr            = np.array(points)                                                                                        # Ensure points are a NumPy array.
        regions_arr           = np.array(regions)                                                                                       # Ensure regions are a NumPy array.
        neighbors_indices_arr = np.array(neighbors_indices)                                                                             # Ensure neighbor indices are an array.
        
        # 1. Setup plot
        plt.figure(figsize=(10, 8))                                                                                                     # Create figure with standard screen DPI.
        ax = plt.gca()                                                                                                                  # Get the active axes.
        ax.set_aspect('equal')                                                                                                          # Set equal scaling.
        
        ax.grid(True, linestyle=':', alpha=0.6)                                                                                         # Configure axes.
        ax.set_xlabel('X Coordinate')                                                                                                   # Set X axis label.
        ax.set_ylabel('Y Coordinate')                                                                                                   # Set Y axis label.
        
        # 2. Define base colors (simplified for lines)
        region_colors = [                                                                                                               # Define palette for region links.
            (51/255, 153/255, 255/255),                                                                                                 # Light Blue
            (77/255, 230/255, 77/255),                                                                                                  # Bright Green
            (255/255, 179/255, 51/255),                                                                                                 # Bright Orange
            (204/255, 77/255, 255/255),                                                                                                 # Bright Purple
            (230/255, 153/255, 51/255),                                                                                                 # Golden Brown
            (255/255, 128/255, 204/255),                                                                                                # Bright Pink
            (179/255, 179/255, 179/255)                                                                                                 # Light Gray
        ]                                                                                                                               # Execute statement.
        
        unique_regions = sorted(list(set(regions)))                                                                                     # Find all unique region IDs.
        
        # 3. Build line segments and scatter points per region
        for region_id in unique_regions:                                                                                                # Iterate through regions.
            region_mask       = (regions_arr == region_id)                                                                              # Create logical mask for the region.
            indices_in_region = np.where(region_mask)[0]                                                                                # Get point indices for this region.
            
            if len(indices_in_region) == 0:                                                                                             # Skip if region is empty.
                continue                                                                                                                # Move to next region.
                
            color_idx  = (region_id - 1) % len(region_colors)                                                                           # Get corresponding color index.
            base_color = region_colors[color_idx]                                                                                       # Assign base color.
            
            segments = []                                                                                                               # List for connection line segments.
            
            for i in indices_in_region:                                                                                                 # Iterate over nodes in region.
                p1 = points[i]                                                                                                          # Get starting node coordinates.
                for j in neighbors_indices[i]:                                                                                          # Iterate over node's neighbors.
                    if j != -1:                                                                                                         # Check if valid neighbor exists.
                        p2 = points[j]                                                                                                  # Get ending node coordinates.
                        segments.append([p1, p2])                                                                                       # Add segment to the list.
            
            if segments:                                                                                                                # If connections exist.
                lc = LineCollection(segments, colors=[base_color], linewidths=0.5, alpha=0.6)                                           # Create collection of lines.
                ax.add_collection(lc)                                                                                                   # Render the lines.
            
            # Plot the nodes on top
            ax.scatter(points[indices_in_region, 0], points[indices_in_region, 1],                                                      # Overlay scatter points.
                       c=[base_color], s=8, marker='.', edgecolors='none', label=f'Region {region_id}', zorder=5)                       # Plot with proper color and legend.
                       
        plt.title(f'Connectivity Graph\n{len(points)} Total Nodes', fontsize=14)                                                        # Set plot title.
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small')                                      # Set legend properties.
        
        # 4. Save output files
        png_file = f"{output_base}.png"                                                                                                 # Save PNG.
        plt.savefig(png_file, format='png', bbox_inches='tight', dpi=300)                                                               # Write PNG image.
        
        svg_file = f"{output_base}.svg"                                                                                                 # Save SVG.
        plt.savefig(svg_file, format='svg', bbox_inches='tight')                                                                        # Write SVG image.
        
        plt.close()                                                                                                                     # Release figure resources.
        
        return True                                                                                                                     # Return success state.
        
    except Exception as e:                                                                                                              # If error occurs during rendering.
        try:                                                                                                                            # Attempt cleanup.
            plt.close()                                                                                                                 # Close the plot.
        except:                                                                                                                         # If close fails.
            pass                                                                                                                        # Ignore.
        return False                                                                                                                    # Return failure state.
