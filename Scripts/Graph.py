"""
Graph — Visualization utilities for point-cloud solutions

Overview:
    Plotting and animation helpers to visualize approximate and reference solutions computed on a
    2D point cloud. Stationary plots show a single snapshot; transient plots show multiple time
    levels and can optionally generate an animation file.

Notes:
    - These functions assume p is an (m, 2) array with columns [x, y].
    - When save=True, images are saved to disk; otherwise figures are displayed interactively.
    - The animation writer defaults to ffmpeg. If ffmpeg is not available, a pillow-based fallback
      writer is used (and the output may be saved as .gif depending on the filename).

Public API:
    Cloud_Stationary          Side-by-side 3D scatter: approximation vs. theoretical solution.
    Cloud_Stationary_1        3D scatter: approximation only.
    Cloud_Transient           Side-by-side 3D scatter over time (interactive stepping or animation).
    Cloud_Transient_Steps     Save side-by-side snapshots at ~3 time levels.
    Cloud_Transient_1         3D scatter over time for approximation only (interactive or animation).
    Cloud_Transient_Steps_1   Save approximation-only snapshots at ~3 time levels.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
        Aula CIMNE-Morelia. México
        SIIIA-MATH: Soluciones de Ingeniería. México

    Date:
        May, 2024.
    Last Modification:
        August, 2026.
"""
## Library importation.
import numpy as np                                                                                      # Array utilities (min/max, linspace, arange).
import matplotlib.pyplot as plt                                                                         # Plotting interface.
from matplotlib import cm                                                                               # Colormaps.
from matplotlib.animation import FuncAnimation                                                          # Animation helper.
import os
import pandas as pd
from scipy.spatial import Delaunay
from Scripts.Neighbors import find_distances
import matplotlib.gridspec as gridspec




def _get_valid_triangulation(p, nom):                                                                   # Internal triangulation helper.
    """
    Compute a valid triangulation for the point cloud p by filtering out Delaunay
    triangles whose baricenters lie outside the water domain (defined in *_contours.csv).
    It processes regions (islands/disconnected bodies) independently to prevent triangles
    from crossing between them, and uses a local cache file to optimize performance.
    """

    if not nom:                                                                                         # Check if a valid filename was provided.
        return None                                                                                     # Return None.

    # Parse dataset details from nom path
    # nom path is typically: .../Results/Equation/LakeName/Scale/Solution
    # So its parent is Scale, grand-parent is LakeName.
    parts = nom.split(os.sep)                                                                           # Split path to parse dataset details.
    dataset = None                                                                                      # Initialize dataset variable.
    scale = None                                                                                        # Initialize scale variable.
    workspace_root = None                                                                               # Initialize workspace root variable.

    if len(parts) >= 4:                                                                                 # Check if path has enough components.
        scale = parts[-2]                                                                               # Extract scale name.
        dataset = parts[-3]                                                                             # Extract dataset name.

        # Climb to find workspace root containing Data/ directory
        temp = nom                                                                                      # Initialize temp path for workspace search.
        for _ in range(10):                                                                             # Climb up directory tree.
            temp = os.path.dirname(temp)                                                                # Go up one directory level.
            if os.path.exists(os.path.join(temp, 'Data')):                                              # Check if Data folder exists.
                workspace_root = temp                                                                   # Set workspace root.
                break                                                                                   # Stop climbing directory tree.

    # Determine triangulation cache path
    cache_path = None                                                                                   # Initialize cache path.
    cloud_path = None                                                                                   # Initialize cloud path.
    if workspace_root and dataset and scale:                                                            # Check if path info was successfully parsed.
        cloud_file = f"{dataset}_cloud.csv"                                                             # Format cloud filename.
        cloud_path = os.path.join(workspace_root, 'Data', dataset, scale, cloud_file)                   # Construct cloud path.
        if not os.path.exists(cloud_path):                                                              # Check if cloud file exists.
            # Fallback search in lake scale folder
            lake_scale_dir = os.path.join(workspace_root, 'Data', dataset, scale)                       # Construct fallback directory.
            if os.path.exists(lake_scale_dir):                                                          # Check if fallback directory exists.
                for f in os.listdir(lake_scale_dir):                                                    # Iterate over fallback directory files.
                    if f.endswith('.csv') and not f.endswith('neighbors.csv') and not f.endswith('triangulation.csv'):  # Find a valid CSV point cloud file.
                        cloud_path = os.path.join(lake_scale_dir, f)                                    # Construct cloud path.
                        break                                                                           # Stop climbing directory tree.
        
        if cloud_path and os.path.exists(cloud_path):                                                   # Check if valid cloud path was found.
            cache_path = cloud_path.replace('.csv', '_triangulation.csv')                               # Set cache path based on cloud path.

    # Note: We do not fallback to Results directory to avoid polluting it with triangulation files.

    # 1. Try to load from cache
    if cache_path and os.path.exists(cache_path):                                                       # Check if cache file exists.
        try:                                                                                            # Attempt operation.
            triangles = np.loadtxt(cache_path, dtype=np.int32, delimiter=',')                           # Load triangulation from cache.
            if triangles.ndim == 1:                                                                     # Check if triangles is 1D.
                triangles = triangles.reshape(1, -1)                                                    # Reshape triangles to 2D.
            return triangles                                                                            # Return cached triangulation.
        except Exception:                                                                               # Ignore exceptions.
            pass                                                                                        # Do nothing.

    if not cloud_path or not os.path.exists(cloud_path):                                                # Check if cloud path is missing.
        return None                                                                                     # Return None.

    try:                                                                                                # Attempt operation.
        data = pd.read_csv(cloud_path)                                                                  # Read point cloud data from CSV.
        x = data['x'].values                                                                            # Extract X coordinates.
        y = data['y'].values                                                                            # Extract Y coordinates.
        if 'region' in data.columns:                                                                    # Check if region data exists.
            regions = data['region'].values                                                             # Extract region information.
        else:                                                                                           # Fallback if no regions.
            regions = np.ones_like(x, dtype=np.int32)                                                   # Assume single region.
    except Exception:                                                                                   # Ignore exceptions.
        return None                                                                                     # Return None.

    # Estimate typical point spacing
    dist = find_distances(np.column_stack([x, y, np.zeros_like(x)]), mode=3)                            # Estimate typical point spacing.
    mean_spacing = np.mean(dist)                                                                        # Calculate average spacing.
    threshold = 2.5 * mean_spacing                                                                      # Max allowed edge length

    # Triangulate the entire cloud
    global_triangles = []                                                                               # Initialize empty list for valid triangles.
    xy = np.column_stack([x, y])                                                                        # Create 2D coordinate array.
    
    try:                                                                                                # Attempt operation.
        tri = Delaunay(xy)                                                                              # Compute Delaunay triangulation.
        simplices = tri.simplices                                                                       # Extract Delaunay simplices.
        
        # Filter triangles by max edge length (Alpha shape concept)
        for t in simplices:                                                                             # Iterate over all triangles.
            pts = xy[t]                                                                                 # Get vertices of current triangle.
            L1 = np.linalg.norm(pts[0] - pts[1])                                                        # Calculate edge length.
            L2 = np.linalg.norm(pts[1] - pts[2])                                                        # Calculate edge length.
            L3 = np.linalg.norm(pts[2] - pts[0])                                                        # Calculate edge length.
            if max(L1, L2, L3) < threshold:                                                             # Filter out large boundary triangles.
                global_triangles.append(t)                                                              # Store valid triangle.
    except Exception:                                                                                   # Ignore exceptions.
        pass                                                                                            # Do nothing.

    if len(global_triangles) > 0:                                                                       # Check if any valid triangles were found.
        global_triangles = np.array(global_triangles, dtype=np.int32)                                   # Convert list to NumPy array.
        if cache_path:                                                                                  # Check if cache path is defined.
            try:                                                                                        # Attempt operation.
                np.savetxt(cache_path, global_triangles, delimiter=',', fmt='%d')                       # Save triangulation to cache.
            except Exception:                                                                           # Ignore exceptions.
                pass                                                                                    # Do nothing.
        return global_triangles                                                                         # Return computed triangulation.

    return None                                                                                         # Return None.


def Cloud_Stationary(p, u_ap, u_ex, save = False, nom = ''):                                            # Function definition.
    """
    Cloud_Stationary
    Render a side-by-side 3D scatter plot of the approximate solution vs the exact solution
    for a stationary PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        u_ex        m               ndarray         Array with the theoretical solution evaluated at all nodes.
        save                        bool            If True, saves the figure to disk instead of displaying.
        nom                         str             Base filename for the saved output (if save=True).
    """

    ## Variable initialization.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    # Calculate absolute error and its base-10 logarithm safely
    err = np.abs(u_ap - u_ex)                                                                           # Calculate absolute error.
    min_err = 0.0                                                                                       # Minimum absolute error is 0.
    max_err = err.max()                                                                                 # Maximum absolute error.
    if max_err == 0:                                                                                    # Check if error is exactly zero.
        max_err = 1.0                                                                                   # Maximum absolute error.

    log_err = np.log10(err + 1e-15)                                                                     # Calculate base-10 logarithm of error.
    min_log = log_err.min()                                                                             # Minimum log error.
    max_log = log_err.max()                                                                             # Maximum log error.

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot on given axes
    def draw_plot(fig_obj, ax1, ax2, ax3, ax4, angle_view=True):                                        # Helper function to render subplots.
        if triangles is not None:                                                                       # Check if triangulation is available.
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err[:], triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err[:], triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
        else:                                                                                           # Fallback if no regions.
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)  # Plot scattered points as fallback.
            s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:], c = u_ex[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)  # Plot scattered points as fallback.
            s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err[:], c = err[:], cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)  # Plot scattered points as fallback.
            s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err[:], c = log_err[:], cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)  # Plot scattered points as fallback.

        for ax in (ax1, ax2, ax3, ax4):                                                                 # Iterate over all subplots.
            ax.set_box_aspect(box_aspect)                                                               # Apply aspect ratio.
            ax.xaxis.pane.fill = False                                                                  # Make X pane transparent.
            ax.yaxis.pane.fill = False                                                                  # Make Y pane transparent.
            ax.zaxis.pane.fill = False                                                                  # Make Z pane transparent.
            ax.xaxis.pane.set_edgecolor('w')                                                            # Set X pane edge color to white.
            ax.yaxis.pane.set_edgecolor('w')                                                            # Set Y pane edge color to white.
            ax.zaxis.pane.set_edgecolor('w')                                                            # Set Z pane edge color to white.
            ax.grid(True, linestyle = '--', alpha = 0.5)                                                # Enable grid lines with transparency.
            ax.set_xlabel('X')                                                                          # Set X-axis label.
            ax.set_ylabel('Y')                                                                          # Set Y-axis label.
            
            # Use tight camera zoom to prevent wasted margins
            if angle_view:                                                                              # Check camera angle setting.
                ax.dist = 6.2                                                                           # Adjust camera zoom.
                ax.view_init(elev = 20, azim = -45)                                                     # Set camera viewing angle.
            else:                                                                                       # Fallback if no regions.
                ax.dist = 5.2                                                                           # Adjust camera zoom.
                ax.set_zticks([])                                                                       # Remove Z-axis ticks for top-down view.
                ax.view_init(elev = 90, azim = 270)                                                     # Set camera viewing angle.

        ax1.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax2.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax3.set_zlim([min_err, max_err])                                                                # Set Z-axis limits.
        ax4.set_zlim([min_log, max_log])                                                                # Set Z-axis limits.

        if angle_view:                                                                                  # Check camera angle setting.
            ax1.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax2.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax3.set_zlabel('Error')                                                                     # Set Z-axis label.
            ax4.set_zlabel('Log10(Error)')                                                              # Set Z-axis label.

        ax1.set_title('Approximation', y = 0.94)                                                        # Set title for subplot.
        ax2.set_title('Theoretical Solution', y = 0.94)                                                 # Set title for subplot.
        ax3.set_title('Absolute Error', y = 0.94)                                                       # Set title for subplot.
        ax4.set_title('Log10(Absolute Error)', y = 0.94)                                                # Set title for subplot.
        
        # Maximize grid layout area and reduce margins based on view angle
        if angle_view:                                                                                  # Check camera angle setting.
            # Take full width since there are no colorbars (Z-axis scale is already visible on the plots)
            fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
        else:                                                                                           # Fallback if no regions.
            # Leave space on the right of each subplot for colorbars in the 2D top view
            fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            
            # Position colorbars at their respective subplots in Top View (right side of each plot)
            fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.
            fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.
            fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.


    if save:                                                                                            # Check if output should be saved.
        # Save perspective view
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])                     # Set up grid layout for subplots.
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)                                             # Call helper to render subplots.
        fig.suptitle('Stationary Solution Comparison (Perspective)', y = 0.98, fontsize=14, fontweight='bold')
        plt.savefig(nom + '.png', bbox_inches = 'tight')                                                # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Save top view
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)                                # Call helper to render subplots.
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y = 0.98, fontsize=14, fontweight='bold')
        plt.savefig(nom + '_top.png', bbox_inches = 'tight')                                            # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to free memory.
    else:                                                                                               # Fallback if no regions.
        # Show both windows interactively
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])                     # Set up grid layout for subplots.
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)                                             # Call helper to render subplots.
        fig.suptitle('Stationary Solution Comparison (Perspective)', y = 0.98, fontsize=14, fontweight='bold')

        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)                                # Call helper to render subplots.
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y = 0.98, fontsize=14, fontweight='bold')

        plt.show()                                                                                      # Display interactive figure.

def Cloud_Stationary_1(p, u_ap, save = False, nom = ''):                                                # Function definition.
    """
    Cloud_Stationary_1
    Render a single 3D scatter plot of the approximate solution for a stationary PDE problem
    where the exact solution is unknown.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        save                        bool            If True, saves the figure to disk instead of displaying.
        nom                         str             Base filename for the saved output (if save=True).
    """

    ## Variable initialization.
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot on a given axis
    def draw_plot(fig_obj, ax1, angle_view=True):                                                       # Helper function to render subplots.
        if triangles is not None:                                                                       # Check if triangulation is available.
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
        else:                                                                                           # Fallback if no regions.
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)  # Plot scattered points as fallback.

        ax1.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax1.set_box_aspect(box_aspect)
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.xaxis.pane.set_edgecolor('w')
        ax1.yaxis.pane.set_edgecolor('w')
        ax1.zaxis.pane.set_edgecolor('w')
        ax1.grid(True, linestyle = '--', alpha = 0.5)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        if angle_view:                                                                                  # Check camera angle setting.
            ax1.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax1.view_init(elev = 20, azim = -45)
        else:                                                                                           # Fallback if no regions.
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)

        ax1.set_title('Approximation')                                                                  # Set title for subplot.
        
        # Colorbar
        fig_obj.subplots_adjust(right = 0.82)
        cbar_ax = fig_obj.add_axes([0.85, 0.15, 0.04, 0.7])
        fig_obj.colorbar(s1, cax = cbar_ax, label = 'Solution Amplitude')                               # Add colorbar to the figure.

    if save:                                                                                            # Check if output should be saved.
        # Save perspective view
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig, ax1, angle_view=True)                                                            # Call helper to render subplots.
        fig.suptitle('Stationary Approximation (Perspective)')
        plt.savefig(nom + '.png', bbox_inches = 'tight')                                                # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Save top view
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)                                                     # Call helper to render subplots.
        fig_top.suptitle('Stationary Approximation (Top View)')
        plt.savefig(nom + '_top.png', bbox_inches = 'tight')                                            # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to free memory.
    else:                                                                                               # Fallback if no regions.
        # Show both windows interactively
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig, ax1, angle_view=True)                                                            # Call helper to render subplots.
        fig.suptitle('Stationary Approximation (Perspective)')

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)                                                     # Call helper to render subplots.
        fig_top.suptitle('Stationary Approximation (Top View)')

        plt.show()                                                                                      # Display interactive figure.


def Cloud_Transient(p, u_ap, u_ex, save = False, nom = ''):                                             # Function definition.
    """
    Cloud_Transient
    Render a side-by-side animated 3D scatter plot of the approximate solution vs the exact solution
    for a transient PDE problem over time. 

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Approximation computed by the routine over all time steps.
        u_ex        m x t           ndarray         Theoretical solution evaluated at all nodes over all time steps.
        save                        bool            If True, saves the animation to disk as a video/gif.
        nom                         str             Base filename for the saved output (if save=True).
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 50)                                                                           # Frame stride for plotting/animation.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for titles).
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    # Compute absolute error globally across all time steps for consistent color scaling
    global_err = np.abs(u_ap - u_ex)
    min_err = 0.0                                                                                       # Minimum absolute error is 0.
    max_err = global_err.max()                                                                          # Maximum absolute error.
    if max_err == 0:                                                                                    # Check if error is exactly zero.
        max_err = 1.0                                                                                   # Maximum absolute error.

    global_log_err = np.log10(global_err + 1e-15)
    min_log = global_log_err.min()
    max_log = global_log_err.max()

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot a single frame on given axes
    def draw_frame(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):                                    # Function definition.
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])
        log_err_k = np.log10(err_k + 1e-15)
        
        if not hasattr(fig_obj, 'surf_artists'):                                                        # Initial setup for first frame.
            fig_obj.surf_artists = {}
            if triangles is not None:                                                                   # Check if triangulation is available.
                s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
                s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
                s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err_k, triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)
                s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err_k, triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)
            else:                                                                                       # Fallback if no regions.
                s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err_k, c = err_k, cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)
                s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err_k, c = log_err_k, cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)
            
            fig_obj.surf_artists['s1'] = s1                                                             # Cache artist.
            fig_obj.surf_artists['s2'] = s2                                                             # Cache artist.
            fig_obj.surf_artists['s3'] = s3                                                             # Cache artist.
            fig_obj.surf_artists['s4'] = s4                                                             # Cache artist.

            for ax in (ax1, ax2, ax3, ax4):                                                             # Iterate over all subplots.
                ax.set_box_aspect(box_aspect)                                                           # Apply aspect ratio.
                ax.xaxis.pane.fill = False                                                              # Make X pane transparent.
                ax.yaxis.pane.fill = False                                                              # Make Y pane transparent.
                ax.zaxis.pane.fill = False                                                              # Make Z pane transparent.
                ax.xaxis.pane.set_edgecolor('w')                                                        # Set X pane edge color to white.
                ax.yaxis.pane.set_edgecolor('w')                                                        # Set Y pane edge color to white.
                ax.zaxis.pane.set_edgecolor('w')                                                        # Set Z pane edge color to white.
                ax.grid(True, linestyle = '--', alpha = 0.5)                                            # Enable grid lines with transparency.
                ax.set_xlabel('X')                                                                      # Set X-axis label.
                ax.set_ylabel('Y')                                                                      # Set Y-axis label.
                
                if angle_view:                                                                          # Check camera angle setting.
                    ax.dist = 6.2                                                                       # Adjust camera zoom.
                    ax.view_init(elev = 20, azim = -45)                                                 # Set camera viewing angle.
                else:                                                                                   # Fallback if no regions.
                    ax.dist = 5.2                                                                       # Adjust camera zoom.
                    ax.set_zticks([])                                                                   # Remove Z-axis ticks for top-down view.
                    ax.view_init(elev = 90, azim = 270)                                                 # Set camera viewing angle.

            ax1.set_zlim([min_val, max_val])                                                            # Set Z-axis limits.
            ax2.set_zlim([min_val, max_val])                                                            # Set Z-axis limits.
            ax3.set_zlim([min_err, max_err])                                                            # Set Z-axis limits.
            ax4.set_zlim([min_log, max_log])                                                            # Set Z-axis limits.

            if angle_view:                                                                              # Check camera angle setting.
                ax1.set_zlabel('U(x, y)')                                                               # Set Z-axis label.
                ax2.set_zlabel('U(x, y)')                                                               # Set Z-axis label.
                ax3.set_zlabel('Error')                                                                 # Set Z-axis label.
                ax4.set_zlabel('Log10(Error)')                                                          # Set Z-axis label.

            ax1.set_title('Approximation', y = 0.94)                                                    # Set title for subplot.
            ax2.set_title('Theoretical Solution', y = 0.94)                                             # Set title for subplot.
            ax3.set_title('Absolute Error', y = 0.94)                                                   # Set title for subplot.
            ax4.set_title('Log10(Absolute Error)', y = 0.94)                                            # Set title for subplot.

            if angle_view:                                                                              # Check camera angle setting.
                fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            else:                                                                                       # Fallback if no regions.
                fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
                fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)                  # Add colorbar to the figure.
                fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)                  # Add colorbar to the figure.
                fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)                  # Add colorbar to the figure.
        else:                                                                                           # Fast path for subsequent frames.
            s1 = fig_obj.surf_artists['s1']                                                             # Retrieve cached artist.
            s2 = fig_obj.surf_artists['s2']                                                             # Retrieve cached artist.
            s3 = fig_obj.surf_artists['s3']                                                             # Retrieve cached artist.
            s4 = fig_obj.surf_artists['s4']                                                             # Retrieve cached artist.
            
            if triangles is not None:                                                                   # Check if triangulation is available.
                for artist, z_data in zip([s1, s2, s3, s4], [u_ap[:, k], u_ex[:, k], err_k, log_err_k]):
                    verts = np.stack((p[triangles, 0], p[triangles, 1], z_data[triangles]), axis=-1)    # Compute new vertices.
                    artist.set_verts(verts)                                                             # Update surface geometry.
                    artist.set_array(np.mean(z_data[triangles], axis=1))                                # Update face colors based on new Z.
            else:                                                                                       # Fallback if no regions.
                for artist, z_data in zip([s1, s2, s3, s4], [u_ap[:, k], u_ex[:, k], err_k, log_err_k]):
                    artist._offsets3d = (p[:, 0], p[:, 1], z_data)                                      # Update scatter coordinates directly.
                    artist.set_array(z_data)                                                            # Update scatter colors.


    if save:                                                                                            # Check if output should be saved.
        # Save perspective animation
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])                     # Set up grid layout for subplots.
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')

        def update_plot_perspective(k):                                                                 # Function definition.
            k = min(k, t - 1)
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
            draw_frame(fig, ax1, ax2, ax3, ax4, k, angle_view=True)
            return fig,
        ani_p = FuncAnimation(fig, update_plot_perspective, frames=np.arange(0, t, step), interval=100)
        
        if nom.endswith('.gif'):
            ani_p.save(nom, writer = 'pillow', fps = 10)
        else:                                                                                           # Fallback if no regions.
            try:                                                                                        # Attempt operation.
                ani_p.save(nom, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")
                nom_gif = nom.replace('.mp4', '.gif') if nom.endswith('.mp4') else nom + '.gif'
                ani_p.save(nom_gif, writer = 'pillow', fps = 10)
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Save top view animation
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')

        def update_plot_top(k):                                                                         # Function definition.
            k = min(k, t - 1)
            tin = float(T[k])
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y = 0.98, fontsize=14, fontweight='bold')
            draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)
            return fig_top,
        ani_t = FuncAnimation(fig_top, update_plot_top, frames=np.arange(0, t, step), interval=100)
        
        # Save top view with suffix '_top'
        nom_top = nom.replace('.mp4', '_top.mp4').replace('.gif', '_top.gif')
        if nom == nom_top:
            nom_top = nom + '_top'
        
        if nom_top.endswith('.gif'):
            ani_t.save(nom_top, writer = 'pillow', fps = 10)
        else:                                                                                           # Fallback if no regions.
            try:                                                                                        # Attempt operation.
                ani_t.save(nom_top, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                nom_gif_top = nom_top.replace('.mp4', '.gif') if nom_top.endswith('.mp4') else nom_top + '.gif'
                ani_t.save(nom_gif_top, writer = 'pillow', fps = 10)
        plt.close(fig_top)                                                                              # Close figure to free memory.

    else:                                                                                               # Fallback if no regions.
        # Show both windows interactively and animate them synchronously
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])                     # Set up grid layout for subplots.
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')

        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        for k in np.arange(0, t, step):
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y = 0.98, fontsize=14, fontweight='bold')
            
            draw_frame(fig, ax1, ax2, ax3, ax4, k, angle_view=True)
            draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)

            plt.pause(0.01)
            
            ax1.clear()                                                                                 # Clear subplot for next frame.
            ax2.clear()                                                                                 # Clear subplot for next frame.
            ax3.clear()                                                                                 # Clear subplot for next frame.
            ax4.clear()                                                                                 # Clear subplot for next frame.
            ax1_t.clear()
            ax2_t.clear()
            ax3_t.clear()
            ax4_t.clear()
            
        # Final step
        tin = float(T[-1])
        fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        
        draw_frame(fig, ax1, ax2, ax3, ax4, t - 1, angle_view=True)
        draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, t - 1, angle_view=False)
        
        plt.pause(0.1)
        plt.close(fig)                                                                                  # Close figure to free memory.
        plt.close(fig_top)                                                                              # Close figure to free memory.


def Cloud_Transient_Steps(p, u_ap, u_ex, nom):                                                          # Function definition.
    """
    Cloud_Transient_Steps
    Render and save a grid of side-by-side 3D scatter plots comparing the approximate solution vs 
    the exact solution at specific snapshot time steps (e.g. t=0, middle, and end).

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Approximation computed by the routine over all time steps.
        u_ex        m x t           ndarray         Theoretical solution evaluated at all nodes over all time steps.
        nom                         str             Base filename for the saved snapshot figures.
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 3)                                                                            # Target about 3 snapshots.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference).
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for filenames).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    # Compute absolute error globally across all time steps for consistent color scaling
    global_err = np.abs(u_ap - u_ex)
    min_err = 0.0                                                                                       # Minimum absolute error is 0.
    max_err = global_err.max()                                                                          # Maximum absolute error.
    if max_err == 0:                                                                                    # Check if error is exactly zero.
        max_err = 1.0                                                                                   # Maximum absolute error.

    global_log_err = np.log10(global_err + 1e-15)
    min_log = global_log_err.min()
    max_log = global_log_err.max()

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot
    def draw_snapshot(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):                                 # Function definition.
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])
        log_err_k = np.log10(err_k + 1e-15)
        
        if triangles is not None:                                                                       # Check if triangulation is available.
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err_k, triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
            s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err_k, triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)  # Plot surface using triangulation.
        else:                                                                                           # Fallback if no regions.
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)  # Plot scattered points as fallback.
            s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)  # Plot scattered points as fallback.
            s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err_k, c = err_k, cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)  # Plot scattered points as fallback.
            s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err_k, c = log_err_k, cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)  # Plot scattered points as fallback.

        for ax in (ax1, ax2, ax3, ax4):                                                                 # Iterate over all subplots.
            ax.set_box_aspect(box_aspect)                                                               # Apply aspect ratio.
            ax.xaxis.pane.fill = False                                                                  # Make X pane transparent.
            ax.yaxis.pane.fill = False                                                                  # Make Y pane transparent.
            ax.zaxis.pane.fill = False                                                                  # Make Z pane transparent.
            ax.xaxis.pane.set_edgecolor('w')                                                            # Set X pane edge color to white.
            ax.yaxis.pane.set_edgecolor('w')                                                            # Set Y pane edge color to white.
            ax.zaxis.pane.set_edgecolor('w')                                                            # Set Z pane edge color to white.
            ax.grid(True, linestyle = '--', alpha = 0.5)                                                # Enable grid lines with transparency.
            ax.set_xlabel('X')                                                                          # Set X-axis label.
            ax.set_ylabel('Y')                                                                          # Set Y-axis label.
            
            if angle_view:                                                                              # Check camera angle setting.
                ax.dist = 6.2                                                                           # Adjust camera zoom.
                ax.view_init(elev = 20, azim = -45)                                                     # Set camera viewing angle.
            else:                                                                                       # Fallback if no regions.
                ax.dist = 5.2                                                                           # Adjust camera zoom.
                ax.set_zticks([])                                                                       # Remove Z-axis ticks for top-down view.
                ax.view_init(elev = 90, azim = 270)                                                     # Set camera viewing angle.

        ax1.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax2.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax3.set_zlim([min_err, max_err])                                                                # Set Z-axis limits.
        ax4.set_zlim([min_log, max_log])                                                                # Set Z-axis limits.

        if angle_view:                                                                                  # Check camera angle setting.
            ax1.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax2.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax3.set_zlabel('Error')                                                                     # Set Z-axis label.
            ax4.set_zlabel('Log10(Error)')                                                              # Set Z-axis label.

        ax1.set_title('Approximation', y = 0.94)                                                        # Set title for subplot.
        ax2.set_title('Theoretical Solution', y = 0.94)                                                 # Set title for subplot.
        ax3.set_title('Absolute Error', y = 0.94)                                                       # Set title for subplot.
        ax4.set_title('Log10(Absolute Error)', y = 0.94)                                                # Set title for subplot.
        
        if angle_view:                                                                                  # Check camera angle setting.
            fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
        else:                                                                                           # Fallback if no regions.
            fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            
            fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.
            fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.
            fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)                      # Add colorbar to the figure.


    ## Create the graphs.
    for k in np.arange(0, t + 1, step):
        if k >= t:
            k = t - 1
        tin = float(T[k])
        nok = nom + '_' + str(format(T[k], '.2f'))

        # Perspective snapshot
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])                     # Set up grid layout for subplots.
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        plt.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        draw_snapshot(fig, ax1, ax2, ax3, ax4, k, angle_view=True)
        plt.savefig(nok + 's.png', bbox_inches = 'tight')                                               # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Top view snapshot
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        plt.suptitle('Solution at t = %1.3f s (Top View)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        draw_snapshot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches = 'tight')                                           # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to free memory.

def Cloud_Transient_1(p, u_ap, save = False, nom = ''):                                                 # Function definition.
    """
    Cloud_Transient_1
    Render an animated 3D scatter plot of the approximate solution for a transient PDE problem
    where the exact solution is unknown.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Approximation computed by the routine over all time steps.
        save                        bool            If True, saves the animation to disk as a video/gif.
        nom                         str             Base filename for the saved output (if save=True).
    """
    ## Variable initialization.
    t       = u_ap.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 50)                                                                           # Frame stride for plotting/animation.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for titles).
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot a single frame
    def draw_frame(ax1, k, angle_view=True):                                                            # Function definition.
        if triangles is not None:                                                                       # Check if triangulation is available.
            ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:                                                                                           # Fallback if no regions.
            ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)

        ax1.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax1.set_box_aspect(box_aspect)
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.xaxis.pane.set_edgecolor('w')
        ax1.yaxis.pane.set_edgecolor('w')
        ax1.zaxis.pane.set_edgecolor('w')
        ax1.grid(True, linestyle = '--', alpha = 0.5)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        if angle_view:                                                                                  # Check camera angle setting.
            ax1.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax1.view_init(elev = 20, azim = -45)
        else:                                                                                           # Fallback if no regions.
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)
        ax1.set_title('Approximation')                                                                  # Set title for subplot.

    if save:                                                                                            # Save an animation to disk.
        # Perspective animation
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        def update_plot_perspective(k):                                                                 # Function definition.
            k = min(k, t - 1)
            ax1.clear()                                                                                 # Clear subplot for next frame.
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
            draw_frame(ax1, k, angle_view=True)
            return fig,
        ani_p = FuncAnimation(fig, update_plot_perspective, frames=np.arange(0, t, step), interval=100)
        
        if nom.endswith('.gif'):
            ani_p.save(nom, writer = 'pillow', fps = 10)
        else:                                                                                           # Fallback if no regions.
            try:                                                                                        # Attempt operation.
                ani_p.save(nom, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")
                nom_gif = nom.replace('.mp4', '.gif') if nom.endswith('.mp4') else nom + '.gif'
                ani_p.save(nom_gif, writer = 'pillow', fps = 10)
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Top view animation
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        def update_plot_top(k):                                                                         # Function definition.
            k = min(k, t - 1)
            ax1_t.clear()
            tin = float(T[k])
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin)
            draw_frame(ax1_t, k, angle_view=False)
            return fig_top,
        ani_t = FuncAnimation(fig_top, update_plot_top, frames=np.arange(0, t, step), interval=100)
        
        # Save top view with suffix '_top'
        nom_top = nom.replace('.mp4', '_top.mp4').replace('.gif', '_top.gif')
        if nom == nom_top:
            nom_top = nom + '_top'
        
        if nom_top.endswith('.gif'):
            ani_t.save(nom_top, writer = 'pillow', fps = 10)
        else:                                                                                           # Fallback if no regions.
            try:                                                                                        # Attempt operation.
                ani_t.save(nom_top, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                nom_gif_top = nom_top.replace('.mp4', '.gif') if nom_top.endswith('.mp4') else nom_top + '.gif'
                ani_t.save(nom_gif_top, writer = 'pillow', fps = 10)
        plt.close(fig_top)                                                                              # Close figure to free memory.

    else:                                                                                               # Fallback if no regions.
        # Show both windows interactively and animate them synchronously
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        
        for k in np.arange(0, t, step):
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin)
            
            draw_frame(ax1, k, angle_view=True)
            draw_frame(ax1_t, k, angle_view=False)

            plt.pause(0.01)
            
            ax1.clear()                                                                                 # Clear subplot for next frame.
            ax1_t.clear()

        # Final step
        tin = float(T[-1])
        fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
        fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin)
        
        draw_frame(ax1, t - 1, angle_view=True)
        draw_frame(ax1_t, t - 1, angle_view=False)
        
        plt.pause(0.1)
        plt.close(fig)                                                                                  # Close figure to free memory.
        plt.close(fig_top)                                                                              # Close figure to free memory.


def Cloud_Transient_Steps_1(p, u_ap, nom = ''):                                                         # Function definition.
    """
    Cloud_Transient_Steps_1
    Render and save a grid of 3D scatter plots showing the approximate solution at specific 
    snapshot time steps for a transient problem where the exact solution is unknown.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Approximation computed by the routine over all time steps.
        nom                         str             Base filename for the saved snapshot figures.
    """

    ## Variable initialization.
    t       = u_ap.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 3)                                                                            # Target about 3 snapshots.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for filenames).
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract Y bounds.
    x_range = x_max - x_min                                                                             # Calculate X range.
    y_range = y_max - y_min                                                                             # Calculate Y range.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Calculate maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio.

    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    # Helper function to plot
    def draw_snapshot(fig_obj, ax1, k, angle_view=True):                                                # Function definition.
        if triangles is not None:                                                                       # Check if triangulation is available.
            ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:                                                                                           # Fallback if no regions.
            ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)

        ax1.set_zlim([min_val, max_val])                                                                # Set Z-axis limits.
        ax1.set_box_aspect(box_aspect)
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.xaxis.pane.set_edgecolor('w')
        ax1.yaxis.pane.set_edgecolor('w')
        ax1.zaxis.pane.set_edgecolor('w')
        ax1.grid(True, linestyle = '--', alpha = 0.5)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        if angle_view:                                                                                  # Check camera angle setting.
            ax1.set_zlabel('U(x, y)')                                                                   # Set Z-axis label.
            ax1.view_init(elev = 20, azim = -45)
        else:                                                                                           # Fallback if no regions.
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)
        ax1.set_title('Approximation')                                                                  # Set title for subplot.

    ## Create the graph.
    for k in np.arange(0, t + 1, step):                                                                 # Iterate selected snapshot indices.
        k = min(k, t - 1)                                                                               # Clamp index to valid range.
        tin = float(T[k])                                                                               # Current normalized time.
        nok = nom + '_' + str(format(T[k], '.2f'))                                                      # Filename tag based on time.
        
        # Perspective snapshot
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        plt.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
        draw_snapshot(fig, ax1, k, angle_view=True)
        plt.savefig(nok + 's.png', bbox_inches = 'tight')                                               # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to free memory.

        # Top view snapshot
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        plt.suptitle('Solution at t = %1.3f s (Top View)' % tin)
        draw_snapshot(fig_top, ax1_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches = 'tight')                                           # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to free memory.
