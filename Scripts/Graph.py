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
    May, 2024.
Last Modification:
    August, 2026.
"""
## Library importation.
import numpy as np                                                                                      # Array utilities (min/max, linspace, arange).
import matplotlib.pyplot as plt                                                                         # Plotting interface.
from matplotlib import cm                                                                               # Colormaps.
from matplotlib import animation
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

    parts = nom.split(os.sep)                                                                           # Split path to parse dataset details.
    dataset = None                                                                                      # Initialize dataset variable.
    scale = None                                                                                        # Initialize scale variable.
    workspace_root = None                                                                               # Initialize workspace root variable.

    if len(parts) >= 4:                                                                                 # Check if path has enough components.
        scale = parts[-2]                                                                               # Extract scale name.
        dataset = parts[-3]                                                                             # Extract dataset name.

        temp = nom                                                                                      # Initialize temp path for workspace search.
        for _ in range(10):                                                                             # Climb up directory tree.
            temp = os.path.dirname(temp)                                                                # Go up one directory level.
            if os.path.exists(os.path.join(temp, 'Data')):                                              # Check if Data folder exists.
                workspace_root = temp                                                                   # Set workspace root.
                break                                                                                   # Stop climbing directory tree.

    cache_path = None                                                                                   # Initialize cache path.
    cloud_path = None                                                                                   # Initialize cloud path.
    if workspace_root and dataset and scale:                                                            # Check if path info was successfully parsed.
        cloud_file = f"{dataset}_cloud.csv"                                                             # Format cloud filename.
        cloud_path = os.path.join(workspace_root, 'Data', dataset, scale, cloud_file)                   # Construct cloud path.
        if not os.path.exists(cloud_path):                                                              # Check if cloud file exists.
            lake_scale_dir = os.path.join(workspace_root, 'Data', dataset, scale)                       # Construct fallback directory.
            if os.path.exists(lake_scale_dir):                                                          # Check if fallback directory exists.
                for f in os.listdir(lake_scale_dir):                                                    # Iterate over fallback directory files.
                    if f.endswith('.csv') and not f.endswith('neighbors.csv') and not f.endswith('triangulation.csv'):
                        cloud_path = os.path.join(lake_scale_dir, f)                                    # Construct cloud path.
                        break                                                                           # Stop climbing directory tree.
        
        if cloud_path and os.path.exists(cloud_path):                                                   # Check if valid cloud path was found.
            cache_path = cloud_path.replace('.csv', '_triangulation.csv')                               # Set cache path based on cloud path.

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
    except Exception:                                                                                   # Ignore exceptions.
        return None                                                                                     # Return None.

    dist = find_distances(np.column_stack([x, y, np.zeros_like(x)]), mode=3)                            # Estimate typical point spacing.
    mean_spacing = np.mean(dist)                                                                        # Calculate average spacing.
    threshold = 2.5 * mean_spacing                                                                      # Max allowed edge length

    global_triangles = []                                                                               # Initialize empty list for valid triangles.
    xy = np.column_stack([x, y])                                                                        # Create 2D coordinate array.
    
    try:                                                                                                # Attempt operation.
        tri = Delaunay(xy)                                                                              # Compute Delaunay triangulation.
        simplices = tri.simplices                                                                       # Extract Delaunay simplices.
        
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


def _setup_3d_axes(ax, angle_view, box_aspect, x_bounds, y_bounds, z_bounds, z_label, title):
    """Helper to standardize 3D axes formatting."""
    ax.set_box_aspect(box_aspect)                                                                       # Apply physical aspect ratio to axes.
    ax.xaxis.pane.fill = False                                                                          # Make X pane transparent.
    ax.yaxis.pane.fill = False                                                                          # Make Y pane transparent.
    ax.zaxis.pane.fill = False                                                                          # Make Z pane transparent.
    ax.xaxis.pane.set_edgecolor('w')                                                                    # Set X pane edge color to white.
    ax.yaxis.pane.set_edgecolor('w')                                                                    # Set Y pane edge color to white.
    ax.zaxis.pane.set_edgecolor('w')                                                                    # Set Z pane edge color to white.
    ax.grid(True, linestyle='--', alpha=0.5)                                                            # Enable grid lines with transparency.
    ax.set_xlabel('X')                                                                                  # Set X-axis label.
    ax.set_ylabel('Y')                                                                                  # Set Y-axis label.
    
    ax.set_xlim(x_bounds)                                                                               # Set X-axis limits.
    ax.set_ylim(y_bounds)                                                                               # Set Y-axis limits.
    ax.set_zlim(z_bounds)                                                                               # Set Z-axis limits.
    ax.autoscale(False)                                                                                 # Disable automatic scaling.

    if angle_view:
        ax.dist = 6.2                                                                                   # Adjust camera zoom for perspective.
        ax.view_init(elev=20, azim=-45)                                                                 # Set perspective viewing angle.
        ax.set_zlabel(z_label)                                                                          # Set Z-axis label.
    else:
        ax.dist = 5.2                                                                                   # Adjust camera zoom for top view.
        ax.set_zticks([])                                                                               # Remove Z-axis ticks for top view.
        ax.view_init(elev=90, azim=270)                                                                 # Set top-down viewing angle.
        
    if title:
        ax.set_title(title, y=0.94)                                                                     # Set title for subplot.


def _render_surface(ax, p, data, triangles, cmap, vmin, vmax):
    """Helper to render trisurf or fallback to scatter."""
    if triangles is not None:                                                                           # Check if triangulation is available.
        return ax.plot_trisurf(p[:, 0], p[:, 1], data, triangles=triangles, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='none', linewidth=0, antialiased=False) # Plot surface using valid triangulation.
    else:
        return ax.scatter(p[:, 0], p[:, 1], zs=data, c=data, cmap=cmap, s=1, vmin=vmin, vmax=vmax)      # Fallback to plotting scattered points.


def _create_colorbars(fig_obj, min_val, max_val, min_err, max_err, min_log, max_log):
    """Helper to attach static colorbars in Top View layout."""
    sm2 = cm.ScalarMappable(cmap=cm.coolwarm, norm=plt.Normalize(vmin=min_val, vmax=max_val))           # Create static ScalarMappable for approximation.
    sm2.set_array([])                                                                                   # Initialize empty array for colorbar.
    cax2 = fig_obj.add_axes([0.90, 0.55, 0.015, 0.35])                                                  # Absolute positioning for approximation colorbar.
    fig_obj.colorbar(sm2, cax=cax2)                                                                     # Add colorbar to fixed axis.

    sm3 = cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=min_err, vmax=max_err))            # Create static ScalarMappable for error.
    sm3.set_array([])                                                                                   # Initialize empty array for colorbar.
    cax3 = fig_obj.add_axes([0.46, 0.05, 0.015, 0.35])                                                  # Absolute positioning for error colorbar.
    fig_obj.colorbar(sm3, cax=cax3)                                                                     # Add colorbar to fixed axis.

    sm4 = cm.ScalarMappable(cmap=cm.plasma, norm=plt.Normalize(vmin=min_log, vmax=max_log))             # Create static ScalarMappable for log error.
    sm4.set_array([])                                                                                   # Initialize empty array for colorbar.
    cax4 = fig_obj.add_axes([0.90, 0.05, 0.015, 0.35])                                                  # Absolute positioning for log error colorbar.
    fig_obj.colorbar(sm4, cax=cax4)                                                                     # Add colorbar to fixed axis.


def _create_colorbar_1(fig_obj, min_val, max_val):
    """Helper to attach a single colorbar for single-plot modes."""
    fig_obj.subplots_adjust(right=0.82)                                                                 # Adjust subplot margins and spacing.
    cbar_ax = fig_obj.add_axes([0.85, 0.15, 0.04, 0.7])
    sm = cm.ScalarMappable(cmap=cm.coolwarm, norm=plt.Normalize(vmin=min_val, vmax=max_val))
    sm.set_array([])
    fig_obj.colorbar(sm, cax=cbar_ax, label='Solution Amplitude')


def _get_aspect_and_bounds(p):
    """Helper to extract physical aspect ratio and bounds."""
    x_min, x_max = p[:, 0].min(), p[:, 0].max()                                                         # Extract global X bounds.
    y_min, y_max = p[:, 1].min(), p[:, 1].max()                                                         # Extract global Y bounds.
    x_range, y_range = x_max - x_min, y_max - y_min                                                     # Calculate physical X and Y ranges.
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                             # Determine maximum range for isometric scaling.
    box_aspect = (x_range, y_range, 0.4 * max_range)                                                    # Set 3D box aspect ratio to preserve geometry.
    return box_aspect, [x_min, x_max], [y_min, y_max]                                                   # Return layout parameters.


def Cloud_Stationary(p, u_ap, u_ex, save=False, nom=''):
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
    min_val, max_val = u_ex.min(), u_ex.max()                                                           # Shared color scale from reference solution.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.
    err = np.abs(u_ap - u_ex)                                                                           # Calculate absolute error.
    min_err, max_err = 0.0, err.max() if err.max() > 0 else 1.0                                         # Set error color scale bounds.
    log_err = np.log10(err + 1e-15)                                                                     # Calculate base-10 logarithm of error.
    min_log, max_log = log_err.min(), log_err.max()                                                     # Set log error color scale bounds.

    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_plot(fig_obj, ax1, ax2, ax3, ax4, angle_view=True):
        _render_surface(ax1, p, u_ap, triangles, cm.coolwarm, min_val, max_val)
        _render_surface(ax2, p, u_ex, triangles, cm.coolwarm, min_val, max_val)
        _render_surface(ax3, p, err, triangles, cm.viridis, min_err, max_err)
        _render_surface(ax4, p, log_err, triangles, cm.plasma, min_log, max_log)

        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
        _setup_3d_axes(ax2, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Theoretical Solution') # Format axes for exact solution.
        _setup_3d_axes(ax3, angle_view, box_aspect, x_bounds, y_bounds, [min_err, max_err], 'Error', 'Absolute Error') # Format axes for error.
        _setup_3d_axes(ax4, angle_view, box_aspect, x_bounds, y_bounds, [min_log, max_log], 'Log10(Error)', 'Log10(Absolute Error)') # Format axes for log error.

        if angle_view:
            fig_obj.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
        else:
            fig_obj.subplots_adjust(left=0.02, right=0.88, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
            _create_colorbars(fig_obj, min_val, max_val, min_err, max_err, min_log, max_log)            # Attach static colorbars.

    if save:
        fig = plt.figure(figsize=(11, 9))                                                               # Create new figure.
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                         # Set up grid layout for subplots.
        ax1, ax2, ax3, ax4 = fig.add_subplot(gs[0, 0], projection='3d'), fig.add_subplot(gs[0, 1], projection='3d'), fig.add_subplot(gs[1, 0], projection='3d'), fig.add_subplot(gs[1, 1], projection='3d') # Initialize 3D subplots.
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)                                             # Render the subplots.
        fig.suptitle('Stationary Solution Comparison (Perspective)', y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.
        plt.savefig(nom + '.png', bbox_inches='tight')                                                  # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top = plt.figure(figsize=(13, 8))                                                           # Create new figure for top view.
        gs_t = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                       # Set up grid layout for top view subplots.
        ax1_t, ax2_t, ax3_t, ax4_t = fig_top.add_subplot(gs_t[0, 0], projection='3d'), fig_top.add_subplot(gs_t[0, 1], projection='3d'), fig_top.add_subplot(gs_t[1, 0], projection='3d'), fig_top.add_subplot(gs_t[1, 1], projection='3d') # Initialize 3D subplots for top view.
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
        plt.savefig(nom + '_top.png', bbox_inches='tight')                                              # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to release memory.
    else:
        fig = plt.figure(figsize=(11, 9))                                                               # Create new figure.
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                         # Set up grid layout for subplots.
        ax1, ax2, ax3, ax4 = fig.add_subplot(gs[0, 0], projection='3d'), fig.add_subplot(gs[0, 1], projection='3d'), fig.add_subplot(gs[1, 0], projection='3d'), fig.add_subplot(gs[1, 1], projection='3d') # Initialize 3D subplots.
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)                                             # Render the subplots.
        fig.suptitle('Stationary Solution Comparison (Perspective)', y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.

        fig_top = plt.figure(figsize=(13, 8))                                                           # Create new figure for top view.
        gs_t = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                       # Set up grid layout for top view subplots.
        ax1_t, ax2_t, ax3_t, ax4_t = fig_top.add_subplot(gs_t[0, 0], projection='3d'), fig_top.add_subplot(gs_t[0, 1], projection='3d'), fig_top.add_subplot(gs_t[1, 0], projection='3d'), fig_top.add_subplot(gs_t[1, 1], projection='3d') # Initialize 3D subplots for top view.
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
        plt.show()                                                                                      # Display the interactive plot window.


def Cloud_Stationary_1(p, u_ap, save=False, nom=''):
    """
    Cloud_Stationary_1
    Render a single 3D scatter plot of the approximate solution for a stationary PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        save                        bool            If True, saves the figure to disk instead of displaying.
        nom                         str             Base filename for the saved output (if save=True).
    """
    min_val, max_val = u_ap.min(), u_ap.max()                                                           # Shared color scale from approximation.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.
    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_plot(fig_obj, ax1, angle_view=True):
        _render_surface(ax1, p, u_ap, triangles, cm.coolwarm, min_val, max_val)
        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
        if not angle_view:
            _create_colorbar_1(fig_obj, min_val, max_val)                                               # Attach single static colorbar.

    if save:
        fig, ax1 = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        draw_plot(fig, ax1, angle_view=True)                                                            # Render the subplots.
        fig.suptitle('Stationary Approximation (Perspective)')                                          # Set main title for the figure.
        plt.savefig(nom + '.png', bbox_inches='tight')                                                  # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)
        fig_top.suptitle('Stationary Approximation (Top View)')                                         # Set main title for the top view figure.
        plt.savefig(nom + '_top.png', bbox_inches='tight')                                              # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to release memory.
    else:
        fig, ax1 = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        draw_plot(fig, ax1, angle_view=True)                                                            # Render the subplots.
        fig.suptitle('Stationary Approximation (Perspective)')                                          # Set main title for the figure.

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)
        fig_top.suptitle('Stationary Approximation (Top View)')                                         # Set main title for the top view figure.
        plt.show()                                                                                      # Display the interactive plot window.


def Cloud_Transient(p, u_ap, u_ex, save=False, nom=''):
    """
    Cloud_Transient
    Render a side-by-side 3D scatter plot of the approximate solution vs the exact solution
    for a transient PDE problem, and optionally save it as an animation.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Matrix with the approximation computed by the routine at each time step.
        u_ex        m x t           ndarray         Matrix with the theoretical solution evaluated at all nodes at each time step.
        save                        bool            If True, saves the animation to disk instead of displaying interactively.
        nom                         str             Base filename for the saved output (if save=True).
    """
    t = u_ex.shape[1]                                                                                   # Get total number of time steps.
    step = max(1, t // 50)                                                                              # Calculate step size to limit frame count.
    T = np.linspace(0, 1, t)                                                                            # Generate array of time values.
    min_val, max_val = u_ex.min(), u_ex.max()                                                           # Shared color scale from reference solution.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.
    global_err = np.abs(u_ap - u_ex)                                                                    # Calculate global absolute error.
    min_err, max_err = 0.0, global_err.max() if global_err.max() > 0 else 1.0                           # Set error color scale bounds.
    global_log_err = np.log10(global_err + 1e-15)                                                       # Calculate base-10 logarithm of global error.
    min_log, max_log = global_log_err.min(), global_log_err.max()                                       # Set log error color scale bounds.

    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_frame(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])                                                         # Calculate error for the current time step.
        log_err_k = np.log10(err_k + 1e-15)                                                             # Calculate log error for the current time step.
        
        if not hasattr(fig_obj, 'surf_artists'):
            s1 = _render_surface(ax1, p, u_ap[:, k], triangles, cm.coolwarm, min_val, max_val)          # Render approximation surface.
            s2 = _render_surface(ax2, p, u_ex[:, k], triangles, cm.coolwarm, min_val, max_val)          # Render exact solution surface.
            s3 = _render_surface(ax3, p, err_k, triangles, cm.viridis, min_err, max_err)                # Render error surface.
            s4 = _render_surface(ax4, p, log_err_k, triangles, cm.plasma, min_log, max_log)             # Render log error surface.
            fig_obj.surf_artists = {'s1': s1, 's2': s2, 's3': s3, 's4': s4}                             # Cache surface artists for animation updates.

            _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
            _setup_3d_axes(ax2, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Theoretical Solution') # Format axes for exact solution.
            _setup_3d_axes(ax3, angle_view, box_aspect, x_bounds, y_bounds, [min_err, max_err], 'Error', 'Absolute Error') # Format axes for error.
            _setup_3d_axes(ax4, angle_view, box_aspect, x_bounds, y_bounds, [min_log, max_log], 'Log10(Error)', 'Log10(Absolute Error)') # Format axes for log error.

            if angle_view:
                fig_obj.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
            else:
                fig_obj.subplots_adjust(left=0.02, right=0.88, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
                _create_colorbars(fig_obj, min_val, max_val, min_err, max_err, min_log, max_log)        # Attach static colorbars.
        else:
            if triangles is not None:                                                                   # Check if triangulation is available.
                for s_key, z_data in zip(['s1', 's2', 's3', 's4'], [u_ap[:, k], u_ex[:, k], err_k, log_err_k]):
                    artist = fig_obj.surf_artists[s_key]                                                # Retrieve cached surface artist.
                    verts = np.stack((p[triangles, 0], p[triangles, 1], z_data[triangles]), axis=-1)    # Compute new 3D vertices for the surface.
                    artist.set_verts(verts)                                                             # Update surface geometry in animation.
                    artist.set_array(np.mean(z_data[triangles], axis=1))                                # Update surface colors based on new Z heights.
            else:
                for s_key, z_data in zip(['s1', 's2', 's3', 's4'], [u_ap[:, k], u_ex[:, k], err_k, log_err_k]):
                    artist = fig_obj.surf_artists[s_key]                                                # Retrieve cached surface artist.
                    artist._offsets3d = (p[:, 0], p[:, 1], z_data)                                      # Update scatter coordinates directly.
                    artist.set_array(z_data)                                                            # Update scatter colors.

    if save:
        fig = plt.figure(figsize=(11, 9))                                                               # Create new figure.
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                         # Set up grid layout for subplots.
        ax1, ax2, ax3, ax4 = fig.add_subplot(gs[0, 0], projection='3d'), fig.add_subplot(gs[0, 1], projection='3d'), fig.add_subplot(gs[1, 0], projection='3d'), fig.add_subplot(gs[1, 1], projection='3d') # Initialize 3D subplots.
        
        def update_p(frame):
            tin = float(T[frame])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.
            draw_frame(fig, ax1, ax2, ax3, ax4, frame, angle_view=True)                                 # Render the current frame.

        ani_p = FuncAnimation(fig, update_p, frames=np.arange(0, t, step), blit=False)                  # Initialize animation object.
        if animation.writers.is_available('ffmpeg'):                                                    # Check if FFmpeg is installed.
            ani_p.save(nom.replace('.gif', '.mp4'), writer='ffmpeg', fps=10)                            # Save animation as high-quality MP4.
        else:                                                                                           # Fallback to Pillow.
            ani_p.save(nom, writer='pillow', fps=10)                                                    # Save animation as GIF.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top = plt.figure(figsize=(13, 8))                                                           # Create new figure for top view.
        gs_t = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                       # Set up grid layout for top view subplots.
        ax1_t, ax2_t, ax3_t, ax4_t = fig_top.add_subplot(gs_t[0, 0], projection='3d'), fig_top.add_subplot(gs_t[0, 1], projection='3d'), fig_top.add_subplot(gs_t[1, 0], projection='3d'), fig_top.add_subplot(gs_t[1, 1], projection='3d') # Initialize 3D subplots for top view.
        
        def update_t(frame):
            tin = float(T[frame])
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
            draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, frame, angle_view=False)

        ani_t = FuncAnimation(fig_top, update_t, frames=np.arange(0, t, step), blit=False)              # Initialize animation object for top view.
        nom_top = nom.replace('.mp4', '_top.mp4').replace('.gif', '_top.gif')
        if nom == nom_top: nom_top = nom + '_top'
        if animation.writers.is_available('ffmpeg'):                                                    # Check if FFmpeg is installed.
            ani_t.save(nom_top.replace('.gif', '.mp4'), writer='ffmpeg', fps=10)                        # Save animation as high-quality MP4.
        else:                                                                                           # Fallback to Pillow.
            ani_t.save(nom_top, writer='pillow', fps=10)                                                # Save animation as GIF.
        plt.close(fig_top)                                                                              # Close figure to release memory.

    else:
        fig = plt.figure(figsize=(11, 9))                                                               # Create new figure.
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                         # Set up grid layout for subplots.
        ax1, ax2, ax3, ax4 = fig.add_subplot(gs[0, 0], projection='3d'), fig.add_subplot(gs[0, 1], projection='3d'), fig.add_subplot(gs[1, 0], projection='3d'), fig.add_subplot(gs[1, 1], projection='3d') # Initialize 3D subplots.

        fig_top = plt.figure(figsize=(13, 8))                                                           # Create new figure for top view.
        gs_t = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                       # Set up grid layout for top view subplots.
        ax1_t, ax2_t, ax3_t, ax4_t = fig_top.add_subplot(gs_t[0, 0], projection='3d'), fig_top.add_subplot(gs_t[0, 1], projection='3d'), fig_top.add_subplot(gs_t[1, 0], projection='3d'), fig_top.add_subplot(gs_t[1, 1], projection='3d') # Initialize 3D subplots for top view.
        
        for k in np.arange(0, t, step):
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.
            fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
            draw_frame(fig, ax1, ax2, ax3, ax4, k, angle_view=True)                                     # Render the current frame.
            draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)
            plt.pause(0.01)                                                                             # Pause to allow UI to update.
            
        tin = float(T[-1])
        fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.
        fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
        draw_frame(fig, ax1, ax2, ax3, ax4, t - 1, angle_view=True)                                     # Render the current frame.
        draw_frame(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, t - 1, angle_view=False)
        plt.pause(0.1)                                                                                  # Pause to allow UI to update.
        plt.close(fig)                                                                                  # Close figure to release memory.
        plt.close(fig_top)                                                                              # Close figure to release memory.


def Cloud_Transient_Steps(p, u_ap, u_ex, nom):
    """
    Cloud_Transient_Steps
    Render and save static side-by-side 3D scatter plots of the approximate solution vs the exact solution
    at a few key time steps during a transient PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Matrix with the approximation computed by the routine at each time step.
        u_ex        m x t           ndarray         Matrix with the theoretical solution evaluated at all nodes at each time step.
        nom                         str             Base filename for the saved output.
    """
    t = u_ex.shape[1]                                                                                   # Get total number of time steps.
    step = max(1, t // 3)                                                                               # Calculate step size to limit frame count.
    T = np.linspace(0, 1, t)                                                                            # Generate array of time values.
    min_val, max_val = u_ex.min(), u_ex.max()                                                           # Shared color scale from reference solution.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.
    global_err = np.abs(u_ap - u_ex)                                                                    # Calculate global absolute error.
    min_err, max_err = 0.0, global_err.max() if global_err.max() > 0 else 1.0                           # Set error color scale bounds.
    global_log_err = np.log10(global_err + 1e-15)                                                       # Calculate base-10 logarithm of global error.
    min_log, max_log = global_log_err.min(), global_log_err.max()                                       # Set log error color scale bounds.

    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_snapshot(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])                                                         # Calculate error for the current time step.
        log_err_k = np.log10(err_k + 1e-15)                                                             # Calculate log error for the current time step.
        
        _render_surface(ax1, p, u_ap[:, k], triangles, cm.coolwarm, min_val, max_val)
        _render_surface(ax2, p, u_ex[:, k], triangles, cm.coolwarm, min_val, max_val)
        _render_surface(ax3, p, err_k, triangles, cm.viridis, min_err, max_err)
        _render_surface(ax4, p, log_err_k, triangles, cm.plasma, min_log, max_log)

        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
        _setup_3d_axes(ax2, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Theoretical Solution') # Format axes for exact solution.
        _setup_3d_axes(ax3, angle_view, box_aspect, x_bounds, y_bounds, [min_err, max_err], 'Error', 'Absolute Error') # Format axes for error.
        _setup_3d_axes(ax4, angle_view, box_aspect, x_bounds, y_bounds, [min_log, max_log], 'Log10(Error)', 'Log10(Absolute Error)') # Format axes for log error.

        if angle_view:
            fig_obj.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
        else:
            fig_obj.subplots_adjust(left=0.02, right=0.88, bottom=0.02, top=0.94, hspace=0.08, wspace=0.08) # Adjust subplot margins and spacing.
            _create_colorbars(fig_obj, min_val, max_val, min_err, max_err, min_log, max_log)            # Attach static colorbars.

    for k in np.arange(0, t + 1, step):
        if k >= t: k = t - 1
        tin = float(T[k])
        nok = nom + '_' + str(format(T[k], '.2f'))

        fig = plt.figure(figsize=(11, 9))                                                               # Create new figure.
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                         # Set up grid layout for subplots.
        ax1, ax2, ax3, ax4 = fig.add_subplot(gs[0, 0], projection='3d'), fig.add_subplot(gs[0, 1], projection='3d'), fig.add_subplot(gs[1, 0], projection='3d'), fig.add_subplot(gs[1, 1], projection='3d') # Initialize 3D subplots.
        fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the figure.
        draw_snapshot(fig, ax1, ax2, ax3, ax4, k, angle_view=True)                                      # Render the current snapshot.
        plt.savefig(nok + 's.png', bbox_inches='tight')                                                 # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top = plt.figure(figsize=(13, 8))                                                           # Create new figure for top view.
        gs_t = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])                       # Set up grid layout for top view subplots.
        ax1_t, ax2_t, ax3_t, ax4_t = fig_top.add_subplot(gs_t[0, 0], projection='3d'), fig_top.add_subplot(gs_t[0, 1], projection='3d'), fig_top.add_subplot(gs_t[1, 0], projection='3d'), fig_top.add_subplot(gs_t[1, 1], projection='3d') # Initialize 3D subplots for top view.
        fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin, y=0.98, fontsize=14, fontweight='bold') # Set main title for the top view figure.
        draw_snapshot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches='tight')                                             # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to release memory.


def Cloud_Transient_1(p, u_ap, save=False, nom=''):
    """
    Cloud_Transient_1
    Render a single 3D scatter plot of the approximate solution for a transient PDE problem,
    and optionally save it as an animation.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Matrix with the approximation computed by the routine at each time step.
        save                        bool            If True, saves the animation to disk instead of displaying interactively.
        nom                         str             Base filename for the saved output (if save=True).
    """
    t = u_ap.shape[1]                                                                                   # Get total number of time steps.
    step = max(1, t // 50)                                                                              # Calculate step size to limit frame count.
    T = np.linspace(0, 1, t)                                                                            # Generate array of time values.
    min_val, max_val = u_ap.min(), u_ap.max()                                                           # Shared color scale from approximation.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.

    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_frame(fig_obj, ax1, k, angle_view=True):
        if not hasattr(fig_obj, 'surf_artists'):
            s1 = _render_surface(ax1, p, u_ap[:, k], triangles, cm.coolwarm, min_val, max_val)          # Render approximation surface.
            fig_obj.surf_artists = {'s1': s1}                                                           # Cache surface artists for animation updates.
            _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
            if not angle_view:
                _create_colorbar_1(fig_obj, min_val, max_val)                                           # Attach single static colorbar.
        else:
            if triangles is not None:                                                                   # Check if triangulation is available.
                artist = fig_obj.surf_artists['s1']                                                     # Retrieve cached surface artist.
                z_data = u_ap[:, k]                                                                     # Extract Z data for current step.
                verts = np.stack((p[triangles, 0], p[triangles, 1], z_data[triangles]), axis=-1)        # Compute new 3D vertices for the surface.
                artist.set_verts(verts)                                                                 # Update surface geometry in animation.
                artist.set_array(np.mean(z_data[triangles], axis=1))                                    # Update surface colors based on new Z heights.
            else:
                artist = fig_obj.surf_artists['s1']                                                     # Retrieve cached surface artist.
                z_data = u_ap[:, k]                                                                     # Extract Z data for current step.
                artist._offsets3d = (p[:, 0], p[:, 1], z_data)                                          # Update scatter coordinates directly.
                artist.set_array(z_data)                                                                # Update scatter colors.

    if save:
        fig, ax1 = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        def update_p(frame):
            tin = float(T[frame])
            fig.suptitle('Approximation at t = %1.3f s (Perspective)' % tin)                            # Set main title for the figure.
            draw_frame(fig, ax1, frame, angle_view=True)                                                # Render the current frame.

        ani_p = FuncAnimation(fig, update_p, frames=np.arange(0, t, step), blit=False)                  # Initialize animation object.
        if animation.writers.is_available('ffmpeg'):                                                    # Check if FFmpeg is installed.
            ani_p.save(nom.replace('.gif', '.mp4'), writer='ffmpeg', fps=10)                            # Save animation as high-quality MP4.
        else:                                                                                           # Fallback to Pillow.
            ani_p.save(nom, writer='pillow', fps=10)                                                    # Save animation as GIF.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        def update_t(frame):
            tin = float(T[frame])
            fig_top.suptitle('Approximation at t = %1.3f s (Top View)' % tin)                           # Set main title for the top view figure.
            draw_frame(fig_top, ax1_t, frame, angle_view=False)

        ani_t = FuncAnimation(fig_top, update_t, frames=np.arange(0, t, step), blit=False)              # Initialize animation object for top view.
        nom_top = nom.replace('.mp4', '_top.mp4').replace('.gif', '_top.gif')
        if nom == nom_top: nom_top = nom + '_top'
        if animation.writers.is_available('ffmpeg'):                                                    # Check if FFmpeg is installed.
            ani_t.save(nom_top.replace('.gif', '.mp4'), writer='ffmpeg', fps=10)                        # Save animation as high-quality MP4.
        else:                                                                                           # Fallback to Pillow.
            ani_t.save(nom_top, writer='pillow', fps=10)                                                # Save animation as GIF.
        plt.close(fig_top)                                                                              # Close figure to release memory.

    else:
        fig, ax1 = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        
        for k in np.arange(0, t, step):
            tin = float(T[k])
            fig.suptitle('Approximation at t = %1.3f s (Perspective)' % tin)                            # Set main title for the figure.
            fig_top.suptitle('Approximation at t = %1.3f s (Top View)' % tin)                           # Set main title for the top view figure.
            draw_frame(fig, ax1, k, angle_view=True)                                                    # Render the current frame.
            draw_frame(fig_top, ax1_t, k, angle_view=False)
            plt.pause(0.01)                                                                             # Pause to allow UI to update.
            
        tin = float(T[-1])
        fig.suptitle('Approximation at t = %1.3f s (Perspective)' % tin)                                # Set main title for the figure.
        fig_top.suptitle('Approximation at t = %1.3f s (Top View)' % tin)                               # Set main title for the top view figure.
        draw_frame(fig, ax1, t - 1, angle_view=True)                                                    # Render the current frame.
        draw_frame(fig_top, ax1_t, t - 1, angle_view=False)
        plt.pause(0.1)                                                                                  # Pause to allow UI to update.
        plt.close(fig)                                                                                  # Close figure to release memory.
        plt.close(fig_top)                                                                              # Close figure to release memory.


def Cloud_Transient_Steps_1(p, u_ap, nom):
    """
    Cloud_Transient_Steps_1
    Render and save static single 3D scatter plots of the approximate solution at a few key time steps
    during a transient PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u_ap        m x t           ndarray         Matrix with the approximation computed by the routine at each time step.
        nom                         str             Base filename for the saved output.
    """
    t = u_ap.shape[1]                                                                                   # Get total number of time steps.
    step = max(1, t // 3)                                                                               # Calculate step size to limit frame count.
    T = np.linspace(0, 1, t)                                                                            # Generate array of time values.
    min_val, max_val = u_ap.min(), u_ap.max()                                                           # Shared color scale from approximation.
    if min_val == max_val: max_val += 1e-6                                                              # Prevent zero-range error in matplotlib.

    box_aspect, x_bounds, y_bounds = _get_aspect_and_bounds(p)
    triangles = _get_valid_triangulation(p, nom)                                                        # Attempt to compute/load valid triangulation.

    def draw_snapshot(fig_obj, ax1, k, angle_view=True):
        _render_surface(ax1, p, u_ap[:, k], triangles, cm.coolwarm, min_val, max_val)
        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)', 'Approximation') # Format axes for approximation.
        if not angle_view:
            _create_colorbar_1(fig_obj, min_val, max_val)                                               # Attach single static colorbar.

    for k in np.arange(0, t + 1, step):
        if k >= t: k = t - 1
        tin = float(T[k])
        nok = nom + '_' + str(format(T[k], '.2f'))

        fig, ax1 = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        fig.suptitle('Approximation at t = %1.3f s (Perspective)' % tin)                                # Set main title for the figure.
        draw_snapshot(fig, ax1, k, angle_view=True)                                                     # Render the current snapshot.
        plt.savefig(nok + 's.png', bbox_inches='tight')                                                 # Save figure to disk.
        plt.close(fig)                                                                                  # Close figure to release memory.

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(7, 6))
        fig_top.suptitle('Approximation at t = %1.3f s (Top View)' % tin)                               # Set main title for the top view figure.
        draw_snapshot(fig_top, ax1_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches='tight')                                             # Save figure to disk.
        plt.close(fig_top)                                                                              # Close figure to release memory.