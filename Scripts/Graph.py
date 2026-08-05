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
        April, 2026.
"""

## Library importation.
import numpy as np                                                                                      # Array utilities (min/max, linspace, arange).
import matplotlib.pyplot as plt                                                                         # Plotting interface.
from matplotlib import cm                                                                               # Colormaps.
from matplotlib.animation import FuncAnimation                                                          # Animation helper.

def _get_valid_triangulation(p, nom):
    """
    Compute a valid triangulation for the point cloud p by filtering out Delaunay
    triangles whose baricenters lie outside the water domain (defined in *_contours.csv).
    It processes regions (islands/disconnected bodies) independently to prevent triangles
    from crossing between them, and uses a local cache file to optimize performance.
    """
    import os
    import pandas as pd
    from scipy.spatial import Delaunay
    import matplotlib.path as mpath

    if not nom:
        return None

    # Parse dataset details from nom path
    # nom path is typically: .../Results/Equation/LakeName/Scale/Variant/Solution
    # So its parent is Variant, grand-parent is Scale, great-grand-parent is LakeName.
    parts = nom.split(os.sep)
    dataset = None
    scale = None
    variant = None
    workspace_root = None

    if len(parts) >= 5:
        variant = parts[-2]
        scale = parts[-3]
        dataset = parts[-4]

        # Climb to find workspace root containing Data/ directory
        temp = nom
        for _ in range(10):
            temp = os.path.dirname(temp)
            if os.path.exists(os.path.join(temp, 'Data')):
                workspace_root = temp
                break

    # Determine triangulation cache path
    cache_path = None
    cloud_path = None
    if workspace_root and dataset and scale and variant:
        cloud_file = f"{dataset}_{variant}.csv"
        cloud_path = os.path.join(workspace_root, 'Data', dataset, scale, cloud_file)
        if not os.path.exists(cloud_path):
            # Fallback search in lake scale folder
            lake_scale_dir = os.path.join(workspace_root, 'Data', dataset, scale)
            if os.path.exists(lake_scale_dir):
                for f in os.listdir(lake_scale_dir):
                    if f.endswith('.csv') and not f.endswith('neighbors.csv') and not f.endswith('triangulation.csv'):
                        cloud_path = os.path.join(lake_scale_dir, f)
                        break
        
        if cloud_path and os.path.exists(cloud_path):
            cache_path = cloud_path.replace('.csv', '_triangulation.csv')

    # Fallback to local Results cache path if Data directory structure is not resolved
    if not cache_path:
        out_dir = os.path.dirname(nom)
        if not out_dir:
            out_dir = '.'
        cache_path = os.path.join(out_dir, 'triangulation.csv')

    # 1. Try to load from cache
    if os.path.exists(cache_path):
        try:
            triangles = np.loadtxt(cache_path, dtype=np.int32, delimiter=',')
            if triangles.ndim == 1:
                triangles = triangles.reshape(1, -1)
            return triangles
        except Exception:
            pass

    if not cloud_path or not os.path.exists(cloud_path):
        return None

    try:
        data = pd.read_csv(cloud_path)
        x = data['x'].values
        y = data['y'].values
        if 'region' in data.columns:
            regions = data['region'].values
        else:
            regions = np.ones_like(x, dtype=np.int32)
    except Exception:
        return None

    # Load contours for parity filtering
    contours_file = f"{dataset}_contours.csv"
    contours_path = os.path.join(workspace_root, 'Data', dataset, contours_file)
    if not os.path.exists(contours_path):
        lake_dir = os.path.join(workspace_root, 'Data', dataset)
        if os.path.exists(lake_dir):
            for f in os.listdir(lake_dir):
                if f.endswith('contours.csv'):
                    contours_path = os.path.join(lake_dir, f)
                    break

    paths = []
    if os.path.exists(contours_path):
        try:
            df_contours = pd.read_csv(contours_path)
            for _, group in df_contours.groupby('region'):
                verts = group[['x', 'y']].values
                if not np.allclose(verts[0], verts[-1]):
                    verts = np.vstack([verts, verts[0]])
                paths.append(mpath.Path(verts))
        except Exception:
            pass

    # 3. Triangulate each region independently
    global_triangles = []
    unique_regions = np.unique(regions)

    for reg_id in unique_regions:
        idx_reg = np.where(regions == reg_id)[0]
        if len(idx_reg) < 3:
            continue
        xy_reg = np.column_stack([x[idx_reg], y[idx_reg]])
        
        try:
            tri = Delaunay(xy_reg)
            t_local = tri.simplices
            # Map local region-based indices to global cloud indices
            t_global = idx_reg[t_local]
            
            # Apply parity filter if contours are available
            if len(paths) > 0:
                for t in t_global:
                    pts = np.column_stack([x[t], y[t]])
                    baricenter = np.mean(pts, axis=0)
                    contain_count = sum(path.contains_point(baricenter) for path in paths)
                    if (contain_count % 2) == 1:
                        global_triangles.append(t)
            else:
                # Fallback to unfiltered Delaunay for this region
                global_triangles.extend(t_global)
        except Exception:
            continue

    if len(global_triangles) > 0:
        global_triangles = np.array(global_triangles, dtype=np.int32)
        try:
            np.savetxt(cache_path, global_triangles, delimiter=',', fmt='%d')
        except Exception:
            pass
        return global_triangles

    return None


def Cloud_Stationary(p, u_ap, u_ex, save = False, nom = ''):
    """
    Cloud_Stationary

    Plot and optionally save the approximated and theoretical stationary solutions side by side,
    along with the absolute error and log10 of absolute error plots in a symmetric 2x2 grid.

    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x 1           ndarray         Array with the computed solution.
        u_ex        m x 1           ndarray         Array with the theoretical solution.
        save                        bool            Save the graphic.
                                                        True: Save the created graphs.
                                                        False: Don't save the created graphs (Default).
        nom                         string          Output file prefix used when save=True (without extension).
        
    Output:
        None
    """

    ## Variable initialization.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    # Calculate absolute error and its base-10 logarithm safely
    err = np.abs(u_ap - u_ex)
    min_err = 0.0
    max_err = err.max()
    if max_err == 0:
        max_err = 1.0

    log_err = np.log10(err + 1e-15)
    min_log = log_err.min()
    max_log = log_err.max()

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot on given axes
    def draw_plot(fig_obj, ax1, ax2, ax3, ax4, angle_view=True):
        if triangles is not None:
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err[:], triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)
            s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err[:], triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:], c = u_ex[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err[:], c = err[:], cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)
            s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err[:], c = log_err[:], cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)

        for ax in (ax1, ax2, ax3, ax4):
            ax.set_box_aspect(box_aspect)
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.xaxis.pane.set_edgecolor('w')
            ax.yaxis.pane.set_edgecolor('w')
            ax.zaxis.pane.set_edgecolor('w')
            ax.grid(True, linestyle = '--', alpha = 0.5)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            # Use tight camera zoom to prevent wasted margins
            if angle_view:
                ax.dist = 6.2
                ax.view_init(elev = 20, azim = -45)
            else:
                ax.dist = 5.2
                ax.set_zticks([])
                ax.view_init(elev = 90, azim = 270)

        ax1.set_zlim([min_val, max_val])
        ax2.set_zlim([min_val, max_val])
        ax3.set_zlim([min_err, max_err])
        ax4.set_zlim([min_log, max_log])

        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax2.set_zlabel('U(x, y)')
            ax3.set_zlabel('Error')
            ax4.set_zlabel('Log10(Error)')

        ax1.set_title('Approximation', y = 0.94)
        ax2.set_title('Theoretical Solution', y = 0.94)
        ax3.set_title('Absolute Error', y = 0.94)
        ax4.set_title('Log10(Absolute Error)', y = 0.94)
        
        # Maximize grid layout area and reduce margins based on view angle
        if angle_view:
            # Take full width since there are no colorbars (Z-axis scale is already visible on the plots)
            fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
        else:
            # Leave space on the right of each subplot for colorbars in the 2D top view
            fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            
            # Position colorbars at their respective subplots in Top View (right side of each plot)
            fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)

    import matplotlib.gridspec as gridspec

    if save:
        # Save perspective view
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)
        fig.suptitle('Stationary Solution Comparison (Perspective)', y = 0.98, fontsize=14, fontweight='bold')
        plt.savefig(nom + '.png', bbox_inches = 'tight')
        plt.savefig(nom + '.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig)

        # Save top view
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y = 0.98, fontsize=14, fontweight='bold')
        plt.savefig(nom + '_top.png', bbox_inches = 'tight')
        plt.savefig(nom + '_top.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig_top)
    else:
        # Show both windows interactively
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        draw_plot(fig, ax1, ax2, ax3, ax4, angle_view=True)
        fig.suptitle('Stationary Solution Comparison (Perspective)', y = 0.98, fontsize=14, fontweight='bold')

        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        draw_plot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, angle_view=False)
        fig_top.suptitle('Stationary Solution Comparison (Top View)', y = 0.98, fontsize=14, fontweight='bold')

        plt.show()

def Cloud_Stationary_1(p, u_ap, save = False, nom = ''):
    """
    Cloud_Stationary_1

    Plot and optionally save the approximated stationary solution (single plot).

    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x 1           ndarray         Array with the computed solution.
        save                        bool            Save the graphic.
                                                        True: Save the created graphs.
                                                        False: Don't save the created graphs (Default).
        nom                         string          Output file prefix used when save=True (without extension).
        
    Output:
        None
    """

    ## Variable initialization.
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot on a given axis
    def draw_plot(fig_obj, ax1, angle_view=True):
        if triangles is not None:
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)

        ax1.set_zlim([min_val, max_val])
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
        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax1.view_init(elev = 20, azim = -45)
        else:
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)

        ax1.set_title('Approximation')
        
        # Colorbar
        fig_obj.subplots_adjust(right = 0.82)
        cbar_ax = fig_obj.add_axes([0.85, 0.15, 0.04, 0.7])
        fig_obj.colorbar(s1, cax = cbar_ax, label = 'Solution Amplitude')

    if save:
        # Save perspective view
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig, ax1, angle_view=True)
        fig.suptitle('Stationary Approximation (Perspective)')
        plt.savefig(nom + '.png', bbox_inches = 'tight')
        plt.savefig(nom + '.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig)

        # Save top view
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)
        fig_top.suptitle('Stationary Approximation (Top View)')
        plt.savefig(nom + '_top.png', bbox_inches = 'tight')
        plt.savefig(nom + '_top.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig_top)
    else:
        # Show both windows interactively
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig, ax1, angle_view=True)
        fig.suptitle('Stationary Approximation (Perspective)')

        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        draw_plot(fig_top, ax1_t, angle_view=False)
        fig_top.suptitle('Stationary Approximation (Top View)')

        plt.show()


def Cloud_Transient(p, u_ap, u_ex, save = False, nom = ''):
    """
    Cloud_Transient

    Plot and optionally save the approximated and theoretical transient solutions over time,
    along with error mapping in a symmetric 2x2 grid.
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 50)                                                                           # Frame stride for plotting/animation.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for titles).
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    # Compute absolute error globally across all time steps for consistent color scaling
    global_err = np.abs(u_ap - u_ex)
    min_err = 0.0
    max_err = global_err.max()
    if max_err == 0:
        max_err = 1.0

    global_log_err = np.log10(global_err + 1e-15)
    min_log = global_log_err.min()
    max_log = global_log_err.max()

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot a single frame on given axes
    def draw_frame(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])
        log_err_k = np.log10(err_k + 1e-15)
        
        if triangles is not None:
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err_k, triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)
            s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err_k, triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err_k, c = err_k, cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)
            s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err_k, c = log_err_k, cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)

        for ax in (ax1, ax2, ax3, ax4):
            ax.set_box_aspect(box_aspect)
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.xaxis.pane.set_edgecolor('w')
            ax.yaxis.pane.set_edgecolor('w')
            ax.zaxis.pane.set_edgecolor('w')
            ax.grid(True, linestyle = '--', alpha = 0.5)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            if angle_view:
                ax.dist = 6.2
                ax.view_init(elev = 20, azim = -45)
            else:
                ax.dist = 5.2
                ax.set_zticks([])
                ax.view_init(elev = 90, azim = 270)

        ax1.set_zlim([min_val, max_val])
        ax2.set_zlim([min_val, max_val])
        ax3.set_zlim([min_err, max_err])
        ax4.set_zlim([min_log, max_log])

        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax2.set_zlabel('U(x, y)')
            ax3.set_zlabel('Error')
            ax4.set_zlabel('Log10(Error)')

        ax1.set_title('Approximation', y = 0.94)
        ax2.set_title('Theoretical Solution', y = 0.94)
        ax3.set_title('Absolute Error', y = 0.94)
        ax4.set_title('Log10(Absolute Error)', y = 0.94)

        if angle_view:
            fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
        else:
            fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            
            fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)

    import matplotlib.gridspec as gridspec

    if save:
        # Save perspective animation
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')

        def update_plot_perspective(k):
            k = min(k, t - 1)
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
            draw_frame(fig, ax1, ax2, ax3, ax4, k, angle_view=True)
            return fig,
        ani_p = FuncAnimation(fig, update_plot_perspective, frames=np.arange(0, t, step), interval=100)
        
        if nom.endswith('.gif'):
            ani_p.save(nom, writer = 'pillow', fps = 10)
        else:
            try:
                ani_p.save(nom, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")
                nom_gif = nom.replace('.mp4', '.gif') if nom.endswith('.mp4') else nom + '.gif'
                ani_p.save(nom_gif, writer = 'pillow', fps = 10)
        plt.close(fig)

        # Save top view animation
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')

        def update_plot_top(k):
            k = min(k, t - 1)
            ax1_t.clear()
            ax2_t.clear()
            ax3_t.clear()
            ax4_t.clear()
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
        else:
            try:
                ani_t.save(nom_top, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                nom_gif_top = nom_top.replace('.mp4', '.gif') if nom_top.endswith('.mp4') else nom_top + '.gif'
                ani_t.save(nom_gif_top, writer = 'pillow', fps = 10)
        plt.close(fig_top)

    else:
        # Show both windows interactively and animate them synchronously
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
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
            
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
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
        plt.close(fig)
        plt.close(fig_top)


def Cloud_Transient_Steps(p, u_ap, u_ex, nom):
    """
    Cloud_Steps

    Save side-by-side snapshots (approximation vs. theoretical) at a few time levels.
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 3)                                                                            # Target about 3 snapshots.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference).
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for filenames).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    # Compute absolute error globally across all time steps for consistent color scaling
    global_err = np.abs(u_ap - u_ex)
    min_err = 0.0
    max_err = global_err.max()
    if max_err == 0:
        max_err = 1.0

    global_log_err = np.log10(global_err + 1e-15)
    min_log = global_log_err.min()
    max_log = global_log_err.max()

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot
    def draw_snapshot(fig_obj, ax1, ax2, ax3, ax4, k, angle_view=True):
        err_k = np.abs(u_ap[:, k] - u_ex[:, k])
        log_err_k = np.log10(err_k + 1e-15)
        
        if triangles is not None:
            s1 = ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s2 = ax2.plot_trisurf(p[:, 0], p[:, 1], u_ex[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
            s3 = ax3.plot_trisurf(p[:, 0], p[:, 1], err_k, triangles = triangles, cmap = cm.viridis, vmin = min_err, vmax = max_err, edgecolors = 'none', linewidth = 0, antialiased = False)
            s4 = ax4.plot_trisurf(p[:, 0], p[:, 1], log_err_k, triangles = triangles, cmap = cm.plasma, vmin = min_log, vmax = max_log, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            s1 = ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s2 = ax2.scatter(p[:, 0], p[:, 1], zs = u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
            s3 = ax3.scatter(p[:, 0], p[:, 1], zs = err_k, c = err_k, cmap = cm.viridis, s = 1, vmin = min_err, vmax = max_err)
            s4 = ax4.scatter(p[:, 0], p[:, 1], zs = log_err_k, c = log_err_k, cmap = cm.plasma, s = 1, vmin = min_log, vmax = max_log)

        for ax in (ax1, ax2, ax3, ax4):
            ax.set_box_aspect(box_aspect)
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.xaxis.pane.set_edgecolor('w')
            ax.yaxis.pane.set_edgecolor('w')
            ax.zaxis.pane.set_edgecolor('w')
            ax.grid(True, linestyle = '--', alpha = 0.5)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            if angle_view:
                ax.dist = 6.2
                ax.view_init(elev = 20, azim = -45)
            else:
                ax.dist = 5.2
                ax.set_zticks([])
                ax.view_init(elev = 90, azim = 270)

        ax1.set_zlim([min_val, max_val])
        ax2.set_zlim([min_val, max_val])
        ax3.set_zlim([min_err, max_err])
        ax4.set_zlim([min_log, max_log])

        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax2.set_zlabel('U(x, y)')
            ax3.set_zlabel('Error')
            ax4.set_zlabel('Log10(Error)')

        ax1.set_title('Approximation', y = 0.94)
        ax2.set_title('Theoretical Solution', y = 0.94)
        ax3.set_title('Absolute Error', y = 0.94)
        ax4.set_title('Log10(Absolute Error)', y = 0.94)
        
        if angle_view:
            fig_obj.subplots_adjust(left = 0.02, right = 0.98, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
        else:
            fig_obj.subplots_adjust(left = 0.02, right = 0.94, bottom = 0.02, top = 0.94, hspace = 0.08, wspace = 0.08)
            
            fig_obj.colorbar(s2, ax = ax2, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s3, ax = ax3, shrink = 0.65, aspect = 15, pad = 0.04)
            fig_obj.colorbar(s4, ax = ax4, shrink = 0.65, aspect = 15, pad = 0.04)

    import matplotlib.gridspec as gridspec

    ## Create the graphs.
    for k in np.arange(0, t + 1, step):
        if k >= t:
            k = t - 1
        tin = float(T[k])
        nok = nom + '_' + str(format(T[k], '.2f'))

        # Perspective snapshot
        fig = plt.figure(figsize = (11, 9))
        gs = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1 = fig.add_subplot(gs[0, 0], projection = '3d')
        ax2 = fig.add_subplot(gs[0, 1], projection = '3d')
        ax3 = fig.add_subplot(gs[1, 0], projection = '3d')
        ax4 = fig.add_subplot(gs[1, 1], projection = '3d')
        
        plt.suptitle('Solution at t = %1.3f s (Perspective)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        draw_snapshot(fig, ax1, ax2, ax3, ax4, k, angle_view=True)
        plt.savefig(nok + 's.png', bbox_inches = 'tight')
        plt.savefig(nok + 's.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig)

        # Top view snapshot
        fig_top = plt.figure(figsize = (13, 8))
        gs_t = gridspec.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
        ax1_t = fig_top.add_subplot(gs_t[0, 0], projection = '3d')
        ax2_t = fig_top.add_subplot(gs_t[0, 1], projection = '3d')
        ax3_t = fig_top.add_subplot(gs_t[1, 0], projection = '3d')
        ax4_t = fig_top.add_subplot(gs_t[1, 1], projection = '3d')
        
        plt.suptitle('Solution at t = %1.3f s (Top View)' % tin, y = 0.98, fontsize=14, fontweight='bold')
        draw_snapshot(fig_top, ax1_t, ax2_t, ax3_t, ax4_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches = 'tight')
        plt.savefig(nok + 's_top.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig_top)

def Cloud_Transient_1(p, u_ap, save = False, nom = ''):
    """
    Cloud_1

    Plot and optionally save the transient approximated solution over time (single plot).

    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x t           ndarray         Array with the computed solution.
        save                        bool            Save the graphic.
                                                        True: Save the created graphs.
                                                        False: Don't save the created graphs (Default).
        nom                         string          Output filename for the animation when save=True.
        
    Output:
        None
    """
    ## Variable initialization.
    t       = u_ap.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 50)                                                                           # Frame stride for plotting/animation.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for titles).
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot a single frame
    def draw_frame(ax1, k, angle_view=True):
        if triangles is not None:
            ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)

        ax1.set_zlim([min_val, max_val])
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
        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax1.view_init(elev = 20, azim = -45)
        else:
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)
        ax1.set_title('Approximation')

    if save:                                                                                            # Save an animation to disk.
        # Perspective animation
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        def update_plot_perspective(k):
            k = min(k, t - 1)
            ax1.clear()
            tin = float(T[k])
            fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
            draw_frame(ax1, k, angle_view=True)
            return fig,
        ani_p = FuncAnimation(fig, update_plot_perspective, frames=np.arange(0, t, step), interval=100)
        
        if nom.endswith('.gif'):
            ani_p.save(nom, writer = 'pillow', fps = 10)
        else:
            try:
                ani_p.save(nom, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")
                nom_gif = nom.replace('.mp4', '.gif') if nom.endswith('.mp4') else nom + '.gif'
                ani_p.save(nom_gif, writer = 'pillow', fps = 10)
        plt.close(fig)

        # Top view animation
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        def update_plot_top(k):
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
        else:
            try:
                ani_t.save(nom_top, writer = 'ffmpeg', fps = 10)
            except Exception as e:
                nom_gif_top = nom_top.replace('.mp4', '.gif') if nom_top.endswith('.mp4') else nom_top + '.gif'
                ani_t.save(nom_gif_top, writer = 'pillow', fps = 10)
        plt.close(fig_top)

    else:
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
            
            ax1.clear()
            ax1_t.clear()

        # Final step
        tin = float(T[-1])
        fig.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
        fig_top.suptitle('Solution at t = %1.3f s (Top View)' % tin)
        
        draw_frame(ax1, t - 1, angle_view=True)
        draw_frame(ax1_t, t - 1, angle_view=False)
        
        plt.pause(0.1)
        plt.close(fig)
        plt.close(fig_top)


def Cloud_Transient_Steps_1(p, u_ap, nom = ''):
    """
    Cloud_Steps_1

    Save approximation-only snapshots at a few time levels.
    
    This helper generates approximately three snapshots across the time interval, saving each one as
    both .png and .svg.
    
    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x t           ndarray         Array with the computed solution.
        nom                         string          Output file prefix used for saved snapshots.
    
    Output:
        None
    """

    ## Variable initialization.
    t       = u_ap.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 3)                                                                            # Target about 3 snapshots.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for filenames).
    min_val = u_ap.min()                                                                                # Color scale minimum (from approximation).
    max_val = u_ap.max()                                                                                # Color scale maximum (from approximation).

    # Calculate physical aspect ratio to keep geographic boundaries unaltered (axis equal mapping)
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0
    box_aspect = (x_range, y_range, 0.4 * max_range)

    triangles = _get_valid_triangulation(p, nom)

    # Helper function to plot
    def draw_snapshot(fig_obj, ax1, k, angle_view=True):
        if triangles is not None:
            ax1.plot_trisurf(p[:, 0], p[:, 1], u_ap[:, k], triangles = triangles, cmap = cm.coolwarm, vmin = min_val, vmax = max_val, edgecolors = 'none', linewidth = 0, antialiased = False)
        else:
            ax1.scatter(p[:, 0], p[:, 1], zs = u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)

        ax1.set_zlim([min_val, max_val])
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
        if angle_view:
            ax1.set_zlabel('U(x, y)')
            ax1.view_init(elev = 20, azim = -45)
        else:
            ax1.set_zticks([])
            ax1.view_init(elev = 90, azim = 270)
        ax1.set_title('Approximation')

    ## Create the graph.
    for k in np.arange(0, t + 1, step):                                                                 # Iterate selected snapshot indices.
        k = min(k, t - 1)                                                                               # Clamp index to valid range.
        tin = float(T[k])                                                                               # Current normalized time.
        nok = nom + '_' + str(format(T[k], '.2f'))                                                      # Filename tag based on time.
        
        # Perspective snapshot
        fig, ax1 = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        plt.suptitle('Solution at t = %1.3f s (Perspective)' % tin)
        draw_snapshot(fig, ax1, k, angle_view=True)
        plt.savefig(nok + 's.png', bbox_inches = 'tight')
        plt.savefig(nok + 's.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig)

        # Top view snapshot
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (7, 6))
        plt.suptitle('Solution at t = %1.3f s (Top View)' % tin)
        draw_snapshot(fig_top, ax1_t, k, angle_view=False)
        plt.savefig(nok + 's_top.png', bbox_inches = 'tight')
        plt.savefig(nok + 's_top.svg', format = 'svg', bbox_inches = 'tight')
        plt.close(fig_top)
