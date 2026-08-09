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
    plot_stationary          3D scatter: approximation vs. optional theoretical solution.
    plot_transient           3D scatter over time (interactive stepping or animation).
    plot_transient_steps     Save snapshots at ~3 time levels.

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
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Core numerical operations.
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib                                                                                                         # Plotting interface.

from matplotlib import animation                                                                                                        # Animation framework.
from matplotlib.animation import FuncAnimation                                                                                          # Animation helper.
from typing import Callable, Optional, Tuple, List, Any                                                                                 # Type hinting.

from mGFD.core.utils import get_valid_triangulation, get_aspect_and_bounds                                                              # Geometry utilities.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def _setup_3d_axes(ax: Any, angle_view: bool, box_aspect: Tuple[float, float, float], x_bounds: List[float], y_bounds: List[float], z_bounds: List[float], z_label: str) -> None:
    """
    _setup_3d_axes
    Helper function to standardize the formatting and limits of 3D axes in Matplotlib.
    
    Input:
        ax                          Any             The 3D axes object to format.
        angle_view                  bool            If True, sets a perspective view; if False, sets a top-down view.
        box_aspect                  Tuple           Physical aspect ratio scaling for the axes box.
        x_bounds                    List            [x_min, x_max] physical bounds.
        y_bounds                    List            [y_min, y_max] physical bounds.
        z_bounds                    List            [z_min, z_max] physical bounds.
        z_label                     str             Label for the Z-axis.
    """
    # 1. Appearance and labels
    ax.set_box_aspect(box_aspect)                                                                                                       # Apply physical aspect ratio to axes.
    ax.xaxis.pane.fill               = False                                                                                            # Make X pane transparent.
    ax.yaxis.pane.fill               = False                                                                                            # Make Y pane transparent.
    ax.zaxis.pane.fill               = False                                                                                            # Make Z pane transparent.
    ax.xaxis.pane.set_edgecolor('w')                                                                                                    # Set X pane edge color to white.
    ax.yaxis.pane.set_edgecolor('w')                                                                                                    # Set Y pane edge color to white.
    ax.zaxis.pane.set_edgecolor('w')                                                                                                    # Set Z pane edge color to white.
    ax.grid(True, linestyle=':', alpha=0.6, color='gray')                                                                               # Enable grid lines with transparency.
    ax.set_xlabel('X', labelpad=10)                                                                                                     # Set X-axis label with padding.
    ax.set_ylabel('Y', labelpad=10)                                                                                                     # Set Y-axis label with padding.
    
    # 2. Bounds
    ax.set_xlim(x_bounds)                                                                                                               # Set X-axis limits.
    ax.set_ylim(y_bounds)                                                                                                               # Set Y-axis limits.
    ax.set_zlim(z_bounds)                                                                                                               # Set Z-axis limits.
    ax.autoscale(False)                                                                                                                 # Disable automatic scaling.

    # 3. Camera angle
    if angle_view:
        ax.dist = 6.2                                                                                                                   # Adjust camera zoom for perspective.
        ax.view_init(elev=20, azim=-45)                                                                                                 # Set perspective viewing angle.
        ax.set_zlabel(z_label, labelpad=10)                                                                                             # Set Z-axis label with padding.
    else:                                                                                                                               # Interactive mode.
        ax.dist = 5.2                                                                                                                   # Adjust camera zoom for top view.
        ax.set_zticks([])                                                                                                               # Remove Z-axis ticks for top view.
        ax.view_init(elev=90, azim=270)                                                                                                 # Set top-down viewing angle.


def _render_surface(ax: Any, p: np.ndarray, data: np.ndarray, triangles: Optional[np.ndarray], cmap: Any, vmin: float, vmax: float) -> Any:
    """
    _render_surface
    Helper function to render a 3D surface plot using a valid Delaunay triangulation, 
    or fallback to a 3D scatter plot if the triangulation is unavailable.
    
    Input:
        ax                          Any             The 3D axes object to draw on.
        p           m x 2           ndarray         Point cloud [x, y] coordinates.
        data        m               ndarray         Z-values (e.g., solution amplitude or error) at each point.
        triangles                   ndarray         Valid Delaunay simplices (or None for fallback).
        cmap                        Any             Matplotlib colormap.
        vmin                        float           Minimum value for color scaling.
        vmax                        float           Maximum value for color scaling.
        
    Output:
        artist                      Any             The created Matplotlib Poly3DCollection or PathCollection.
    """
    # 1. Rendering
    if triangles is not None:                                                                                                           # Check if triangulation is available.
        return ax.plot_trisurf(p[:, 0], p[:, 1], data, triangles=triangles, cmap=cmap,
                                vmin=vmin, vmax=vmax, edgecolors='none',
                                linewidth=0, antialiased=False)                                                                         # Plot surface using valid triangulation.
    else:                                                                                                                               # Interactive mode.
        return ax.scatter(p[:, 0], p[:, 1], zs=data, c=data, cmap=cmap, s=1, vmin=vmin, vmax=vmax)                                      # Fallback to plotting scattered points.


def _prepare_plot_data(p: np.ndarray, u: np.ndarray, nom: str) -> Tuple[float, float, Tuple[float, float, float], List[float], List[float], Optional[np.ndarray]]:
    """
    _prepare_plot_data
    Helper to compute common plotting bounds, color scales, and triangulation.
    
    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u           m               ndarray         Array or Matrix with the solution evaluated at all nodes.
        nom                         str             Base filename for the saved output (used to load cached triangulation).
        
    Output:
        min_val                     float           Minimum value of the solution array.
        max_val                     float           Maximum value of the solution array.
        box_aspect                  Tuple           Physical aspect ratio scaling for the axes box.
        x_bounds                    List            [x_min, x_max] physical bounds.
        y_bounds                    List            [y_min, y_max] physical bounds.
        triangles                   ndarray         Valid Delaunay simplices (or None for fallback).
    """
    # 1. Calculate extrema
    min_val, max_val               = float(u.min()), float(u.max())                                                                     # Shared color scale from the data.
    
    if min_val == max_val: max_val += 1e-6                                                                                              # Prevent zero-range error in matplotlib.
    
    # 2. Fetch triangulation and bounds
    box_aspect, x_bounds, y_bounds = get_aspect_and_bounds(p)                                                                           # Compute bounds and aspect ratio.
    triangles                      = get_valid_triangulation(p, nom)                                                                    # Attempt to compute/load valid triangulation.
    
    return min_val, max_val, box_aspect, x_bounds, y_bounds, triangles                                                                  # Return preparation components.


def _generate_static_views(render_func: Callable, title: str, save_path: Optional[str] = None, show: bool = False, verbose: bool = True) -> None:
    """
    _generate_static_views
    Helper to generate both perspective and top views for a static snapshot.
    
    Input:
        render_func                 Callable        Callback function to render the specific snapshot on a given axis.
        title                       str             Title for the figure.
        save_path                   str             Base filename for the saved output (None if interactive).
        show                        bool            If True, calls plt.show() after rendering.
        verbose                     bool            If True, prints confirmation of saved figures.
    """
    # 1. View definition
    views = [(True, 'Perspective', ''), (False, 'Top View', '_top')]                                                                    # Define angle, name, and suffix for both views.
    
    # 2. View iteration
    for angle_view, view_name, suffix in views:                                                                                         # Iterate through views.
        fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(10, 8))                                                  # Create figure and 3D axis with larger size.
        fig.suptitle(f'{title} ({view_name})', fontsize=16, fontweight='bold', y=0.95)                                                  # Set main title for the figure centered and bold.
        render_func(fig, ax, angle_view)                                                                                                # Execute the rendering callback.
        
        # 3. Save handling
        if save_path:                                                                                                                   # Check if a save path was provided.
            plt.savefig(f'{save_path}{suffix}.png', bbox_inches='tight')                                                                # Save figure to disk.
            
            if verbose:                                                                                                                 # Check verbosity flag.
                logger.info(f'\tSaved figure to {save_path}{suffix}.png')                                                               # Print confirmation of save.
            
            plt.close(fig)                                                                                                              # Close figure to release memory.
            
    # 4. Interactive display
    if show:                                                                                                                            # Check if figures should be shown.
        plt.show()                                                                                                                      # Display the interactive plot window.


def plot_stationary(p: np.ndarray, u: np.ndarray, save: bool = False, nom: str = '', title: str = 'Solution', verbose: bool = True) -> None:
    """
    plot_stationary
    Render a single 3D scatter plot of the solution for a stationary PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u           m               ndarray         Array with the solution evaluated at all nodes.
        save                        bool            If True, saves the figure to disk instead of displaying.
        nom                         str             Base filename for the saved output (if save=True).
        title                       str             Title for the figure and Z-axis label.
        verbose                     bool            If True, prints progress and errors to console.
    """
    # 1. Plot data preparation
    min_val, max_val, box_aspect, x_bounds, y_bounds, triangles = _prepare_plot_data(p, u, nom)

    # 2. Helper function definition
    def draw_plot(fig_obj: Any, ax1: Any, angle_view: bool = True) -> None:
        """
        draw_plot
        Helper function to render a single stationary plot.
        
        Input:
            fig_obj                 Any             Matplotlib figure object.
            ax1                     Any             The 3D axes object to draw on.
            angle_view              bool            If True, sets a perspective view; if False, sets a top-down view.
        """
        surf = _render_surface(ax1, p, u, triangles, matplotlib.colormaps['coolwarm'], min_val, max_val)                                # Render solution surface.
        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)')                                  # Format axes for the solution.
        
        if not angle_view and not hasattr(fig_obj, 'colorbar_added'):                                                                   # Check if a colorbar is needed.
            fig_obj.colorbar(surf, ax=ax1, fraction=0.046, pad=0.04)                                                                    # Add colorbar to the plot.
            fig_obj.colorbar_added = True                                                                                               # Mark colorbar as added.

    # 3. Render execution
    if save:                                                                                                                            # Check if save flag is enabled.
        _generate_static_views(draw_plot, title, save_path=nom, verbose=verbose)                                                        # Generate and save static views.
    else:                                                                                                                               # Interactive mode.
        _generate_static_views(draw_plot, title, show=True, verbose=verbose)                                                            # Generate and show interactive views.


def plot_transient(p: np.ndarray, u: np.ndarray, save: bool = False, nom: str = '', title: str = 'Solution', verbose: bool = True) -> None:
    """
    plot_transient
    Render a single 3D scatter plot of the solution for a transient PDE problem,
    and optionally save it as an animation.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u           m x t           ndarray         Matrix with the solution evaluated at all nodes at each time step.
        save                        bool            If True, saves the animation to disk instead of displaying interactively.
        nom                         str             Base filename for the saved output (if save=True).
        title                       str             Title for the figure and Z-axis label.
        verbose                     bool            If True, prints progress and errors to console.
    """
    # 1. Variable initialization
    t                                                           = int(u.shape[1])                                                       # Get total number of time steps.
    step                                                        = max(1, t // 50)                                                       # Calculate step size to limit frame count.
    T                                                           = np.linspace(0, 1, t)                                                  # Generate array of time values.
    min_val, max_val, box_aspect, x_bounds, y_bounds, triangles = _prepare_plot_data(p, u, nom)

    # 2. Helper function definition
    def draw_frame(fig_obj: Any, ax1: Any, k: int, angle_view: bool = True) -> None:
        """
        draw_frame
        Helper function to render a single animation frame.
        
        Input:
            fig_obj                 Any             Matplotlib figure object.
            ax1                     Any             The 3D axes object to draw on.
            k                       int             Current time step index.
            angle_view              bool            If True, sets a perspective view; if False, sets a top-down view.
        """
        if not hasattr(fig_obj, 'surf_artists'):
            s1 = _render_surface(ax1, p, u[:, k], triangles, matplotlib.colormaps['coolwarm'], min_val, max_val)                        # Render solution surface.
            fig_obj.surf_artists = {'s1': s1}                                                                                           # Cache surface artists for animation updates.
            _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)')                              # Format axes for the solution.
            
            if not angle_view and not hasattr(fig_obj, 'colorbar_added'):
                fig_obj.colorbar(s1, ax=ax1, fraction=0.046, pad=0.04)                                                                  # Add colorbar to the plot.
                fig_obj.colorbar_added = True                                                                                           # Mark colorbar as added.
        else:
            if triangles is not None:                                                                                                   # Check if triangulation is available.
                artist = fig_obj.surf_artists['s1']                                                                                     # Retrieve cached surface artist.
                z_data = u[:, k]                                                                                                        # Extract Z data for current step.
                verts  = np.stack((p[triangles, 0], p[triangles, 1], z_data[triangles]), axis=-1)                                       # Compute new 3D vertices for the surface.
                artist.set_verts(verts)                                                                                                 # Update surface geometry in animation.
                artist.set_array(np.mean(z_data[triangles], axis=1))                                                                    # Update surface colors based on new Z heights.
            else:
                artist            = fig_obj.surf_artists['s1']                                                                          # Retrieve cached surface artist.
                z_data            = u[:, k]                                                                                             # Extract Z data for current step.
                artist._offsets3d = (p[:, 0], p[:, 1], z_data)                                                                          # Update scatter coordinates directly.
                artist.set_array(z_data)                                                                                                # Update scatter colors.

    # 3. Render execution
    if save:                                                                                                                            # Check if save flag is enabled.
        def save_animation(angle_view: bool, view_name: str, suffix: str, verbose: bool = verbose) -> None:
            """
            save_animation
            Helper function to initialize and save an animation for a specific view.
            
            Input:
                angle_view          bool            If True, sets a perspective view; if False, sets a top-down view.
                view_name           str             Name of the view for the title (e.g. 'Perspective').
                suffix              str             Filename suffix for saving (e.g. '_top').
                verbose             bool            If True, prints confirmation of saved animations.
            """
            fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(10, 8))                                              # Create figure and 3D axis for animation.
            
            def update(frame: int) -> None:
                tin = float(T[frame])                                                                                                   # Get the current time step value.
                fig.suptitle(f'{title} at t = {tin:1.3f} s ({view_name})', fontsize=16, fontweight='bold', y=0.95)                      # Set main title for the figure.
                draw_frame(fig, ax, frame, angle_view=angle_view)                                                                       # Render the current frame.

            ani     = FuncAnimation(fig, update, frames=np.arange(0, t, step), blit=False)  # type: ignore                              # Initialize animation object.
            out_nom = nom                                                                                                               # Set initial output filename.
            
            if suffix:                                                                                                                  # Check if suffix is provided.
                out_nom = nom.replace('.mp4', f'{suffix}.mp4').replace('.gif', f'{suffix}.gif')                                         # Append suffix before extension.
                if out_nom == nom: out_nom = nom + suffix                                                                               # Fallback suffix append.
                
            save_nom = out_nom if out_nom.endswith('.gif') else out_nom + '.gif'
            ani.save(save_nom, writer='pillow', fps=10)                                                                                 # Save animation as GIF.
            
            if verbose:
                logger.info(f'\tSaved animation to {save_nom}')                                                                         # Print confirmation of save.
                
            plt.close(fig)                                                                                                              # Close figure to release memory.
            
        save_animation(True, 'Perspective', '', verbose=verbose)                                                                        # Save perspective animation.
        save_animation(False, 'Top View', '_top', verbose=verbose)                                                                      # Save top view animation.

    # 4. Interactive animation
    else:                                                                                                                               # Interactive mode.
        fig, ax1       = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(10, 8))                                           # Create figure for perspective view.
        fig_top, ax1_t = plt.subplots(1, 1, subplot_kw={"projection": "3d"}, figsize=(10, 8))                                           # Create figure for top view.
        
        for k in range(0, t, step):
            tin = float(T[k])                                                                                                           # Get physical time value.
            fig.suptitle(f'{title} at t = {tin:1.3f} s (Perspective)', fontsize=16, fontweight='bold', y=0.95)                          # Set main title for the figure.
            fig_top.suptitle(f'{title} at t = {tin:1.3f} s (Top View)', fontsize=16, fontweight='bold', y=0.95)                         # Set main title for the top view figure.
            draw_frame(fig, ax1, int(k), angle_view=True)                                                                               # Render the current frame.
            draw_frame(fig_top, ax1_t, int(k), angle_view=False)                                                                        # Render top view frame.
            plt.pause(0.01)                                                                                                             # Pause to allow UI to update.
            
        tin = float(T[-1])                                                                                                              # Get final time value.
        fig.suptitle(f'{title} at t = {tin:1.3f} s (Perspective)', fontsize=16, fontweight='bold', y=0.95)                              # Set main title for the figure.
        fig_top.suptitle(f'{title} at t = {tin:1.3f} s (Top View)', fontsize=16, fontweight='bold', y=0.95)                             # Set main title for the top view figure.
        draw_frame(fig, ax1, t - 1, angle_view=True)                                                                                    # Render the current frame.
        draw_frame(fig_top, ax1_t, t - 1, angle_view=False)                                                                             # Render final top view frame.
        plt.pause(0.1)                                                                                                                  # Pause to allow UI to update.
        plt.close(fig)                                                                                                                  # Close figure to release memory.
        plt.close(fig_top)                                                                                                              # Close figure to release memory.


def plot_transient_steps(p: np.ndarray, u: np.ndarray, nom: str, title: str = 'Solution', verbose: bool = True) -> None:
    """
    plot_transient_steps
    Render and save static single 3D scatter plots of the solution at a few key time steps
    during a transient PDE problem.

    Input:
        p           m x 4           ndarray         Point cloud with columns [x, y, flag, region].
        u           m x t           ndarray         Matrix with the solution evaluated at all nodes at each time step.
        nom                         str             Base filename for the saved output.
        title                       str             Title for the figure and Z-axis label.
        verbose                     bool            If True, prints confirmation of saved figures.
    """
    # 1. Variable initialization
    t                                                           = int(u.shape[1])                                                       # Get total number of time steps.
    step                                                        = max(1, t // 3)                                                        # Calculate step size to limit frame count.
    T                                                           = np.linspace(0, 1, t)                                                  # Generate array of time values.
    min_val, max_val, box_aspect, x_bounds, y_bounds, triangles = _prepare_plot_data(p, u, nom)

    # 2. Helper function definition
    def draw_snapshot(fig_obj: Any, ax1: Any, k: int, angle_view: bool = True) -> None:
        """
        draw_snapshot
        Helper function to render a specific snapshot step.
        
        Input:
            fig_obj                 Any             Matplotlib figure object.
            ax1                     Any             The 3D axes object to draw on.
            k                       int             Current time step index.
            angle_view              bool            If True, sets a perspective view; if False, sets a top-down view.
        """
        surf = _render_surface(ax1, p, u[:, k], triangles, matplotlib.colormaps['coolwarm'], min_val, max_val)                          # Render solution surface.
        _setup_3d_axes(ax1, angle_view, box_aspect, x_bounds, y_bounds, [min_val, max_val], 'U(x, y)')                                  # Format axes for the solution.
        if not angle_view and not hasattr(fig_obj, 'colorbar_added'):                                                                   # Check if a colorbar is needed.
            fig_obj.colorbar(surf, ax=ax1, fraction=0.046, pad=0.04)                                                                    # Add colorbar to the plot.
            fig_obj.colorbar_added = True                                                                                               # Mark colorbar as added.

    # 3. Snapshot iteration
    for k in range(0, t + 1, step):                                                                                                     # Iterate through selected snapshots.
        if k >= t: k = int(t - 1)                                                                                                       # Ensure index does not exceed limits.
        tin = float(T[k])
        nok = nom + '_' + str(format(T[k], '.2f'))                                                                                      # Format snapshot filename.

        def render_cb(fig: Any, ax: Any, angle_view: bool) -> None:
            """
            render_cb
            Callback function to execute the snapshot rendering with the fixed time step k.
            
            Input:
                fig                 Any             Matplotlib figure object.
                ax                  Any             The 3D axes object to draw on.
                angle_view          bool            If True, sets a perspective view; if False, sets a top-down view.
            """
            draw_snapshot(fig, ax, int(k), angle_view)                                                                                  # Execute draw_snapshot with the captured step k.
            
        _generate_static_views(render_cb, f'{title} at t = {tin:1.3f} s', save_path=nok + 's', verbose=verbose)                         # Generate and save views.
