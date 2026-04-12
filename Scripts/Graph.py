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

def Cloud_Stationary(p, u_ap, u_ex, save = False, nom = ''):
    """
    Cloud_Stationary

    Plot and optionally save the approximated and theoretical stationary solutions side by side.

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
    
    Notes:
        The color scale (vmin/vmax) is taken from u_ex to make both plots comparable.
    """

    ## Variable initialization.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw = {"projection": "3d"}, figsize = (10, 5))          # Create side-by-side 3D axes.
    
    ## Plot approximated solution.
    ax1.scatter(p[:, 0], p[:, 1], u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # 3D scatter colored by u_ap.
    ax1.set_zlim([min_val, max_val])                                                                    # Keep z-limits consistent across axes.
    ax1.view_init(elev = 9, azim = -50)                                                                 # Choose a consistent camera angle.
    ax1.set_title('Approximation')                                                                      # Label left plot.
    
    ## Plot theoretical/reference solution.
    ax2.scatter(p[:, 0], p[:, 1], u_ex[:], c = u_ex[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # 3D scatter colored by u_ex.
    ax2.set_zlim([min_val, max_val])                                                                    # Keep z-limits consistent across axes.
    ax2.view_init(elev = 9, azim = -50)                                                                 # Choose a consistent camera angle.
    ax2.set_title('Theoretical Solution')                                                               # Label right plot.
    
    fig.suptitle('Stationary Solution Comparison')                                                      # Shared figure title.
    
    if save:                                                                                            # Save images to disk.
        plt.savefig(nom + '.png')                                                                       # Save raster image.
        plt.savefig(nom + '.svg', format = 'svg')                                                       # Save vector image.
    else:                                                                                               # Show interactively.
        plt.show()                                                                                      # Display the figure.
    
    plt.close()                                                                                         # Always close to free memory/resources.

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

    fig, (ax1) = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (5, 5))                # Create a single 3D axis.
    
    ## Plot approximated solution.
    ax1.scatter(p[:, 0], p[:, 1], u_ap[:], c = u_ap[:], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # 3D scatter colored by u_ap.
    ax1.set_zlim([min_val, max_val])                                                                    # Set z-limits to match color scale.
    ax1.view_init(elev = 9, azim = -50)                                                                 # Choose a consistent camera angle.
        
    fig.suptitle('Stationary Approximation')                                                            # Figure title.
    
    if save:                                                                                            # Save images to disk.
        plt.savefig(nom + '.png')                                                                       # Save raster image.
        plt.savefig(nom + '.svg', format = 'svg')                                                       # Save vector image.
    else:                                                                                               # Show interactively.
        plt.show()                                                                                      # Display the figure.
    
    plt.close()                                                                                         # Always close to free memory/resources.


def Cloud_Transient(p, u_ap, u_ex, save = False, nom = ''):
    """
    Cloud_Transient

    Plot and optionally save the approximated and theoretical transient solutions over time.
    
    When save=False, this function iterates over a subset of time steps and updates the plots interactively.
    When save=True, this function builds a Matplotlib animation and attempts to save it using ffmpeg
    (falling back to the pillow writer if ffmpeg is not available).

    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x t           ndarray         Array with the computed solution.
        u_ex        m x t           ndarray         Array with the theoretical solution.
        save                        bool            Save the graphic.
                                                        True: Save the created graphs.
                                                        False: Don't save the created graphs (Default).
        nom                         string          Output filename for the animation when save=True.
        
    Output:
        None
    
    Notes:
        The time grid shown in titles is normalized as T = linspace(0, 1, t), matching the convention used
        by the transient solvers in this repository.
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 50)                                                                           # Frame stride for plotting/animation.
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for titles).
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference solution).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference solution).

    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw = {"projection": "3d"}, figsize = (10, 5))          # Create side-by-side 3D axes.
    
    if save:                                                                                            # Save an animation to disk.
        def update_plot(k):                                                                             # Frame update callback for FuncAnimation.
            k = min(k, t - 1)                                                                           # Clamp frame index.
            ax1.clear()                                                                                 # Clear left axis.
            ax2.clear()                                                                                 # Clear right axis.
            tin = float(T[k])                                                                           # Current normalized time.
            fig.suptitle('Solution at t = %1.3f s.' % tin)                                              # Update figure title.
            
            ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot approximation at frame k.
            ax1.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax1.set_title('Approximation')                                                              # Left title.
            
            ax2.scatter(p[:, 0], p[:, 1], u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot reference at frame k.
            ax2.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax2.set_title('Theoretical Solution')                                                       # Right title.
            
            return fig,                                                                                 # Return artists (blit=False, but keep signature).
    
        ani = FuncAnimation(fig, update_plot, frames = np.arange(0, t + 1, step), blit = False)         # Build animation over selected frames.
        try:                                                                                            # Prefer ffmpeg when available.
            ani.save(nom, writer = 'ffmpeg', fps = 10)                                                  # Save video (e.g., .mp4).
        except Exception as e:                                                                          # Fall back when ffmpeg is not available.
            print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")            # Inform the user about the fallback.
            if nom.endswith('.mp4'):                                                                    # If filename suggests mp4, produce a gif instead.
                nom_gif = nom.replace('.mp4', '.gif')                                                   # Replace extension for pillow writer.
                ani.save(nom_gif, writer = 'pillow', fps = 10)                                          # Save as GIF.
                print(f"Video saved as {nom_gif} instead of {nom}")                                     # Inform about the actual saved file.
            else:                                                                                       # Keep the requested name/extension.
                ani.save(nom, writer = 'pillow', fps = 10)                                              # Save using pillow writer.
        plt.close()                                                                                     # Close figure after saving.

    else:
        for k in np.arange(0, t, step):                                                                 # Iterate frames for interactive preview.
            tin = float(T[k])                                                                           # Current normalized time.
            fig.suptitle('Solution at t = %1.3f s.' % tin)                                              # Update figure title.

            ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot approximation at time k.
            ax1.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax1.set_title('Approximation')                                                              # Left title.
            
            ax2.scatter(p[:, 0], p[:, 1], u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot reference at time k.
            ax2.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax2.set_title('Theoretical Solution')                                                       # Right title.

            plt.pause(0.01)                                                                             # Allow GUI event loop to update.
            ax1.clear()                                                                                 # Clear for the next frame.
            ax2.clear()                                                                                 # Clear for the next frame.

        tin = float(T[-1])                                                                              # Final time for display.
        fig.suptitle('Solution at t = %1.3f s.' % tin)                                                  # Final figure title.

        ax1.scatter(p[:, 0], p[:, 1], u_ap[:, -1], c = u_ap[:, -1], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot final approximation.
        ax1.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax1.set_title('Approximation')                                                                  # Left title.
        
        ax2.scatter(p[:, 0], p[:, 1], u_ex[:, -1], c = u_ex[:, -1], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot final reference.
        ax2.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax2.set_title('Theoretical Solution')                                                           # Right title.

        plt.pause(0.1)                                                                                  # Keep final frame visible briefly.
        plt.close()                                                                                     # Close figure.


def Cloud_Transient_Steps(p, u_ap, u_ex, nom):
    """
    Cloud_Steps

    Save side-by-side snapshots (approximation vs. theoretical) at a few time levels.
    
    This helper generates approximately three snapshots across the time interval, saving each one as
    both .png and .svg.

    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        u_ap        m x t           ndarray         Array with the computed solution.
        u_ex        m x t           ndarray         Array with the theoretical solution.
        nom                         string          Output file prefix used for saved snapshots.
    
    Output:
        None
    """

    ## Variable initialization.
    t       = u_ex.shape[1]                                                                             # Number of time steps.
    step    = max(1, t // 3)                                                                            # Target about 3 snapshots.
    min_val = u_ex.min()                                                                                # Shared color scale minimum (from reference).
    max_val = u_ex.max()                                                                                # Shared color scale maximum (from reference).
    T       = np.linspace(0, 1, t)                                                                      # Normalized time grid (for filenames).

    ## Create the graphs.
    for k in np.arange(0, t + 1, step):                                                                 # Iterate selected snapshot indices.
        if k >= t:                                                                                      # Clamp last index to valid range.
            k = t - 1                                                                                   # Use the final time step.
        fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw = {"projection": "3d"}, figsize = (10, 5))      # Create side-by-side 3D axes.
        tin = float(T[k])                                                                               # Current normalized time.
        plt.suptitle('Solution at t = %1.3f s.' % tin)                                                  # Figure title.
        ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Approximation snapshot.
        ax1.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax1.set_title('Approximation')                                                                  # Left title.
        ax2.scatter(p[:, 0], p[:, 1], u_ex[:, k], c = u_ex[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Reference snapshot.
        ax2.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax2.set_title('Theoretical Solution')                                                           # Right title.
        nok = nom + '_' + str(format(T[k], '.2f'))                                                      # Filename tag based on time.
        plt.savefig(nok + 's.png')                                                                      # Save raster snapshot.
        plt.savefig(nok + 's.svg', format = 'svg')                                                      # Save vector snapshot.
        plt.close()                                                                                     # Close to free resources.


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

    fig, (ax1) = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (5, 5))                # Create a single 3D axis.

    if save:                                                                                            # Save an animation to disk.
        def update_plot(k):                                                                             # Frame update callback for FuncAnimation.
            k = min(k, t - 1)                                                                           # Clamp frame index.
            ax1.clear()                                                                                 # Clear axis for redraw.
            tin = float(T[k])                                                                           # Current normalized time.

            fig.suptitle('Solution at t = %1.3f s.' % tin)                                              # Update figure title.
            ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot approximation at frame k.
            ax1.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax1.set_title('Approximation')                                                              # Axis title.
            ax1.view_init(90, 270)                                                                      # Top-down view.
            ax1.set_zticks([])                                                                          # Hide z ticks for a cleaner top-down view.

            return fig,                                                                                 # Return artists (blit=False, but keep signature).
    
        ani = FuncAnimation(fig, update_plot, frames = np.arange(0, t + 1, step), blit = False)         # Build animation over selected frames.
        try:                                                                                            # Prefer ffmpeg when available.
            ani.save(nom, writer = 'ffmpeg', fps = 10)                                                  # Save video (e.g., .mp4).
        except Exception as e:                                                                          # Fall back when ffmpeg is not available.
            print(f"Warning: ffmpeg not available, using pillow writer instead. Error: {e}")            # Inform the user about the fallback.
            if nom.endswith('.mp4'):                                                                    # If filename suggests mp4, produce a gif instead.
                nom_gif = nom.replace('.mp4', '.gif')                                                   # Replace extension for pillow writer.
                ani.save(nom_gif, writer = 'pillow', fps = 10)                                          # Save as GIF.
                print(f"Video saved as {nom_gif} instead of {nom}")                                     # Inform about the actual saved file.
            else:                                                                                       # Keep the requested name/extension.
                ani.save(nom, writer = 'pillow', fps = 10)                                              # Save using pillow writer.
        plt.close()                                                                                     # Close figure after saving.

    else:
        for k in np.arange(0, t, step):                                                                 # Iterate frames for interactive preview.
            tin = float(T[k])                                                                           # Current normalized time.
            fig.suptitle('Solution at t = %1.3f s.' % tin)                                              # Update figure title.
            
            ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot approximation at time k.
            ax1.set_zlim([min_val, max_val])                                                            # Keep z-limits consistent.
            ax1.set_title('Approximation')                                                              # Axis title.
            ax1.view_init(90, 270)                                                                      # Top-down view.
            ax1.set_zticks([])                                                                          # Hide z ticks for a cleaner top-down view.

            plt.pause(0.1)                                                                              # Allow GUI event loop to update.
            ax1.clear()                                                                                 # Clear axis for the next frame.

        tin = float(T[-1])                                                                              # Final time for display.
        fig.suptitle('Solution at t = %1.3f s.' % tin)                                                  # Final figure title.
        
        ax1.scatter(p[:, 0], p[:, 1], u_ap[:, -1], c = u_ap[:, -1], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Plot final approximation.
        ax1.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax1.set_title('Approximation')                                                                  # Axis title.
        ax1.view_init(90, 270)                                                                          # Top-down view.
        ax1.set_zticks([])                                                                              # Hide z ticks.
        
        plt.pause(0.1)                                                                                  # Keep final frame visible briefly.
        plt.close()                                                                                     # Close figure.


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

    ## Create the graph.
    for k in np.arange(0, t + 1, step):                                                                 # Iterate selected snapshot indices.
        k = min(k, t - 1)                                                                               # Clamp index to valid range.
        fig, (ax1) = plt.subplots(1, 1, subplot_kw = {"projection": "3d"}, figsize = (5, 5))            # Create a single 3D axis.
        tin = float(T[k])                                                                               # Current normalized time.
        
        plt.suptitle('Solution at t = %1.3f s.' % tin)                                                  # Figure title.
        ax1.scatter(p[:, 0], p[:, 1], u_ap[:, k], c = u_ap[:, k], cmap = cm.coolwarm, s = 1, vmin = min_val, vmax = max_val)
                                                                                                        # Approximation snapshot.
        ax1.set_zlim([min_val, max_val])                                                                # Keep z-limits consistent.
        ax1.set_title('Approximation')                                                                  # Axis title.
        ax1.view_init(90, 270)                                                                          # Top-down view.
        ax1.set_zticks([])                                                                              # Hide z ticks.
        
        nok = nom + '_' + str(format(T[k], '.2f'))                                                      # Filename tag based on time.
        plt.savefig(nok + 's.png')                                                                      # Save raster snapshot.
        plt.savefig(nok + 's.svg', format = 'svg')                                                      # Save vector snapshot.
        plt.close()                                                                                     # Close to free resources.
