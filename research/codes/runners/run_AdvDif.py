"""
run_AdvDif — Reference batch for the 2D Advection–Diffusion equation

Overview:
    This script runs the first-order transient Advection–Diffusion reference problem on all available
    point clouds under Data/, using the meshless mGFD solver.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with TimeDerivative1 (implicit scheme by default)
    - Compute error metrics and save outputs to Results/
    - Plot/save static step snapshots (optional) and a transient animation (always)

Public API:
    process_cloud       Process a single point cloud for the Advection-Diffusion equation.
    main                Batch runner entry point for the Advection-Diffusion problem.

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
    September, 2026.
"""

## Library importation.
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.path manipulation.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Numerical arrays and math.

from typing import List, Any, Optional                                                                                                  # Type hinting.

from mGFD import Cloud, Dirichlet, PDE, Solver                                                                                          # OOP Architecture classes.
from mGFD.solvers.results import SolverResult                                                                                           # Standardized SolverResult dataclass.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_transient                                                                                               # Plotting utilities for the results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory (for local utils).
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for stationary/transient runs.

from utils.batch_utils import run_batch_suite, save_metrics                                                                             # Batch suite runner and metrics persistence.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.

## Problem parameters.
v: float = 0.1                                                                                                                          # Diffusion coefficient.
a: float = 0.3                                                                                                                          # Transport velocity on the x direction.
b: float = 0.2                                                                                                                          # Transport velocity on the y direction.

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the transient Advection-Diffusion problem on a single point cloud file.

    Input:
        dataset         1           str             Dataset folder name under Data/ (e.g., 'Catemaco', 'Chapala').
        scale           1           str             Cloud scale folder (e.g., '1', '2').
        cloud_path      1           str             Path to input CSV with point cloud.
        results_path    1           str             Base output directory (typically <repo>/Results).
        save            1           bool            Whether to save the solution arrays.
        verbose         1           bool            If True, prints progress and errors to console.
        **kwargs                    Any             Configuration values from main orchestrator.

    Output:
        None
    """
    # 0. Input validation
    if not isinstance(dataset, str):                                                                                                    # Validate dataset argument.
        raise TypeError("Dataset name must be a string.")                                                                               # Raise explicit error on bad input.
    if not isinstance(scale, str):                                                                                                      # Validate scale argument.
        raise TypeError("Scale must be a string.")                                                                                      # Raise explicit error on bad input.
    if not isinstance(cloud_path, str) or not os.path.exists(cloud_path):                                                               # Validate cloud path.
        raise ValueError("Cloud path must be a valid existing file path.")                                                              # Raise explicit error on bad input.

    # 1. Variable initialization
    region_id = f'{dataset}/{scale}'                                                                                                    # Region identifier.
    out_dir   = os.path.join(results_path, 'Advection-Diffusion', dataset)                                                              # Output directory for this region.
    os.makedirs(out_dir, exist_ok=True)                                                                                                 # Ensure output directory exists (even if save=False).

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                                  # Progress message for the batch run.

    # 2. Data Loading & Neighbor Cache
    p = load_points(cloud_path, verbose=False)                                                                                          # Load point cloud into (m, 3) array [x, y, flag].

    interior = p[p[:, 2] == 0]                                                                                                          # Extract interior nodes.

    if len(interior) > 0:                                                                                                               # Calculate pulse center dynamically inside the interior region.
        x_mean, y_mean = np.mean(p[:, 0]), np.mean(p[:, 1])                                                                             # Compute arithmetic mean.
        dists          = (interior[:, 0] - x_mean)**2 + (interior[:, 1] - y_mean)**2                                                    # Compute distances to mean.
        idx            = np.argmin(dists)                                                                                               # Find closest interior point.
        xc, yc         = interior[idx, 0], interior[idx, 1]                                                                             # Assign dynamic center coordinates.
    else:                                                                                                                               # Fallback for degenerate clouds.
        xc, yc         = np.mean(p[:, 0]), np.mean(p[:, 1])                                                                             # Assign mean center coordinates.

    def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                                 # Define dynamic analytical solution with xc, yc.
        """
        f
        Advection-Diffusion analytical solution.

        Input:
            x               m           ndarray         x coordinates.
            y               m           ndarray         y coordinates.
            t_val                       float           Time step value.
            coef            3           List[float]     Coefficients [v, a, b].

        Output:
            f_val           m           ndarray         Evaluated analytical solution.
        """
        denom = coef[0] * (4 * t_val + 1)                                                                                               # Compute diffusion denominator.
        decay = 1 / (4 * t_val + 1)                                                                                                     # Compute temporal amplitude factor.
        dist2 = (x - coef[1] * t_val - xc)**2 + (y - coef[2] * t_val - yc)**2                                                           # Quadratic advective distance.
        return decay * np.exp(-dist2 / denom)                                                                                           # Return evaluated analytical solution.

    nvec            = kwargs.get('nvec', 12)                                                                                            # Extract neighbor count from config, default 12.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Extract verbose flag.
    upwind          = kwargs.get('upwind', True)                                                                                        # Extract upwind flag, default True for AdvDif.
    device          = kwargs.get('device', 'cpu')                                                                                       # Extract device backend, default cpu.
    t_param         = kwargs.get('t', None)                                                                                             # Extract optional explicit t.
    cfl_param       = kwargs.get('cfl', 0.5)                                                                                            # Extract optional target CFL.
    config_id       = f'nvec_{nvec}_{device}_upwind_{upwind}'                                                                           # Create unique config identifier for the sweep.

    # 2. Data Loading
    # (p loaded in block 2 above)

    # 3. Solver Execution via OOP Interface
    cloud           = Cloud.from_array(p)                                                                                               # Instantiate Cloud.
    ic_func         = lambda x, y: f(x, y, 0.0, [v, a, b])                                                                              # Analytical initial condition u(x, y, 0).
    domain          = cloud.set_boundary(Dirichlet(f))                                                                                  # Set Dirichlet boundary.
    pde             = PDE(operator=[-a, -b, 2 * v, 0, 2 * v, 0], ic=ic_func, order=1, coef=[v, a, b])                                   # Formulate Advection-Diffusion PDE.
    solver = Solver(domain, pde, device=device, cfl=cfl_param, nvec=nvec, upwind=upwind, verbose=verbose_solvers, implicit=True)        # Create Solver.
    dt_input        = (1.0 / int(t_param)) if t_param is not None else None                                                             # Explicit dt or None for adaptive CFL.
    result          = solver.solve(dt=dt_input)                                                                                         # Solve PDE with adaptive CFL or custom step.
    u_ap, vec       = result.solution, result.neighbors                                                                                 # Unpack solution and neighbors.
    comp_time       = result.compute_time                                                                                               # Get compute time.

    # 4. Exact Solution and Metrics
    t_used          = result.t_steps if result.t_steps is not None else (int(t_param) if t_param is not None else 2)                    # Retrieve actual simulated time steps from solver.
    T_arr           = np.linspace(0, 1, t_used)                                                                                         # Reconstruct time vector.
    X               = p[:, 0, None]                                                                                                     # Node X coordinates broadcast shape (m, 1)
    Y               = p[:, 1, None]                                                                                                     # Node Y coordinates broadcast shape (m, 1).
    T               = T_arr[None, :]                                                                                                    # Time steps broadcast shape (1, t).
    denom_t         = v * (4 * T + 1)                                                                                                   # Compute spatio-temporal denominator.
    decay_t         = 1 / (4 * T + 1)                                                                                                   # Compute spatio-temporal amplitude factor.
    dist2_t         = (X - a * T - xc)**2 + (Y - b * T - yc)**2                                                                         # Compute spatio-temporal distance.
    u_ex            = decay_t * np.exp(-dist2_t / denom_t)                                                                              # Vectorized computation of analytical matrix.
    metrics         = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)                                      # Compute comprehensive transient error metrics.

    metrics['CFL']        = result.cfl                                                                                                  # Store CFL number.
    metrics['dt']         = result.dt                                                                                                   # Store time step size.
    metrics['Time_Steps'] = result.t_steps                                                                                              # Store total time step count.

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                                                # Print average error for quick inspection.

    # 5. Output persistence
    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        plot_scales = kwargs.get('plot_scales', ['2'])                                                                                  # Extract target plotting scales.
        if plot_scales is None or 'all' in plot_scales or scale in [str(s) for s in plot_scales]:                                       # Filter rendering for selected scale(s).
            cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                               # Extract clean cloud name.
            appx_nom   = os.path.join(out_dir, f'Appx_{config_id}_scale_{scale}_{cloud_name}')                                          # Define numerical approximation plot path.
            plot_transient(p, u_ap, save=True, nom=appx_nom, title='Transient Approximation', verbose=verbose)                          # Save transient animation of approximation.

            exact_nom = os.path.join(out_dir, f'Exact_scale_{scale}_{cloud_name}')                                                      # Define theoretical exact solution plot path.
            if not (os.path.exists(exact_nom + '.mp4') or os.path.exists(exact_nom + '.gif')):                                          # Avoid regenerating exact animation if already present.
                plot_transient(p, u_ex, save=True, nom=exact_nom, title='Theoretical Solution', verbose=verbose)                        # Save exact transient animation.

def main(**kwargs: Any) -> None:
    """
    main
    Entry point for the Advection-Diffusion batch script.

    Input:
        **kwargs                    Any             Configuration values from main orchestrator.

    Output:
        None
    """
    run_batch_suite(process_cloud, DATA_ROOT, RESULTS_ROOT, **kwargs)                                                                   # Execute universal batch orchestrator.

if __name__ == "__main__":                                                                                                              # Evaluate condition.
    main()                                                                                                                              # Execute statement.
