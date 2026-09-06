"""
run_Perturbation — Stationary advection-dominated perturbation reference

Overview:
    This batch solves a stationary, advection-dominated perturbation problem (a Poisson-like operator
    with advection terms and very small diffusion) on each point cloud under Data/. It uses the meshless
    mGFD stationary solver with Adv=True.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Solve the stationary problem with Stationary(..., upwind=True)
    - Compute error metrics and write outputs to Results/

Public API:
    process_cloud       Process a single point cloud for the Perturbation problem.
    main                Batch runner entry point for the Perturbation problem.

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
from typing import Any                                                                                                                  # Type hinting.

from mGFD import Cloud, Dirichlet, PDE, Solver                                                                                          # OOP Architecture classes.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_stationary                                                                                              # Plotting utilities for the results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory (for local utils).
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for stationary/transient runs.

from utils.batch_utils import run_batch_suite, save_metrics                                                                             # Batch suite runner and metrics persistence.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.

## Variables for the problem.
vx: float = 1.0                                                                                                                         # Advection velocity in x.
vy: float = 0.0                                                                                                                         # Advection velocity in y.
D: float  = 1e-5                                                                                                                        # Diffusion coefficient (very small -> boundary layers).

def phi(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    phi
    Boundary condition for the problem.

    Input:
        x               m           ndarray         x coordinates.
        y               m           ndarray         y coordinates.

    Output:
        phi_val         m           ndarray         Evaluated boundary condition.
    """
    return (x**2 - np.exp(-(1 - x) / D)) * y * (1 - y)                                                                                  # Return output values.

def f(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    f
    Right-hand side forcing term.

    Input:
        x               m           ndarray         x coordinates.
        y               m           ndarray         y coordinates.

    Output:
        f_val           m           ndarray         Evaluated forcing term.
    """
    exp_term = np.exp((x - 1) / D)                                                                                                      # Boundary layer exponential factor.
    term1    = -D * (2 * exp_term - 2 * x**2 + y * (exp_term / D**2 - 2) * (y - 1))                                                     # Diffusion component of forcing.
    term2    = -y * (2 * x - exp_term / D) * (y - 1)                                                                                    # Advective component of forcing.
    return term1 + term2                                                                                                                # Return combined forcing term.

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the stationary perturbation problem (upwind=True) on a single point cloud file.

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
    out_dir   = os.path.join(results_path, 'Perturbation', dataset)                                                                     # Output directory for this region.
    os.makedirs(out_dir, exist_ok=True)                                                                                                 # Ensure output directory exists.

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                                  # Progress message for the batch run.

    nvec            = kwargs.get('nvec', 12)                                                                                            # Extract neighbor count from config, default 12.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Extract verbose flag.
    upwind          = kwargs.get('upwind', True)                                                                                        # Extract upwind flag, default True for Perturbation.
    device          = kwargs.get('device', 'cpu')                                                                                       # Extract device backend, default cpu.
    config_id       = f'nvec_{nvec}_{device}_upwind_{upwind}'                                                                           # Create unique config identifier for the sweep.

    # 2. Data Loading
    p         = load_points(cloud_path, verbose=False)                                                                                  # Load point cloud into (m, 3) array [x, y, flag].
    L         = np.vstack([[vx], [vy], [-D], [0], [-D], [0]])                                                                           # Operator coefficients for PDE.

    # 3. Solver Execution via OOP Interface
    cloud     = Cloud.from_array(p)                                                                                                     # Instantiate Cloud.
    domain    = cloud.set_boundary(Dirichlet(phi))                                                                                      # Set Dirichlet boundary.
    pde       = PDE(operator=L, source=f, order=0)                                                                                      # Formulate Stationary PDE with L.
    solver    = Solver(domain, pde, device=device, nvec=nvec, upwind=upwind, verbose=verbose_solvers, implicit=True)                    # Create Solver.
    result    = solver.solve()                                                                                                          # Solve PDE.
    u_ap, vec = result.solution, result.neighbors                                                                                       # Unpack solution and neighbors.
    comp_time = result.compute_time                                                                                                     # Get compute time.

    # 4. Exact Solution and Metrics
    u_ex    = phi(p[:, 0], p[:, 1])                                                                                                     # Compute exact theoretical solution locally.
    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time=comp_time)                                             # Compute comprehensive stationary error metrics.

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (RMSE): {metrics["RMSE"]}')                                                                               # Print RMSE error for quick inspection.

    # 5. Output persistence
    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        plot_scales = kwargs.get('plot_scales', ['2'])                                                                                  # Extract target plotting scales.
        if plot_scales is None or 'all' in plot_scales or scale in [str(s) for s in plot_scales]:                                       # Filter rendering for selected scale(s).
            cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                               # Extract clean cloud name.
            appx_nom   = os.path.join(out_dir, f'Appx_{config_id}_scale_{scale}_{cloud_name}')                                          # Define numerical approximation plot path.
            plot_stationary(p, u_ap, save=True, nom=appx_nom, title='Stationary Appx', verbose=verbose)                                 # Save stationary approximation plot.

            exact_nom = os.path.join(out_dir, f'Exact_scale_{scale}_{cloud_name}')                                                      # Define theoretical exact solution plot path.
            if not os.path.exists(exact_nom + '.png'):                                                                                  # Avoid regenerating exact plot if already present.
                plot_stationary(p, u_ex, save=True, nom=exact_nom, title='Theoretical Solution', verbose=verbose)                       # Save exact stationary plot.

def main(**kwargs: Any) -> None:
    """
    main
    Entry point for the Perturbation batch script.

    Input:
        **kwargs                    Any             Configuration values from main orchestrator.

    Output:
        None
    """
    run_batch_suite(process_cloud, DATA_ROOT, RESULTS_ROOT, **kwargs)                                                                   # Execute universal batch orchestrator.

if __name__ == "__main__":                                                                                                              # Evaluate condition.
    main()                                                                                                                              # Execute statement.
