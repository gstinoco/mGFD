"""
run_Poisson — Reference batch for the 2D Poisson equation

Overview:
    This script runs the stationary Poisson reference problem on all available point clouds under Data/
    (both Clouds and Holes datasets), using the meshless mGFD stationary solver.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv (or any configured scales)
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with Stationary
    - Compute error metrics and save outputs to Results/
    - Plot/save a stationary comparison figure

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
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.path manipulation.
import time                                                                                                                             # Time tracking for execution performance.
import json                                                                                                                             # JSON serialization for metrics.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Numerical arrays and math.
import pandas as pd                                                                                                                     # Dataframes and series for new v0.10.0 interface.
from typing import Optional, List, Callable, Any                                                                                        # Type hinting.

import mGFD.io.export_vtk as ExportVTK                                                                                                  # VTK export utilities for ParaView.

from mGFD import Stationary                                                                                                             # First-order transient solver to run the reference case.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_stationary                                                                                              # Plotting utilities for the results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory (for local utils).
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for stationary/transient runs.

from utils.batch_utils import iter_clouds, load_neighbors, save_neighbors, run_batch_suite, save_metrics                                # Dataset loading + neighbor cache helpers.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.
SCALES: tuple     = ('1', '2', '3')                                                                                                     # Scales to process under each dataset.

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
    return 2 * np.exp(2 * x + y)

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
    return 10 * np.exp(2 * x + y)

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the stationary Poisson-type problem on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1', '2').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays.
        verbose                     bool            If True, prints progress and errors to console.
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
    out_dir   = os.path.join(results_path, 'Poisson', dataset)                                                                          # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                                                               # Ensure output directory exists.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                                  # Progress message for the batch run.

    nvec            = kwargs.get('nvec', 12)                                                                                            # Extract neighbor count from config, default 12.
    linear_solver   = kwargs.get('linear_solver', 'spsolve')                                                                            # Extract solver backend, default spsolve.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Extract verbose flag.
    device          = kwargs.get('device', 'cpu')                                                                                       # Extract device backend, default cpu.
    matrix_free     = kwargs.get('matrix_free', False)                                                                                  # Extract matrix_free flag, default False.
    preconditioner  = kwargs.get('preconditioner', None)                                                                                # Extract preconditioner, default None.
    config_id       = f'nvec_{nvec}_{linear_solver}_{device}_mf{matrix_free}_pre{preconditioner}'                                                                                    # Create unique config identifier for the sweep.

    # 2. Data Loading & Neighbor Cache
    p    = load_points(cloud_path, verbose=False)                                                                                       # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, nvec)                                                                                             # Load cached neighbor list if present.
    L    = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                                    # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.
    
    # 3. Solver Execution
    # --- A. Using Callable ---
    res_call = Stationary(                                                                                                              # Solve the stationary Poisson problem (Callable).
        p, phi, f, operator = L, vec = vec0, nvec = nvec, linear_solver = linear_solver, device = device, matrix_free = matrix_free, preconditioner = preconditioner, verbose = verbose_solvers                      # Execute with dynamic config.
    )                                                                                                                                   # Extract solver result object.
    
    # --- B. Using Numpy Arrays ---
    phi_arr = phi(p[:, 0], p[:, 1])                                                                                                     # Precompute boundary array.
    f_arr   = f(p[:, 0], p[:, 1])                                                                                                       # Precompute forcing term array.
    res_arr = Stationary(                                                                                                               # Solve using array inputs.
        p, phi_arr, f_arr, operator = L, vec = vec0, nvec = nvec, linear_solver = linear_solver, device = device, matrix_free = matrix_free, preconditioner = preconditioner, verbose = False                        # Silent execution for array test.
    )
    assert np.allclose(res_call.solution, res_arr.solution), "Mismatch between Callable and Array solver outputs."                      # Validate output equivalence.
    
    # --- C. Using Pandas DataFrames/Series ---
    phi_pd = pd.Series(phi_arr.tolist())                                                                                                # Wrap array in Pandas Series.
    f_pd   = pd.Series(f_arr.tolist())                                                                                                  # Wrap array in Pandas Series.
    res_pd = Stationary(                                                                                                                # Solve using Pandas inputs.
        p, phi_pd, f_pd, operator = L, vec = vec0, nvec = nvec, linear_solver = linear_solver, device = device, matrix_free = matrix_free, preconditioner = preconditioner, verbose = False                          # Silent execution for Pandas test.
    )
    assert np.allclose(res_call.solution, res_pd.solution), "Mismatch between Callable and Pandas solver outputs."                      # Validate output equivalence.
    
    u_ap, vec  = res_call.solution, res_call.neighbors                                                                                  # Unpack approximate solution and neighbor list.
    comp_time  = res_call.compute_time                                                                                                  # Get solver execution time from v0.10.0 dataclass.
    
    # 4. Exact Solution and Metrics
    u_ex    = phi_arr                                                                                                                   # Compute exact theoretical solution locally (already computed as phi_arr).
    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time = comp_time)                                           # Compute comprehensive stationary error metrics.
    
    # Track extra compute times for numpy/pandas to demonstrate overhead
    metrics['Time_Array']  = res_arr.compute_time                                                                                       # Array execution time.
    metrics['Time_Pandas'] = res_pd.compute_time                                                                                        # Pandas execution time.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (RMSE): {metrics["RMSE"]}')                                                                               # Print RMSE error for quick inspection.

    # 5. Output persistence
    if vec0 is None:                                                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, nvec, vec)                                                                                           # Save vec to the canonical cache file.

    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        if scale == '3':                                                                                                                # Only for scale 3.
            if config_id.startswith('nvec_16_spsolve'):                                                                                 # Only plot baseline config.
                plot_stationary(p, u_ap, save=True, nom=os.path.join(out_dir, f'Approximation_{config_id}'),
                            title='Stationary Appx', verbose=verbose)                                                                       # Save 3D scatter image.
        
            exact_nom = os.path.join(out_dir, 'Exact')                                                                                  # Define exact solution filename.
            if not os.path.exists(exact_nom + '.png'):                                                                                  # Avoid regenerating the exact solution.
                plot_stationary(p, u_ex, save=True, nom=exact_nom,
                            title='Theoretical Solution', verbose=verbose)                                                              # Create independent plot of exact solution.

def main(**kwargs: Any) -> None:
    """
    main
    Entry point for the Poisson batch script.

    Input:
        **kwargs                    Any             Configuration values from main orchestrator.

    Output:
        None
    """
    run_batch_suite(process_cloud, DATA_ROOT, RESULTS_ROOT, SCALES, **kwargs)                                                           # Execute universal batch orchestrator.

if __name__ == "__main__":
    main()
