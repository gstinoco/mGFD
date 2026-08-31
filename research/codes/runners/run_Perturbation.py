"""
run_Perturbation — Stationary advection-dominated perturbation reference

Overview:
    This batch solves a stationary, advection-dominated perturbation problem (a Poisson-like operator
    with advection terms and very small diffusion) on each point cloud under Data/. It uses the meshless
    mGFD stationary solver with Adv=True.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load tagged neighbor caches (the cache depends on the advection parameters)
    - Validate the neighbor cache for interior nodes and recompute if insufficient neighbors exist
    - Solve the stationary problem with Stationary(..., upwind=True)
    - Compute error metrics and write outputs to Results/

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

## Variables for the problem.
vx: float  = 1.0                                                                                                                        # Advection velocity in x.
vy: float  = 0.0                                                                                                                        # Advection velocity in y.
D: float   = 1e-5                                                                                                                       # Diffusion coefficient (very small -> boundary layers).

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
    return (x**2 - np.exp(-(1 - x) / D)) * y * (1 - y)

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
    return -(1e-5) * (2 * np.exp((x - 1) / (1e-5)) - 2 * x**2 + y * (np.exp((x - 1) / (1e-5)) / (1e-5)**2 - 2) * (y - 1)) - y * (2 * x - np.exp((x - 1) / (1e-5)) / (1e-5)) * (y - 1)

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the stationary perturbation problem (upwind=True) on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Catemaco', 'Chapala').
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
    out_dir   = os.path.join(results_path, 'Perturbation', dataset)                                                                     # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                                                               # Ensure output directory exists.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                                  # Progress message for the batch run.

    nvec            = kwargs.get('nvec', 12)                                                                                            # Extract neighbor count from config, default 12.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Extract verbose flag.
    upwind          = kwargs.get('upwind', True)                                                                                        # Extract upwind flag, default True for Perturbation.
    device          = kwargs.get('device', 'cpu')                                                                                       # Extract device backend, default cpu.
    input_types     = kwargs.get('input_types', ['callable', 'array', 'pandas'])                                                        # Extract input_types, default all.
    config_id       = f'nvec_{nvec}_{device}_upwind_{upwind}'                                                                           # Create unique config identifier for the sweep.

    # 2. Data Loading & Neighbor Cache
    p    = load_points(cloud_path, verbose=False)                                                                                       # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, nvec)                                                                                             # Load cached neighbor list if present.
    L    = np.vstack([[vx], [vy], [-D], [0], [-D], [0]])                                                                                # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.
    
    # 3. Solver Execution
    u_ap, vec = None, vec0                                                                                                              # Initialize solution and neighbors.
    comp_time = 0.0                                                                                                                     # Initialize compute time.
    phi_arr = phi(p[:, 0], p[:, 1])                                                                                                     # Precompute boundary array.
    f_arr   = f(p[:, 0], p[:, 1])                                                                                                       # Precompute forcing term array.

    # --- A. Using Callable ---
    if 'callable' in input_types:                                                                                                       # Check if callable test is enabled.
        res_call = Stationary(                                                                                                          # Solve the stationary perturbation problem (Callable).
            p, phi, f, operator = L, upwind = upwind, vec = vec0, nvec = nvec, device = device, verbose = verbose_solvers               # Execute with dynamic config.
        )                                                                                                                               # Extract solver result object.
        u_ap, vec  = res_call.solution, res_call.neighbors                                                                              # Unpack approximate solution and neighbor list.
        comp_time  = res_call.compute_time                                                                                              # Get solver execution time from v0.10.0 dataclass.
    
    # --- B. Using Numpy Arrays ---
    if 'array' in input_types:                                                                                                          # Check if array test is enabled.
        res_arr = Stationary(                                                                                                           # Solve using array inputs.
            p, phi_arr, f_arr, operator = L, upwind = upwind, vec = vec0, nvec = nvec, device = device, verbose = False                 # Silent execution for array test.
        )
        if u_ap is not None:                                                                                                            # If previous result exists, validate.
            assert np.allclose(u_ap, res_arr.solution), "Mismatch between Callable and Array solver outputs."                           # Validate output equivalence.
        else:                                                                                                                           # If callable was skipped.
            u_ap, vec  = res_arr.solution, res_arr.neighbors                                                                            # Unpack approximate solution and neighbor list.
            comp_time  = res_arr.compute_time                                                                                           # Get solver execution time from v0.10.0 dataclass.
    
    # --- C. Using Pandas DataFrames/Series ---
    if 'pandas' in input_types:                                                                                                         # Check if pandas test is enabled.
        phi_pd = pd.Series(phi_arr.tolist())                                                                                            # Wrap array in Pandas Series.
        f_pd   = pd.Series(f_arr.tolist())                                                                                              # Wrap array in Pandas Series.
        res_pd = Stationary(                                                                                                            # Solve using Pandas inputs.
            p, phi_pd, f_pd, operator = L, upwind = upwind, vec = vec0, nvec = nvec, device = device, verbose = False                 # Silent execution for Pandas test.
        )
        if u_ap is not None:                                                                                                            # If previous result exists, validate.
            assert np.allclose(u_ap, res_pd.solution), "Mismatch between Array/Callable and Pandas solver outputs."                     # Validate output equivalence.
        else:                                                                                                                           # If no previous result exists.
            u_ap, vec  = res_pd.solution, res_pd.neighbors                                                                              # Unpack approximate solution and neighbor list.
            comp_time  = res_pd.compute_time                                                                                            # Get solver execution time from v0.10.0 dataclass.
            
    if u_ap is None: raise ValueError("No valid input_types were specified.")                                                           # Safety fallback.
    
    # 4. Exact Solution and Metrics
    u_ex    = phi_arr                                                                                                                   # Compute exact theoretical solution locally (already computed as phi_arr).
    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time=comp_time)                                             # Compute comprehensive stationary error metrics.
    
    # Track extra compute times for numpy/pandas to demonstrate overhead
    if 'array' in input_types: metrics['Time_Array']  = res_arr.compute_time                                                            # Array execution time if computed.
    if 'pandas' in input_types: metrics['Time_Pandas'] = res_pd.compute_time                                                            # Pandas execution time if computed.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (RMSE): {metrics["RMSE"]}')                                                                               # Print RMSE error for quick inspection.

    # 5. Output persistence
    if vec0 is None:                                                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, nvec, vec)                                                                                           # Save vec to the canonical cache file.

    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                                   # Extract clean cloud name.
        if scale == '3':                                                                                                                # Only for scale 3.
            if config_id.startswith('nvec_20_spsolve') or kwargs.get('plot_approximations', False):                                     # Only plot baseline config or if explicitly requested.
                plot_stationary(p, u_ap, save=True, nom=os.path.join(out_dir, f'Appx_{config_id}_{cloud_name}'),
                            title='Stationary Appx', verbose=verbose)                                                                   # Save 3D scatter image.
        
            exact_nom = os.path.join(out_dir, f'Exact_{cloud_name}')                                                                    # Define exact solution filename linked to the cloud.
            if not os.path.exists(exact_nom + '.png'):                                                                                  # Check for PNG rendering.
                plot_stationary(p, u_ex, save=True, nom=exact_nom,
                            title='Theoretical Solution', verbose=verbose)                                                              # Create independent plot of exact solution.

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

if __name__ == "__main__":
    main()
