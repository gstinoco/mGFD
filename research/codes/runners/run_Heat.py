"""
run_Heat — Reference batch for the 2D Heat equation

Overview:
    This script runs the first-order transient Heat reference problem on all available point clouds under
    Data/, using the meshless mGFD solver.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with TimeDerivative1 (implicit scheme by default)
    - Compute error metrics and save outputs to Results/
    - Plot/save static step snapshots (optional) and a transient animation (always)

Public API:
    process_cloud       Process a single point cloud for the Heat equation.
    main                Batch runner entry point for the Heat problem.

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
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Numerical arrays and math.
import pandas as pd                                                                                                                     # Dataframes and series for new v0.10.0 interface.
from typing import List, Any                                                                                                            # Type hinting.

from mGFD import TimeDerivative1                                                                                                        # First-order transient solver to run the reference case.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_transient                                                                                               # Plotting utilities for the results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory (for local utils).
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for stationary/transient runs.

from utils.batch_utils import load_neighbors, save_neighbors, run_batch_suite, save_metrics                                             # Dataset loading + neighbor cache helpers.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.

## Problem parameters.
v: float = 0.2                                                                                                                          # Diffusion coefficient.
t: int   = 2000                                                                                                                         # Number of time steps.

def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Heat analytical solution.
    """
    return np.exp(-2 * np.pi**2 * coef[0] * t_val) * np.sin(np.pi * x) * np.sin(np.pi * y)                                              # Return output values.

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the transient Heat problem on a single point cloud file.

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
    out_dir   = os.path.join(results_path, 'Heat', dataset)                                                                             # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                                                               # Ensure output directory exists (even if save=False).
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                                  # Progress message for the batch run.

    nvec            = kwargs.get('nvec', 12)                                                                                            # Extract neighbor count from config, default 12.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Extract verbose flag.
    device          = kwargs.get('device', 'cpu')                                                                                       # Extract device backend, default cpu.
    input_types     = kwargs.get('input_types', ['callable'])                                                                           # Extract input_types, default callable.
    config_id       = f'nvec_{nvec}_{device}'                                                                                           # Create unique config identifier for the sweep.

    # 2. Data Loading & Neighbor Cache
    p    = load_points(cloud_path, verbose=False)                                                                                       # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, nvec)                                                                                             # Load cached neighbor list if present.
    L    = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                            # Operator coefficients for the Laplacian with diffusion.
    
    # 3. Solver Execution
    u_ap, vec = None, vec0                                                                                                              # Initialize solution and neighbors.
    comp_time = 0.0                                                                                                                     # Initialize compute time.
    t_param   = kwargs.get('t', None)                                                                                                   # Extract optional explicit t.
    cfl_param = kwargs.get('cfl', 0.5)                                                                                                  # Extract optional target CFL.
    res_main  = None                                                                                                                    # Main solver result reference.
        
    # --- A. Using Callable ---
    if 'callable' in input_types:                                                                                                       # Check if callable test is enabled.
        res_call = TimeDerivative1(                                                                                                     # Solve the transient heat problem (Callable).
            p, f, t=t_param, cfl=cfl_param, coef=[v], operator=L, vec=vec0, nvec=nvec, implicit=True, lam=0.5, device=device, verbose=verbose_solvers # Assign p, f, t.
        )                                                                                                                               # Extract solver result object.
        u_ap, vec  = res_call.solution, res_call.neighbors                                                                              # Unpack approximate solution and neighbor list.
        comp_time  = res_call.compute_time                                                                                              # Get solver execution time.
        res_main   = res_call                                                                                                           # Track main result object.
        
    t_used = res_main.t_steps if res_main is not None and res_main.t_steps is not None else (t_param if t_param is not None else 2000)  # Determine actual time steps executed.
    T_arr = np.linspace(0, 1, t_used)                                                                                                   # Reconstruct time vector.
    spatial_part  = np.sin(np.pi * p[:, 0]) * np.sin(np.pi * p[:, 1])                                                                   # Spatial part of solution.
    temporal_part = np.exp(-2 * v * np.pi**2 * T_arr)                                                                                   # Temporal decay part.
    f_arr         = np.outer(spatial_part, temporal_part)                                                                               # Vectorized outer product for exact solution.
        
    # --- B. Using Numpy Arrays ---
    if 'array' in input_types:                                                                                                          # Check if array test is enabled.
        res_arr = TimeDerivative1(                                                                                                      # Solve using array inputs.
            p, f_arr, t=t_used, coef=[v], operator=L, vec=vec0, nvec=nvec, implicit=True, lam=0.5, device=device, verbose=False         # Assign p, f_arr, t.
        )                                                                                                                               # Execute statement.
        if u_ap is not None:                                                                                                            # If previous result exists, validate.
            assert np.allclose(u_ap, res_arr.solution), "Mismatch between Callable and Array solver outputs."                           # Validate output equivalence.
        else:                                                                                                                           # If callable was skipped.
            u_ap, vec  = res_arr.solution, res_arr.neighbors                                                                            # Unpack approximate solution and neighbor list.
            comp_time  = res_arr.compute_time                                                                                           # Get solver execution time.
            res_main   = res_arr                                                                                                        # Track main result object.
    
    # --- C. Using Pandas DataFrames/Series ---
    if 'pandas' in input_types:                                                                                                         # Check if pandas test is enabled.
        f_pd = pd.DataFrame(f_arr)                                                                                                      # Wrap spatiotemporal array in Pandas DataFrame.
        res_pd = TimeDerivative1(                                                                                                       # Solve using Pandas inputs.
            p, f_pd, t=t_used, coef=[v], operator=L, vec=vec0, nvec=nvec, implicit=True, lam=0.5, device=device, verbose=False          # Assign p, f_pd, t.
        )                                                                                                                               # Execute statement.
        if u_ap is not None:                                                                                                            # If previous result exists, validate.
            assert np.allclose(u_ap, res_pd.solution), "Mismatch between Array/Callable and Pandas solver outputs."                     # Validate output equivalence.
        else:                                                                                                                           # If no previous result exists.
            u_ap, vec  = res_pd.solution, res_pd.neighbors                                                                              # Unpack approximate solution and neighbor list.
            comp_time  = res_pd.compute_time                                                                                            # Get solver execution time.
            res_main   = res_pd                                                                                                         # Track main result object.
            
    if u_ap is None: raise ValueError("No valid input_types were specified.")                                                           # Safety fallback.
    
    # 4. Exact Solution and Metrics
    u_ex    = f_arr                                                                                                                     # Compute exact theoretical solution locally.
    metrics = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)                                              # Compute comprehensive transient error metrics.
    
    if res_main is not None:                                                                                                            # Track CFL metrics.
        metrics['CFL']        = res_main.cfl                                                                                            # Store CFL number.
        metrics['dt']         = res_main.dt                                                                                             # Store time step size.
        metrics['Time_Steps'] = res_main.t_steps                                                                                        # Store total time step count.
    
    # Track extra compute times for numpy/pandas to demonstrate overhead
    if 'array' in input_types: metrics['Time_Array']  = res_arr.compute_time                                                            # Array execution time.
    if 'pandas' in input_types: metrics['Time_Pandas'] = res_pd.compute_time                                                            # Pandas execution time.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                                                # Print average error for quick inspection.

    # 5. Output persistence
    if vec0 is None:                                                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, nvec, vec)                                                                                           # Save vec to the canonical cache file.

    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                                   # Extract clean cloud name.
        if scale == '3':                                                                                                                # Only for scale 3.
            if config_id.startswith('nvec_20_spsolve') or kwargs.get('plot_approximations', False):                                     # Only plot baseline config or if explicitly requested.
                plot_transient(p, u_ap, save=True, nom=os.path.join(out_dir, f'Appx_{config_id}_{cloud_name}'),                         # Assign plot_transient(p, u_ap, save.
                                            title='Transient Approximation', verbose=verbose)                                           # Save transient animation.
            exact_nom = os.path.join(out_dir, f'Exact_{cloud_name}')                                                                    # Define exact solution filename linked to the cloud.
            if not (os.path.exists(exact_nom + '.mp4') or os.path.exists(exact_nom + '.gif')):                                          # Avoid regenerating the exact solution.
                plot_transient(p, u_ex, save=True, nom=exact_nom,                                                                       # Assign plot_transient(p, u_ex, save.
                                            title='Theoretical Solution', verbose=verbose)                                              # Save exact transient animation.

def main(**kwargs: Any) -> None:
    """
    main
    Entry point for the Heat batch script.

    Input:
        **kwargs                    Any             Configuration values from main orchestrator.

    Output:
        None
    """
    run_batch_suite(process_cloud, DATA_ROOT, RESULTS_ROOT, **kwargs)                                                                   # Execute universal batch orchestrator.

if __name__ == "__main__":                                                                                                              # Evaluate condition.
    main()                                                                                                                              # Execute statement.
