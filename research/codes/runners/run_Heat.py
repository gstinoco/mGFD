"""
run_Heat — Reference batch for the 2D Heat equation

Overview:
    This script runs the first-order transient Heat reference problem on all available point clouds under
    Data/ (both Clouds and Holes datasets), using the meshless mGFD solver.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with TimeDerivative1 (implicit scheme by default)
    - Compute error metrics and save outputs to Results/
    - Plot/save static step snapshots (optional) and a transient animation (always)

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
import os                                                                                                                       # Filesystem and path utilities.
import sys                                                                                                                      # sys.path manipulation so this script can import project modules.
import time                                                                                                                     # Time tracking for execution performance.
import json                                                                                                                     # JSON serialization for metrics.
import logging                                                                                                                  # Standard logging module.
import numpy as np                                                                                                              # Numerical arrays and math.
from typing import Optional, List                                                                                               # Type hinting.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                     # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                                                       # Enable imports like "from mGFD import TimeDerivative1".

import utils.metrics as Errors                                                                                                        # Error metrics for stationary/transient runs.
import mGFD.io.export_vtk as ExportVTK                                                                                          # VTK export utilities for ParaView.
from mGFD import TimeDerivative1                                                                                                # First-order transient solver to run the reference case.
from mGFD.viz.graph import plot_transient                                                                                       # Plotting utilities for the results.
from mGFD.io.io import load_points                                                                                              # Point cloud loading utility.
from utils.batch_utils import iter_clouds, load_neighbors, save_neighbors                                                             # Dataset loading + neighbor cache helpers.

logger = logging.getLogger(__name__)                                                                                            # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                   # Basic logger configuration.

DATA_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                                 # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                           # Output results root directory.
SCALES: tuple = ('1', '2', '3', '4', '5')                                                                                       # Scales to process under each dataset.
NVEC: int = 12                                                                                                                  # Neighbor count used by the solver.

## Problem parameters.
v: float = 0.2                                                                                                                  # Diffusion coefficient.
t: int = 2000                                                                                                                   # Number of time steps.

def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Heat analytical solution.
    """
    return np.exp(-2 * np.pi**2 * coef[0] * t_val) * np.cos(np.pi * x) * np.cos(np.pi * y)

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True) -> None:      # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Heat benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1', '2').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays and step plots.
        verbose                     bool            If True, prints progress and errors to console.

    Output:
        None
    """
    # 0. Input validation
    if not isinstance(dataset, str):                                                                                            # Validate dataset argument.
        raise TypeError("Dataset name must be a string.")                                                                       # Raise explicit error on bad input.
    if not isinstance(scale, str):                                                                                              # Validate scale argument.
        raise TypeError("Scale must be a string.")                                                                              # Raise explicit error on bad input.
    if not isinstance(cloud_path, str) or not os.path.exists(cloud_path):                                                       # Validate cloud path.
        raise ValueError("Cloud path must be a valid existing file path.")                                                      # Raise explicit error on bad input.

    # 1. Variable initialization
    region_id = f'{dataset}/{scale}'                                                                                            # Region identifier.
    out_dir = os.path.join(results_path, 'Heat', dataset, scale)                                                                # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                                                       # Ensure output directory exists (even if save=False).
    
    if verbose:                                                                                                                 # Check if verbosity is enabled.
        logger.info(f'Working on region: {region_id}')                                                                          # Progress message for the batch run.

    # 2. Data Loading & Neighbor Cache
    p = load_points(cloud_path, verbose=False)                                                                                  # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, NVEC)                                                                                     # Load cached neighbor list if present.
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                       # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.
    
    # 3. Solver Execution
    start_time = time.time()                                                                                                    # Start execution timer.
    u_ap, vec = TimeDerivative1(                                                                                                # Solve Heat with implicit time stepping.
        p, f, t, [v], operator = L, implicit = True, lam = 0.5, vec = vec0, nvec = NVEC, verbose = False                        # Solve with cached neighbors when available.
    )                                                                                                                           # Unpack approximate solution and neighbor list.
    comp_time = time.time() - start_time                                                                                        # Compute execution duration.
    
    # 4. Exact Solution and Metrics
    T_arr = np.linspace(0, 1, t)                                                                                                # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                                                # Initialize exact solution matrix.
    for k in range(t):                                                                                                          # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T_arr[k], [v])                                                                         # Compute exact theoretical solution.
        
    metrics = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)                                      # Compute comprehensive transient error metrics.
    
    if verbose:                                                                                                                 # Check if verbosity is enabled.
        logger.info(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                                        # Print average error for quick inspection.

    # 5. Output persistence
    if vec0 is None:                                                                                                            # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                                                   # Save vec to the canonical cache file.

    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                                        # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                                                       # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                                                      # Write structured metrics as JSON.

    # 6. VTK and Graphical rendering
    if save:                                                                                                                    # Save solution to VTK format if requested.
        ExportVTK.export_transient_vtk(p, u_ap, u_ex, t, T_arr, out_dir, basename="Heat_Solution", cloud_path=cloud_path, verbose=verbose) # Save VTK time series to disk.
        plot_transient(p, u_ap, save=True, nom=os.path.join(out_dir, 'Heat_Approximation'), title='Transient Approximation', verbose=verbose) # Save transient animation.

def main() -> None:
    """
    main
    Entry point for the Heat batch script.
    """
    Save: bool = True                                                                                                           # Choose whether VTK outputs must be saved.
    Verbose: bool = True                                                                                                        # Choose whether prints should be visible.

    if Verbose:                                                                                                                 # Check if verbosity is enabled.
        logger.info(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                             # Print batch discovery info.
        
    found: int = 0                                                                                                              # Counter to detect empty runs.
    for dataset, scale, cloud_path in iter_clouds(DATA_ROOT, SCALES):                                                           # Iterate all discovered cloud CSVs.
        found += 1                                                                                                              # Count discovered inputs.
        process_cloud(dataset, scale, cloud_path, RESULTS_ROOT, Save, verbose=Verbose)                                          # Run one case and write outputs.
        
    if found == 0:                                                                                                              # Provide a clear message when no inputs are found.
        if Verbose:                                                                                                             # Check if verbosity is enabled.
            logger.warning(f'No point clouds found under {DATA_ROOT} for scales={SCALES}.')                                     # Report empty discovery outcome.

if __name__ == "__main__":
    main()
