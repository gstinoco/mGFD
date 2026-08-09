"""
run_Heat — Reference batch for the 2D Heat equation

Overview:
    This script runs the first-order transient Heat reference problem on all available point clouds under
    Data/ (both Clouds and Holes datasets), using the meshless mGFD solver.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with TimeDerivative1 (explicit scheme by default)
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
import os                                                                                               # Filesystem and path utilities.
import sys                                                                                              # sys.path manipulation so this script can import project modules.
import time                                                                                             # Time tracking for execution performance.
import numpy as np                                                                                      # Numerical arrays and math.
import json                                                                                             # JSON serialization for metrics.

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import TimeDerivative1".

import metrics as Errors                                                                         # Error metrics for stationary/transient runs.
import mGFD.io.export_vtk as ExportVTK                                                                   # VTK export utilities for ParaView.
from mGFD import TimeDerivative1                                                                        # First-order transient solver to run the reference case.
from mGFD.viz.graph import plot_transient
from mGFD.io.io import load_points
from batch_utils import iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.


def process_cloud(dataset, scale, cloud_path, results_path, save, verbose=True):                      # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Heat benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1.00x', '2.00x').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays and step plots.
        verbose                     bool            If True, prints progress and errors to console.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}'                                                                    # Region identifier.
    if verbose:
        print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path, verbose=verbose)                                                                         # Load point cloud into (m, 3) array [x, y, flag].

    vec0 = load_neighbors(cloud_path, NVEC)                                                             # Load cached neighbor list if present.
    
    start_time = time.time()                                                                            # Start execution timer.
    u_ap, vec = TimeDerivative1(                                                                        # Solve Heat with implicit time stepping.
        p, f, t, [v], operator = L, implicit = True, lam = 0.5, vec = vec0, nvec = NVEC, verbose = verbose              # Solve with cached neighbors when available.
    )                                                                                                   # Unpack approximate solution and neighbor list.
    comp_time = time.time() - start_time                                                                # Compute execution duration.
    
    T = np.linspace(0, 1, t)                                                                            # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                        # Initialize exact solution matrix.
    for k in range(t):                                                                                  # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [v])                                                     # Compute exact theoretical solution.
    
    if vec0 is None:                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save vec to the canonical cache file.

    metrics = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)              # Compute comprehensive transient error metrics.
    if verbose:
        print(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                          # Print average error for quick inspection.

    out_dir = os.path.join(results_path, 'Heat', dataset, scale)                                        # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists (even if save=False).
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                               # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                              # Write structured metrics as JSON.

    if save:                                                                                            # Save solution to VTK format if requested.
        T = np.linspace(0, 1, t)                                                                        # Reconstruct time vector.
        ExportVTK.export_transient_vtk(p, u_ap, u_ex, t, T, out_dir, basename="Heat_Solution", cloud_path=cloud_path, verbose=verbose) # Save VTK time series to disk.
        
        plot_transient(p, u_ap, save=True, nom=os.path.join(out_dir, 'Heat_Approximation'), title='Transient Approximation', verbose=verbose)


DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('1', '2', '3', '4', '5')                                                           # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Problem parameters.
v = 0.2                                                                                                 # Diffusion coefficient.
t = 2000                                                                                                # Number of time steps.

## Functions for the problem.
f = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)   # Heat analytical solution.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                   # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = True                                                                                             # Choose whether VTK outputs must be saved.
Verbose = True                                                                                          # Choose whether prints should be visible.

## Solve the problem using a meshless Generalized Finite Difference approach.
if Verbose:
    print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
found = 0                                                                                               # Counter to detect empty runs.
for dataset, scale, cloud_path in iter_clouds(DATA_ROOT, SCALES):                              # Iterate all discovered cloud CSVs.
    found += 1
    process_cloud(dataset, scale, cloud_path, RESULTS_ROOT, Save, verbose=Verbose)                      # Run one case and write outputs.
if found == 0:
    if Verbose:
        print(f'No point clouds found under {DATA_ROOT} for scales={SCALES}.')
