"""
run_AdvDif — Reference batch for the 2D Advection–Diffusion equation

Overview:
    This script runs the first-order transient Advection–Diffusion reference problem on all available
    point clouds under Data/ (both Clouds and Holes datasets), using the meshless mGFD solver.

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


import mGFD.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
import mGFD.ExportVTK as ExportVTK                                                                   # VTK export utilities for ParaView.
from mGFD import TimeDerivative1                                                                        # First-order transient solver to run the reference case.
from mGFD.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.


def process_cloud(dataset, scale, cloud_path, results_path, save):                                      # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Advection–Diffusion benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1.00x', '2.00x').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays and step plots.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}'                                                                    # Region identifier.
    print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path)                                                                         # Load point cloud into (m, 3) array [x, y, flag].

    vec0 = load_neighbors(cloud_path, NVEC)                                                             # Load cached neighbor list if present.
    
    start_time = time.time()                                                                            # Start execution timer.
    u_ap, u_ex, vec = TimeDerivative1(                                                                  # Solve Advection-Diffusion with implicit time stepping.
        p, f, t, [v, a, b], operator = L, implicit = True, lam = 0.5, upwind = True, vec = vec0, nvec = NVEC
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    comp_time = time.time() - start_time                                                                # Compute execution duration.
    
    if vec0 is None:                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save vec to the canonical cache file.

    metrics = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)              # Compute comprehensive transient error metrics.
    print(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                          # Print average error for quick inspection.

    out_dir = os.path.join(results_path, 'Advection-Diffusion', dataset, scale)                         # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists (even if save=False).
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                               # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                              # Write structured metrics as JSON.

    if save:                                                                                            # Save solution to VTK format if requested.
        T = np.linspace(0, 1, t)                                                                        # Reconstruct time vector.
        ExportVTK.export_transient_vtk(p, u_ap, u_ex, t, T, out_dir, basename="AdvDif_Solution", cloud_path=cloud_path) # Save VTK time series to disk.


DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('0.25x', '0.50x', '1.00x', '1.50x')                                                           # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Problem parameters.
v = 0.1                                                                                                 # Diffusion coefficient.
a = 0.3                                                                                                 # Transport velocity on the x direction.
b = 0.2                                                                                                 # Transport velocity on the y direction.
t = 2000                                                                                                # Number of time steps.

## Functions for the problem.
f = lambda x, y, t, coef: (1 / (4 * t + 1)) * np.exp(                                                   
    - (x - coef[1] * t - 0.5)**2 / (coef[0] * (4 * t + 1)) - (y - coef[2] * t - 0.5)**2 / (coef[0] * (4 * t + 1))
)                                                                                                       # Advection–diffusion analytical solution.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[-a], [-b], [2 * v], [0], [2 * v], [0]])                                                 # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = True                                                                                             # Choose whether VTK outputs must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
for dataset, scale, cloud_path in iter_clouds(DATA_ROOT, scales = SCALES):                              # Iterate all discovered cloud CSVs.
    process_cloud(dataset, scale, cloud_path, RESULTS_ROOT, Save)                                       # Run one case and write outputs.
