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

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import Stationary".

import json                                                                                             # JSON serialization for metrics.

import mGFD.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
import mGFD.ExportVTK as ExportVTK                                                                   # VTK export utilities for ParaView.
from mGFD import Stationary                                                                             # Core solver to run the reference case.Poisson problem.
from mGFD.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.

def process_cloud(dataset, scale, cloud_path, results_path, save):                                      # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Poisson benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1.00x', '2.00x').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}'                                                                    # Region identifier.
    print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path)                                                                         # Load point cloud into (m, 3) array [x, y, flag].

    vec0 = load_neighbors(cloud_path, NVEC)                                                             # Load cached neighbor list if present.
    
    start_time = time.time()                                                                            # Start execution timer.
    u_ap, u_ex, vec = Stationary(                                                                       # Solve the stationary Poisson problem.
        p, phi, f, operator = L, vec = vec0, nvec = NVEC                                                # Use cached neighbors when available.
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    comp_time = time.time() - start_time                                                                # Compute execution duration.
    
    if vec0 is None:                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save vec to the canonical cache file.

    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time=comp_time)             # Compute comprehensive stationary error metrics.
    print(f'\tError (RMSE): {metrics["RMSE"]}')                                                         # Print RMSE error for quick inspection.

    out_dir = os.path.join(results_path, 'Poisson', dataset, scale)                                     # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists.
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                               # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                              # Write structured metrics as JSON.

    if save:                                                                                            # Save solution to VTK format if requested.
        ExportVTK.export_stationary_vtk(p, u_ap, u_ex, out_dir, basename="Poisson_Solution", cloud_path=cloud_path) # Save VTK data to disk.

DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('0.25x', '0.50x', '1.00x', '1.50x', '2.00x', '3.00x')                                         # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Functions for the problem.
phi = lambda x, y: 2 * np.exp(2 * x + y)                                                                # Boundary condition for the problem.
f = lambda x, y: 10 * np.exp(2 * x + y)                                                                 # Right-hand side forcing term.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[0], [0], [2], [0], [2], [0]])                                                           # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = True                                                                                             # Choose whether VTK outputs must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
found = 0                                                                                               # Counter to detect empty runs.
for dataset, scale, cloud_path in iter_clouds(DATA_ROOT, SCALES):                                       # Iterate all discovered cloud CSVs.
    found += 1                                                                                          # Count discovered inputs.
    process_cloud(dataset, scale, cloud_path, RESULTS_ROOT, Save)                                       # Run one case and write outputs.
if found == 0:                                                                                          # Provide a clear message when no inputs are found.
    print(f'No point clouds found under {DATA_ROOT} for scales={SCALES}.')                              # Report empty discovery outcome.
