"""
run_Perturbation — Stationary advection-dominated perturbation reference

Overview:
    This batch solves a stationary, advection-dominated perturbation problem (a Poisson-like operator
    with advection terms and very small diffusion) on each point cloud under Data/ (both Clouds and
    Holes datasets). It uses the meshless mGFD stationary solver with Adv=True.

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load tagged neighbor caches (the cache depends on the advection parameters)
    - Validate the neighbor cache for interior nodes and recompute if insufficient neighbors exist
    - Solve the stationary problem with Stationary(..., Adv=True)
    - Compute error metrics and write outputs to Results/

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
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import Stationary".

import metrics as Errors                                                                         # Error metrics for stationary/transient runs.
import mGFD.io.export_vtk as ExportVTK                                                                   # VTK export utilities for ParaView.
from mGFD import Stationary                                                                             # Core solver to run the reference case.rturbation problem.
from mGFD.viz.graph import plot_stationary
from mGFD.io.io import load_points
from batch_utils import iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.


def process_cloud(dataset, scale, cloud_path, results_path, save, verbose=True):                      # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the perturbation stationary problem (upwind=True) on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1.00x', '2.00x').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays.
        verbose                     bool            If True, prints progress and errors to console.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}'                                                                    # Region identifier.
    if verbose:
        print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path, verbose=verbose)                                                                         # Load point cloud into (m, 3) array [x, y, flag].

    vec0 = load_neighbors(cloud_path, NVEC, tag = TAG)                                                  # Load tagged cached neighbor list if present.
    
    start_time = time.time()                                                                            # Start execution timer.
    u_ap, vec = Stationary(                                                                             # Solve the stationary advection-dominated problem.
        p, phi, f, operator = L, upwind = True, vec = vec0, nvec = NVEC                                 # Use cached neighbors and Upwind scheme.
    )                                                                                                   # Unpack approximate solution and neighbor list.
    comp_time = time.time() - start_time                                                                # Compute execution duration.
    u_ex = phi(p[:, 0], p[:, 1])                                                                        # Compute exact theoretical solution.
    
    if vec0 is None:                                                                                    # If we recomputed neighbors, persist them to disk.
        save_neighbors(cloud_path, NVEC, vec, tag = TAG)                                                # Save tagged neighbor cache for future runs.

    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time=comp_time)             # Compute comprehensive stationary error metrics.
    if verbose:
        print(f'\tError (RMSE): {metrics["RMSE"]}')                                                         # Print RMSE error for quick inspection.

    out_dir = os.path.join(results_path, 'Perturbation', dataset, scale)                                # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists.
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                               # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                              # Write structured metrics as JSON.

    if save:                                                                                            # Save solution to VTK format if requested.
        ExportVTK.export_stationary_vtk(p, u_ap, u_ex, out_dir, basename="Perturbation_Solution", cloud_path=cloud_path, verbose=verbose) # Save VTK data to disk.
        
        plot_stationary(p, u_ap, save=True, nom=os.path.join(out_dir, 'Perturbation_Approximation'), title='Stationary Approximation', verbose=verbose)
        
        if scale == '5':                                                                                # Only for scale 5, create an independent plot of the exact solution.
            plot_stationary(p, u_ex, save=True, nom=os.path.join(out_dir, 'Perturbation_Exact'), title='Theoretical Solution', verbose=verbose)


DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('1', '2', '3', '4', '5')                                         # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Variables for the problem.
vx = 1                                                                                                  # Advection velocity in x.
vy = 0                                                                                                  # Advection velocity in y.
D = 1e-5                                                                                                # Diffusion coefficient (very small -> boundary layers).
TAG = f'adv_vx{vx:g}_vy{vy:g}'                                                                          # Neighbor-cache tag tied to (vx, vy).

## Functions for the problem.
phi = lambda x, y: (x**2 - np.exp(-(1 - x) / D)) * y * (1 - y)                                          # Boundary condition for the problem.
f = lambda x, y: -(1e-5) * (2 * np.exp((x - 1) / (1e-5)) - 2 * x**2 + y * (np.exp((x - 1) / (1e-5)) / (1e-5)**2 - 2) * (y - 1)) - y * (2 * x - np.exp((x - 1) / (1e-5)) / (1e-5)) * (y - 1)
                                                                                                        # Right-hand side forcing term.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[vx], [vy], [-D], [0], [-D], [0]])                                                       # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

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
