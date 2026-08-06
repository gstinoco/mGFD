"""
run_Perturbation2 — Stationary perturbation reference (non-advection mode)

Overview:
    This batch solves the same perturbation-style stationary reference used in run_Perturbation.py,
    but configures the stationary solver with Adv=False. It runs on each point cloud under Data/
    (both Clouds and Holes datasets).

Workflow:
    - Discover input clouds under Data/*/(0.50x|1.00x|1.50x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load neighbor caches when available (or compute + save them)
    - Validate the neighbor cache for interior nodes and recompute if insufficient neighbors exist
    - Solve the stationary problem with Stationary(..., Adv=False)
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
import numpy as np                                                                                      # Numerical arrays and math.
import json                                                                                             # JSON serialization for metrics.

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import Stationary".

import Scripts.Graph as Graph                                                                           # Plotting helpers for stationary/transient solutions.
import Scripts.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
import Scripts.ExportVTK as ExportVTK                                                                   # VTK export utilities for ParaView.
from mGFD import Stationary                                                                             # Core solver to run the reference case.
from Scripts.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.

def process_cloud(dataset, scale, cloud_path, results_path, save):                             # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the perturbation stationary problem (Adv=False) on a single point cloud file.

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
    if vec0 is not None:                                                                                # Validate cache consistency for interior nodes.
        inne_n = p[:, 2] == 0                                                                           # Interior node mask (flag==0).
        if np.any(inne_n):                                                                              # Only validate when interior nodes exist.
            counts = np.sum(vec0[inne_n] != -1, axis = 1)                                               # Count valid neighbors per interior node.
            if int(np.min(counts)) < int(NVEC):                                                         # If any interior node has insufficient neighbors, recompute.
                vec0 = None                                                                             # Force neighbor recomputation in the solver.
    u_ap, u_ex, vec = Stationary(                                                                       # Solve the stationary perturbation problem.
        p, phi, f, operator = L, Adv = False, vec = vec0, nvec = NVEC                                   # Use cached neighbors when valid.
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    if vec0 is None:                                                                                    # If we recomputed neighbors, persist them to disk.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save neighbor cache for future runs.

    metrics = Errors.Compute_Metrics_Stationary(p, vec, u_ap, u_ex)                                     # Compute comprehensive stationary error metrics.
    print(f'\tError (RMSE): {metrics["RMSE"]}')                                                         # Print RMSE error for quick inspection.

    out_dir = os.path.join(results_path, 'Perturbation2', dataset, scale)                               # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists.
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                # Output path for JSON metrics report.
    with open(metrics_path, 'w') as file:                                                               # Open metrics report file.
        json.dump(metrics, file, indent=4)                                                              # Write structured metrics as JSON.

    if save:                                                                                            # Save solution to VTK format if requested.
        ExportVTK.export_stationary_vtk(p, u_ap, u_ex, out_dir, basename="Perturbation2_Solution")      # Save VTK data to disk.

DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('0.25x', '0.50x', '1.00x', '1.50x', '2.00x', '3.00x')                                         # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Variables for the problem.
vx = 1                                                                                                  # Advection velocity in x (used inside operator L only).
vy = 0                                                                                                  # Advection velocity in y (used inside operator L only).
D = 1e-5                                                                                                # Diffusion coefficient (very small -> boundary layers).

## Functions for the problem.
phi = lambda x, y: (x**2 - np.exp(-(1 - x) / D)) * y * (1 - y)                                          # Boundary condition for the problem.
f = lambda x, y: -(1e-5) * (2 * np.exp((x - 1) / (1e-5)) - 2 * x**2 + y * (np.exp((x - 1) / (1e-5)) / (1e-5)**2 - 2) * (y - 1)) - y * (2 * x - np.exp((x - 1) / (1e-5)) / (1e-5)) * (y - 1)
                                                                                                        # Right-hand side forcing term.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[vx], [vy], [-D], [0], [-D], [0]])                                                       # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = True                                                                                             # Choose whether VTK outputs must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
for dataset, scale, cloud_path in iter_clouds(DATA_ROOT, scales = SCALES):                              # Iterate all discovered cloud CSVs.
    process_cloud(dataset, scale, cloud_path, RESULTS_ROOT, Save)                                       # Run one case and write outputs.
