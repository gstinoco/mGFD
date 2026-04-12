"""
run_Wave — Reference batch for the 2D Wave equation

Overview:
    This script runs the second-order transient Wave reference problem on all available point clouds under
    Data/ (both Clouds and Holes datasets), using the meshless mGFD solver.

Workflow:
    - Discover input clouds under Data/*/(2x|3x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load cached neighbor lists when available (or compute + save them)
    - Solve the PDE with TimeDerivative2 (explicit scheme by default)
    - Compute error metrics and save outputs to Results/
    - Plot/save static step snapshots (optional) and a transient animation (always)

All the codes presented below were developed by:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

With the funding of:
    Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
    Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
    Aula CIMNE-Morelia. México
    SIIIA-MATH: Soluciones de Ingeniería. México

Date:
    May, 2024.

Last Modification:
    April, 2026.
"""

## Library importation.
import os                                                                                               # Filesystem and path utilities.
import sys                                                                                              # sys.path manipulation so this script can import project modules.
import numpy as np                                                                                      # Numerical arrays and math.

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import TimeDerivative2".

import Scripts.Graph as Graph                                                                           # Plotting helpers for stationary/transient solutions.
import Scripts.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
from mGFD import TimeDerivative2                                                                        # Second-order transient solver to run the reference case.
from Scripts.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.

def process_cloud(dataset, scale, variant, cloud_path, results_path, save):                             # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Wave benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '2x', '3x').
        variant                     str             Variant label emitted by iter_clouds (e.g., 'cloud', 'cloud_exterior').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays and step plots.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}/{variant}'                                                          # Human-readable region identifier.
    print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path)                                                                         # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, NVEC)                                                             # Load cached neighbor list if present.
    u_ap, u_ex, vec = TimeDerivative2(                                                                  # Solve Wave with explicit time stepping (2nd order).
        p, f, g, t, [c], operator = L, implicit = False, lam = 1, vec = vec0, nvec = NVEC               # Solve with cached neighbors when available.
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    if vec0 is None:                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save vec to the canonical cache file.

    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)                                                     # Compute per-node transient error metric.
    print(f'\tError: {np.mean(er)}')                                                                    # Print average error for quick inspection.

    out_dir = os.path.join(results_path, 'Wave', dataset, scale, variant)                               # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists (even if save=False).
    if save:                                                                                            # Save solution CSVs and step visualization if requested.
        computed_solution_path = os.path.join(out_dir, 'Computed Solution.csv')                         # Output path for numerical solution.
        np.savetxt(computed_solution_path, u_ap, delimiter = ',', fmt = '%.8f')                         # Save computed solution (all time steps).

        theoretical_solution_path = os.path.join(out_dir, 'Theoretical Solution.csv')                   # Output path for exact/theoretical solution.
        np.savetxt(theoretical_solution_path, u_ex, delimiter = ',', fmt = '%.8f')                      # Save exact solution (all time steps).

        plot_path = os.path.join(out_dir, 'Solution')                                                   # Base path used by plotting helper.
        Graph.Cloud_Transient_Steps(p, u_ap, u_ex, nom = plot_path)                                     # Save a subset of transient snapshots.

    error_path = os.path.join(out_dir, 'Error.txt')                                                     # Output path for scalar error report.
    with open(error_path, 'w') as file:                                                                 # Open error report file.
        file.write(str(np.mean(er)))                                                                    # Write mean error as text.

    plot_path = os.path.join(out_dir, 'Solution.mp4')                                                   # Output path for transient animation.
    Graph.Cloud_Transient(p, u_ap, u_ex, save = True, nom = plot_path)                                  # Save transient animation (mp4/gif depending on system).

DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('2x', '3x', '4x')                                                                             # Scales to process under each dataset.
NVEC = 8                                                                                                # Neighbor count used by the solver.

## Problem parameters.
c = float(np.sqrt(1 / 2))                                                                               # Wave propagation coefficient.
t = 2000                                                                                                # Number of time steps.

## Functions for the problem.
f = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                                   # Initial displacement / forcing generator.
g = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))                          # Initial velocity / forcing generator.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                             # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = False                                                                                            # Choose whether CSVs and step plots must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
for dataset, scale, variant, cloud_path in iter_clouds(DATA_ROOT, scales = SCALES):                     # Iterate all discovered cloud CSVs.
    process_cloud(dataset, scale, variant, cloud_path, RESULTS_ROOT, Save)                              # Run one case and write outputs.
