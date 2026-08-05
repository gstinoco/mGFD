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
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import Stationary".

import Scripts.Graph as Graph                                                                           # Plotting helpers for stationary/transient solutions.
import Scripts.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
from mGFD import Stationary                                                                             # Stationary solver for the reference Poisson problem.
from Scripts.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.

def process_cloud(dataset, scale, variant, cloud_path, results_path, save):                             # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the Poisson benchmark on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '1.00x', '2.00x').
        variant                     str             Variant label emitted by iter_clouds (e.g., 'cloud', 'cloud_exterior').
        cloud_path                  str             Path to input CSV with point cloud.
        results_path                str             Base output directory (typically <repo>/Results).
        save                        bool            Whether to save the solution arrays.

    Output:
        None
    """
    region_id = f'{dataset}/{scale}/{variant}'                                                          # Human-readable region identifier.
    print(f'Working on region: {region_id}')                                                            # Progress message for the batch run.

    p = load_points(cloud_path)                                                                         # Load point cloud into (m, 3) array [x, y, flag].
    vec0 = load_neighbors(cloud_path, NVEC)                                                             # Load cached neighbor list if present.
    u_ap, u_ex, vec = Stationary(                                                                       # Solve the stationary Poisson problem.
        p, phi, f, operator = L, vec = vec0, nvec = NVEC                                                # Use cached neighbors when available.
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    if vec0 is None:                                                                                    # If there was no cache, persist computed neighbors.
        save_neighbors(cloud_path, NVEC, vec)                                                           # Save vec to the canonical cache file.

    er = Errors.Cloud_Stationary(p, vec, u_ap, u_ex)                                                    # Compute per-node stationary error metric.
    print(f'\tError: {np.mean(er)}')                                                                    # Print average error for quick inspection.

    out_dir = os.path.join(results_path, 'Poisson', dataset, scale, variant)                            # Output directory for this region.
    os.makedirs(out_dir, exist_ok = True)                                                               # Ensure output directory exists.
    if save:                                                                                            # Save solution arrays if requested.
        computed_solution_path = os.path.join(out_dir, 'Computed Solution.csv')                         # Output path for numerical solution.
        np.savetxt(computed_solution_path, u_ap, delimiter = ',', fmt = '%.8f')                         # Save computed solution.

        theoretical_solution_path = os.path.join(out_dir, 'Theoretical Solution.csv')                   # Output path for exact/theoretical solution.
        np.savetxt(theoretical_solution_path, u_ex, delimiter = ',', fmt = '%.8f')                      # Save exact solution.
    
    error_path = os.path.join(out_dir, 'Error.txt')                                                     # Output path for scalar error report.
    with open(error_path, 'w') as file:                                                                 # Open error report file.
        file.write(str(np.mean(er)))                                                                    # Write mean error as text.

    plot_path = os.path.join(out_dir, 'Solution')                                                       # Base path used by plotting helper.
    Graph.Cloud_Stationary(p, u_ap, u_ex, save = True, nom = plot_path)                                 # Save stationary comparison plot(s).

DATA_ROOT = os.path.join(BASE_DIR, 'Data')                                                              # Input dataset root directory.
RESULTS_ROOT = os.path.join(BASE_DIR, 'Results')                                                        # Output results root directory.
SCALES = ('0.50x', '1.00x', '1.50x')                                                                    # Scales to process under each dataset.
NVEC = 12                                                                                               # Neighbor count used by the solver.


## Functions for the problem.
phi = lambda x, y: 2 * np.exp(2 * x + y)                                                                # Boundary condition for the problem.
f = lambda x, y: 10 * np.exp(2 * x + y)                                                                 # Right-hand side forcing term.

## Operator L = [D, E, A, B, C, F].
L = np.vstack([[0], [0], [2], [0], [2], [0]])                                                           # Operator coefficients for Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + Fu.

## Save policy.
Save = False                                                                                            # Choose whether CSVs must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
found = 0                                                                                               # Counter to detect empty runs.
for dataset, scale, variant, cloud_path in iter_clouds(DATA_ROOT, SCALES):                              # Iterate all discovered cloud CSVs.
    found += 1                                                                                          # Count discovered inputs.
    process_cloud(dataset, scale, variant, cloud_path, RESULTS_ROOT, Save)                              # Run one case and write outputs.
if found == 0:                                                                                          # Provide a clear message when no inputs are found.
    print(f'No point clouds found under {DATA_ROOT} for scales={SCALES}.')                              # Report empty discovery outcome.
