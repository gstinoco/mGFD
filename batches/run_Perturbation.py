"""
run_Perturbation — Stationary advection-dominated perturbation reference

Overview:
    This batch solves a stationary, advection-dominated perturbation problem (a Poisson-like operator
    with advection terms and very small diffusion) on each point cloud under Data/ (both Clouds and
    Holes datasets). It uses the meshless mGFD stationary solver with Adv=True.

Workflow:
    - Discover input clouds under Data/*/(2x|3x)/*.csv
    - Load point clouds into the (m, 3) format [x, y, flag]
    - Load tagged neighbor caches (the cache depends on the advection parameters)
    - Validate the neighbor cache for interior nodes and recompute if insufficient neighbors exist
    - Solve the stationary problem with Stationary(..., Adv=True)
    - Compute error metrics and write outputs to Results/

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
from mGFD import Stationary                                                                             # Stationary solver used for the perturbation problem.
from Scripts.IO import load_points, iter_clouds, load_neighbors, save_neighbors                         # Dataset loading + neighbor cache helpers.


def process_cloud(dataset, scale, variant, cloud_path, results_path, save):                             # Run one cloud case and write outputs to Results/.
    """
    process_cloud
    Run the perturbation (advection-dominated) stationary problem on a single point cloud file.

    Input:
        dataset                     str             Dataset folder name under Data/ (e.g., 'Clouds', 'Holes').
        scale                       str             Cloud scale folder (e.g., '2x', '3x').
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
    vec0 = load_neighbors(cloud_path, NVEC, tag = TAG)                                                  # Load tagged cached neighbor list if present.
    if vec0 is not None:                                                                                # Validate cache consistency for interior nodes.
        inne_n = p[:, 2] == 0                                                                           # Interior node mask (flag==0).
        if np.any(inne_n):                                                                              # Only validate when interior nodes exist.
            counts = np.sum(vec0[inne_n] != -1, axis = 1)                                               # Count valid neighbors per interior node.
            if int(np.min(counts)) < int(NVEC):                                                         # If any interior node has insufficient neighbors, recompute.
                vec0 = None                                                                             # Force neighbor recomputation in the solver.
    u_ap, u_ex, vec = Stationary(                                                                       # Solve the stationary advection-dominated problem.
        p, phi, f, operator = L, Adv = True, vec = vec0, nvec = NVEC                                    # Use cached neighbors when valid.
    )                                                                                                   # Unpack approximate/exact solutions and neighbor list.
    if vec0 is None:                                                                                    # If we recomputed neighbors, persist them to disk.
        save_neighbors(cloud_path, NVEC, vec, tag = TAG)                                                # Save tagged neighbor cache for future runs.

    er = Errors.Cloud_Stationary(p, vec, u_ap, u_ex)                                                    # Compute per-node stationary error metric.
    print(f'\tError: {np.mean(er)}')                                                                    # Print average error for quick inspection.

    out_dir = os.path.join(results_path, 'Perturbation', dataset, scale, variant)                       # Output directory for this region.
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
SCALES = ('2x', '3x', '4x')                                                                             # Scales to process under each dataset.
NVEC = 8                                                                                                # Neighbor count used by the solver.

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
Save = False                                                                                            # Choose whether CSVs must be saved.

## Solve the problem using a meshless Generalized Finite Difference approach.
print(f'Processing point clouds from {DATA_ROOT} (scales={SCALES}).')                                   # Print batch discovery info.
for dataset, scale, variant, cloud_path in iter_clouds(DATA_ROOT, scales = SCALES):                     # Iterate all discovered cloud CSVs.
    process_cloud(dataset, scale, variant, cloud_path, RESULTS_ROOT, Save)                              # Run one case and write outputs.
