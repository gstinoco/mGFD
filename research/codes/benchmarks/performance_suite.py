"""
performance_suite.py

Overview:
    Benchmarks the mGFD PDE solvers against different algebraic backends:
    - SciPy spsolve (Direct LU)
    - SciPy bicgstab (Iterative Krylov)
    - SciPy gmres (Iterative Krylov)
    
    Exploits the `compute_time` metric built natively into the v0.10.0 `SolverResult`.

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
    August, 2026.
Last Modification:
    August, 2026.
"""

import os                                                                                                                               # OS module for path manipulation.
import sys                                                                                                                              # Sys module for path injection.
import numpy as np                                                                                                                      # Numerical arrays and math.
import pandas as pd                                                                                                                     # Dataframes handling.

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))                                            # Ensure mGFD is loadable from src.
from mGFD import Stationary                                                                                                             # Stationary PDE solver.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.

def phi_func(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Boundary conditions (Dirichlet).
    """
    return np.sin(np.pi * x) * np.sin(np.pi * y)                                                                                        # Exact analytical value for the solution.

def f_func(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Right-hand side forcing term.
    """
    return -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)                                                                      # Evaluated RHS forcing term.

def main() -> None:
    """
    Main execution loop for benchmarking performance.
    """
    print("=" * 60)                                                                                                                     # Print visual header separator.
    print(" mGFD Performance Suite: Linear Solvers Benchmark ")                                                                         # Print script title.
    print("=" * 60)                                                                                                                     # Print visual header separator.
    
    densities = ['0.50x', '1.00x']                                                                                                      # We could test more dense clouds later.
    solvers   = ['spsolve', 'bicgstab', 'gmres']                                                                                        # Backends to evaluate in the loop.
    operator  = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                               # Laplacian matrix format.
    
    for density in densities:                                                                                                           # Iterate over the available point cloud densities.
        # 1. Load a test point cloud
        cloud_path = os.path.abspath(os.path.join(os.path.dirname(__file__), f'../../Data/Chapala/1/Chapala_cloud.csv'))                # Resolve relative path to Chapala cloud.
        if not os.path.exists(cloud_path):                                                                                              # Check if the cloud exists.
            print(f"Skipping {density} (file not found: {cloud_path})")                                                                 # Alert user if cloud missing.
            continue                                                                                                                    # Jump to the next density profile.
            
        p = load_points(cloud_path)                                                                                                     # Load spatial node coordinates and flags.
        print(f"\n[Cloud: Chapala {density}] (Nodes: {p.shape[0]})")                                                                    # Log cloud metadata.
        print("-" * 50)                                                                                                                 # Sub-header separation.
        
        # 2. Precompute arrays for maximum speed
        phi_arr = phi_func(p[:, 0], p[:, 1])                                                                                            # Precompute Dirichlet boundaries numerically.
        f_arr   = f_func(p[:, 0], p[:, 1])                                                                                              # Precompute RHS forcing numerically.
        
        # 3. Warm-up / Pre-cache neighbors
        print("  Warming up neighbors...")                                                                                              # Log initialization phase.
        res_warm  = Stationary(p, phi_arr, f_arr, operator=operator, verbose=False)                                                     # Execute dummy pass to generate kdtree.
        vec_cache = res_warm.neighbors                                                                                                  # Extract neighbors cache from dataclass.
        
        # 4. Benchmark loop
        for solver in solvers:                                                                                                          # Iterate over the different algebraic backends.
            try:                                                                                                                        # Use a safe wrapper to catch algebraic decomposition errors.
                res = Stationary(p, phi_arr, f_arr, operator=operator, vec=vec_cache, linear_solver=solver, verbose=False)              # Execute stationary step.
                if res.converged:                                                                                                       # Check if convergence flag is True.
                    print(f"  > {solver.ljust(10)}: {res.compute_time:.4f} seconds")                                                    # Print internal timing metric.
                else:                                                                                                                   # If divergence or error occurred.
                    print(f"  > {solver.ljust(10)}: FAILED to converge")                                                                # Alert solver crash.
            except Exception as e:                                                                                                      # Catch any lower-level SciPy exceptions.
                print(f"  > {solver.ljust(10)}: ERROR ({e})")                                                                           # Dump exception to console.

if __name__ == "__main__":
    main()                                                                                                                              # Run main flow.
