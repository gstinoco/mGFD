"""
Stationary — CPU Backend for Stationary PDEs

Overview:
    CPU implementation for solving stationary PDEs using SciPy sparse solvers (SuperLU direct solver).

Public API:
    solve_cpu                   Core CPU execution routine for the stationary solver.

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
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Core numerical operations.

from scipy.sparse.linalg import spsolve                                                                                                 # Direct sparse linear solver.
from typing import Callable, Optional, Tuple, Union                                                                                     # Type hinting.

import mGFD.spatial.gammas as Gammas                                                                                                    # Gammas calculation and sparse matrix builder.
import mGFD.spatial.neighbors as Neighbors                                                                                              # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def solve_cpu(p: np.ndarray,                                                                                                            # Function definition.
              phi: Union[Callable, np.ndarray, float, int],
              f: Union[Callable, np.ndarray, float, int],
              operator: np.ndarray,
              upwind: bool,
              vec: Optional[np.ndarray],
              nvec: int,
              verbose: bool) -> Tuple[np.ndarray, np.ndarray, bool]:
    """CPU backend for Stationary solver using direct sparse LU factorization."""
    
    m = len(p[:, 0])                                                                                                                    # Total nodes.
    if verbose:                                                                                                                         # Verbosity.
        logger.info(f"Solving Stationary problem for {m} nodes on CPU...")                                                              # Progress info.
    
    u_ap = np.zeros([m])                                                                                                                # u_ap init.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Inner nodes.

    if upwind:                                                                                                                          # Upwind check.
        a = operator[0][0] if operator.ndim == 2 else operator[0]                                                                       # Vel x.
        b = operator[1][0] if operator.ndim == 2 else operator[1]                                                                       # Vel y.

    if callable(phi):                                                                                                                   # If boundary condition is a function.
        u_ap[boun_n] = np.asarray(phi(p[boun_n, 0], p[boun_n, 1]))                                                                      # Apply function.
    elif isinstance(phi, np.ndarray):                                                                                                   # Array.
        u_ap[boun_n] = phi[boun_n]                                                                                                      # Apply array.
    elif isinstance(phi, (int, float)):                                                                                                 # Constant.
        u_ap[boun_n] = phi                                                                                                              # Apply constant.

    if vec is None:                                                                                                                     # No vec provided.
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                              # Upwind neighbors.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Central neighbors.

    L = operator[:-1]                                                                                                                   # Extracted operator.
    K = Gammas.compute_sparse_matrix(p, vec, L)                                                                                         # K sparse matrix construction.
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                                                                           # Right-hand side vector.
    
    un = spsolve(K, R)                                                                                                                  # SciPy direct spsolve.
    u_ap[inne_n] = un[inne_n]                                                                                                           # Unpack to interior.
    
    if verbose: logger.info("\tCPU Solver finished successfully.")                                                                      # Success.
    
    return u_ap, vec, True                                                                                                              # Return core data.
