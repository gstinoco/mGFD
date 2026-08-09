"""
Stationary — Stationary PDEs solver

Overview:
    Numerical solver for stationary PDEs (no time derivatives) using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    Stationary          Main solver function for stationary problems.

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
from typing import Callable, Optional, Tuple                                                                                            # Type hinting.

import mGFD.core.gammas as Gammas                                                                                                       # Gammas calculation and sparse matrix builder.
import mGFD.core.neighbors as Neighbors                                                                                                 # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def Stationary(p: np.ndarray, phi: Callable, f: Callable, operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2]]), upwind: bool = False, vec: Optional[np.ndarray] = None, nvec: int = 12, verbose: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Numerical solution of partial differential equations with no time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu = -f(x,y)
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        phi                         Callable        Function declared with the boundary condition.
        f                           Callable        Function declared with the right-hand side of the equation.
        operator    6 x 1           ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
        verbose                     bool            If True, prints solver progress.
    
    Output:
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        vec         m x nvec        ndarray         Array with the correspondence of the nvec neighbors of each node.  
    """
    # 0. Input validation
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise ValueError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                       # Raise explicit error on bad input.
    if not callable(phi):                                                                                                               # Validate boundary condition function.
        raise TypeError("Boundary condition 'phi' must be a callable function.")                                                        # Raise explicit error on bad input.
    if not callable(f):                                                                                                                 # Validate RHS function.
        raise TypeError("Right-hand side 'f' must be a callable function.")                                                             # Raise explicit error on bad input.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise ValueError("Operator must be a numpy array with at least 5 coefficients.")                                                # Raise explicit error on bad input.

    # 1. Variable initialization
    m      = len(p[:, 0])                                                                                                               # The total number of nodes is calculated.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f"Solving Stationary problem for {m} nodes...")                                                                     # Print solver progress.
        
    u_ap   = np.zeros([m])                                                                                                              # u_ap initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Save the inner nodes.

    # 2. Extract advection velocities for Upwind scheme.
    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = operator[0][0] if operator.ndim == 2 else operator[0]                                                                       # Value of the velocity on x.
        b = operator[1][0] if operator.ndim == 2 else operator[1]                                                                       # Value of the velocity on y.

    # 3. Apply Boundary Conditions
    u_ap[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                                                                      # The boundary condition is assigned.
    
    # 4. Neighbor search for all the nodes.
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind:                                                                                                                      # If an Upwind stencil is requested.
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                                     # Neighbor search with the proper routine.
        else:                                                                                                                           # All the other cases.
            vec = Neighbors.compute_neighbors(p, nvec)                                                                                  # Neighbor search with the proper routine.

    # 5. Computation of Gamma values
    L = operator[:-1]                                                                                                                   # The values of the differential operator are assigned.
    K = Gammas.compute_sparse_matrix(p, vec, L)                                                                                         # K computation with the required Gammas (Sparse).
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                                                                           # Right-hand side of the equation.
    
    # 6. Solution of the linear system (Generalized Finite Differences)
    un           = spsolve(K, R)                                                                                                        # Direct sparse solve of the linear system.
    u_ap[inne_n] = un[inne_n]                                                                                                           # Save the computed solution to the interior nodes.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info("\tSolver finished successfully.")                                                                                  # Print completion message.
        
    return u_ap, vec                                                                                                                    # Return computed values.
