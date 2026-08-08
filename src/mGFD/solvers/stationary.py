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
import numpy as np                                                                  # Numerical arrays and math.
from scipy.sparse.linalg import spsolve                                             # Direct sparse linear solver.
from typing import Callable, Optional, Tuple

import mGFD.core.gammas as Gammas                                                     # Gammas calculation and sparse matrix builder.
import mGFD.core.neighbors as Neighbors                                               # Neighbor search routines.

def Stationary(p: np.ndarray, phi: Callable, f: Callable, operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2]]), upwind: bool = False, vec: Optional[np.ndarray] = None, nvec: int = 12) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Numerical solution of partial differential equations with no time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu = -f(x,y)
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        phi                         function        Function declared with the boundary condition.
        f                           function        Function declared with the right-hand side of the equation.
        operator                    ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
    
    Output:
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        u_ex        m               ndarray         Array with the theoretical solution.
        vec         m x o           ndarray         Array with the correspondence of the o neighbors of each node.  
    """
    # Variable initialization
    m      = len(p[:, 0])                                                           # The total number of nodes is calculated.
    u_ap   = np.zeros([m])                                                          # u_ap initialization with zeros.
    u_ex   = np.zeros([m])                                                          # u_ex initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                        # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                           # Save the inner nodes.

    # Values for the velocities form the operator
    if upwind:                                                                      # If an Upwind stencil is requested.
        a = operator[0][0] if operator.ndim == 2 else operator[0]                   # Value of the velocity on x.
        b = operator[1][0] if operator.ndim == 2 else operator[1]                   # Value of the velocity on y.

    # Boundary conditions
    u_ap[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                  # The boundary condition is assigned.
    
    ## Neighbor search for all the nodes.
    if vec is None:
        if upwind:                                                                  # If an Upwind stencil is requested.
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                 # Neighbor search with the proper routine.
        else:                                                                       # All the other cases.
            vec = Neighbors.compute_neighbors(p, nvec)                              # Neighbor search with the proper routine.

    # Computation of Gamma values
    L = operator[:-1]                                                               # The values of the differential operator are assigned.
    K = Gammas.compute_sparse_matrix(p, vec, L)                                     # K computation with the required Gammas (Sparse).
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                       # Right-hand side of the equation.
    
    # A Generalized Finite Differences Method
    un           = spsolve(K, R)                                                    # Direct sparse solve of the linear system.
    u_ap[inne_n] = un[inne_n]                                                       # Save the computed solution to the interior nodes.
    
    # Theoretical Solution
    u_ex = phi(p[:,0], p[:,1])                                                      # The theoretical solution is computed.

    return u_ap, u_ex, vec