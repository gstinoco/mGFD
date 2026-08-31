"""
TimeDerivative2 — CPU Backend for Second-order Transient PDEs

Overview:
    CPU implementation for solving hyperbolic PDEs (like the Wave equation) using pre-factorized direct sparse solvers.

Public API:
    solve_cpu                   Core CPU execution routine for the hyperbolic solver.

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

from scipy.sparse import diags, eye                                                                                                     # SciPy sparse matrices.
from scipy.sparse.linalg import factorized                                                                                              # Direct pre-factorized sparse LU solver.
from typing import Callable, Optional, Tuple, List, Union                                                                               # Type hinting.

import mGFD.spatial.gammas as Gammas                                                                                                    # Gammas calculation and sparse matrix builder.
import mGFD.spatial.neighbors as Neighbors                                                                                              # Neighbor search routines.

from mGFD.exceptions import DimensionMismatchError                                                                                      # Custom exceptions.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def solve_cpu(p: np.ndarray,                                                                                                            # Function definition.
              f: Union[Callable, np.ndarray, float, int],
              g: Union[Callable, np.ndarray, float, int],
              t: int,
              coef: List[float],
              operator: np.ndarray,
              upwind: bool,
              vec: Optional[np.ndarray],
              nvec: int,
              implicit: bool,
              lam: float,
              verbose: bool) -> Tuple[np.ndarray, np.ndarray, bool]:
    """CPU backend for TimeDerivative2 using direct pre-factorized sparse LU solver."""

    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    if verbose: logger.info(f"Solving Wave/Hyperbolic problem ({t} steps) for {m} nodes on CPU...")
        
    T      = np.linspace(0, 1, t)                                                                                                       # Time discretization array.
    dt     = T[1] - T[0]                                                                                                                # Time step size.
    u_ap   = np.zeros([m, t])                                                                                                           # Numerical approximation matrix.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Boolean mask for boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Boolean mask for interior nodes.

    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]                                                                     # X-velocity (D coefficient).
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]                                                                     # Y-velocity (E coefficient).
    
    if callable(f):                                                                                                                     # If data is a function.
        for k in range(t):                                                                                                              # Loop over time.
            u_ap[boun_n, k] = np.asarray(f(p[boun_n, 0], p[boun_n, 1], T[k], coef))                                                     # Boundary condition (Dirichlet).
        u_ap[:, 0] = np.asarray(f(p[:, 0], p[:, 1], T[0], coef))                                                                        # Initial condition across all nodes.
    elif isinstance(f, np.ndarray):                                                                                                     # If data is an array.
        if f.ndim == 2 and f.shape == (m, t):                                                                                           # Spatiotemporal data array.
            u_ap[boun_n, :] = f[boun_n, :]                                                                                              # Spatiotemporal boundary conditions.
            u_ap[:, 0] = f[:, 0]                                                                                                        # Initial conditions.
        elif f.ndim == 1 and f.shape[0] == m:                                                                                           # Spatial constant data array.
            for k in range(t): u_ap[boun_n, k] = f[boun_n]
            u_ap[:, 0] = f                                                                                                              # Constant initial condition.
        if isinstance(f, np.ndarray) and f.shape not in [(m, t), (m,)]:                                                                 # Check numeric matrix shape matches dimensions.
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")                                       # Raise explicit error.
    elif isinstance(f, (int, float)):                                                                                                   # If data is a constant scalar.
        u_ap[boun_n, :] = f                                                                                                             # Constant boundary condition.
        u_ap[:, 0] = f                                                                                                                  # Constant initial condition.
    
    if callable(g):                                                                                                                     # If initial velocity is a function.
        g_eval = np.asarray(g(p[:, 0], p[:, 1], T[0], coef))                                                                            # Evaluate initial velocity function.
    elif isinstance(g, np.ndarray):                                                                                                     # If data is an array.
        if g.ndim == 1 and g.shape[0] == m:                                                                                             # Check correct spatial dimensions.
            g_eval = g                                                                                                                  # Use provided initial velocity array.
        else:
            raise DimensionMismatchError(f"Data array 'g' must have shape ({m},).")                                                     # Raise explicit error.
    else:
        g_eval = np.full(m, float(g))                                                                                                   # Use scalar value.
    
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                              # Upwind-biased neighbor selection.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Standard distance-based neighbors.

    L         = operator[:-1]                                                                                                           # Original operator weights.
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                                 # Build sparse spatial differentiation matrix.
    K         = (dt**2) * K_spatial                                                                                                     # Scale by time step squared.
    
    if not implicit:                                                                                                                    # Explicit scheme.
        K2 = eye(m) + (1/2)*K                                                                                                           # Explicit matrix for k = 1.
        K4 = 2*eye(m) + K                                                                                                               # Explicit matrix for k = 2, ..., t.
        
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step (k=1).
                g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                      # Evaluate velocity for first step.
                un              = K2.dot(u_ap[:, k - 1]) + dt*g_val                                                                     # New time level (k=1).
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
            else:
                un              = K4.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]                                                               # New time level (k>=2).
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
    else:                                                                                                                               # Implicit scheme.
        Id_inner = diags(inne_n.astype(float))                                                                                          # Diagonal mask for inner nodes.
        Id_bound = diags(boun_n.astype(float))                                                                                          # Diagonal mask for boundary nodes.
        
        A1       = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                             # LHS Matrix for k=1.
        B1       = Id_inner @ (eye(m) + (0.5 - lam) * K)                                                                                # RHS Matrix for k=1.
        
        A2       = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                             # LHS Matrix for k>=2.
        A2       = A2.tocsc()                                                                                                           # CSC format for SuperLU factorization.
        B2       = Id_inner @ (2*eye(m) + (1 - 2*lam) * K)                                                                              # RHS Matrix for k>=2.
        C2       = Id_inner @ (eye(m) - lam * K)                                                                                        # RHS Matrix for k-2 term.
        
        solve1   = factorized(A1.tocsc())                                                                                               # Pre-factorize LHS for k=1.
        solve2   = factorized(A2)                                                                                                       # Pre-factorize LHS for k>=2.
            
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step logic (k=1).
                g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                      # Evaluate velocity for first step.
                RHS             = B1.dot(u_ap[:, k - 1]) + dt*g_val                                                                     # Right-hand side from previous step + du/dt.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary condition for time k>=2.
                un              = solve1(RHS)                                                                                           # Fast pre-factorized solve for k=1.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
            else:
                RHS             = B2.dot(u_ap[:, k - 1]) - C2.dot(u_ap[:, k - 2])                                                       # Right-hand side from previous steps.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary condition for time k>=2.
                un              = solve2(RHS)                                                                                           # Fast pre-factorized solve for k>=2.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
                    
    if verbose: logger.info("\tCPU Solver finished successfully.")
    
    return u_ap, vec, True
