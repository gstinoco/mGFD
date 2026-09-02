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
              f: Optional[Union[Callable, np.ndarray, float, int]],                                                                     # Execute statement.
              g: Union[Callable, np.ndarray, float, int],                                                                               # Execute statement.
              t: int,                                                                                                                   # Execute statement.
              coef: List[float],                                                                                                        # Execute statement.
              operator: np.ndarray,                                                                                                     # Execute statement.
              upwind: bool,                                                                                                             # Execute statement.
              vec: Optional[np.ndarray],                                                                                                # Execute statement.
              nvec: int,                                                                                                                # Execute statement.
              implicit: bool,                                                                                                           # Execute statement.
              lam: float,                                                                                                               # Execute statement.
              verbose: bool = True,                                                                                                     # Assign verbose: bool.
              ic: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                             # Assign ic: Optional[Union[Callable, np.ndarray, float, int]].
              bc: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                             # Assign bc: Optional[Union[Callable, np.ndarray, float, int]].
              source: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                         # Assign source: Optional[Union[Callable, np.ndarray, float, int]].
              t_span: Tuple[float, float] = (0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray, bool]:                                         # Assign t_span: Tuple[float, float].
    """
    solve_cpu
    CPU backend execution routine for hyperbolic/second-order transient PDEs using pre-factorized sparse LU solver.

    Input:
        p           m x 3           ndarray         Point cloud array [x, y, flag].
        f                           Optional        Legacy/unified initial condition function or array.
        g                           Union           Initial velocity function, array, or scalar.
        t                           int             Total number of time steps.
        coef                        List            Physical coefficients for problem formulation.
        operator    6 x 1           ndarray         Differential operator coefficients [D, E, A, B, C, F_react].
        upwind                      bool            If True, uses upwind neighbor stencil.
        vec         m x nvec        ndarray         Cached neighbor matrix (optional).
        nvec                        int             Maximum number of neighbors per node.
        implicit                    bool            If True, uses implicit time-stepping scheme.
        lam                         float           Theta parameter for implicit scheme.
        verbose                     bool            If True, prints execution logs.
        ic                          Optional        Initial condition function, array, or scalar.
        bc                          Optional        Boundary condition function, array, or scalar.
        source                      Optional        Independent forcing source term F(x, y, t).
        t_span      Tuple[float]                    Physical time domain tuple (t_start, t_end).

    Output:
        u_ap        m x t           ndarray         Approximate solution matrix over time.
        vec         m x nvec        ndarray         Computed/cached neighbor matrix.
        converged                   bool            Solver convergence status flag.
    """

    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    if verbose: logger.info(f"Solving Hyperbolic problem ({t} steps) for {m} nodes on CPU...")                                          # Evaluate condition.
        
    T      = np.linspace(t_span[0], t_span[1], t)                                                                                       # Physical time discretization array across t_span.
    dt     = T[1] - T[0]                                                                                                                # Time step size.
    u_ap   = np.zeros([m, t])                                                                                                           # Numerical approximation matrix.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Boolean mask for boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Boolean mask for interior nodes.

    inne_idx = np.where(inne_n)[0]                                                                                                      # Integer index array for interior nodes.
    boun_idx = np.where(boun_n)[0]                                                                                                      # Integer index array for boundary nodes.

    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]                                                                     # X-velocity (D coefficient).
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]                                                                     # Y-velocity (E coefficient).
    
    ic_use = f if ic is None else ic                                                                                                    # Initial condition fallback to f.
    bc_use = f if bc is None else bc                                                                                                    # Boundary condition fallback to f.

    # 1. Evaluate Initial Condition across all nodes
    if callable(ic_use):                                                                                                                # Evaluate condition.
        try:                                                                                                                            # Execute statement.
            u_ap[:, 0] = np.asarray(ic_use(p[:, 0], p[:, 1], T[0], coef))                                                               # Assign u_ap[:, 0].
        except TypeError:                                                                                                               # Execute statement.
            u_ap[:, 0] = np.asarray(ic_use(p[:, 0], p[:, 1]))                                                                           # Fallback to 2-arg lambda.
    elif isinstance(ic_use, np.ndarray):                                                                                                # Evaluate condition.
        if ic_use.ndim == 2 and ic_use.shape[0] == m: u_ap[:, 0] = ic_use[:, 0]                                                         # Evaluate condition.
        elif ic_use.ndim == 1 and ic_use.shape[0] == m: u_ap[:, 0] = ic_use                                                             # Evaluate condition.
    elif isinstance(ic_use, (int, float)):                                                                                              # Evaluate condition.
        u_ap[:, 0] = ic_use                                                                                                             # Assign u_ap[:, 0].

    # 2. Evaluate Boundary Conditions across boundary nodes over time (t > 0)
    if callable(bc_use):                                                                                                                # Evaluate condition.
        try:                                                                                                                            # Execute statement.
            u_ap[boun_n, 1:] = np.asarray(bc_use(p[boun_n, 0, None], p[boun_n, 1, None], T[None, 1:], coef))                            # Vectorized 2D space-time evaluation for t > 0.
        except Exception:                                                                                                               # Execute statement.
            try:                                                                                                                        # Execute statement.
                bc_nodes = np.asarray(bc_use(p[boun_n, 0], p[boun_n, 1], T[0], coef))                                                   # 1D node vectorization attempt.
                if bc_nodes.ndim == 1 and bc_nodes.shape[0] == len(boun_idx):                                                           # Evaluate condition.
                    u_ap[boun_n, 1:] = bc_nodes[:, None]                                                                                # Broadcast static spatial boundary across time.
                else:                                                                                                                   # Execute fallback branch.
                    for k in range(1, t): u_ap[boun_n, k] = np.asarray(bc_use(p[boun_n, 0], p[boun_n, 1], T[k], coef))                  # Fallback loop.
            except TypeError:                                                                                                           # Execute statement.
                bc_nodes = np.asarray(bc_use(p[boun_n, 0], p[boun_n, 1]))                                                               # 2-arg lambda attempt.
                u_ap[boun_n, 1:] = bc_nodes[:, None]                                                                                    # Broadcast across time.
            except Exception:                                                                                                           # Execute statement.
                for k in range(1, t): u_ap[boun_n, k] = np.asarray(bc_use(p[boun_n, 0], p[boun_n, 1], T[k], coef))                      # Fallback loop.
    elif isinstance(bc_use, np.ndarray):                                                                                                # Evaluate condition.
        if bc_use.ndim == 2 and bc_use.shape == (m, t): u_ap[boun_n, 1:] = bc_use[boun_n, 1:]                                           # Evaluate condition.
        elif bc_use.ndim == 1 and bc_use.shape[0] == m:                                                                                 # Evaluate condition.
            for k in range(1, t): u_ap[boun_n, k] = bc_use[boun_n]                                                                      # Iterate over collection.
        if isinstance(bc_use, np.ndarray) and bc_use.shape not in [(m, t), (m,)]:                                                       # Evaluate condition.
            raise DimensionMismatchError(f"Data array 'bc' must have shape ({m}, {t}) or ({m},).")                                      # Raise exception.
    elif isinstance(bc_use, (int, float)):                                                                                              # Evaluate condition.
        u_ap[boun_n, 1:] = bc_use                                                                                                       # Assign u_ap[boun_n, 1:].

    # 3. Evaluate Initial Velocity Condition g(x, y)
    if callable(g):                                                                                                                     # Evaluate condition.
        try:                                                                                                                            # Execute statement.
            g_eval = np.asarray(g(p[:, 0], p[:, 1], T[0], coef))                                                                        # Assign g_eval.
        except TypeError:                                                                                                               # Execute statement.
            g_eval = np.asarray(g(p[:, 0], p[:, 1]))                                                                                    # Fallback to 2-arg lambda.
    elif isinstance(g, np.ndarray):                                                                                                     # Evaluate condition.
        if g.ndim == 1 and g.shape[0] == m: g_eval = g                                                                                  # Evaluate condition.
        else: raise DimensionMismatchError(f"Data array 'g' must have shape ({m},).")                                                   # Raise exception.
    else:                                                                                                                               # Execute fallback branch.
        g_eval = np.full(m, float(g))                                                                                                   # Assign g_eval.

    # 4. Evaluate Independent Source/Forcing Term F(x, y, t)
    F_mat = None                                                                                                                        # Assign F_mat.
    if source is not None:                                                                                                              # Evaluate condition.
        F_mat = np.zeros((m, t), dtype=float)                                                                                           # Assign F_mat.
        if callable(source):                                                                                                            # Evaluate condition.
            try:                                                                                                                        # Execute statement.
                F_mat = np.asarray(source(p[:, 0, None], p[:, 1, None], T[None, :], coef))                                              # Assign F_mat.
            except Exception:                                                                                                           # Execute statement.
                try:                                                                                                                    # Execute statement.
                    src_nodes = np.asarray(source(p[:, 0], p[:, 1], T[0], coef))                                                        # 1D node vectorization attempt.
                    if src_nodes.ndim == 1 and src_nodes.shape[0] == m:                                                                 # Evaluate condition.
                        F_mat[:, :] = src_nodes[:, None]                                                                                # Broadcast static spatial source across time.
                    else:                                                                                                               # Execute fallback branch.
                        for k in range(t): F_mat[:, k] = np.asarray(source(p[:, 0], p[:, 1], T[k], coef))                               # Fallback loop.
                except TypeError:                                                                                                       # Execute statement.
                    src_nodes = np.asarray(source(p[:, 0], p[:, 1]))                                                                    # 2-arg lambda attempt.
                    F_mat[:, :] = src_nodes[:, None]                                                                                    # Broadcast across time.
                except Exception:                                                                                                       # Execute statement.
                    for k in range(t): F_mat[:, k] = np.asarray(source(p[:, 0], p[:, 1], T[k], coef))                                   # Fallback loop.
        elif isinstance(source, np.ndarray):                                                                                            # Evaluate condition.
            if source.ndim == 2 and source.shape == (m, t): F_mat = source                                                              # Evaluate condition.
            elif source.ndim == 1 and source.shape[0] == m:                                                                             # Evaluate condition.
                for k in range(t): F_mat[:, k] = source                                                                                 # Iterate over collection.
        elif isinstance(source, (int, float)):                                                                                          # Evaluate condition.
            F_mat[:, :] = float(source)                                                                                                 # Assign F_mat[:, :].
    
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                              # Upwind-biased neighbor selection.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Standard distance-based neighbors.

    L         = operator.flatten()                                                                                                      # Flatten operator weights (5 or 6 elements).
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                                 # Build sparse spatial differentiation matrix (includes F_react).
    K         = (dt**2) * K_spatial                                                                                                     # Scale by time step squared.
    
    if not implicit:                                                                                                                    # Explicit scheme.
        K2 = eye(m) + (1/2)*K                                                                                                           # Explicit matrix for k = 1.
        K4 = 2*eye(m) + K                                                                                                               # Explicit matrix for k = 2, ..., t.
        
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step (k=1).
                g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                      # Evaluate velocity for first step.
                un              = K2.dot(u_ap[:, k - 1]) + dt*g_val                                                                     # New time level (k=1).
                if F_mat is not None:                                                                                                   # Evaluate condition.
                    un[inne_idx] += (dt**2 / 2) * F_mat[inne_idx, k-1]                                                                  # Source term for k=1.
                u_ap[inne_idx, k] = un[inne_idx]                                                                                        # Update interior nodes.
            else:                                                                                                                       # Execute fallback branch.
                un              = K4.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]                                                               # New time level (k>=2).
                if F_mat is not None:                                                                                                   # Evaluate condition.
                    un[inne_idx] += (dt**2) * F_mat[inne_idx, k-1]                                                                      # Source term for k>=2.
                u_ap[inne_idx, k] = un[inne_idx]                                                                                        # Update interior nodes.
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
        RHS      = np.empty(m, dtype=float)                                                                                             # Pre-allocated RHS buffer.
            
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step logic (k=1).
                g_val  = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                     # Evaluate velocity for first step.
                RHS[:] = B1.dot(u_ap[:, k - 1])                                                                                         # Right-hand side from previous step.
                RHS   += dt*g_val                                                                                                       # Add initial velocity term.
                if F_mat is not None:                                                                                                   # Evaluate condition.
                    RHS[inne_idx] += (dt**2) * (lam * F_mat[inne_idx, k] + (0.5 - lam) * F_mat[inne_idx, k-1])                          # Source term for k=1.
                RHS[boun_idx] = u_ap[boun_idx, k]                                                                                       # Inject exact boundary condition.
                u_ap[:, k]    = solve1(RHS)                                                                                             # Direct update across all nodes.
            else:                                                                                                                       # Execute fallback branch.
                RHS[:]  = B2.dot(u_ap[:, k - 1])                                                                                        # Multiply B2 by u^(k-1).
                RHS    -= C2.dot(u_ap[:, k - 2])                                                                                        # Subtract C2 by u^(k-2).
                if F_mat is not None:                                                                                                   # Evaluate condition.
                    F_k = lam * F_mat[inne_idx, k] + (1 - 2*lam) * F_mat[inne_idx, k-1] + lam * F_mat[inne_idx, k-2]                    # Spatiotemporal source term.
                    RHS[inne_idx] += (dt**2) * F_k                                                                                      # Source term for k>=2.
                RHS[boun_idx] = u_ap[boun_idx, k]                                                                                       # Inject exact boundary condition.
                u_ap[:, k]    = solve2(RHS)                                                                                             # Direct update across all nodes.
                    
    if verbose: logger.info("\tCPU Solver finished successfully.")                                                                      # Evaluate condition.
    
    return u_ap, vec, True                                                                                                              # Return output values.
