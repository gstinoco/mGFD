"""
TimeDerivative1 — CPU Backend for First-order Transient PDEs

Overview:
    CPU implementation for solving parabolic PDEs (like the Heat equation) via implicit time-stepping.

Public API:
    solve_cpu                   Core CPU execution routine for the parabolic solver.

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
from typing import Callable, Optional, Tuple, List, Union                                                                               # Type hinting.
from scipy.sparse.linalg import factorized, bicgstab, gmres, spsolve                                                                    # SciPy sparse linear solvers.

import mGFD.spatial.gammas as Gammas                                                                                                    # Gammas calculation and sparse matrix builder.
import mGFD.spatial.neighbors as Neighbors                                                                                              # Neighbor search routines.

from mGFD.exceptions import ParameterError, DimensionMismatchError                                                                      # Custom exceptions.
from mGFD.solvers._backends.cpu.preconditioners import compute_preconditioner                                                           # Backend imports.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def solve_cpu(p: np.ndarray,                                                                                                            # Function definition.
              f: Union[Callable, np.ndarray, float, int],
              t: int,
              coef: List[float],
              operator: np.ndarray,
              upwind: bool,
              vec: Optional[np.ndarray],
              nvec: int,
              implicit: bool,
              lam: float,
              linear_solver: str,
              preconditioner: Optional[str],
              matrix_free: bool,
              verbose: bool) -> Tuple[np.ndarray, np.ndarray, bool]:
    """CPU backend for TimeDerivative1."""

    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    if verbose: logger.info(f"Solving Transient problem ({t} steps) for {m} nodes on CPU...")
        
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
        if isinstance(f, np.ndarray) and f.shape not in [(m, t), (m,)]:                                                                 # Check if numeric matrix shape matches dimensions.
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")                                       # Raise explicit error.
    elif isinstance(f, (int, float)):                                                                                                   # If data is a constant scalar.
        u_ap[boun_n, :] = f                                                                                                             # Constant boundary condition.
        u_ap[:, 0] = f                                                                                                                  # Constant initial condition.
    
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                              # Upwind-biased neighbor selection.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Standard distance-based neighbors.

    L = operator[:-1]
    if matrix_free:
        from scipy.sparse.linalg import LinearOperator                                                                                  # SciPy sparse linear solvers.
        K_matvec = Gammas.compute_K_matvec(p, vec, L)                                                                                   # Generate the on-the-fly matrix-vector closure.
    else:                                                                                                                               # CPU Iterative Solver.
        K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                             # Build sparse spatial differentiation matrix.
        K         = dt * K_spatial                                                                                                      # Scale by time step.
    
    converged = True                                                                                                                    # Assume convergence by default.
    if matrix_free and preconditioner is not None:
        raise ParameterError("Preconditioners are not currently supported in matrix_free=True mode.")                                   # Preconditioners need explicit matrix elements.
    
    if not implicit:
        if matrix_free:
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                un              = u_ap[:, k-1] + dt * K_matvec(u_ap[:, k-1])                                                            # Explicit matrix-free multiplication.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.
        else:                                                                                                                           # CPU Iterative Solver.
            K2 = eye(m) + K                                                                                                             # LHS Explicit Matrix.
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                un              = K2.dot(u_ap[:, k-1])                                                                                  # Explicit matrix-vector multiplication.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.
    else:                                                                                                                               # CPU Iterative Solver.
        if matrix_free:
            if linear_solver == "spsolve":                                                                                              # Direct pre-factorized solver.
                raise ParameterError("Direct solver 'spsolve' is incompatible with matrix_free=True.")                                  # Direct solvers need explicit matrix.
                
            def A_matvec(x):                                                                                                            # Function definition.
                """
                A_matvec
                
                Closure for the matrix-free application of (I - lam * dt * K) * x.
                
                Input:
                    x           m               ndarray         Input vector to be multiplied.
                    
                Output:
                    y           m               ndarray         Resulting vector.
                """
                y = x - lam * dt * K_matvec(x)                                                                                          # Theta parameter applied to inner nodes.
                y[boun_n] = x[boun_n]                                                                                                   # Identity applied to boundary nodes.
                return y                                                                                                                # Return the result.
            
            A = LinearOperator(shape=(m, m), matvec=A_matvec, dtype=np.float64)                                                         # type: ignore
            M = None                                                                                                                    # No preconditioner in matrix-free.
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                u_prev = u_ap[:, k-1]                                                                                                   # Previous step solution.
                RHS = u_prev + (1 - lam) * dt * K_matvec(u_prev)                                                                        # Matrix-free RHS calculation.
                RHS[boun_n] = u_ap[boun_n, k]                                                                                           # Inject exact boundary conditions.
                if linear_solver == "bicgstab":                                                                                         # Iterative solver (BiCGStab).
                    un, info = bicgstab(A, RHS, x0=u_prev, M=M)                                                                         # Solve with previous step as initial guess.
                elif linear_solver == "gmres":                                                                                          # Iterative solver (GMRES).
                    un, info = gmres(A, RHS, x0=u_prev, M=M)                                                                            # Solve with previous step as initial guess.
                if info != 0:                                                                                                           # If not strictly converged.
                    converged = False                                                                                                   # Mark convergence failure.
                    if verbose: logger.warning(f"{linear_solver.upper()} did not converge perfectly (code {info}) at time {k}.")        # Warn on convergence issues.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.
        else:                                                                                                                           # CPU Iterative Solver.
            Id_inner = diags(inne_n.astype(float))                                                                                      # Diagonal mask for inner nodes.
            Id_bound = diags(boun_n.astype(float))                                                                                      # Diagonal mask for boundary nodes.
            
            A        = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                         # LHS Matrix: Theta parameter applied to inner, Identity to boundary.
            A        = A.tocsc()                                                                                                        # Convert to CSC format for efficient SuperLU factorization.
            B        = Id_inner @ (eye(m) + (1 - lam) * K)                                                                              # RHS Matrix: Zeros for boundaries, explicit part for inner.
            
            M        = compute_preconditioner(A, preconditioner) if preconditioner else None                                            # type: ignore
            
            if linear_solver == "spsolve":                                                                                              # Direct pre-factorized solver.
                solve = factorized(A)                                                                                                   # Pre-factorize LHS for fast repeated solves.
            elif linear_solver not in ["bicgstab", "gmres"]:                                                                            # Invalid iterative solver choice.
                raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")       # Raise explicit error.
                
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                RHS             = B.dot(u_ap[:, k-1])                                                                                   # Right-hand side from previous step.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary conditions.
                
                if linear_solver == "spsolve":                                                                                          # Direct pre-factorized solver.
                    un       = solve(RHS)                                                                                               # Solve global system for time level k.
                elif linear_solver == "bicgstab":                                                                                       # Iterative solver (BiCGStab).
                    un, info = bicgstab(A, RHS, x0=u_ap[:, k-1], M=M)                                                                   # Solve with previous step as initial guess.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")                   # Warn on convergence issues.
                elif linear_solver == "gmres":                                                                                          # Iterative solver (GMRES).
                    un, info = gmres(A, RHS, x0=u_ap[:, k-1], M=M)                                                                      # Solve with previous step as initial guess.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")                      # Warn on convergence issues.
                    
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.
                
    if verbose: logger.info("\tCPU Solver finished successfully.")
    
    return u_ap, vec, converged
