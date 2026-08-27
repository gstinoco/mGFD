import logging
import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import factorized, bicgstab, gmres, spsolve
from typing import Callable, Optional, Tuple, List, Union

from mGFD.exceptions import ParameterError, DimensionMismatchError
import mGFD.spatial.gammas as Gammas
from mGFD.solvers._backends.cpu.preconditioners import compute_preconditioner
import mGFD.spatial.neighbors as Neighbors

logger = logging.getLogger(__name__)

def solve_cpu(p: np.ndarray, 
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
              linear_solver: str,
              preconditioner: Optional[str],
              matrix_free: bool,
              verbose: bool) -> Tuple[np.ndarray, np.ndarray, bool]:
    """CPU backend for TimeDerivative2."""

    m      = p.shape[0]
    if verbose: logger.info(f"Solving Transient problem ({t} steps) for {m} nodes on CPU...")
        
    T      = np.linspace(0, 1, t)
    dt     = T[1] - T[0]
    u_ap   = np.zeros([m, t])
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)
    inne_n = p[:, 2] == 0

    if upwind:
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]

    if callable(f):
        for k in range(t):
            u_ap[boun_n, k] = np.asarray(f(p[boun_n, 0], p[boun_n, 1], T[k], coef))
        u_ap[:, 0] = np.asarray(f(p[:, 0], p[:, 1], T[0], coef))
    elif isinstance(f, np.ndarray):
        if f.ndim == 2 and f.shape == (m, t):
            u_ap[boun_n, :] = f[boun_n, :]
            u_ap[:, 0] = f[:, 0]
        elif f.ndim == 1 and f.shape[0] == m:
            for k in range(t):
                u_ap[boun_n, k] = f[boun_n]
            u_ap[:, 0] = f
        else:
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")
    elif isinstance(f, (int, float)):
        u_ap[boun_n, :] = f
        u_ap[:, 0]      = f

    if callable(g):
        g_eval = np.asarray(g(p[:, 0], p[:, 1], T[0], coef))
    elif isinstance(g, np.ndarray):
        if g.ndim == 1 and g.shape[0] == m:
            g_eval = g
        else:
            raise DimensionMismatchError(f"Data array 'g' must have shape ({m},).")
    else:
        g_eval = np.full(m, float(g))
    
    if vec is None:
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)
        else: vec = Neighbors.compute_neighbors(p, nvec)

    L         = operator[:-1]
    if matrix_free:
        from scipy.sparse.linalg import LinearOperator
        K_matvec = Gammas.compute_K_matvec(p, vec, L)
    else:
        K_spatial = Gammas.compute_sparse_matrix(p, vec, L)
        K         = (dt**2) * K_spatial
    
    converged = True
    if matrix_free and preconditioner is not None:
        raise ParameterError("Preconditioners are not currently supported in matrix_free=True mode.")
    
    if not implicit:
        if matrix_free:
            for k in range(1, t):
                if k == 1:
                    g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval
                    un              = u_ap[:, k - 1] + 0.5 * (dt**2) * K_matvec(u_ap[:, k - 1]) + dt*g_val
                    u_ap[inne_n, k] = un[inne_n]
                else:
                    un              = 2 * u_ap[:, k - 1] + (dt**2) * K_matvec(u_ap[:, k - 1]) - u_ap[:, k - 2]
                    u_ap[inne_n, k] = un[inne_n]
        else:
            K2 = eye(m) + (1/2)*K
            K4 = 2*eye(m) + K
            
            for k in range(1, t):
                if k == 1:
                    g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval
                    un              = K2.dot(u_ap[:, k - 1]) + dt*g_val
                    u_ap[inne_n, k] = un[inne_n]
                else:
                    un              = K4.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]
                    u_ap[inne_n, k] = un[inne_n]
    else:
        if matrix_free:
            if linear_solver == "spsolve":
                raise ParameterError("Direct solver 'spsolve' is incompatible with matrix_free=True.")
            
            def A1_matvec(x):
                """
                A1_matvec
                
                Closure for the matrix-free application of (I - lam * 0.5 * dt**2 * K) * x.
                
                Input:
                    x           m               ndarray         Input vector to be multiplied.
                    
                Output:
                    y           m               ndarray         Resulting vector.
                """
                y = x - lam * (1/2) * (dt**2) * K_matvec(x)
                y[boun_n] = x[boun_n]
                return y
                
            def A2_matvec(x):
                """
                A2_matvec
                
                Closure for the matrix-free application of (I - lam * dt**2 * K) * x.
                
                Input:
                    x           m               ndarray         Input vector to be multiplied.
                    
                Output:
                    y           m               ndarray         Resulting vector.
                """
                y = x - lam * (dt**2) * K_matvec(x)
                y[boun_n] = x[boun_n]
                return y
                
            A1 = LinearOperator(shape=(m, m), matvec=A1_matvec, dtype=np.float64)                                                       # type: ignore
            A2 = LinearOperator(shape=(m, m), matvec=A2_matvec, dtype=np.float64)                                                       # type: ignore
            M1 = None
            M2 = None
            
            for k in range(1, t):
                if k == 1:
                    g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval
                    u_prev = u_ap[:, k - 1]
                    RHS = u_prev + (1 - lam) * (1/2) * (dt**2) * K_matvec(u_prev) + dt*g_val
                    RHS[boun_n] = u_ap[boun_n, k]
                    if linear_solver == "bicgstab":
                        un, info = bicgstab(A1, RHS, M=M1)
                    elif linear_solver == "gmres":
                        un, info = gmres(A1, RHS, M=M1)
                    if info != 0:
                        converged = False
                        if verbose: logger.warning(f"{linear_solver.upper()} did not converge perfectly (code {info}) at time {k}.")
                    u_ap[inne_n, k] = un[inne_n]
                else:
                    u_prev1 = u_ap[:, k - 1]
                    u_prev2 = u_ap[:, k - 2]
                    RHS = 2 * u_prev1 + (1 - lam) * (dt**2) * K_matvec(u_prev1) - u_prev2
                    RHS[boun_n] = u_ap[boun_n, k]
                    if linear_solver == "bicgstab":
                        un, info = bicgstab(A2, RHS, x0=u_prev1, M=M2)
                    elif linear_solver == "gmres":
                        un, info = gmres(A2, RHS, x0=u_prev1, M=M2)
                    if info != 0:
                        converged = False
                        if verbose: logger.warning(f"{linear_solver.upper()} did not converge perfectly (code {info}) at time {k}.")
                    u_ap[inne_n, k] = un[inne_n]
        else:
            Id_inner = diags(inne_n.astype(float))
            Id_bound = diags(boun_n.astype(float))
            
            A1       = Id_inner @ (eye(m) - lam * (1/2) * K) + Id_bound
            B1       = Id_inner @ (eye(m) + (1 - lam) * (1/2) * K)
            
            A2       = Id_inner @ (eye(m) - lam * K) + Id_bound
            A2       = A2.tocsc()
            B2       = Id_inner @ (2*eye(m) + (1 - lam) * K)
            
            M1       = compute_preconditioner(A1.tocsc(), preconditioner) if preconditioner else None                                   # type: ignore
            M2       = compute_preconditioner(A2, preconditioner) if preconditioner else None                                           # type: ignore
            
            if linear_solver == "spsolve":
                solve1   = factorized(A1.tocsc())
                solve2   = factorized(A2)
            elif linear_solver not in ["bicgstab", "gmres"]:
                raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")
                
            for k in range(1, t):
                if k == 1:
                    g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval
                    RHS             = B1.dot(u_ap[:, k - 1]) + dt*g_val
                    RHS[boun_n]     = u_ap[boun_n, k]
                    
                    if linear_solver == "spsolve":
                        un       = solve1(RHS)
                    elif linear_solver == "bicgstab":
                        un, info = bicgstab(A1, RHS, M=M1)
                        if info != 0:
                            converged = False
                            if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")
                    elif linear_solver == "gmres":
                        un, info = gmres(A1, RHS, M=M1)
                        if info != 0:
                            converged = False
                            if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")
                    u_ap[inne_n, k] = un[inne_n]
                else:
                    RHS             = B2.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]
                    RHS[boun_n]     = u_ap[boun_n, k]
                    
                    if linear_solver == "spsolve":
                        un = solve2(RHS)
                    elif linear_solver == "bicgstab":
                        un, info = bicgstab(A2, RHS, x0=u_ap[:, k-1], M=M2)
                        if info != 0:
                            converged = False
                            if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")
                    elif linear_solver == "gmres":
                        un, info = gmres(A2, RHS, x0=u_ap[:, k-1], M=M2)
                        if info != 0:
                            converged = False
                            if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")
                        
                    u_ap[inne_n, k] = un[inne_n]
                    
    if verbose: logger.info("\tCPU Solver finished successfully.")
    
    return u_ap, vec, converged
