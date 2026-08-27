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
            for k in range(t): u_ap[boun_n, k] = f[boun_n]
            u_ap[:, 0] = f
        if isinstance(f, np.ndarray) and f.shape not in [(m, t), (m,)]:
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")
    elif isinstance(f, (int, float)):
        u_ap[boun_n, :] = f
        u_ap[:, 0] = f
    
    if vec is None:
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)
        else: vec = Neighbors.compute_neighbors(p, nvec)

    L = operator[:-1]
    if matrix_free:
        from scipy.sparse.linalg import LinearOperator
        K_matvec = Gammas.compute_K_matvec(p, vec, L)
    else:
        K_spatial = Gammas.compute_sparse_matrix(p, vec, L)
        K         = dt * K_spatial
    
    converged = True
    if matrix_free and preconditioner is not None:
        raise ParameterError("Preconditioners are not currently supported in matrix_free=True mode.")
    
    if not implicit:
        if matrix_free:
            for k in range(1, t):
                un              = u_ap[:, k-1] + dt * K_matvec(u_ap[:, k-1])
                u_ap[inne_n, k] = un[inne_n]
        else:
            K2 = eye(m) + K
            for k in range(1, t):
                un              = K2.dot(u_ap[:, k-1])
                u_ap[inne_n, k] = un[inne_n]
    else:
        if matrix_free:
            if linear_solver == "spsolve":
                raise ParameterError("Direct solver 'spsolve' is incompatible with matrix_free=True.")
                
            def A_matvec(x):
                """
                A_matvec
                
                Closure for the matrix-free application of (I - lam * dt * K) * x.
                
                Input:
                    x           m               ndarray         Input vector to be multiplied.
                    
                Output:
                    y           m               ndarray         Resulting vector.
                """
                y = x - lam * dt * K_matvec(x)
                y[boun_n] = x[boun_n]
                return y
            
            A = LinearOperator(shape=(m, m), matvec=A_matvec, dtype=np.float64)                                                         # type: ignore
            M = None
            for k in range(1, t):
                u_prev = u_ap[:, k-1]
                RHS = u_prev + (1 - lam) * dt * K_matvec(u_prev)
                RHS[boun_n] = u_ap[boun_n, k]
                if linear_solver == "bicgstab":
                    un, info = bicgstab(A, RHS, x0=u_prev, M=M)
                elif linear_solver == "gmres":
                    un, info = gmres(A, RHS, x0=u_prev, M=M)
                if info != 0:
                    converged = False
                    if verbose: logger.warning(f"{linear_solver.upper()} did not converge perfectly (code {info}) at time {k}.")
                u_ap[inne_n, k] = un[inne_n]
        else:
            Id_inner = diags(inne_n.astype(float))
            Id_bound = diags(boun_n.astype(float))
            
            A        = Id_inner @ (eye(m) - lam * K) + Id_bound
            A        = A.tocsc()
            B        = Id_inner @ (eye(m) + (1 - lam) * K)
            
            M        = compute_preconditioner(A, preconditioner) if preconditioner else None                                            # type: ignore
            
            if linear_solver == "spsolve":
                solve = factorized(A)
            elif linear_solver not in ["bicgstab", "gmres"]:
                raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")
                
            for k in range(1, t):
                RHS             = B.dot(u_ap[:, k-1])
                RHS[boun_n]     = u_ap[boun_n, k]
                
                if linear_solver == "spsolve":
                    un       = solve(RHS)
                elif linear_solver == "bicgstab":
                    un, info = bicgstab(A, RHS, x0=u_ap[:, k-1], M=M)
                    if info != 0:
                        converged = False
                        if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")
                elif linear_solver == "gmres":
                    un, info = gmres(A, RHS, x0=u_ap[:, k-1], M=M)
                    if info != 0:
                        converged = False
                        if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")
                    
                u_ap[inne_n, k] = un[inne_n]     
                
    if verbose: logger.info("\tCPU Solver finished successfully.")
    
    return u_ap, vec, converged
