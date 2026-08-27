import logging
import numpy as np
from scipy.sparse import diags, eye
from typing import Callable, Optional, Tuple, List, Union

from mGFD.exceptions import ParameterError, DimensionMismatchError
import mGFD.spatial.gammas as Gammas
from mGFD.solvers._backends.cuda.preconditioners import compute_preconditioner
import mGFD.spatial.neighbors as Neighbors

logger = logging.getLogger(__name__)

def solve_cuda(p: np.ndarray, 
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
    """CUDA backend for TimeDerivative1."""

    try:
        import cupy as cp
        from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix, csc_matrix as cp_csc_matrix
        from cupyx.scipy.sparse.linalg import factorized as cp_factorized, bicgstab as cp_bicgstab, gmres as cp_gmres, spsolve as cp_spsolve
    except ImportError:
        raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]'.")

    if matrix_free:
        raise ParameterError("matrix_free=True is incompatible with device='cuda'.")

    m      = p.shape[0]
    if verbose: logger.info(f"Solving Transient problem ({t} steps) for {m} nodes on CUDA...")
        
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
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)
    K         = dt * K_spatial
    
    converged = True
    
    if not implicit:
        K2 = eye(m) + K
        K2 = cp_csr_matrix(K2)
        u_ap = cp.asarray(u_ap)
        inne_n = cp.asarray(inne_n)
        for k in range(1, t):
            un              = K2.dot(u_ap[:, k-1])
            u_ap[inne_n, k] = un[inne_n]
    else:
        Id_inner = diags(inne_n.astype(float))
        Id_bound = diags(boun_n.astype(float))
        
        A        = Id_inner @ (eye(m) - lam * K) + Id_bound
        A        = A.tocsc()
        B        = Id_inner @ (eye(m) + (1 - lam) * K)
        
        A = cp_csc_matrix(A)
        B = cp_csr_matrix(B)
        u_ap = cp.asarray(u_ap)
        boun_n = cp.asarray(boun_n)
        inne_n = cp.asarray(inne_n)
        
        # We compute the preconditioner on CPU sparse matrices if needed, since cupy lacks an ILU preconditioner directly accessible sometimes,
        # but the local preconditioner builder handles it returning a cp.LinearOperator if cuda is available.
        # We pass the CuPy matrix A to compute_preconditioner.
        M = compute_preconditioner(A, preconditioner) if preconditioner else None
        
        if linear_solver == "spsolve":
            solve = cp_factorized(A)
        elif linear_solver not in ["bicgstab", "gmres"]:
            raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")
            
        for k in range(1, t):
            RHS             = B.dot(u_ap[:, k-1])
            RHS[boun_n]     = u_ap[boun_n, k]
            
            if linear_solver == "spsolve":
                un       = solve(RHS)
            elif linear_solver == "bicgstab":
                un, info = cp_bicgstab(A, RHS, x0=u_ap[:, k-1], M=M)
                if info != 0:
                    converged = False
                    if verbose: logger.warning(f"CUDA BiCGStab did not converge perfectly (code {info}) at time {k}.")
            elif linear_solver == "gmres":
                un, info = cp_gmres(A, RHS, x0=u_ap[:, k-1], M=M)
                if info != 0:
                    converged = False
                    if verbose: logger.warning(f"CUDA GMRES did not converge perfectly (code {info}) at time {k}.")
                
            u_ap[inne_n, k] = un[inne_n]     
            
    u_ap = u_ap.get()
    
    if verbose: logger.info("\tCUDA Solver finished successfully.")
    
    return u_ap, vec, converged
