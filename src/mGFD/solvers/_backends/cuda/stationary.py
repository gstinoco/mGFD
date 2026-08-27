import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Callable, Optional, Tuple, Union                                                                                     # Type hinting.

from mGFD.exceptions import ParameterError                                                                                              # Custom exceptions.
import mGFD.spatial.gammas as Gammas                                                                                                    # Gammas calculation and sparse matrix builder.
from mGFD.solvers._backends.cuda.preconditioners import compute_preconditioner                                                          # CUDA Preconditioners.
import mGFD.spatial.neighbors as Neighbors                                                                                              # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def solve_cuda(p: np.ndarray, 
               phi: Union[Callable, np.ndarray, float, int],
               f: Union[Callable, np.ndarray, float, int],
               operator: np.ndarray,
               upwind: bool,
               vec: Optional[np.ndarray],
               nvec: int,
               linear_solver: str,
               preconditioner: Optional[str],
               matrix_free: bool,
               verbose: bool) -> Tuple[np.ndarray, np.ndarray, bool]:
    """CUDA backend for Stationary solver."""
    
    try:                                                                                                                                # Import CuPy libraries.
        import cupy as cp                                                                                                               # Core CuPy arrays.
        from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix                                                                      # Sparse CuPy matrix.
        from cupyx.scipy.sparse.linalg import spsolve as cp_spsolve                                                                     # Direct solver.
        from cupyx.scipy.sparse.linalg import bicgstab as cp_bicgstab                                                                   # Iterative BiCGStab.
        from cupyx.scipy.sparse.linalg import gmres as cp_gmres                                                                         # Iterative GMRES.
    except ImportError:                                                                                                                 # Catch missing dependency.
        raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]' to use device='cuda'.")                # Raise clear error.

    m = len(p[:, 0])                                                                                                                    # Total nodes.
    if verbose:                                                                                                                         # Verbosity.
        logger.info(f"Solving Stationary problem for {m} nodes on CUDA...")                                                             # Progress info.
    
    u_ap = np.zeros([m])                                                                                                                # u_ap init (CPU).
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
    
    if matrix_free:                                                                                                                     # Matrix-free check.
        raise ParameterError("matrix_free=True is incompatible with device='cuda'.")                                                    # CuPy requires explicit sparse matrices.
        
    K = Gammas.compute_sparse_matrix(p, vec, L)                                                                                         # K sparse generation on CPU.
    K = cp_csr_matrix(K)                                                                                                                # Transfer matrix to GPU.
    
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                                                                           # RHS on CPU.
    R = cp.asarray(R)                                                                                                                   # Transfer vector to GPU.
    
    converged = True                                                                                                                    # Default to true.
    M = compute_preconditioner(K, preconditioner) if preconditioner else None                                                           # Optional precond.

    if linear_solver == "spsolve":                                                                                                      # Direct solver.
        un = cp_spsolve(K, R)                                                                                                           # CuPy spsolve.
    elif linear_solver == "bicgstab":                                                                                                   # Iterative BiCGStab.
        un, info = cp_bicgstab(K, R, M=M)                                                                                               # CuPy bicgstab.
        if info != 0:                                                                                                                   # If fail.
            converged = False                                                                                                           # Mark fail.
            if verbose: logger.warning(f"CUDA BiCGStab did not converge perfectly (code {info}).")                                      # Log fail.
    elif linear_solver == "gmres":                                                                                                      # Iterative GMRES.
        un, info = cp_gmres(K, R, M=M)                                                                                                  # CuPy gmres.
        if info != 0:                                                                                                                   # If fail.
            converged = False                                                                                                           # Mark fail.
            if verbose: logger.warning(f"CUDA GMRES did not converge perfectly (code {info}).")                                         # Log fail.
    else:                                                                                                                               # Unknown.
        raise ParameterError(f"Unsupported linear_solver '{linear_solver}'.")                                                           # Raise error.

    un = un.get()                                                                                                                       # Pull from GPU.
    u_ap[inne_n] = un[inne_n]                                                                                                           # Unpack to interior.
    if verbose: logger.info("\tCUDA Solver finished successfully.")                                                                     # Success.
    
    return u_ap, vec, converged                                                                                                         # Return core data.
