"""
TimeDerivative1 — CUDA Backend for First-order Transient PDEs

Overview:
    CUDA implementation for solving parabolic PDEs (like the Heat equation) via pre-factorized direct cuSOLVER sparse solvers.

Public API:
    solve_cuda                  Core CUDA execution routine for the parabolic solver.

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

import mGFD.spatial.gammas as Gammas                                                                                                    # Gammas calculation and sparse matrix builder.
import mGFD.spatial.neighbors as Neighbors                                                                                              # Neighbor search routines.

from mGFD.exceptions import DimensionMismatchError                                                                                      # Custom exceptions.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def solve_cuda(p: np.ndarray,                                                                                                           # Function definition.
               f: Optional[Union[Callable, np.ndarray, float, int]],                                                                    # Execute statement.
               t: int,                                                                                                                  # Execute statement.
               coef: List[float],                                                                                                       # Execute statement.
               operator: np.ndarray,                                                                                                    # Execute statement.
               upwind: bool,                                                                                                            # Execute statement.
               vec: Optional[np.ndarray],                                                                                               # Execute statement.
               nvec: int,                                                                                                               # Execute statement.
               implicit: bool,                                                                                                          # Execute statement.
               lam: float,                                                                                                              # Execute statement.
               verbose: bool = True,                                                                                                    # Assign verbose: bool.
               ic: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                            # Assign ic: Optional[Union[Callable, np.ndarray, float, int]].
               bc: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                            # Assign bc: Optional[Union[Callable, np.ndarray, float, int]].
               source: Optional[Union[Callable, np.ndarray, float, int]] = None,                                                        # Assign source: Optional[Union[Callable, np.ndarray, float, int]].
               t_span: Tuple[float, float] = (0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray, bool]:                                        # Assign t_span: Tuple[float, float].
    """CUDA backend for TimeDerivative1 using CuPy pre-factorized direct sparse LU solver."""

    try:                                                                                                                                # Attempt CuPy import.
        import cupy as cp                                                                                                               # CuPy GPU array operations.
        from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix, csc_matrix as cp_csc_matrix                                         # CuPy sparse matrices.
        from cupyx.scipy.sparse.linalg import factorized as cp_factorized                                                               # CuPy direct sparse factorized solver.
    except ImportError:                                                                                                                 # Catch missing CuPy.
        raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]'.")                                     # Friendly error message.

    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    if verbose: logger.info(f"Solving Transient problem ({t} steps) for {m} nodes on CUDA...")                                          # Evaluate condition.
        
    T      = np.linspace(t_span[0], t_span[1], t)                                                                                       # Physical time discretization array across t_span.
    dt     = T[1] - T[0]                                                                                                                # Time step size.
    u_ap   = np.zeros([m, t])                                                                                                           # Numerical approximation matrix.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Boolean mask for boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Boolean mask for interior nodes.

    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]                                                                     # X-velocity (D coefficient).
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]                                                                     # Y-velocity (E coefficient).

    ic_use = f if ic is None else ic                                                                                                    # Initial condition fallback to f.
    bc_use = f if bc is None else bc                                                                                                    # Boundary condition fallback to f.

    # 1. Evaluate Initial Condition across all nodes
    if callable(ic_use):                                                                                                                # Evaluate condition.
        u_ap[:, 0] = np.asarray(ic_use(p[:, 0], p[:, 1], T[0], coef))                                                                   # Assign u_ap[:, 0].
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
            for k in range(1, t): u_ap[boun_n, k] = np.asarray(bc_use(p[boun_n, 0], p[boun_n, 1], T[k], coef))                          # Fallback loop.
    elif isinstance(bc_use, np.ndarray):                                                                                                # Evaluate condition.
        if bc_use.ndim == 2 and bc_use.shape == (m, t): u_ap[boun_n, 1:] = bc_use[boun_n, 1:]                                           # Evaluate condition.
        elif bc_use.ndim == 1 and bc_use.shape[0] == m:                                                                                 # Evaluate condition.
            for k in range(1, t): u_ap[boun_n, k] = bc_use[boun_n]                                                                      # Iterate over collection.
        if isinstance(bc_use, np.ndarray) and bc_use.shape not in [(m, t), (m,)]:                                                       # Evaluate condition.
            raise DimensionMismatchError(f"Data array 'bc' must have shape ({m}, {t}) or ({m},).")                                      # Raise exception.
    elif isinstance(bc_use, (int, float)):                                                                                              # Evaluate condition.
        u_ap[boun_n, 1:] = bc_use                                                                                                       # Assign u_ap[boun_n, 1:].

    # 3. Evaluate Independent Source/Forcing Term F(x, y, t)
    F_mat = None                                                                                                                        # Assign F_mat.
    if source is not None:                                                                                                              # Evaluate condition.
        F_mat = np.zeros((m, t), dtype=float)                                                                                           # Assign F_mat.
        if callable(source):                                                                                                            # Evaluate condition.
            try:                                                                                                                        # Execute statement.
                F_mat = np.asarray(source(p[:, 0, None], p[:, 1, None], T[None, :], coef))                                              # Assign F_mat.
            except Exception:                                                                                                           # Execute statement.
                for k in range(t): F_mat[:, k] = np.asarray(source(p[:, 0], p[:, 1], T[k], coef))                                       # Iterate over collection.
        elif isinstance(source, np.ndarray):                                                                                            # Evaluate condition.
            if source.ndim == 2 and source.shape == (m, t): F_mat = source                                                              # Evaluate condition.
            elif source.ndim == 1 and source.shape[0] == m:                                                                             # Evaluate condition.
                for k in range(t): F_mat[:, k] = source                                                                                 # Iterate over collection.
        elif isinstance(source, (int, float)):                                                                                          # Evaluate condition.
            F_mat[:, :] = float(source)                                                                                                 # Assign F_mat[:, :].
    
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind: vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                              # Upwind-biased neighbor selection.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Standard distance-based neighbors.

    L = operator.flatten()                                                                                                              # Flatten operator coefficients (5 or 6 elements).
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                                 # Build sparse spatial differentiation matrix (includes F_react).
    K         = dt * K_spatial                                                                                                          # Scale by time step.
    
    u_ap_gpu     = cp.asarray(u_ap)                                                                                                     # Transfer solution matrix to VRAM.
    F_mat_gpu    = cp.asarray(F_mat) if F_mat is not None else None                                                                     # Transfer source matrix to VRAM.
    boun_idx_gpu = cp.where(cp.asarray(boun_n))[0]                                                                                      # Integer index array for boundary nodes on GPU.
    inne_idx_gpu = cp.where(cp.asarray(inne_n))[0]                                                                                      # Integer index array for interior nodes on GPU.
    
    if not implicit:                                                                                                                    # Evaluate condition.
        K2_gpu = cp_csr_matrix(eye(m) + K)                                                                                              # Transfer LHS explicit matrix to GPU.
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            un = K2_gpu.dot(u_ap_gpu[:, k-1])                                                                                           # Explicit matrix-vector multiplication in VRAM.
            if F_mat_gpu is not None:                                                                                                   # Evaluate condition.
                un[inne_idx_gpu] += dt * F_mat_gpu[inne_idx_gpu, k-1]                                                                   # Inject GPU source term F^(k-1).
            u_ap_gpu[inne_idx_gpu, k] = un[inne_idx_gpu]                                                                                # Update interior nodes on GPU.
    else:                                                                                                                               # Implicit Direct cuSOLVER Solver.
        Id_inner = diags(inne_n.astype(float))                                                                                          # Diagonal mask for inner nodes.
        Id_bound = diags(boun_n.astype(float))                                                                                          # Diagonal mask for boundary nodes.
        
        A        = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                             # LHS Matrix.
        A_gpu    = cp_csc_matrix(A)                                                                                                     # Transfer LHS matrix to GPU CSC.
        B_gpu    = cp_csr_matrix(Id_inner @ (eye(m) + (1 - lam) * K))                                                                   # Transfer RHS matrix to GPU CSR.
        
        solve_gpu = cp_factorized(A_gpu)                                                                                                # Pre-factorize LHS on GPU ONCE.
            
        for k in range(1, t):                                                                                                           # Loop over all time steps in VRAM.
            RHS                      = B_gpu.dot(u_ap_gpu[:, k-1])                                                                      # RHS product in VRAM.
            if F_mat_gpu is not None:                                                                                                   # Evaluate condition.
                RHS[inne_idx_gpu]   += dt * (lam * F_mat_gpu[inne_idx_gpu, k] + (1 - lam) * F_mat_gpu[inne_idx_gpu, k-1])               # Inject theta-weighted GPU source term.
            RHS[boun_idx_gpu]        = u_ap_gpu[boun_idx_gpu, k]                                                                        # Inject boundary condition.
            u_ap_gpu[:, k]           = solve_gpu(RHS[:, None]).ravel()                                                                  # GPU pre-factorized solve.
            
    u_ap = u_ap_gpu.get()                                                                                                               # Pull final result array to CPU RAM.
    cp.get_default_memory_pool().free_all_blocks()                                                                                      # Free CuPy VRAM blocks to prevent pool fragmentation across dataset sweeps.
    
    if verbose: logger.info("\tCUDA Solver finished successfully.")                                                                     # Evaluate condition.
    
    return u_ap, vec, True                                                                                                              # Return output values.
