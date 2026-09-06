"""
TimeDerivative2 — CUDA Backend for Second-order Transient PDEs

Overview:
    CUDA implementation for solving hyperbolic PDEs (like the Wave equation) using pre-factorized direct cuSOLVER sparse solvers.

Public API:
    solve_cuda                  Core CUDA execution routine for the hyperbolic solver.

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
               g: Union[Callable, np.ndarray, float, int],                                                                              # Execute statement.
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
               t_span: Tuple[float, float] = (0.0, 1.0),                                                                                # Assign t_span: Tuple[float, float].
               damping: float = 0.0,                                                                                                    # Assign damping: float (velocity damping eta).
               alpha: float = 0.0) -> Tuple[np.ndarray, np.ndarray, bool]:                                                              # Assign alpha: float (HHT-alpha high-frequency dissipation).
    """CUDA backend for TimeDerivative2 using CuPy pre-factorized direct sparse LU solver."""

    try:                                                                                                                                # Attempt CuPy import.
        import cupy as cp                                                                                                               # CuPy GPU array operations.
        from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix, csc_matrix as cp_csc_matrix                                         # CuPy sparse matrices.
        from cupyx.scipy.sparse.linalg import factorized as cp_factorized                                                               # CuPy direct sparse factorized solver.
    except ImportError:                                                                                                                 # Catch missing CuPy.
        raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]'.")                                     # Friendly error message.

    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    if verbose: logger.info(f"Solving Hyperbolic problem ({t} steps) for {m} nodes on CUDA...")                                         # Evaluate condition.
        
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
            res = np.asarray(bc_use(p[boun_n, 0, None], p[boun_n, 1, None], T[None, 1:], coef))                                         # Vectorized 2D space-time evaluation for t > 0.
            if res.shape == (len(boun_idx), t - 1):                                                                                     # Exact shape match.
                u_ap[boun_n, 1:] = res                                                                                                  # Assign 2D matrix.
            elif res.shape == (len(boun_idx), 1):                                                                                       # Static spatial boundary column.
                u_ap[boun_n, 1:] = np.broadcast_to(res, (len(boun_idx), t - 1))                                                         # Broadcast across time steps.
            else:                                                                                                                       # Fallback on shape mismatch.
                raise ValueError("Shape mismatch")                                                                                      # Force fallback execution.
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
                res = np.asarray(source(p[:, 0, None], p[:, 1, None], T[None, :], coef))                                                # 2D space-time evaluation attempt.
                if res.shape == (m, t):                                                                                                 # Exact 2D matrix shape match.
                    F_mat = res                                                                                                         # Assign matrix.
                elif res.shape == (m, 1):                                                                                               # Single spatial column returned.
                    F_mat = np.broadcast_to(res, (m, t)).copy()                                                                         # Broadcast across time steps.
                elif res.shape == (m,):                                                                                                 # 1D node vector returned.
                    F_mat = np.broadcast_to(res[:, None], (m, t)).copy()                                                                # Broadcast across time steps.
                else:                                                                                                                   # Fallback on shape mismatch.
                    raise ValueError("Shape mismatch")                                                                                  # Force fallback execution.
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
        if upwind and (a != 0.0 or b != 0.0): vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                   # Upwind-biased neighbor selection.
        else: vec = Neighbors.compute_neighbors(p, nvec)                                                                                # Standard distance-based neighbors.

    L         = operator.flatten()                                                                                                      # Flatten operator weights (5 or 6 elements).
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                                 # Build sparse spatial differentiation matrix (includes F_react).
    K         = (dt**2) * K_spatial                                                                                                     # Scale by time step squared.
    
    u_ap_gpu     = cp.asarray(u_ap)                                                                                                     # Transfer solution matrix to VRAM.
    F_mat_gpu    = cp.asarray(F_mat) if F_mat is not None else None                                                                     # Transfer source matrix to VRAM.
    boun_idx_gpu = cp.where(cp.asarray(boun_n))[0]                                                                                      # Integer index array for boundary nodes on GPU.
    inne_idx_gpu = cp.where(cp.asarray(inne_n))[0]                                                                                      # Integer index array for interior nodes on GPU.
    
    # Calculate HHT-alpha and Newmark parameters
    if alpha != 0.0:                                                                                                                    # HHT-alpha specified.
        gamma_n = 0.5 - alpha                                                                                                           # Compute HHT gamma.
        beta_n  = 0.25 * ((1.0 - alpha) ** 2)                                                                                           # Compute HHT beta.
    else:                                                                                                                               # Standard Newmark-beta.
        gamma_n = 0.5                                                                                                                   # Standard gamma.
        beta_n  = lam                                                                                                                   # Standard beta (lam).

    damp_scale = damping * dt                                                                                                           # Dimensionless velocity damping scale.

    if not implicit:                                                                                                                    # Explicit scheme on GPU.
        K2_gpu     = cp_csr_matrix(eye(m) + (1/2)*K)                                                                                    # Transfer explicit matrix k=1.
        K4_gpu     = cp_csr_matrix((2.0 - damp_scale)*eye(m) + K)                                                                       # Transfer explicit matrix k>=2.
        K_prev_gpu = cp_csr_matrix((damp_scale - 1.0)*eye(m))                                                                           # Previous step matrix on GPU.
        
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step (k=1).
                g_val     = g_eval                                                                                                      # Use pre-evaluated initial velocity vector.
                g_val_gpu = cp.asarray(g_val)                                                                                           # Transfer velocity array to VRAM.
                un        = K2_gpu.dot(u_ap_gpu[:, k - 1]) + dt*g_val_gpu                                                               # New time level (k=1).
                if F_mat_gpu is not None:                                                                                               # Evaluate condition.
                    un[inne_idx_gpu]   += (dt**2 / 2) * F_mat_gpu[inne_idx_gpu, k-1]                                                    # GPU source term for k=1.
                u_ap_gpu[inne_idx_gpu, k] = un[inne_idx_gpu]                                                                            # Update interior nodes.
            else:                                                                                                                       # Execute fallback branch.
                un                      = K4_gpu.dot(u_ap_gpu[:, k - 1]) + K_prev_gpu.dot(u_ap_gpu[:, k - 2])                           # New time level (k>=2) with damping.
                if F_mat_gpu is not None:                                                                                               # Evaluate condition.
                    un[inne_idx_gpu]   += (dt**2) * F_mat_gpu[inne_idx_gpu, k-1]                                                        # GPU source term for k>=2.
                u_ap_gpu[inne_idx_gpu, k] = un[inne_idx_gpu]                                                                            # Update interior nodes.
    else:                                                                                                                               # Implicit GPU Iterative Solver with Warm Start.
        from cupyx.scipy.sparse.linalg import bicgstab as cp_bicgstab                                                                   # GPU iterative BiCGSTAB solver.

        Id_inner = diags(inne_n.astype(float))                                                                                          # Diagonal mask for inner nodes.
        Id_bound = diags(boun_n.astype(float))                                                                                          # Diagonal mask for boundary nodes.
        
        coef_lhs = (1.0 + gamma_n * damp_scale)                                                                                         # LHS Identity coefficient.
        coef_b1  = (1.0 + (gamma_n - 1.0) * damp_scale)                                                                                 # B1 Identity coefficient.
        coef_b2  = (2.0 + (2.0 * gamma_n - 1.0) * damp_scale)                                                                           # B2 Identity coefficient.
        coef_c2  = (1.0 + (gamma_n - 1.0) * damp_scale)                                                                                 # C2 Identity coefficient.
        
        beta_alpha_factor = beta_n * (1.0 + alpha)                                                                                      # Beta-alpha factor.

        A1_gpu   = cp_csr_matrix(Id_inner @ (coef_lhs * eye(m) - beta_alpha_factor * K) + Id_bound)                                     # LHS Matrix k=1 CSR.
        B1_gpu   = cp_csr_matrix(Id_inner @ (coef_b1 * eye(m) + (0.5 - beta_alpha_factor) * K))                                         # RHS Matrix k=1 CSR.
        
        A2_gpu   = cp_csr_matrix(Id_inner @ (coef_lhs * eye(m) - beta_alpha_factor * K) + Id_bound)                                     # LHS Matrix k>=2 CSR.
        B2_gpu   = cp_csr_matrix(Id_inner @ (coef_b2 * eye(m) + (1.0 + alpha - 2.0 * beta_alpha_factor) * K))                           # RHS Matrix k>=2 CSR.
        C2_gpu   = cp_csr_matrix(Id_inner @ (coef_c2 * eye(m) - beta_alpha_factor * K))                                                 # RHS Matrix k-2 CSR.
        
        x_iter   = cp.copy(u_ap_gpu[:, 0])                                                                                              # Initialize warm start solution vector.
            
        for k in range(1, t):                                                                                                           # Loop over all time steps in VRAM.
            if k == 1:                                                                                                                  # Initial time step (k=1).
                g_val     = g_eval                                                                                                      # Use pre-evaluated initial velocity vector.
                g_val_gpu = cp.asarray(g_val)                                                                                           # Transfer velocity array to VRAM.
                RHS       = B1_gpu.dot(u_ap_gpu[:, k - 1]) + dt*g_val_gpu                                                               # Right-hand side for k=1.
                if F_mat_gpu is not None:                                                                                               # Evaluate condition.
                    baf       = beta_alpha_factor                                                                                       # Alias beta-alpha factor.
                    term1_gpu = baf * F_mat_gpu[inne_idx_gpu, k] + (0.5 - baf) * F_mat_gpu[inne_idx_gpu, k-1]                           # GPU forcing term k=1.
                    RHS[inne_idx_gpu]   += (dt**2) * term1_gpu                                                                          # GPU source term for k=1.
                
                RHS[boun_idx_gpu]        = u_ap_gpu[boun_idx_gpu, k]                                                                    # Inject exact boundary condition.
                x_iter, _                = cp_bicgstab(A1_gpu, RHS, x0=x_iter, rtol=1e-10, atol=1e-10, maxiter=50)                      # GPU BiCGSTAB iterative solve for k=1.
                u_ap_gpu[:, k]           = x_iter                                                                                       # Direct update across all nodes.
            else:                                                                                                                       # Execute fallback branch.
                RHS                      = B2_gpu.dot(u_ap_gpu[:, k - 1]) - C2_gpu.dot(u_ap_gpu[:, k - 2])                              # Right-hand side for k>=2.
                if F_mat_gpu is not None:                                                                                               # Evaluate condition.
                    F_k_gpu = (beta_alpha_factor * F_mat_gpu[inne_idx_gpu, k] +                                                         # GPU source component at step k.
                               (1.0 + alpha - 2.0 * beta_alpha_factor) * F_mat_gpu[inne_idx_gpu, k-1] +                                 # GPU source component at step k-1.
                               beta_alpha_factor * F_mat_gpu[inne_idx_gpu, k-2])                                                        # GPU source component at step k-2.
                    RHS[inne_idx_gpu]    += (dt**2) * F_k_gpu                                                                           # GPU source term for k>=2.
                
                RHS[boun_idx_gpu] = u_ap_gpu[boun_idx_gpu, k]                                                                           # Inject exact boundary condition.
                x_iter, _         = cp_bicgstab(A2_gpu, RHS, x0=x_iter, rtol=1e-10, atol=1e-10, maxiter=50)                             # GPU BiCGSTAB iterative solve for k>=2.
                u_ap_gpu[:, k]    = x_iter                                                                                              # Direct update across all nodes.
                    
    u_ap = u_ap_gpu.get()                                                                                                               # Pull final solution array to CPU RAM.
    cp.get_default_memory_pool().free_all_blocks()                                                                                      # Free CuPy VRAM blocks to prevent pool fragmentation across dataset sweeps.
    
    if verbose: logger.info("\tCUDA Solver finished successfully.")                                                                     # Evaluate condition.
    
    return u_ap, vec, True                                                                                                              # Return output values.
