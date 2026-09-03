"""
Gammas — Stencil weights and linear-algebra helpers for meshless GFD

Overview:
    Core numerical routines used by the mGFD solvers to:
    - Build generalized finite-difference (GFD) stencil weights on a 2D point cloud
    - Assemble right-hand sides for Dirichlet problems
    - Apply a precomputed stencil to a vector (matrix-free operator application)
    - Solve linear systems with a BiCGStab Krylov method (SciPy-accelerated when available)

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.
    L       (5,) or (5, 1) array-like
            Differential operator coefficients for:
                D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy
            provided as [D, E, A, B, C].

Public API:
    Cloud               Build a full dense stencil matrix K (mainly for small problems / debugging).
    RHS                 Build RHS vector with Dirichlet boundary values and interior forcing.
    CloudStencil        Build stencil in a sparse-like form: diagonal + neighbor weights.
    ApplyCloudStencil   Apply a CloudStencil (diag, w) to a vector u with vec neighbor indexing.
    BiCGStab            Matrix-free BiCGStab solver with optional SciPy acceleration.

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
import numpy as np                                                                                                                      # Core numerical operations.
import numba as nb                                                                                                                      # JIT compiler.

from scipy.sparse import csr_matrix                                                                                                     # Core sparse matrix representations.
from typing import Callable, Optional, Tuple, Union                                                                                     # Type hinting.
from scipy.sparse.linalg import LinearOperator, bicgstab                                                                                # SciPy iterative solver interface.

@nb.njit(cache=True, fastmath=True)                                                                                                     # Assign @nb.njit(cache.
def _nnls_numba(A: np.ndarray, b: np.ndarray, max_iter: int = 100) -> np.ndarray:
    """
    _nnls_numba
    Non-Negative Least Squares solver using Lawson-Hanson algorithm, compiled with Numba JIT.
    
    Input:
        A           m x n           ndarray             Design matrix.
        b           m               ndarray             Right-hand side vector.
        max_iter                    int                 Maximum number of iterations.
    
    Output:
        x           n               ndarray             Solution vector with non-negative constraints.
    """
    # 1. Variable initialization
    m, n       = A.shape                                                                                                                # Matrix dimensions.
    P          = np.zeros(n, dtype=np.bool_)                                                                                            # Active set mask (positive).
    Z          = np.ones(n, dtype=np.bool_)                                                                                             # Inactive set mask (zero).
    x          = np.zeros(n, dtype=np.float64)                                                                                          # Solution vector.
    w          = np.dot(A.T, (b - np.dot(A, x)))                                                                                        # Dual vector (gradient).
    iter_count = 0                                                                                                                      # Iteration counter.
    
    # 2. Main iterative loop
    while np.any(Z) and iter_count < max_iter:                                                                                          # While there are inactive variables.
        w_z = w.copy()                                                                                                                  # Copy gradient.
        
        for i in range(n):                                                                                                              # Iterate over all variables.
            if not Z[i]:                                                                                                                # If variable is active.
                w_z[i] = -1e9                                                                                                           # Penalize active variables to ignore them.
        
        max_w = np.max(w_z)                                                                                                             # Find maximum gradient in inactive set.
        
        if max_w <= 1e-10:                                                                                                              # If optimal condition met.
            break                                                                                                                       # Exit loop.
        
        j    = np.argmax(w_z)                                                                                                           # Index of maximum gradient.
        P[j] = True                                                                                                                     # Move to active set.
        Z[j] = False                                                                                                                    # Remove from inactive set.
        
        while True:                                                                                                                     # Inner loop for active set optimization.
            P_count = 0                                                                                                                 # Count active variables.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i]:                                                                                                                # If active.
                    P_count += 1                                                                                                        # Increment count.
            
            A_P     = np.zeros((m, P_count), dtype=np.float64)                                                                          # Allocate submatrix.
            col_idx = 0                                                                                                                 # Column index.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i]:                                                                                                                # If active.
                    A_P[:, col_idx] = A[:, i]                                                                                           # Copy column.
                    col_idx        += 1                                                                                                 # Increment column index.
            
            if A_P.shape[1] == 0:                                                                                                       # Degenerate case.
                break                                                                                                                   # Exit inner loop.
            
            z_P     = np.dot(np.linalg.pinv(A_P), b)                                                                                    # Solve least squares for active set.
            z       = np.zeros(n, dtype=np.float64)                                                                                     # Allocate full vector.
            col_idx = 0                                                                                                                 # Column index.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i]:                                                                                                                # If active.
                    z[i]     = z_P[col_idx]                                                                                             # Copy solution.
                    col_idx += 1                                                                                                        # Increment column index.
            
            all_pos = True                                                                                                              # Check if all active variables are positive.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i] and z[i] <= 0:                                                                                                  # If active and negative.
                    all_pos = False                                                                                                     # Flag as not all positive.
                    break                                                                                                               # Stop checking.
            
            if all_pos:                                                                                                                 # If feasible.
                x = z                                                                                                                   # Update solution.
                break                                                                                                                   # Exit inner loop.
            
            alpha = 1.0                                                                                                                 # Initialization of step size.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i] and z[i] <= 0:                                                                                                  # If active and negative.
                    alpha_i = x[i] / (x[i] - z[i])                                                                                      # Compute max step size.
                    
                    if alpha_i < alpha:                                                                                                 # If smaller than current step size.
                        alpha = alpha_i                                                                                                 # Update step size.
            
            x = x + alpha * (z - x)                                                                                                     # Apply step.
            
            for i in range(n):                                                                                                          # Iterate over all variables.
                if P[i] and abs(x[i]) < 1e-12:                                                                                          # If active and zero.
                    P[i] = False                                                                                                        # Remove from active set.
                    Z[i] = True                                                                                                         # Move to inactive set.
        
        w = np.dot(A.T, (b - np.dot(A, x)))                                                                                             # Update gradient.
        iter_count += 1                                                                                                                 # Increment iteration counter.
    
    return x                                                                                                                            # Return solution.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Assign @nb.njit(cache.
def _compute_cloud_dense_jit(p: np.ndarray, vec: np.ndarray, L: np.ndarray) -> np.ndarray:
    """
    _compute_cloud_dense_jit
    Assemble a full dense stencil matrix K for a 2D point cloud, compiled with Numba JIT.
    
    Input:
        p           m x 3           ndarray             Point cloud [x, y, flag].
        vec         m x nvec        ndarray[int]        Neighbor indices per node (padded with -1).
        L           5               ndarray             Operator coefficients [D, E, A, B, C].
     
    Output:
        K           m x m           ndarray             Dense matrix with stencil weights.
    """
    # 1. Variable initialization
    m    = p.shape[0]                                                                                                                   # Total number of nodes in the cloud.
    nvec = vec.shape[1]                                                                                                                 # The maximum number of neighbors (vec width).
    K    = np.zeros((m, m), dtype=np.float64)                                                                                           # Dense stencil matrix initialization.
    
    # 2. Dense Gammas computation
    for i in nb.prange(m):                                                                                                              # type: ignore
        if p[i, 2] == 0:                                                                                                                # Interior node: compute a stencil row.
            nvec_i = 0                                                                                                                  # Neighbor count initialization.
            
            for j in range(nvec):                                                                                                       # Iterate over neighbor slots.
                if vec[i, j] != -1:                                                                                                     # If valid neighbor.
                    nvec_i += 1                                                                                                         # Increment neighbor count.
            
            if nvec_i == 0:                                                                                                             # If no neighbors available.
                K[i, i] = 1.0                                                                                                           # Fallback to identity.
            else:                                                                                                                       # Execute fallback branch.
                dx            = np.zeros(nvec_i, dtype=np.float64)                                                                      # Neighbor x-offsets (relative to node i).
                dy            = np.zeros(nvec_i, dtype=np.float64)                                                                      # Neighbor y-offsets (relative to node i).
                neigh_indices = np.zeros(nvec_i, dtype=np.int64)                                                                        # Neighbor indices array.
                
                idx           = 0                                                                                                       # Local index.
                
                for j in range(nvec):                                                                                                   # Iterate over neighbor slots.
                    if vec[i, j] != -1:                                                                                                 # If valid neighbor.
                        vec1               = vec[i, j]                                                                                  # Fetch neighbor global index.
                        dx[idx]            = p[vec1, 0] - p[i, 0]                                                                       # Compute x-offset.
                        dy[idx]            = p[vec1, 1] - p[i, 1]                                                                       # Compute y-offset.
                        neigh_indices[idx] = vec1                                                                                       # Store neighbor index.
                        idx               += 1                                                                                          # Increment local index.
                
                M = np.zeros((5, nvec_i), dtype=np.float64)                                                                             # Reconstruction matrix initialization.
                
                for j in range(nvec_i):                                                                                                 # Populate reconstruction matrix.
                    M[0, j] = dx[j]                                                                                                     # dx term.
                    M[1, j] = dy[j]                                                                                                     # dy term.
                    M[2, j] = dx[j]**2                                                                                                  # dx^2 term.
                    M[3, j] = dx[j] * dy[j]                                                                                             # dx*dy term.
                    M[4, j] = dy[j]**2                                                                                                  # dy^2 term.
                
                L_derivatives = L[:5]                                                                                                   # Spatial derivative coefficients [D, E, A, B, C].
                F_react       = L[5] if L.shape[0] > 5 else 0.0                                                                         # Reaction term coefficient F(p_i) per Eq. 190 of paper.
                M_pinv        = np.linalg.pinv(M)                                                                                       # Compute pseudoinverse of reconstruction matrix.
                YY            = np.dot(M_pinv, L_derivatives)                                                                           # Compute neighbor weights.
                K[i, i]       = -np.sum(YY) + F_react                                                                                   # Central weight includes reaction term per Eq. 190.
                
                for j in range(nvec_i):                                                                                                 # Assign neighbor weights.
                    K[i, neigh_indices[j]] = YY[j]                                                                                      # Set weight for neighbor j.
        
        elif p[i, 2] == 1 or p[i, 2] == 2:                                                                                              # Boundary node: Dirichlet identity row.
            K[i, i] = 1.0                                                                                                               # Enforce u_i = boundary value.
            
    return K                                                                                                                            # Return dense stencil matrix.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Assign @nb.njit(cache.
def _compute_cloud_stencil_jit(p: np.ndarray, vec: np.ndarray, L: np.ndarray, reg_factor: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    _compute_cloud_stencil_jit
    Compute a stencil representation (diag, w) for a 2D point cloud, compiled with Numba JIT.
    
    Input:
        p           m x 3           ndarray             Point cloud [x, y, flag].
        vec         m x nvec        ndarray[int]        Neighbor indices per node (padded with -1).
        L           5 or 6          ndarray             Operator coefficients [D, E, A, B, C] or [D, E, A, B, C, F_react].
        reg_factor                  float               Regularization factor.
    
    Output:
        diag        m               ndarray             Central weights.
        w           m x nvec        ndarray             Neighbor weights aligned with vec indices.
    """
    # 1. Variable initialization
    m             = p.shape[0]                                                                                                          # Number of nodes.
    nvec          = vec.shape[1]                                                                                                        # Neighbor slots per node.
    diag          = np.zeros(m, dtype=np.float64)                                                                                       # Diagonal weights (central coefficient).
    w             = np.zeros((m, nvec), dtype=np.float64)                                                                               # Neighbor weights for each node.
    L_derivatives = L[:5]                                                                                                               # Spatial derivative coefficients [D, E, A, B, C].
    F_react       = L[5] if L.shape[0] > 5 else 0.0                                                                                     # Reaction term coefficient F(p_i) per Eq. 190 of paper.
    
    # 2. Stencil computation
    for i in nb.prange(m):                                                                                                              # type: ignore
        if p[i, 2] == 0:                                                                                                                # Interior node: compute stencil weights.
            nvec_i = 0                                                                                                                  # Count valid neighbors.
            
            for j in range(nvec):                                                                                                       # Iterate over neighbor slots.
                if vec[i, j] != -1:                                                                                                     # If valid neighbor.
                    nvec_i += 1                                                                                                         # Increment neighbor count.
            
            if nvec_i == 0:                                                                                                             # Degenerate case: no neighbors available.
                diag[i] = 1.0                                                                                                           # Fall back to identity row.
            else:                                                                                                                       # Execute fallback branch.
                dx  = np.empty(nvec_i, dtype=np.float64)                                                                                # Neighbor x-offsets.
                dy  = np.empty(nvec_i, dtype=np.float64)                                                                                # Neighbor y-offsets.
                idx = 0                                                                                                                 # Local index.
                
                for j in range(nvec):                                                                                                   # Iterate over neighbor slots.
                    if vec[i, j] != -1:                                                                                                 # If valid neighbor.
                        vec1    = vec[i, j]                                                                                             # Fetch neighbor global index.
                        dx[idx] = p[vec1, 0] - p[i, 0]                                                                                  # Compute dx offset.
                        dy[idx] = p[vec1, 1] - p[i, 1]                                                                                  # Compute dy offset.
                        idx    += 1                                                                                                     # Increment local index.
                
                h_i = 0.0                                                                                                               # Compute local stencil scale h_i.
                for j in range(nvec_i):                                                                                                 # Find maximum neighbor distance.
                    d_j = np.sqrt(dx[j]*dx[j] + dy[j]*dy[j])                                                                            # Distance to neighbor j.
                    if d_j > h_i: h_i = d_j                                                                                             # Update h_i.
                if h_i < 1e-12: h_i = 1.0                                                                                               # Guard against division by zero.

                dx_hat = dx / h_i                                                                                                       # Non-dimensionalized x-offset.
                dy_hat = dy / h_i                                                                                                       # Non-dimensionalized y-offset.

                M = np.empty((5, nvec_i), dtype=np.float64)                                                                             # Scaled reconstruction matrix.
                for j in range(nvec_i):                                                                                                 # Populate scaled matrix.
                    M[0, j] = dx_hat[j]                                                                                                 # dx_hat term.
                    M[1, j] = dy_hat[j]                                                                                                 # dy_hat term.
                    M[2, j] = dx_hat[j]**2                                                                                              # dx_hat^2 term.
                    M[3, j] = dx_hat[j] * dy_hat[j]                                                                                     # dx_hat*dy_hat term.
                    M[4, j] = dy_hat[j]**2                                                                                              # dy_hat^2 term.

                L_scaled = np.empty(5, dtype=np.float64)                                                                                # Scaled operator coefficients.
                L_scaled[0] = L_derivatives[0] / h_i                                                                                    # D / h_i.
                L_scaled[1] = L_derivatives[1] / h_i                                                                                    # E / h_i.
                L_scaled[2] = L_derivatives[2] / (h_i**2)                                                                               # A / h_i^2.
                L_scaled[3] = L_derivatives[3] / (h_i**2)                                                                               # B / h_i^2.
                L_scaled[4] = L_derivatives[4] / (h_i**2)                                                                               # C / h_i^2.

                YY       = np.zeros(nvec_i, dtype=np.float64)                                                                           # Placeholder for neighbor weights.
                valid_YY = False                                                                                                        # Validity flag for neighbor weights.
                w_dist   = np.empty(nvec_i, dtype=np.float64)                                                                           # Distance weighting array.
                for j in range(nvec_i):                                                                                                 # Compute quadratic distance weights.
                    r2_j = dx_hat[j]**2 + dy_hat[j]**2                                                                                  # Squared distance to neighbor j.
                    w_dist[j] = 1.0 / (r2_j + 1e-4)                                                                                     # Inverse distance squared weighting (1 / r^2).

                M_W    = np.empty((5, nvec_i), dtype=np.float64)                                                                        # Distance-weighted reconstruction matrix.
                for j in range(nvec_i):                                                                                                 # Populate distance-weighted columns.
                    M_W[0, j] = M[0, j] * w_dist[j]                                                                                     # Weighted dx.
                    M_W[1, j] = M[1, j] * w_dist[j]                                                                                     # Weighted dy.
                    M_W[2, j] = M[2, j] * w_dist[j]                                                                                     # Weighted dx^2.
                    M_W[3, j] = M[3, j] * w_dist[j]                                                                                     # Weighted dx*dy.
                    M_W[4, j] = M[4, j] * w_dist[j]                                                                                     # Weighted dy^2.

                G_weighted = np.dot(M_W, M.T)                                                                                           # Distance-weighted normal matrix.
                tr         = np.trace(G_weighted)                                                                                       # Trace of weighted normal matrix.
                reg        = (1e-12 + reg_factor * tr) if tr > 0.0 else 1e-12                                                           # Adaptive Tikhonov regularization.
                
                c_coeffs   = np.linalg.solve(G_weighted + reg * np.eye(5), L_scaled)                                                    # Solve for coefficient vector c.
                YY         = w_dist * np.dot(M.T, c_coeffs)                                                                             # Compute distance-weighted stencil weights.
                valid_YY   = True                                                                                                       # Mark solution as valid.

                for j in range(nvec_i):                                                                                                 # Check solution for NaN/Inf.
                    if not np.isfinite(YY[j]):                                                                                          # If invalid element found.
                        valid_YY = False                                                                                                # Mark invalid.
                        break                                                                                                           # Stop checking.

                if not valid_YY:                                                                                                        # Fallback if matrix solve failed.
                    diag[i] = 1.0                                                                                                       # Identity diagonal.
                else:                                                                                                                   # Execute valid branch.
                    # 4. Result assignment with M-matrix non-negativity projection
                    for j in range(nvec_i):                                                                                             # Iterate over neighbor slots.
                        if YY[j] < 0.0:                                                                                                 # If negative off-diagonal weight (anti-diffusion mode).
                            YY[j] = 0.0                                                                                                 # Clamp to 0.0 to guarantee negative semi-definite spectrum.

                    diag[i] = -np.sum(YY) + F_react                                                                                     # Central weight includes reaction term per Eq. 190.
                    idx     = 0                                                                                                         # Local index.
                    
                    for j in range(nvec):                                                                                               # Iterate over neighbor slots.
                        if vec[i, j] != -1:                                                                                             # If valid neighbor.
                            w[i, j] = YY[idx]                                                                                           # Store weight aligned to vec order.
                            idx    += 1                                                                                                 # Increment local index.
                    
        elif p[i, 2] == 1 or p[i, 2] == 2:                                                                                              # Boundary node: identity row.
            diag[i] = 1.0                                                                                                               # Enforce u_i = boundary value.
            
    return diag, w                                                                                                                      # Return stencil representation.

def Cloud(p: np.ndarray, vec: np.ndarray, L: np.ndarray) -> np.ndarray:
    """
    Cloud
    Assemble a full dense stencil matrix K for a 2D point cloud.
    
    For each interior node i, this computes neighbor offsets (dx, dy), constructs the least-squares
    reconstruction matrix:
        M = [dx, dy, dx^2, dx*dy, dy^2]^T
    and solves for weights that approximate the differential operator specified by L. The resulting
    row i of K contains the stencil weights for node i and its neighbors.
    
    Boundary nodes (flag 1 or 2) use a Dirichlet row: K[i, i] = 1, and all neighbor weights are 0.
    
    Input:
        p           m x 3           array-like          Point cloud [x, y, flag].
        vec         m x nvec        array-like[int]     Neighbor indices per node (padded with -1).
        L           5 x 1           array-like          Operator coefficients [D, E, A, B, C] (shape (5,) or (5,1)).
     
    Output:
        K           m x m           ndarray             Dense matrix with stencil weights.
    
    Notes:
        This routine allocates an (m, m) dense matrix and is therefore O(m^2) in memory. For typical
        mGFD runs, prefer CloudStencil() + ApplyCloudStencil() to avoid dense assembly.
    """
    # 1. Variable initialization and type enforcement
    p_np   = np.asarray(p, dtype=np.float64)                                                                                            # Enforce float64 type for points.
    vec_np = np.asarray(vec, dtype=np.int64)                                                                                            # Enforce int64 type for neighbor indices.
    L_np   = np.asarray(L, dtype=np.float64).flatten()                                                                                  # Enforce float64 type and flatten operator.
    
    # 2. Gammas computation and matrix assembly via JIT
    K = _compute_cloud_dense_jit(p_np, vec_np, L_np)                                                                                    # Delegate to JIT compiled dense solver.
    
    return K                                                                                                                            # Return dense stencil matrix.

def RHS(p: np.ndarray, boun_n: np.ndarray, inne_n: np.ndarray, phi: Union[Callable, np.ndarray, float, int], f: Union[Callable, np.ndarray, float, int]) -> np.ndarray:
    """
    RHS
    Build the right-hand-side vector for a Dirichlet problem on a point cloud.
    
    Interior entries are filled from the forcing f(x, y). Boundary entries are filled from phi(x, y).
    
    Input:
        p           m x 3           ndarray             Point cloud [x, y, flag].
        boun_n      m               ndarray             Boundary mask (True on boundary nodes).
        inne_n      m               ndarray             Interior mask (True on interior nodes).
        phi                         Union               Boundary condition phi(x, y) or constant/array data.
        f                           Union               Forcing function f(x, y) or constant/array data.
    
    Output:
        R           m               ndarray             RHS vector consistent with Dirichlet enforcement.
    """
    # 1. Variable initialization
    m         = len(p[:, 0])                                                                                                            # Total number of nodes.
    R         = np.zeros([m])                                                                                                           # RHS initialization.

    # 2. Vector assignment
    if callable(f):                                                                                                                     # If forcing is a function.
        R[inne_n] = f(p[inne_n, 0], p[inne_n, 1])                                                                                       # Interior forcing term via function.
    elif isinstance(f, np.ndarray):                                                                                                     # If forcing is data array.
        R[inne_n] = f[inne_n]                                                                                                           # Interior forcing term via data.
    elif isinstance(f, (int, float)):                                                                                                   # If forcing is constant.
        R[inne_n] = f                                                                                                                   # Interior forcing term via constant.

    if callable(phi):                                                                                                                   # If boundary condition is a function.
        R[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                                                                     # Boundary Dirichlet values via function.
    elif isinstance(phi, np.ndarray):                                                                                                   # If boundary condition is data array.
        R[boun_n] = phi[boun_n]                                                                                                         # Boundary Dirichlet values via data.
    elif isinstance(phi, (int, float)):                                                                                                 # If boundary condition is constant.
        R[boun_n] = phi                                                                                                                 # Boundary Dirichlet values via constant.

    return R                                                                                                                            # Return assembled RHS.

def CloudStencil(p: np.ndarray, vec: np.ndarray, L: np.ndarray, reg_factor: float = 1e-4) -> Tuple[np.ndarray, np.ndarray]:
    """
    CloudStencil
    Compute a stencil representation (diag, w) for a 2D point cloud without assembling a dense matrix.
    
    For each interior node i, this returns:
        diag[i]          the central weight
        w[i, :nvec_i]    the weights for each neighbor in vec[i, :nvec_i]
    
    Boundary nodes (flag 1 or 2) are set to an identity row: diag[i] = 1 and all neighbor weights 0.
    
    Input:
        p           m x 3           ndarray             Point cloud [x, y, flag].
        vec         m x nvec        ndarray             Neighbor indices per node (padded with -1).
        L           5 x 1           ndarray             Operator coefficients [D, E, A, B, C].
        reg_factor                  float               Regularization factor (default: 1e-4).
    
    Output:
        diag        m               ndarray             Central weights.
        w           m x nvec        ndarray             Neighbor weights aligned with vec indices.
    
    Notes:
        The interior weights are computed by solving a regularized normal equation on M @ M.T
        (a 5x5 system) when possible, and falling back to a pseudoinverse if needed.
    """
    # 1. Variable initialization and type enforcement
    p_np       = np.asarray(p, dtype=np.float64)                                                                                        # Enforce float64 type for points.
    vec_np     = np.asarray(vec, dtype=np.int64)                                                                                        # Enforce int64 type for neighbor indices.
    L_np       = np.asarray(L, dtype=np.float64).flatten()                                                                              # Enforce float64 type and flatten operator.
    reg_factor = reg_factor                                                                                                             # Enforce float type for regularization parameter.

    # 2. CloudStencil computation via JIT
    diag, w = _compute_cloud_stencil_jit(p_np, vec_np, L_np, reg_factor)                                                                # Delegate to JIT compiled stencil solver.

    return diag, w                                                                                                                      # Return stencil representation (central + neighbor weights).

def ApplyCloudStencil(u: np.ndarray, vec: np.ndarray, diag: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    ApplyCloudStencil
    Apply a precomputed cloud stencil (diag, w) to a vector u.
    
    The neighbor indexing vec uses -1 to indicate unused slots. This implementation maps all negative
    indices to an extra padding entry u_pad[m] = 0, so that missing neighbor slots contribute zero.
    
    Input:
        u           m               ndarray             Values at nodes.
        vec         m x nvec        ndarray             Neighbor indices per node (padded with -1).
        diag        m               ndarray             Central weights.
        w           m x nvec        ndarray             Neighbor weights aligned with vec.
    
    Output:
        Lu          m               ndarray             Stencil application result.
    """
    # 1. Variable initialization
    m            = int(u.shape[0])                                                                                                      # Number of nodes.
    idx          = vec.copy()                                                                                                           # Copy indices to safely rewrite padding entries.
    idx[idx < 0] = m                                                                                                                    # Map padding (-1) to the extra padded index m.
    u_pad        = np.concatenate([u, np.zeros(1, dtype = u.dtype)])                                                                    # Append u_pad[m]=0 for missing neighbors.
    neigh_vals   = u_pad[idx]                                                                                                           # Gather neighbor values with padding handled.
    
    # 2. Stencil application
    return diag * u + np.sum(w * neigh_vals, axis = 1)                                                                                  # Apply diag + weighted neighbor sum.

def BiCGStab(matvec: Callable, b: np.ndarray, x0: Optional[np.ndarray] = None, tol: float = 1e-10, max_iter: int = 2000) -> np.ndarray:
    """
    BiCGStab
    Solve the linear system A x = b using the BiCGStab Krylov method (matrix-free).
    
    Input:
        matvec                      Callable            Function that returns A @ x for a given x.
        b           m               ndarray             Right-hand side vector.
        x0          m               ndarray             Optional initial guess.
        tol                         float               Relative tolerance used by stopping criteria.
        max_iter                    int                 Maximum number of iterations.
    
    Output:
        x                           ndarray             Approximate solution vector.
    
    Notes:
        This implementation first tries to use SciPy's bicgstab through a LinearOperator. If SciPy
        is unavailable or the SciPy call fails, it falls back to a pure NumPy implementation.
    """
    # 1. Try SciPy Implementation
    try:                                                                                                                                # Prefer SciPy when available (faster/robust).
        b_np = np.asarray(b, dtype = float)                                                                                             # Normalize RHS to float array.
        n    = int(b_np.shape[0])                                                                                                       # System size.
        
        class MatrixFreeOp(LinearOperator):                                                                                             # Define a custom LinearOperator subclass.
            """
            MatrixFreeOp
            Custom LinearOperator subclass for SciPy iterative solvers.
            """
            def __init__(self):
                """
                __init__
                Initialize the LinearOperator required by scipy.sparse.linalg.lsqr.
                Defines the shape and datatype of the implicit operator matrix.
                """
                self.shape = (n, n)                                                                                                     # Shape of the matrix.
                self.dtype = np.dtype(float)                                                                                            # Set operator data type.
            
            def _matvec(self, x):                                                                       
                """
                _matvec
                Matrix-vector multiplication method.
                
                Input:
                    x           (n,)            ndarray             Vector to multiply.
                    
                Output:
                                (n,)            ndarray             Result of A @ x.
                """
                return matvec(x)                                                                                                        # Delegate to the closure function.

        A    = MatrixFreeOp()                                                                                                           # Matrix-free linear operator.
        
        if x0 is None:                                                                                                                  # Optional initial guess handling.
            x0_np   = None                                                                                                              # No initial guess passed to SciPy.
        else:                                                                                                                           # Initial guess provided by caller.
            x0_np   = np.asarray(x0, dtype = float)                                                                                     # Normalize initial guess to float array.

        try:                                                                                                                            # Old SciPy API uses tol.
            x, info = bicgstab(A, b_np, x0 = x0_np, atol = tol, maxiter = max_iter)                                                     # Attempt SciPy BiCGStab (old signature).
        except TypeError:                                                                                                               # New SciPy API uses rtol/atol.
            x, info = bicgstab(A, b_np, x0 = x0_np, rtol = tol, atol = 0.0, maxiter = max_iter)                                         # Attempt SciPy BiCGStab (new signature).

        if info == 0:                                                                                                                   # Converged solve.
            return np.asarray(x, dtype = float)                                                                                         # Return SciPy solution.
    
    except Exception:                                                                                                                   # SciPy missing or solver failed unexpectedly.
        pass                                                                                                                            # Fall back to pure NumPy implementation below.

    # 2. Pure NumPy Implementation initialization
    if x0 is None:                                                                                                                      # Initialize x for pure NumPy BiCGStab.
        x = np.zeros_like(b, dtype = float)                                                                                             # Default initial guess: zeros.
    else:                                                                                                                               # Use provided initial guess as starting iterate.
        x = np.asarray(x0, dtype = float).copy()                                                                                        # Use provided initial guess.

    b     = np.asarray(b, dtype = float)                                                                                                # Normalize RHS to float array.
    r     = b - matvec(x)                                                                                                               # Initial residual.
    r0    = r.copy()                                                                                                                    # Shadow residual (fixed).
    rho   = 1.0                                                                                                                         # BiCGStab scalar rho.
    alpha = 1.0                                                                                                                         # BiCGStab scalar alpha.
    omega = 1.0                                                                                                                         # BiCGStab scalar omega.
    v     = np.zeros_like(b, dtype = float)                                                                                             # v vector.
    p     = np.zeros_like(b, dtype = float)                                                                                             # p search direction.

    # 3. Check for trivial solution
    b_norm = np.linalg.norm(b)                                                                                                          # Norm of RHS for relative stopping.
    if b_norm == 0:                                                                                                                     # Trivial system with zero RHS.
        return x                                                                                                                        # x=0 is a valid solution.

    # 4. Iterate BiCGStab steps
    for _ in range(max_iter):                                                                                                           # Iterate BiCGStab steps.
        rho_new = float(np.dot(r0, r))                                                                                                  # rho_new = <r0, r>.
        
        if rho_new == 0:                                                                                                                # Breakdown condition.
            break                                                                                                                       # Exit iterations.
        
        beta  = (rho_new / rho) * (alpha / omega)                                                                                       # Update coefficient beta.
        p     = r + beta * (p - omega * v)                                                                                              # Update search direction.
        v     = matvec(p)                                                                                                               # Apply A to p.
        denom = float(np.dot(r0, v))                                                                                                    # denom = <r0, v>.
        
        if denom == 0:                                                                                                                  # Breakdown condition.
            break                                                                                                                       # Exit iterations.
        
        alpha = rho_new / denom                                                                                                         # Update alpha.
        s     = r - alpha * v                                                                                                           # Intermediate residual.
        
        if np.linalg.norm(s) / b_norm < tol:                                                                                            # Converged after the alpha step.
            x = x + alpha * p                                                                                                           # Update solution.
            break                                                                                                                       # Exit iterations.
        
        t  = matvec(s)                                                                                                                  # Apply A to s.
        tt = float(np.dot(t, t))                                                                                                        # tt = <t, t>.
        
        if tt == 0:                                                                                                                     # Degenerate update (t is zero).
            break                                                                                                                       # Exit iterations.
        
        omega = float(np.dot(t, s)) / tt                                                                                                # Update omega.
        x     = x + alpha * p + omega * s                                                                                               # Update solution.
        r     = s - omega * t                                                                                                           # Update residual.
        
        if np.linalg.norm(r) / b_norm < tol:                                                                                            # Check convergence.
            break                                                                                                                       # Exit iterations.
        
        rho = rho_new                                                                                                                   # Carry rho forward.

    return x                                                                                                                            # Return final iterate.

@nb.njit(parallel=True, fastmath=True)                                                                                                  # Assign @nb.njit(parallel.
def _compute_K_matvec_numba(x: np.ndarray, diag: np.ndarray, w: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """
    _compute_K_matvec_numba
    Numba JIT-compiled closure helper to compute the matrix-vector multiplication K * x in parallel.
    
    Input:
        x           m               ndarray         Input vector to multiply.
        diag        m               ndarray         Main diagonal elements of the implicit sparse matrix.
        w           m x nvec        ndarray         Stencil weight elements (off-diagonal).
        vec         m x nvec        ndarray         Neighbor list mapping for the off-diagonal elements.
        
    Output:
        y           m               ndarray         Resulting vector from K * x.
    """
    # 1. Variable initialization
    m    = x.shape[0]                                                                                                                   # Total number of nodes.
    nvec = vec.shape[1]                                                                                                                 # Maximum number of neighbors.
    y    = np.zeros(m, dtype=np.float64)                                                                                                # Output vector initialization.
    
    # 2. Parallel Matrix-Vector Multiplication
    for i in nb.prange(m):                                                                                                              # type: ignore # Iterate over all nodes.
        val = diag[i] * x[i]                                                                                                            # Initialize with the main diagonal contribution.
        for j in range(nvec):                                                                                                           # Iterate over the node's stencil.
            col = vec[i, j]                                                                                                             # Fetch target neighbor index.
            if col != -1:                                                                                                               # If valid neighbor slot.
                val += w[i, j] * x[col]                                                                                                 # Accumulate neighbor contribution.
        y[i] = val                                                                                                                      # Save computed row value.
        
    return y                                                                                                                            # Return resulting vector.

def compute_K_matvec(p: np.ndarray, vec: np.ndarray, L: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    compute_K_matvec
    Generates an on-the-fly matrix-vector product closure for K * x without building K explicitly.

    Input:
        p           m x 3           ndarray             Array with the coordinates of the nodes and the boundary flag.
        vec         m x nvec        ndarray             Neighbor indices per node (padded with -1).
        L           5               ndarray             Array with the weights for the operator.

    Output:
        matvec                      Callable            A function that accepts an array x and returns K * x.
    """
    diag, w = CloudStencil(p, vec, L, reg_factor=1e-8)                                                                                  # Retrieve the raw stencil weight values.
    
    def matvec(x: np.ndarray) -> np.ndarray:                                                                                            # Define the closure over x.
        """
        matvec
        
        Closure function for the matrix-vector multiplication K * x without explicitly building K.
        
        Input:
            x           m               ndarray         Input vector to be multiplied.
            
        Output:
            y           m               ndarray         Resulting vector from K * x.
        """
        return _compute_K_matvec_numba(x, diag, w, vec)                                                                                 # Delegate to Numba JIT loop.
        
    return matvec                                                                                                                       # Return the callable operator.

def compute_sparse_matrix(p: np.ndarray, vec: np.ndarray, L: np.ndarray) -> csr_matrix:
    """
    compute_sparse_matrix
    Builds a sparse csr_matrix for the stencil K directly, avoiding dense allocation.

    This function utilizes the CloudStencil function to retrieve the raw stencil weights
    and constructs the scipy.sparse.csr_matrix efficiently without looping linearly over every single index.

    Input:
        p           m x 3           ndarray             Array with the coordinates of the nodes and the boundary flag.
        vec         m x nvec        ndarray             Neighbor indices per node (padded with -1).
        L           5               ndarray             Array with the weights for the operator.

    Output:
        K           m x m           csr_matrix          The highly-efficient Sparse operator representation.
    """
    # 1. Stencil calculation
    diag, w = CloudStencil(p, vec, L, reg_factor=1e-8)                                                                                  # Retrieve the raw stencil weight values.
    m       = p.shape[0]                                                                                                                # Obtain the total number of nodes.

    # 2. Masking valid neighbors
    valid   = vec != -1                                                                                                                 # Create a boolean mask of valid neighbors.
    rows    = np.repeat(np.arange(m), vec.shape[1]).reshape(m, -1)[valid]                                                               # Map every valid neighbor to its origin row index.
    cols    = vec[valid]                                                                                                                # Obtain the target column indices corresponding to the neighbors.
    data    = w[valid]                                                                                                                  # Obtain the neighbor-specific weight data values.

    # 3. Diagonal concatenation
    rows    = np.concatenate([rows, np.arange(m)])                                                                                      # Append the central row index positions for the main diagonal.
    cols    = np.concatenate([cols, np.arange(m)])                                                                                      # Append the central column index positions for the main diagonal.
    data    = np.concatenate([data, diag])                                                                                              # Append the main diagonal (central) weights.

    return csr_matrix((data, (rows, cols)), shape=(m, m))                                                                               # Construct and return the compressed sparse row matrix.


