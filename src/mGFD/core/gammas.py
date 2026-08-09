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

from scipy.optimize import nnls                                                                                                         # Core numerical operations.
from scipy.sparse import csr_matrix                                                                                                     # Core sparse matrix representations.
from typing import Callable, Optional, Tuple, List                                                                                      # Type hinting.
from scipy.sparse.linalg import LinearOperator, bicgstab                                                                                # SciPy iterative solver interface.

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
    # 1. Variable initialization.
    nvec = int(len(vec[0, :]))                                                                                                          # The maximum number of neighbors (vec width).
    m    = int(len(p[:, 0]))                                                                                                            # Total number of nodes in the cloud.
    K    = np.zeros([m, m])                                                                                                             # Dense stencil matrix initialization.
    
    # 2. Gammas computation and matrix assembly.
    for i in np.arange(m):                                                                                                              # For each node in the cloud.
        if p[i, 2] == 0:                                                                                                                # Interior node: compute a stencil row.
            nvec_i = int(np.sum(vec[i, :] != -1))                                                                                       # Count valid neighbors (skip padding -1).
            dx     = np.zeros([nvec_i])                                                                                                 # Neighbor x-offsets (relative to node i).
            dy     = np.zeros([nvec_i])                                                                                                 # Neighbor y-offsets (relative to node i).
            
            for j in np.arange(nvec_i):                                                                                                 # For each neighbor of node i.
                vec1  = int(vec[i, j])                                                                                                  # Neighbor global index.
                dx[j] = p[vec1, 0] - p[i, 0]                                                                                            # Compute dx to neighbor.
                dy[j] = p[vec1, 1] - p[i, 1]                                                                                            # Compute dy to neighbor.
            
            M       = np.vstack([[dx], [dy], [dx**2], [dx * dy], [dy**2]])                                                              # Reconstruction matrix (5, nvec_i).
            M       = np.linalg.pinv(M)                                                                                                 # Pseudoinverse to solve least squares.
            YY      = M @ L                                                                                                             # Weights for neighbor contributions (nvec_i,).
            Gamma   = np.vstack([-np.sum(YY), YY]).transpose()                                                                          # Include central weight as -sum(neighbors).
            K[i, i] = Gamma[0, 0]                                                                                                       # Set diagonal (central) weight.
            
            for j in np.arange(nvec_i):                                                                                                 # Fill neighbor weights in row i.
                K[i, int(vec[i, j])] = Gamma[0, j + 1]                                                                                  # Set weight for neighbor j.
            
        if p[i, 2] == 1 or p[i, 2] == 2:                                                                                                # Boundary node: Dirichlet identity row.
            K[i, i] = 1                                                                                                                 # Enforce u_i = boundary value.
            
            for j in np.arange(nvec):                                                                                                   # For each neighbor slot.
                if int(vec[i, j]) != -1:                                                                                                # Skip padding indices.
                    K[i, int(vec[i, j])] = 0                                                                                            # Zero out neighbor weights.
    
    return K                                                                                                                            # Return dense stencil matrix.

def RHS(p: np.ndarray, boun_n: np.ndarray, inne_n: np.ndarray, phi: Callable, f: Callable) -> np.ndarray:
    """
    RHS
    Build the right-hand-side vector for a Dirichlet problem on a point cloud.
    
    Interior entries are filled from the forcing f(x, y). Boundary entries are filled from phi(x, y).
    
    Input:
        p           m x 3           ndarray             Point cloud [x, y, flag].
        boun_n      m               ndarray             Boundary mask (True on boundary nodes).
        inne_n      m               ndarray             Interior mask (True on interior nodes).
        phi                         Callable            Boundary function phi(x, y).
        f                           Callable            Forcing function f(x, y).
    
    Output:
        R           m               ndarray             RHS vector consistent with Dirichlet enforcement.
    """
    # 1. Variable initialization
    m         = int(len(p[:, 0]))                                                                                                       # Total number of nodes.
    R         = np.zeros([m])                                                                                                           # RHS initialization.

    # 2. Vector assignment
    R[inne_n] = f(p[inne_n, 0], p[inne_n, 1])                                                                                           # Interior forcing term.
    R[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                                                                         # Boundary Dirichlet values.

    return R                                                                                                                            # Return assembled RHS.

def CloudStencil(p: np.ndarray, vec: np.ndarray, L: np.ndarray, reg_factor: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
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
        reg_factor                  float               Regularization factor (default: 1e-8).
    
    Output:
        diag        m               ndarray             Central weights.
        w           m x nvec        ndarray             Neighbor weights aligned with vec indices.
    
    Notes:
        The interior weights are computed by solving a regularized normal equation on M @ M.T
        (a 5x5 system) when possible, and falling back to a pseudoinverse if needed.
    """
    # 1. Variable initialization
    nvec = int(vec.shape[1])                                                                                                            # Neighbor slots per node.
    m    = int(p.shape[0])                                                                                                              # Number of nodes.
    diag = np.zeros(m, dtype = float)                                                                                                   # Diagonal weights (central coefficient).
    w    = np.zeros((m, nvec), dtype = float)                                                                                           # Neighbor weights for each node.

    # 2. CloudStencil computation
    for i in np.arange(m):                                                                                                              # Loop over nodes.
        if p[i, 2] == 0:                                                                                                                # Interior node: compute stencil weights.
            nvec_i = int(np.sum(vec[i, :] != -1))                                                                                       # Number of valid neighbors for node i.
            
            if nvec_i == 0:                                                                                                             # Degenerate case: no neighbors available.
                diag[i] = 1.0                                                                                                           # Fall back to identity row.
                continue                                                                                                                # Continue to next node.
            
            neigh = vec[i, :nvec_i].astype(int)                                                                                         # Neighbor indices (trim padding).
            dx    = p[neigh, 0] - p[i, 0]                                                                                               # Neighbor dx offsets.
            dy    = p[neigh, 1] - p[i, 1]                                                                                               # Neighbor dy offsets.
            M     = np.vstack([dx, dy, dx**2, dx * dy, dy**2])                                                                          # Reconstruction matrix (5, nvec_i).
            YY    = None                                                                                                                # Placeholder for neighbor weights.
            
            try:                                                                                                                        # 1. Try NNLS for strictly positive weights
                G_trace = float(np.trace(M @ M.T))                                                                                      # Trace of the M @ M.T matrix.
                lam     = (1e-12 + reg_factor * G_trace) if G_trace > 0.0 else 1e-12                                                    # Regularization parameter.
                A_aug   = np.vstack([M, np.sqrt(lam) * np.eye(nvec_i)])                                                                 # Augmented matrix for NNLS.
                b_aug   = np.concatenate([np.array(L).reshape(-1), np.zeros(nvec_i)])                                                   # Augmented RHS vector for NNLS.
                YY, _   = nnls(A_aug, b_aug)                                                                                            # Non-negative least squares solve.
                
                if np.sum(YY) < 1e-10:                                                                                                  # Pure advection failure fallback
                    YY = None
            
            except Exception as e:                                                                                                      # Catch any exception during NNLS.
                YY = None                                                                                                               # Nullify weights on failure.
            
            if YY is None or (not np.all(np.isfinite(YY))):                                                                             # 2. Fallback to standard Mínimos Cuadrados
                try:                                                                                                                    # Try regularized least squares.
                    G   = M @ M.T                                                                                                       # Normal matrix.
                    tr  = float(np.trace(G))                                                                                            # Trace of the normal matrix.
                    reg = (1e-12 + reg_factor * tr) if tr > 0.0 else 1e-12                                                              # Adaptive regularization parameter.
                    c   = np.linalg.solve(G + reg * np.eye(5), L)                                                                       # Solve for regularized coefficients.
                    YY  = (M.T @ c).reshape(-1)                                                                                         # Compute final neighbor weights.
                except Exception:                                                                                                       # Catch linear algebra errors.
                    YY = None                                                                                                           # Nullify weights on failure.

            if YY is None or (not np.all(np.isfinite(YY))):                                                                             # 3. Absolute fallback to pseudoinverse
                YY = (np.linalg.pinv(M) @ L).reshape(-1)                                                                                # Direct pseudoinverse solve.
                
            if not np.all(np.isfinite(YY)):                                                                                             # If still invalid, use identity row.
                diag[i] = 1.0                                                                                                           # Identity diagonal.
                continue                                                                                                                # Skip neighbor weights.
            diag[i]        = -float(np.sum(YY))                                                                                         # Central weight ensures sum-to-zero structure.
            w[i, :nvec_i]  = YY                                                                                                         # Store neighbor weights aligned to vec order.
        
        elif p[i, 2] == 1 or p[i, 2] == 2:                                                                                              # Boundary node: identity row.
            diag[i] = 1.0                                                                                                               # Enforce u_i = boundary value.

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
    for _ in range(int(max_iter)):                                                                                                      # Iterate BiCGStab steps.
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
