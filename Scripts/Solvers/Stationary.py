"""
Stationary — Stationary PDEs solver

Overview:
    Numerical solver for stationary PDEs (no time derivatives) using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    Stationary          Main solver function for stationary problems.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
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
import numpy as np

from Scripts import Gammas
from Scripts.Solvers.Utils import *

try:
    from scipy.sparse import coo_matrix, eye
    from scipy.sparse.linalg import LinearOperator, bicgstab, spsolve, splu
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def Stationary(p, phi, f, operator = np.vstack([[0], [0], [2], [0], [2]]), Adv = False, vec = None, nvec = 8):
    """
    Stationary
    Numerical solution of stationary PDEs (no time derivatives) using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

    The problem solved is:
        D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy = -f(x, y)
    with Dirichlet boundary conditions:
        u(x, y) = phi(x, y) on boundary nodes.

    Input:
        p           m x 3           ndarray         Point cloud with columns [x, y, flag].
        phi                         function        Dirichlet boundary function phi(x, y).
        f                           function        Right-hand side function f(x, y).
        operator                    ndarray         Operator coefficients. Use a 6-entry vector [D, E, A, B, C, F].
                                                    Only [D, E, A, B, C] are used to build the stencil.
        Adv                         bool            If True, uses upwind-biased neighbor selection based on (D, E).
        vec                         ndarray|None    Optional neighbor list (m x nvec, padded with -1).
        nvec                        int             Maximum number of neighbors used by the stencil.

    Output:
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        u_ex        m               ndarray         Array with the theoretical solution (phi evaluated at all nodes).
        vec         m x o           ndarray         Neighbor list used by the solver.

    Notes:
        If SciPy is available, the routine attempts a matrix-free BiCGStab solve first; otherwise it
        falls back to a custom BiCGStab implementation.
    """
    ## Variable initialization.
    m              = int(p.shape[0])                                                                    # Total number of nodes.
    nvec           = int(nvec)                                                                          # Requested maximum neighbors.
    u_ap           = np.zeros([m])                                                                      # Approximate solution.
    u_ex           = np.zeros([m])                                                                      # Theoretical solution.
    boun_n, inne_n = _node_masks(p)                                                                     # Boundary and interior masks.

    ## Boundary conditions.
    u_ap[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                                      # Impose Dirichlet boundary values.

    ## Neighbor computation / validation.
    vec, nvec    = _prepare_neighbors(p, vec, nvec, Adv, operator)                                      # Build/validate neighbor list.

    ## Gamma (stencil) computation.
    L            = operator[:-1]                                                                        # Operator coefficients without reaction term in RHS builder.
    diag, w      = Gammas.CloudStencil(p, vec, L)                                                       # Stencil diagonal and weights.

    ## Right-hand side.
    R            = Gammas.RHS(p, boun_n, inne_n, phi, f)                                                # RHS vector with boundary values included.
    
    u_ap[boun_n] = R[boun_n]                                                                            # Re-impose boundary values.
    if not np.any(inne_n):                                                                              # Early exit when there are no interior unknowns.
        u_ex = phi(p[:, 0], p[:, 1])                                                                    # Theoretical solution on all nodes.
        return u_ap, u_ex, vec                                                                          # Nothing to solve if no interior nodes.

    ## Linear system on interior nodes only.
    inn_idx    = np.where(inne_n)[0]                                                                    # Global indices of interior nodes.
    u0         = u_ap.copy()                                                                            # Baseline vector with boundary values fixed.
    stencil_u0 = Gammas.ApplyCloudStencil(u0, vec, diag, w)[inn_idx]                                    # Constant part from boundary contributions.
    rhs        = R[inn_idx] - stencil_u0                                                                # Effective RHS for interior unknowns.

    def matvec_int(x):                                                                                  # Matrix-vector product restricted to interior nodes.
        u          = u0.copy()                                                                          # Start from boundary-fixed vector.
        u[inn_idx] = x                                                                                  # Insert interior unknowns.
        return Gammas.ApplyCloudStencil(u, vec, diag, w)[inn_idx] - stencil_u0                          # Apply operator and remove constant term.

    x = None                                                                                            # Placeholder for interior solution.
    try:                                                                                                # Prefer SciPy sparse solvers when available.
        if not HAS_SCIPY:
            raise ImportError("SciPy is required for this operation")                                   # Load SciPy operators/solvers lazily.
        
        class MatrixFreeOp(LinearOperator):                                                             # Define a custom LinearOperator subclass.
            def __init__(self):                                                                         # Initialize operator.
                self.shape = (inn_idx.size, inn_idx.size)                                               # Set operator shape.
                self.dtype = np.dtype(float)                                                            # Set operator data type.
            def _matvec(self, x):                                                                       # Matrix-vector multiplication method.
                return matvec_int(x)                                                                    # Delegate to the closure function.
        
        Aop = MatrixFreeOp()                                                                            # Matrix-free linear operator on interior nodes.
        try:                                                                                            # Attempt solve with SciPy BiCGStab.
            x_it, info = bicgstab(Aop, rhs, tol=1e-10, maxiter=5000)                                    # Older SciPy API.
        except TypeError:                                                                               # Newer SciPy API uses rtol/atol.
            x_it, info = bicgstab(Aop, rhs, rtol=1e-10, atol=0.0, maxiter=5000)                         # Newer SciPy API.

        if info == 0:                                                                                   # Converged iterative solve.
            x = np.asarray(x_it, dtype=float)                                                           # Converged iterative solution.
        else:                                                                                           # Assemble sparse matrix and solve directly.
            inn_map          = np.full(m, -1, dtype=np.int32)                                           # Global -> local mapping for interior nodes.
            inn_map[inn_idx] = np.arange(inn_idx.size, dtype=np.int32)                                  # Fill mapping.
            rows = []                                                                                   # Sparse COO rows.
            cols = []                                                                                   # Sparse COO cols.
            data = []                                                                                   # Sparse COO data.
            for row_k, gi in enumerate(inn_idx):                                                        # For each interior equation row.
                rows.append(row_k)                                                                      # Diagonal entry row.
                cols.append(row_k)                                                                      # Diagonal entry col.
                data.append(diag[gi])                                                                   # Diagonal weight.
                for t, gj in enumerate(vec[gi]):                                                        # Iterate neighbors of this node.
                    if gj == -1:                                                                        # Missing neighbor slot.
                        continue                                                                        # Skip missing neighbors.
                    j     = int(gj)                                                                     # Neighbor global index.
                    col_k = int(inn_map[j])                                                             # Neighbor local index (if interior).
                    if col_k != -1:                                                                     # Keep only interior-interior couplings.
                        rows.append(row_k)                                                              # Off-diagonal entry row.
                        cols.append(col_k)                                                              # Off-diagonal entry col.
                        data.append(w[gi, t])                                                           # Neighbor weight.
            A = coo_matrix((data, (rows, cols)), shape=(inn_idx.size, inn_idx.size)).tocsr()            # Assemble sparse matrix for interior system.
            x = spsolve(A, rhs)                                                                         # Direct sparse solve.
    except ImportError:                                                                                 # SciPy not available (module can still solve via custom BiCGStab).
        pass                                                                                            # Fall back to custom BiCGStab below.
    except Exception:                                                                                   # SciPy solver failed unexpectedly.
        pass                                                                                            # Fall back to custom BiCGStab below.

    if x is None:                                                                                       # If SciPy path failed, run the custom solver.
        x = Gammas.BiCGStab(matvec_int, rhs, tol=1e-10, max_iter=5000)                                  # Custom BiCGStab implementation.
    u_ap[inn_idx] = x                                                                                   # Insert interior solution.
    
    ## Theoretical solution.
    u_ex = phi(p[:, 0], p[:, 1])                                                                        # Evaluate exact/boundary function.

    return u_ap, u_ex, vec                                                                              # Return approximate/exact solutions and neighbor list.