"""
TimeDerivative2 — Second-order transient PDEs solver

Overview:
    Numerical solver for PDEs with a second-order time derivative using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    TimeDerivative2     Main solver function for second-order transient problems.

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
    from scipy.sparse.linalg import LinearOperator, bicgstab, spsolve, splu
    from scipy.sparse import coo_matrix, eye
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def TimeDerivative2(p, f, g, t, coef, operator = np.vstack([[0], [0], [2], [0], [2]]), implicit = False, lam = 0.5, Adv = False, vec = None, nvec = 8):
    """
    TimeDerivative2
    Numerical solution of PDEs with a second-order time derivative using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

    The problem solved is:
        u_tt = D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy

    Input:
        p           m x 3           ndarray         Point cloud with columns [x, y, flag].
        f                           function        Reference solution / boundary function f(x, y, t, coef).
                                                    Boundary and initial displacement are taken from this function.
        g                           function        Initial velocity function g(x, y, t, coef) for u_t.
        t                           int             Number of time steps to be considered.
        coef                        ndarray         Coefficients of the specific problem.
        operator                    ndarray         Operator coefficients. Use a 6-entry vector [D, E, A, B, C, F].
                                                    Only [D, E, A, B, C] are used to build the stencil.
        implicit                    bool            Select whether or not use an implicit scheme.
                                                        True: Implicit scheme used.
                                                        False: Explicit scheme used (Default).
        lam                         float           Lambda parameter for the implicit scheme.
                                                        lam = 0.0 gives a fully implicit formulation.
                                                        0 < lam < 1 gives a semi-implicit formulation.
        Adv                         bool            If True, uses upwind-biased neighbor selection based on (D, E).
        vec                         ndarray|None    Optional neighbor list (m x nvec, padded with -1).
        nvec                        int             Maximum number of neighbors used by the stencil.

    Output:
        u_ap        m x t           ndarray         Array with the approximation computed by the routine.
        u_ex        m x t           ndarray         Array with the theoretical solution f(x, y, t, coef).
        vec         m x o           ndarray         Neighbor list used by the solver.

    Notes:
        The explicit wave update can be unstable depending on geometry/stencil. The routine monitors
        growth and switches to an implicit solve (optionally with stencil expansion) when needed.
    """
    
    ## Variable initialization.
    m              = int(p.shape[0])                                                                    # Total number of nodes.
    nvec           = int(nvec)                                                                          # Requested maximum neighbors.
    T              = np.linspace(0, 1, t)                                                               # Time grid on [0, 1].
    dt             = T[1] - T[0]                                                                        # Time step size.
    u_ap           = np.zeros([m, t])                                                                   # Approximate solution in time.
    u_ex           = np.zeros([m, t])                                                                   # Theoretical solution in time.
    boun_n, inne_n = _node_masks(p)                                                                     # Boundary and interior masks.

    ## Boundary conditions over time.
    for k in np.arange(t):                                                                              # For each time step.
        u_ap[boun_n, k] = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                                     # Impose boundary values at time T[k].

    ## Initial condition.
    u_ap[:, 0] = f(p[:, 0], p[:, 1], T[0], coef)                                                        # Initial state at t = 0.

    ## Neighbor computation / validation.
    vec, nvec = _prepare_neighbors(p, vec, nvec, Adv, operator)                                         # Build/validate neighbor list.

    ## Stencil for wave operator scaled by dt^2.
    L       = (dt**2) * operator[:-1]                                                                   # Scale spatial operator by dt^2.
    diag, w = Gammas.CloudStencil(p, vec, L)                                                            # Diagonal and neighbor weights.

    def unstable_state(arr, ref_scale):                                                                 # Detect instability for the explicit wave update.
        if not np.all(np.isfinite(arr)):                                                                # NaN/Inf check.
            return True                                                                                 # Unstable by definition.
        amax = float(np.max(np.abs(arr))) if arr.size else 0.0                                          # Max absolute amplitude.
        return (not np.isfinite(amax)) or (amax > 1e3 * max(1.0, ref_scale))                            # Growth relative to boundary/initial scale.

    # The explicit scheme for wave equation can be unstable depending on stencil/geometry.
    # Try explicit first (when requested), but fall back to implicit if growth is detected.
    fallback_to_implicit = bool(implicit)                                                               # Start in implicit if user requested it.
    if not fallback_to_implicit:                                                                        # Attempt explicit scheme first.
        g1              = g(p[:, 0], p[:, 1], T[1], coef)                                               # Initial velocity term contribution.
        Ku0             = Gammas.ApplyCloudStencil(u_ap[:, 0], vec, diag, w)                            # Operator applied to u at t=0.
        un1             = u_ap[:, 0] + 0.5 * Ku0 + dt * g1                                              # Second-order start step (k=1).
        u_ap[inne_n, 1] = un1[inne_n]                                                                   # Save interior values.

        ref_scale = float(np.max(np.abs(u_ap[boun_n, 1]))) if np.any(boun_n) else float(np.max(np.abs(u_ap[:, 0])))
                                                                                                        # Reference scale after first step.
        if unstable_state(u_ap[inne_n, 1], ref_scale):                                                  # Check stability after the first step.
            fallback_to_implicit = True                                                                 # Switch to implicit if unstable.
        else:                                                                                           # Proceed with explicit update.
            for k in np.arange(2, t):                                                                   # Time loop for explicit update.
                Ku              = Gammas.ApplyCloudStencil(u_ap[:, k - 1], vec, diag, w)                # Operator applied to previous time state.
                un              = 2 * u_ap[:, k - 1] + Ku - u_ap[:, k - 2]                              # Central difference in time.
                u_ap[inne_n, k] = un[inne_n]                                                            # Save interior update.
                ref_scale       = float(np.max(np.abs(u_ap[boun_n, k]))) if np.any(boun_n) else float(np.max(np.abs(u_ap[:, k - 1])))
                                                                                                        # Reference scale at time k for stability check.
                if unstable_state(u_ap[inne_n, k], ref_scale):                                          # Detect blow-up / NaNs.
                    fallback_to_implicit = True                                                         # Switch to implicit solver.
                    break                                                                               # Stop explicit loop.

    if fallback_to_implicit:                                                                            # Implicit solver path.
        try:                                                                                            # Import SciPy lazily only when implicit is taken.
            if not HAS_SCIPY:
                raise ImportError("SciPy is required for implicit schemes")                             # SciPy sparse matrices and solvers. for implicit wave step.
        except ImportError as e:                                                                        # If SciPy is missing, implicit stepping is not available.
            raise ImportError('TimeDerivative2: implicit scheme requires SciPy (pip install scipy)') from e
                                                                                                        # Explain dependency to the user.

        ## Lambda handling for implicit family (clamped for robustness).
        if not bool(implicit):                                                                          # If we are here due to fallback, force fully implicit.
            lam_eff = 0.0                                                                               # Fallback mode forces fully implicit.
        else:                                                                                           # Otherwise use user-provided lambda (clamped above).
            lam_eff = float(lam)                                                                        # User-provided lambda.
            if not np.isfinite(lam_eff):                                                                # Guard against NaN/inf input.
                lam_eff = 0.0                                                                           # Guard against NaN/inf.
            if lam_eff < 0.0:                                                                           # Clamp to [0, 1].
                lam_eff = 0.0                                                                           # Clamp to [0, 1].
            if lam_eff > 1.0:                                                                           # Clamp to [0, 1].
                lam_eff = 1.0                                                                           # Clamp to [0, 1].
            if lam_eff >= 0.99:                                                                         # Avoid near-singular implicit formulation.
                lam_eff = 0.5                                                                           # Avoid near-singular formulation.

        inn_idx        = np.flatnonzero(inne_n)                                                         # Interior node global indices.
        ninn           = int(inn_idx.size)                                                              # Number of interior nodes.
        local          = -np.ones(m, dtype = np.int32)                                                  # Global -> local map (default -1).
        local[inn_idx] = np.arange(ninn, dtype = np.int32)                                              # Fill mapping.

        vec_inn  = vec[inn_idx]                                                                         # Neighbor indices for interior rows.
        w_inn    = w[inn_idx]                                                                           # Neighbor weights for interior rows.
        bcols    = -np.ones_like(vec_inn, dtype = np.int32)                                             # Boundary neighbor indices (per row).
        bweights = np.zeros_like(w_inn, dtype = float)                                                  # Boundary neighbor weights (per row).

        ## Matrix builder: assembles the interior-interior block and stores boundary couplings for RHS terms.
        def build_matrix(alpha):                                                                        # Assemble (I - alpha * L) on interior nodes.
            rows = []                                                                                   # Sparse COO rows.
            cols = []                                                                                   # Sparse COO cols.
            data = []                                                                                   # Sparse COO data.
            for r, i in enumerate(inn_idx):                                                             # Loop over interior rows for assembly.
                rows.append(r)                                                                          # Row index in interior system.
                cols.append(r)                                                                          # Column index for diagonal entry.
                data.append(1.0 - float(alpha) * float(diag[i]))                                        # Diagonal: I - alpha * L.
                for kk in range(nvec):                                                                  # Loop over neighbor slots for this row.
                    j = int(vec_inn[r, kk])                                                             # Global index of neighbor.
                    if j < 0:                                                                           # Skip missing neighbors (padding).
                        continue                                                                        # Skip missing neighbors.
                    wij = float(w_inn[r, kk])                                                           # Weight for neighbor contribution.
                    if inne_n[j]:                                                                       # Interior neighbor contributes to matrix.
                        rows.append(r)                                                                  # Interior-interior coupling.
                        cols.append(int(local[j]))                                                      # Map global neighbor to local column.
                        data.append(-float(alpha) * wij)                                                # Off-diagonal: -alpha * w_ij.
                    else:                                                                               # Boundary neighbor contributes to RHS.
                        bcols[r, kk] = j                                                                # Store boundary neighbor index.
                        bweights[r, kk] = wij                                                           # Store boundary weight.
            return coo_matrix((data, (rows, cols)), shape = (ninn, ninn)).tocsr()                       # Return CSR matrix.

        alpha1  = 0.5 * (1.0 - lam_eff)                                                                 # Alpha for k=1 system.
        alpha2  = (1.0 - lam_eff)                                                                       # Alpha for k>=2 system.
        A1      = build_matrix(alpha1)                                                                  # System matrix for the k=1 solve.
        A2      = build_matrix(alpha2)                                                                  # System matrix for the k>=2 solves.

        lu1     = None                                                                                  # LU factorization handler for A1.
        lu2     = None                                                                                  # LU factorization handler for A2.
        lu1_reg = None                                                                                  # Regularized LU for A1 (built on demand).
        lu2_reg = None                                                                                  # Regularized LU for A2 (built on demand).

        A1_csc  = A1.tocsc()                                                                            # Convert A1 to CSC for LU.
        A2_csc  = A2.tocsc()                                                                            # Convert A2 to CSC for LU.
        try:                                                                                            # LU factorization for the k=1 system.
            lu1 = splu(A1_csc)                                                                          # LU factorization for k=1 system.
        except Exception:                                                                               # LU may fail for singular/ill-conditioned matrices.
            lu1 = None                                                                                  # Keep lu1=None to use iterative/direct fallback.
        try:                                                                                            # LU factorization for the k>=2 system.
            lu2 = splu(A2_csc)                                                                          # LU factorization for k>=2 system.
        except Exception:                                                                               # LU may fail for singular/ill-conditioned matrices.
            lu2 = None                                                                                  # Keep lu2=None to use iterative/direct fallback.

        ## Boundary contribution for implicit systems: bc = sum(w_ij * u_b(j)).
        def bc_term(u_b_full):                                                                          # Compute RHS contribution from boundary neighbors.
            bc = np.zeros(ninn)                                                                         # Initialize RHS boundary contribution.
            for kk in range(nvec):                                                                      # Accumulate over neighbor slots.
                j = bcols[:, kk]                                                                        # Boundary neighbor indices for slot kk.
                mask = j >= 0                                                                           # Rows where this slot is a boundary neighbor.
                if np.any(mask):                                                                        # Only add when there are boundary couplings.
                    bc[mask] += bweights[mask, kk] * u_b_full[j[mask]]                                  # Add boundary contribution for active rows.
            return bc                                                                                   # Return bc vector for interior equations.

        ## Linear solver wrapper: LU first (fast), then iterative, then direct sparse solve as fallback.
        def solve_with_lu(A, A_csc, lu, lu_reg, rhs):                                                   # Solve A x = rhs with LU/iterative fallbacks.
            if lu is not None:                                                                          # Primary path: use LU factorization.
                x = lu.solve(rhs)                                                                       # Solve A x = rhs.
                if not np.all(np.isfinite(x)):                                                          # Detect invalid numbers.
                    x = spsolve(A, rhs)                                                                 # Fallback to direct sparse solve.
                else:                                                                                   # Optional sanity checks on magnitude.
                    scale_rhs = float(np.max(np.abs(rhs))) if rhs.size else 1.0                         # Scale reference for rhs.
                    maxabs_x  = float(np.max(np.abs(x))) if x.size else 0.0                             # Max magnitude of solution.
                    if (not np.isfinite(maxabs_x)) or (maxabs_x > 1e12 * max(1.0, scale_rhs)):          # Detect extreme amplification.
                        ## Regularize only when needed: helps with near-singular matrices in some clouds.
                        if lu_reg is None:                                                              # Build regularized LU only once if needed.
                            try:                                                                        # Attempt LU on A + eps*I.
                                lu_reg = splu(A_csc + (1e-8 * eye(A.shape[0], format = 'csc')))         # Regularized LU (A + eps I).
                            except Exception:                                                           # Regularized LU may fail as well.
                                lu_reg = None                                                           # Keep lu_reg=None if regularization fails.
                        if lu_reg is not None:                                                          # Re-solve if regularized LU succeeded.
                            x = lu_reg.solve(rhs)                                                       # Re-solve using the regularized LU.
            else:                                                                                       # Secondary path: use iterative + direct fallback.
                try:                                                                                    # Attempt BiCGStab (new SciPy API).
                    x, info = bicgstab(A, rhs, rtol = 1e-10, atol = 0.0, maxiter = 5000)                # Iterative solve (new SciPy API).
                except TypeError:                                                                       # Attempt BiCGStab (old SciPy API).
                    x, info = bicgstab(A, rhs, tol = 1e-10, maxiter = 5000)                             # Iterative solve (old SciPy API).
                if info != 0 or (not np.all(np.isfinite(x))):                                           # Fallback on non-convergence or invalid numbers.
                    x = spsolve(A, rhs)                                                                 # Direct sparse solve fallback.
            return x, lu_reg                                                                            # Return solution and updated regularized LU.

        ## First time level (k = 1).
        u_b_next_full         = np.zeros(m)                                                             # Boundary-only vector.
        u_b_next_full[boun_n] = u_ap[boun_n, 1]                                                         # Boundary values already filled.
        Ku0                   = Gammas.ApplyCloudStencil(u_ap[:, 0], vec, diag, w)                      # Operator applied to u at k=0.
        rhs1                  = u_ap[inn_idx, 0] + (0.5 * lam_eff) * Ku0[inn_idx] + dt * g(p[inn_idx, 0], p[inn_idx, 1], T[1], coef)
                                                                                                        # RHS for k=1 interior system.
        rhs1                  = rhs1 + float(alpha1) * bc_term(u_b_next_full)                           # Add boundary contributions.
        x1, lu1_reg           = solve_with_lu(A1, A1_csc, lu1, lu1_reg, rhs1)                           # Solve linear system.
        u_ap[inn_idx, 1]      = x1                                                                      # Save interior values.

        expand_neighbors = False                                                                        # Flag to request neighbor expansion.
        nvec_retry       = None                                                                         # Target nvec for retry.
        for k in np.arange(2, t):                                                                       # Time loop for implicit wave update (k >= 2).
            u_b_next_full         = np.zeros(m)                                                         # Boundary-only vector.
            u_b_next_full[boun_n] = u_ap[boun_n, k]                                                     # Boundary at time k.
            Ku_prev               = Gammas.ApplyCloudStencil(u_ap[:, k - 1], vec, diag, w)              # Operator applied to u at k-1.
            rhs                   = 2.0 * u_ap[inn_idx, k - 1] - u_ap[inn_idx, k - 2] + lam_eff * Ku_prev[inn_idx]
                                                                                                        # RHS for k>=2 interior system.
            rhs                   = rhs + float(alpha2) * bc_term(u_b_next_full)                        # Add boundary contributions.
            x, lu2_reg            = solve_with_lu(A2, A2_csc, lu2, lu2_reg, rhs)                        # Solve linear system.

            ref_scale = float(np.max(np.abs(u_ap[boun_n, k]))) if np.any(boun_n) else float(np.max(np.abs(u_ap[inn_idx, k - 1]))) # Reference scale.
            if unstable_state(x, ref_scale):                                                            # Detect instability of the implicit solve result.
                nvec_retry       = _next_nvec(nvec)                                                     # Request a larger stencil size.
                if nvec_retry is not None:
                    expand_neighbors = True                                                             # Retry only if there is a larger stencil available.
                    break                                                                               # Exit time loop to retry with a larger stencil.
                else:
                    raise FloatingPointError(f'TimeDerivative2 became unstable and no neighbor expansion was possible (max nvec={nvec} reached)')

            u_ap[inn_idx, k] = x                                                                        # Commit interior update at time k.

        if expand_neighbors:                                                                            # Re-run with expanded stencil if instability detected.
            if nvec_retry is None or nvec_retry <= nvec:                                                # Sanity check on requested expansion.
                raise FloatingPointError('TimeDerivative2 became unstable and no neighbor expansion was possible')
                                                                                                        # Hard failure when no expansion remains.
            return TimeDerivative2(p, f, g, t, coef, operator = operator, implicit = True, lam = 0.0, Adv = Adv, vec = None, nvec = nvec_retry)
                                                                                                        # Retry fully implicit with more neighbors.

    ## Theoretical solution.
    for k in np.arange(t):                                                                              # For all the time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], coef)                                                    # Evaluate exact solution at time T[k].

    return u_ap, u_ex, vec                                                                              # Return approximate/exact solutions and neighbor list.

