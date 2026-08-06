"""
TimeDerivative1 — First-order transient PDEs solver

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

def TimeDerivative1(p, f, t, coef, operator = np.vstack([[0], [0], [2], [0], [2]]), implicit = False, lam = 0.5, Adv = False, vec = None, nvec = 8):
    """
    TimeDerivative1
    Numerical solution of PDEs with a first-order time derivative using a Meshless Generalized Finite
    Difference scheme on a 2D cloud of points.

    The problem solved is:
        u_t = D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy

    Input:
        p           m x 3           ndarray         Point cloud with columns [x, y, flag].
        f                           function        Reference solution / boundary function f(x, y, t, coef).
                                                    Boundary and initial conditions are taken from this function.
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
                                                        Values close to 1 are avoided because they weaken the implicit part.
        Adv                         bool            If True, uses upwind-biased neighbor selection based on (D, E).
        vec                         ndarray|None    Optional neighbor list (m x nvec, padded with -1).
        nvec                        int             Maximum number of neighbors used by the stencil.

    Output:
        u_ap        m x t           ndarray         Array with the approximation computed by the routine.
        u_ex        m x t           ndarray         Array with the theoretical solution f(x, y, t, coef).
        vec         m x o           ndarray         Neighbor list used by the solver.

    Notes:
        When the explicit scheme becomes unstable (NaNs/infs or rapid growth), the routine switches to
        a fully implicit retry with an expanded neighbor stencil when available.
    """

    ## Variable initialization.
    m              = int(p.shape[0])                                                                    # Total number of nodes.
    nvec           = int(nvec)                                                                          # Requested maximum neighbors.
    T              = np.linspace(0, 1, t)                                                               # Time grid on [0, 1].
    dt             = T[1] - T[0]                                                                        # Time step size.
    u_ap           = np.zeros([m, t])                                                                   # Approximate solution in time.
    u_ex           = np.zeros([m, t])                                                                   # Theoretical solution in time.
    boun_n, inne_n = _node_masks(p)                                                                     # Boundary and interior masks.
    
    # Boundary conditions.
    for k in np.arange(t):                                                                              # For each time step.
        u_ap[boun_n, k] = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                                     # Impose boundary values at time T[k].
  
    # Initial condition
    u_ap[:, 0]          = f(p[:, 0], p[:, 1], T[0], coef)                                               # Initial condition at t = 0.
    
    ## Neighbor computation / validation.
    vec, nvec           = _prepare_neighbors(p, vec, nvec, Adv, operator)                               # Build/validate neighbor list.

    ## Stencil for spatial operator.
    L_op                = operator[:-1]                                                                 # Operator coefficients (spatial part).
    diag_op, w_op       = Gammas.CloudStencil(p, vec, L_op)                                             # Diagonal and neighbor weights.
    
    ## Stability heuristic for explicit stepping.
    row_sum             = np.abs(diag_op[inne_n]) + np.sum(np.abs(w_op[inne_n]), axis = 1)              # Row sum magnitude heuristic for stability/conditioning.
    max_row_sum         = float(np.max(row_sum)) if row_sum.size else 0.0                               # Worst-case row sum magnitude.
    
    if (nvec < 12) and (row_sum.size > 0) and (not np.isfinite(max_row_sum) or (max_row_sum > 5e5)):    # Expand stencil when row sums look ill-conditioned.
                                                                                                        # Expand stencil when row sums look ill-conditioned.
        nvec          = 12                                                                              # Expand stencil to reduce conditioning issues.
        vec, nvec     = _prepare_neighbors(p, None, nvec, Adv, operator)                                # Recompute neighbors for new nvec.
        diag_op, w_op = Gammas.CloudStencil(p, vec, L_op)                                               # Recompute stencil for expanded neighbors.
        row_sum       = np.abs(diag_op[inne_n]) + np.sum(np.abs(w_op[inne_n]), axis = 1)                # Recompute stability heuristic.
        max_row_sum   = float(np.max(row_sum)) if row_sum.size else 0.0                                 # Updated worst-case row sum magnitude.
    dt_stable         = (0.2 / max_row_sum) if (max_row_sum > 0.0 and np.isfinite(max_row_sum)) else dt # Heuristic stable dt for explicit scheme.
    nsub              = int(np.ceil(dt / dt_stable)) if dt_stable > 0.0 else 1                          # Number of explicit substeps.
    if nsub < 1:                                                                                        # Guard against invalid substep count.
        nsub = 1                                                                                        # Guard against invalid values.

    ## If explicit requires too many sub-steps, switch to implicit.
    stiff_fallback = (not bool(implicit)) and ((not np.isfinite(max_row_sum)) or (nsub > 250))          # Trigger implicit fallback when explicit is too stiff.
    lam_eff        = 0.0 if stiff_fallback else float(lam)                                              # Effective lambda for implicit scheme.
    use_implicit   = bool(implicit) or stiff_fallback                                                   # Choose implicit path.

    if use_implicit:                                                                                    # Implicit path (stable for stiff operators).
        try:                                                                                            # Import SciPy only when the implicit branch is taken.
            if not HAS_SCIPY:
                raise ImportError("SciPy is required for implicit schemes")                         # Sparse assembly + solvers for implicit step.
        except ImportError as e:                                                                        # If SciPy is missing, implicit stepping is not available.
            raise ImportError('TimeDerivative1: implicit scheme requires SciPy (pip install scipy)') from e
                                                                                                        # Explain dependency to the user.

        inn_idx        = np.flatnonzero(inne_n)                                                         # Interior node global indices.
        ninn           = int(inn_idx.size)                                                              # Number of interior nodes.
        local          = -np.ones(m, dtype = np.int32)                                                  # Global -> local map (default -1).
        local[inn_idx] = np.arange(ninn, dtype = np.int32)                                              # Fill mapping.

        ## Sparse matrix assembly: we build only the interior-interior block, and store boundary couplings separately.
        rows = []                                                                                       # Sparse COO rows.
        cols = []                                                                                       # Sparse COO cols.
        data = []                                                                                       # Sparse COO data.

        vec_inn  = vec[inn_idx]                                                                         # Neighbor indices for interior rows.
        w_inn    = w_op[inn_idx]                                                                        # Neighbor weights for interior rows.
        bcols    = -np.ones_like(vec_inn, dtype = np.int32)                                             # Boundary neighbor indices (per row).
        bweights = np.zeros_like(w_inn, dtype = float)                                                  # Boundary neighbor weights (per row).

        alpha    = dt * (1.0 - lam_eff)                                                                 # Implicit weight for operator contribution.
        for r, i in enumerate(inn_idx):                                                                 # Loop over interior rows for assembly.
            rows.append(r)                                                                              # Row index in interior system.
            cols.append(r)                                                                              # Column index for diagonal entry.
            data.append(1.0 - alpha * float(diag_op[i]))                                                # Diagonal: I - alpha * L.
            for kk in range(nvec):                                                                      # Loop over neighbor slots for this row.
                j = int(vec_inn[r, kk])                                                                 # Global index of neighbor.
                if j < 0:                                                                               # Skip missing neighbors (padding).
                    continue                                                                            # Skip missing neighbors.
                wij = float(w_inn[r, kk])                                                               # Weight for neighbor contribution.
                if inne_n[j]:                                                                           # Interior neighbor contributes to matrix.
                    rows.append(r)                                                                      # Interior-interior coupling.
                    cols.append(int(local[j]))                                                          # Map global neighbor to local column.
                    data.append(-alpha * wij)                                                           # Off-diagonal: -alpha * w_ij.
                else:                                                                                   # Boundary neighbor contributes to RHS.
                    bcols[r, kk] = j                                                                    # Store boundary neighbor index.
                    bweights[r, kk] = wij                                                               # Store boundary weight (applied to known u_b).

        A = coo_matrix((data, (rows, cols)), shape = (ninn, ninn)).tocsr()                              # System matrix for interior unknowns.
        A_csc  = A.tocsc()                                                                              # CSC format for sparse LU.
        lu     = None                                                                                   # LU factorization handler (built once if possible).
        lu_reg = None                                                                                   # Regularized LU handler (built on demand).
        try:                                                                                            # Try LU factorization once (reused each step).
            lu = splu(A_csc)                                                                            # LU factorization of the constant system matrix.
        except Exception:                                                                               # LU may fail for singular/ill-conditioned matrices.
            lu = None                                                                                   # Keep lu=None to use iterative/direct fallback.

        expand_neighbors = False                                                                        # Flag to request neighbor expansion.
        nvec_retry       = None                                                                         # Target nvec for retry.
        for k in np.arange(1, t):                                                                       # Time loop (implicit update).
            u_curr                = u_ap[:, k - 1].copy()                                               # Previous time state.
            u_b_next              = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                           # Boundary at next time.
            u_curr[boun_n]        = u_b_next                                                            # Insert boundary into working vector.
            u_b_next_full         = np.zeros(m)                                                         # Boundary-only vector for BC term.
            u_b_next_full[boun_n] = u_b_next                                                            # Fill boundary entries.

            if lam_eff != 0.0:                                                                          # Semi-implicit family uses lam term on RHS.
                Lu  = Gammas.ApplyCloudStencil(u_ap[:, k - 1], vec, diag_op, w_op)                      # L(u^n) on all nodes.
                rhs = u_ap[inn_idx, k - 1] + dt * lam_eff * Lu[inn_idx]                                 # Semi-implicit RHS: u^n + dt*lam*L(u^n).
            else:                                                                                       # Fully implicit fallback uses RHS = u^n.
                rhs = u_ap[inn_idx, k - 1].copy()                                                       # Fully implicit (lam=0): RHS is u^n.

            ## Boundary term: move known boundary neighbor contributions to the RHS.
            bc = np.zeros(ninn)                                                                         # Boundary contribution vector.
            for kk in range(nvec):                                                                      # Accumulate boundary contributions from stored couplings.
                j    = bcols[:, kk]                                                                     # Boundary neighbor indices for slot kk.
                mask = j >= 0                                                                           # Rows that have a boundary neighbor in this slot.
                if np.any(mask):                                                                        # Only process rows that have a boundary neighbor in this slot.
                    bc[mask] += bweights[mask, kk] * u_b_next_full[j[mask]]                             # Add contribution from known boundary values.
            rhs = rhs + alpha * bc                                                                      # Add alpha * sum(w_ij * u_b) to RHS.

            if lu is not None:                                                                          # Preferred fast path: LU solve.
                x         = lu.solve(rhs)                                                               # Fast solve using sparse LU factorization.
                scale_rhs = float(np.max(np.abs(rhs))) if rhs.size else 1.0                             # Reference scale for sanity checks.
                maxabs_x  = float(np.max(np.abs(x))) if x.size else 0.0                                 # Solution magnitude check.
                if (not np.all(np.isfinite(x))) or (maxabs_x > 1e12 * max(1.0, scale_rhs)):             # Detect invalid numbers or extreme amplification.
                    ## If LU gives NaNs/Infs or extreme growth, retry with a tiny diagonal regularization.
                    if lu_reg is None:                                                                  # Build regularized LU only once if needed.
                        try:                                                                            # Attempt LU on A + eps*I to regularize.
                            lu_reg = splu(A_csc + (1e-8 * eye(ninn, format = 'csc')))                   # Regularized LU (A + eps I).
                        except Exception:                                                               # Regularized LU may also fail.
                            lu_reg = None                                                               # Keep lu_reg=None when regularization fails.
                    if lu_reg is not None:                                                              # Re-solve using regularized LU.
                        x = lu_reg.solve(rhs)                                                           # Re-solve using the regularized LU.
                if not np.all(np.isfinite(x)):                                                          # If still invalid, use direct solve.
                    x = spsolve(A, rhs)                                                                 # As a last resort, use direct sparse solve.
            else:                                                                                       # If LU is unavailable, try iterative solve.
                try:                                                                                    # Attempt BiCGStab.
                    x, info = bicgstab(A, rhs, rtol = 1e-10, atol = 0.0, maxiter = 5000)                # Iterative solve (new SciPy API).
                except TypeError:                                                                       # Old SciPy API uses tol instead of rtol/atol.
                    x, info = bicgstab(A, rhs, tol = 1e-10, maxiter = 5000)                             # Iterative solve (old SciPy API).
                if info != 0 or (not np.all(np.isfinite(x))):                                           # Fallback when iterative solve fails.
                    x = spsolve(A, rhs)                                                                 # Fallback if iterative solve fails.

            ## Growth/NaN detection: if solution blows up relative to boundary, retry with more neighbors.
            u_scale = float(np.max(np.abs(u_b_next))) if u_b_next.size else 1.0                         # Boundary scale used as a reference.
            x_max   = float(np.max(np.abs(x))) if x.size else 0.0                                       # Interior scale used to detect blow-up.
            if (not np.all(np.isfinite(x))) or (not np.isfinite(x_max)) or (x_max > 1e6 * max(1.0, u_scale)):
                                                                                                        # Detect runaway growth in the implicit solution.
                                                                                                        # Detect runaway growth.
                nvec_retry = _next_nvec(nvec)                                                           # Request a larger stencil size.
                if nvec_retry is not None:                                                              # Retry only if there is a larger stencil available.
                    expand_neighbors = True                                                             # Trigger a full retry with the larger stencil.
                    break                                                                               # Exit time loop to rebuild with more neighbors.
            u_ap[inn_idx, k] = x                                                                        # Save interior update at time k.
            u_ap[boun_n, k]  = u_b_next                                                                 # Save boundary update at time k.
        if expand_neighbors:                                                                            # Re-run with expanded stencil if instability detected.
            if nvec_retry is None or nvec_retry <= nvec:                                                # Sanity check on requested expansion.
                raise FloatingPointError('TimeDerivative1 became unstable and no neighbor expansion was possible')
                                                                                                        # No larger stencil available.
                                                                                                        # Hard failure when no expansion remains.
            return TimeDerivative1(p, f, t, coef, operator = operator, implicit = True, lam = 0.0, Adv = Adv, vec = None, nvec = nvec_retry) # Retry fully implicit with more neighbors.
    else:                                                                                               # Explicit path (with optional sub-stepping).
        dt_sub = dt / nsub                                                                              # Explicit sub-step size.
        for k in np.arange(1, t):                                                                       # Time loop (explicit update).
            u_curr = u_ap[:, k - 1].copy()                                                              # Start from previous time level.
            for s in range(nsub):                                                                       # Sub-steps to satisfy stability estimate.
                Ku              = Gammas.ApplyCloudStencil(u_curr, vec, diag_op, w_op)                  # L(u) at current substep.
                u_curr[inne_n]  = u_curr[inne_n] + dt_sub * Ku[inne_n]                                  # Explicit Euler update on interior nodes.
                t_new           = T[k - 1] + (s + 1) * dt_sub                                           # Substep time for boundary evaluation.
                u_curr[boun_n]  = f(p[boun_n, 0], p[boun_n, 1], t_new, coef)                            # Re-impose boundary at substep time.
                if not np.all(np.isfinite(u_curr[inne_n])):                                             # Detect NaNs/Infs during explicit stepping.
                    use_implicit = True                                                                 # Detect numerical blow-up and switch.
                    break                                                                               # Abort sub-stepping.
            if use_implicit:                                                                            # If explicit blew up, leave time loop to retry implicit.
                break                                                                                   # Break time loop.
            u_ap[:, k] = u_curr                                                                         # Commit explicit update at time k.

        if use_implicit:                                                                                # Implicit fallback for explicit instability.
            nvec_next = 12 if nvec < 12 else nvec                                                       # Minimum nvec used for implicit retry.
            return TimeDerivative1(p, f, t, coef, operator = operator, implicit = True, lam = 0.0, Adv = Adv, vec = None, nvec = nvec_next)
                                                                                                        # Retry fully implicit.
        
    ## Theoretical solution.
    for k in np.arange(t):                                                                              # For all the time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], coef)                                                    # Evaluate exact solution at time T[k].

    return u_ap, u_ex, vec                                                                              # Return approximate/exact solutions and neighbor list.

