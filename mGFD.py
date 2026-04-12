"""
mGFD — Meshless Generalized Finite Differences

Overview:
    Meshless solvers for 2D stationary and transient PDEs using generalized finite differences (GFD).
    Spatial derivatives are approximated using local reconstructions over a node neighborhood (nvec).

Public API:
    Stationary          Stationary problems with Dirichlet boundary conditions.
    TimeDerivative1     First-order-in-time problems (heat / advection-diffusion family).
    TimeDerivative2     Second-order-in-time problems (wave family).

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior nodes; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.
            If vec is not provided, it is computed from p using Neighbors.Cloud / Neighbors.CloudAdv.

Operator conventions:
    operator                array-like
            A 6-coefficient vector [D, E, A, B, C, F] (shape (6,) or (6, 1)).
            The spatial stencil is built with L = operator[:5] = [D, E, A, B, C], interpreted as:
                D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy
            The reaction term F*u is reserved in the last coefficient, but it is not applied by the
            current implementation. For the Laplacian, use [0, 0, 2, 0, 2, 0].
            When Adv=True, neighbor selection is upwind-biased using velocities (D, E).

Notes:
    NumPy is required. SciPy is optional, but implicit schemes require SciPy for sparse linear algebra.
    Transient solvers use a normalized time grid T = linspace(0, 1, t).
    When instability is detected, the solver may retry with expanded neighborhoods (8→12→16→20→30).

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx
    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
        Aula CIMNE-Morelia. México
        SIIIA-MATH: Soluciones de Ingeniería. México
    Date:
        May, 2024.
    Last Modification:
        April, 2026.
"""

## Library importation.
import numpy as np                                                                                      # Core numerical operations.
import Scripts.Gammas as Gammas                                                                         # Stencil (gamma) computation and application.
import Scripts.Neighbors as Neighbors                                                                   # Neighbor search utilities.

def _import_scipy_stationary():
    '''
    _import_scipy_stationary
    Lazy import of the SciPy sparse utilities used by Stationary().
    
    Notes:
        This keeps SciPy as an optional dependency for the module. If SciPy is not installed,
        this helper raises ImportError and the caller may fall back to non-SciPy paths.
    
    Output:
        LinearOperator            SciPy LinearOperator class.
        bicgstab                  SciPy BiCGStab iterative solver.
        spsolve                   SciPy direct sparse solver.
        coo_matrix                SciPy COO sparse matrix constructor.
    '''
    try:                                                                                                # Import SciPy lazily to keep SciPy optional for the module.
        from scipy.sparse.linalg import LinearOperator, bicgstab, spsolve                               # Stationary: matrix-free operator + iterative/direct solvers.
        from scipy.sparse import coo_matrix                                                             # Stationary fallback: sparse matrix assembly.
    except ImportError as e:                                                                            # SciPy is not installed in the environment.
        raise ImportError('SciPy is required for this operation') from e                                # Provide a clear dependency error.
    return LinearOperator, bicgstab, spsolve, coo_matrix                                                # Return imported symbols to the caller.

def _import_scipy_implicit():
    '''
    _import_scipy_implicit
    Lazy import of the SciPy sparse utilities used by implicit schemes in TimeDerivative1/2.
    
    Notes:
        This keeps explicit schemes usable without SciPy. If SciPy is not installed, this helper
        raises ImportError and the implicit branches should surface a clear dependency message.
    
    Output:
        coo_matrix                SciPy COO sparse matrix constructor.
        eye                       SciPy sparse identity constructor.
        bicgstab                  SciPy BiCGStab iterative solver.
        spsolve                   SciPy direct sparse solver.
        splu                      SciPy sparse LU factorization.
    '''
    try:                                                                                                # Import SciPy lazily to keep explicit schemes usable without SciPy.
        from scipy.sparse import coo_matrix, eye                                                        # Sparse matrix utilities (assembly + identity).
        from scipy.sparse.linalg import bicgstab, spsolve, splu                                         # Sparse solvers (iterative/direct/LU).
    except ImportError as e:                                                                            # SciPy is not installed in the environment.
        raise ImportError('SciPy is required for implicit schemes') from e                              # Provide a clear dependency error.
    return coo_matrix, eye, bicgstab, spsolve, splu                                                     # Return imported symbols to the caller.

def _node_masks(p):
    '''
    _node_masks
    Build boolean masks for boundary and interior nodes.
    
    Input:
        p           m x 3           ndarray         Array with node coordinates and a node flag in column 2.
                                                    flag = 0 for interior, flag = 1/2 for boundary.
    
    Output:
        boun_n      m               ndarray         Boolean mask: True on boundary nodes (flag 1 or 2).
        inne_n      m               ndarray         Boolean mask: True on interior nodes (flag 0).
    '''
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                            # Boundary nodes mask.
    inne_n = p[:, 2] == 0                                                                               # Interior nodes mask.
    return boun_n, inne_n                                                                               # Return both masks.

def _prepare_neighbors(p, vec, nvec, Adv, operator):
    '''
    _prepare_neighbors
    Validate/load/build the neighbor array used by the stencil.
    
    Input:
        p           m x 3           ndarray         Array with node coordinates and type flag.
        vec                         ndarray|None    Optional precomputed neighbors (m x nvec, padded with -1).
        nvec                        int             Requested maximum number of neighbors.
        Adv                         bool            If True, use upwind neighbor selection (CloudAdv).
        operator                    ndarray         Operator coefficients (used to extract velocities when Adv=True).
                                                    It assumes operator[0] and operator[1] correspond to D and E.
    
    Output:
        vec         m x nvec        ndarray         Neighbor indices per node (int32).
        nvec                        int             Effective nvec (taken from vec if vec provided).
    '''
    m    = int(p.shape[0])                                                                              # Number of nodes.
    nvec = int(nvec)                                                                                    # Ensure integer nvec.

    if vec is None:                                                                                     # Compute neighbors when not provided.
        if Adv:                                                                                         # Upwind neighbor selection (advection problems).
            a   = float(operator[0][0])                                                                 # Velocity in x (for upwind).
            b   = float(operator[1][0])                                                                 # Velocity in y (for upwind).
            vec = Neighbors.CloudAdv(p, a, b, nvec)                                                     # Compute neighbors with upwind preference.
        else:                                                                                           # Standard neighbor selection (no upwind bias).
            vec = Neighbors.Cloud(p, nvec)                                                              # Compute standard neighbors.
    else:                                                                                               # Validate and normalize a provided vec.
        vec = np.asarray(vec)                                                                           # Normalize input to ndarray.
        if vec.ndim == 1:                                                                               # Allow passing a single row as 1D.
            vec = vec.reshape(1, -1)                                                                    # Support single-row vec passed as 1D.
        if vec.shape[0] != m:                                                                           # Ensure one neighbor row per node.
            raise ValueError('vec has incorrect number of rows')                                        # Shape mismatch.
        nvec = int(vec.shape[1])                                                                        # Trust nvec from provided vec.

    vec = vec.astype(np.int32, copy = False)                                                            # Ensure fast indexing dtype.
    return vec, nvec                                                                                    # Return neighbors and effective nvec.

def _next_nvec(nvec):
    '''
    _next_nvec
    Choose the next neighbor count when expanding stencils.
    
    Input:
        nvec                        int             Current nvec.
    
    Output:
        nvec_next                   int|None        Next nvec in {12,16,20,30} or None if no expansion left.
    '''
    nvec = int(nvec)                                                                                    # Ensure integer nvec.
    if nvec < 12:                                                                                       # Promote to the first expansion size.
        return 12                                                                                       # First expansion.
    if nvec < 16:                                                                                       # Promote to the second expansion size.
        return 16                                                                                       # Second expansion.
    if nvec < 20:                                                                                       # Promote to the third expansion size.
        return 20                                                                                       # Third expansion.
    if nvec < 30:                                                                                       # Promote to the largest supported expansion size.
        return 30                                                                                       # Max expansion used in batches.
    return None                                                                                         # No further expansion possible.

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
        LinearOperator, bicgstab, spsolve, coo_matrix = _import_scipy_stationary()                      # Load SciPy operators/solvers lazily.
        Aop = LinearOperator((inn_idx.size, inn_idx.size), matvec=matvec_int, dtype=float)              # Matrix-free linear operator on interior nodes.
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
            for row_k, gi in enumerate(inn_idx):                                                   # For each interior equation row.
                rows.append(row_k)                                                                      # Diagonal entry row.
                cols.append(row_k)                                                                      # Diagonal entry col.
                data.append(diag[gi])                                                                   # Diagonal weight.
                for t, gj in enumerate(vec[gi]):                                                   # Iterate neighbors of this node.
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
            coo_matrix, eye, bicgstab, spsolve, splu = _import_scipy_implicit()                         # Sparse assembly + solvers for implicit step.
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
        for r, i in enumerate(inn_idx):                                                            # Loop over interior rows for assembly.
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
        try:                                                                                            # Import SciPy only when the implicit branch is taken.
            coo_matrix, eye, bicgstab, spsolve, splu = _import_scipy_implicit()                         # Sparse assembly + solvers for implicit wave step.
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
            for r, i in enumerate(inn_idx):                                                        # Loop over interior rows for assembly.
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
                expand_neighbors = nvec_retry is not None                                               # Retry only if there is a larger stencil available.
                break                                                                                   # Exit time loop to retry with a larger stencil.

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
