"""
TimeDerivative1 — First-order transient PDEs solver

Overview:
    Numerical solver for PDEs with a first-order time derivative using a Meshless Generalized Finite
    Difference (mGFD) scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    TimeDerivative1     Main solver function for first-order transient problems.

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
import time                                                                                                                             # Timing for SolverResult.

from scipy.sparse import diags, eye                                                                                                     # Sparse matrix operations.
from scipy.sparse.linalg import factorized, bicgstab, gmres, spsolve                                                                    # Direct and iterative sparse linear solvers.
from typing import Callable, Optional, Tuple, List, Union, Any                                                                          # Type hinting.

from mGFD.exceptions import CloudShapeError, InputTypeError, DimensionMismatchError, OperatorFormatError, ParameterError                # Custom exceptions.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.core.adapters import extract_cloud, repack_solution                                                                           # Pandas/Xarray adapters.
import mGFD.core.gammas as Gammas                                                                                                       # Gammas calculation and sparse matrix builder.
import mGFD.core.neighbors as Neighbors                                                                                                 # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def TimeDerivative1(p: Union[np.ndarray, Any], 
                    f: Union[Callable, np.ndarray, float, int, Any],
                    t: int,
                    coef: List[float],
                    operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2], [0]]),
                    upwind: bool = False,
                    vec: Optional[np.ndarray] = None,
                    nvec: int = 20,
                    implicit: bool = False,
                    lam: float = 0.5,
                    linear_solver: str = "spsolve",
                    preconditioner: Optional[str] = None,
                    matrix_free: bool = False,
                    device: str = "cpu",
                    verbose: bool = True) -> SolverResult:
    """
    Numerical solution of partial differential equations with first-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \\frac{\\partial u}{\\partial t} = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        f                           Union           Boundary/Initial condition f(x, y, t) or data array.
        t                           int             Number of time steps to compute.
        coef                        List            Physical coefficients for the problem formulation.
        operator    6 x 1           ndarray         Array with the weights for the operator.
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
        implicit                    bool            If True, uses an implicit scheme (Theta Method).
        lam                         float           Theta parameter for implicit method (0.5 = Crank-Nicolson, 1.0 = Backward Euler).
        linear_solver               str             Algebraic backend: 'spsolve', 'bicgstab', or 'gmres'.
        preconditioner              Optional[str]   Preconditioner method: 'ilu', 'jacobi', or None.
        matrix_free                 bool            If True, uses on-the-fly matrix-vector multiplication (requires iterative solver).
        device                      str             'cpu' or 'cuda'. If 'cuda', offloads solvers to the GPU (requires CuPy).
        verbose                     bool            If True, prints solver progress.
    
    Output:
        SolverResult                Returns a structured result containing the approximation and the neighbor list.
    """
    
    start_time = time.perf_counter()                                                                                                    # Start solver timer.
    
    # 0. Input validation
    p_orig = p                                                                                                                          # Preserve original format.
    p = extract_cloud(p)                                                                                                                # Extract NumPy array if Pandas/Xarray.
    f = extract_cloud(f) if not callable(f) and not isinstance(f, (float, int)) else f                                                  # Extract array from Pandas/Xarray for f.
    
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                  # Raise explicit error on bad input.
    if not (callable(f) or isinstance(f, (np.ndarray, float, int))):                                                                    # Validate RHS type.
        raise InputTypeError("Forcing function 'f' must be a callable, ndarray, or numeric constant.")                                  # Raise explicit error on bad input.
    if not isinstance(t, int) or t <= 0:                                                                                                # Validate time steps.
        raise ParameterError("Number of time steps 't' must be a positive integer.")                                                    # Raise explicit error on bad input.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients.")                                       # Raise explicit error on bad input.
    if device not in ["cpu", "cuda"]:                                                                                                   # Validate device choice.
        raise ParameterError("Unsupported device. Choose 'cpu' or 'cuda'.")                                                             # Raise explicit error.
    if device == "cuda" and matrix_free:                                                                                                # Validate matrix_free and cuda.
        raise ParameterError("matrix_free=True is incompatible with device='cuda'.")                                                    # CuPy requires explicit sparse matrices.

    # 1. Variable initialization
    m      = p.shape[0]                                                                                                                 # Total number of nodes.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f"Solving Transient problem ({t} steps) for {m} nodes...")                                                          # Print solver progress.
        
    T      = np.linspace(0, 1, t)                                                                                                       # Time discretization array.
    dt     = T[1] - T[0]                                                                                                                # Time step size.
    u_ap   = np.zeros([m, t])                                                                                                           # Numerical approximation matrix.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Boolean mask for boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Boolean mask for interior nodes.

    # 2. Extract advection velocities for Upwind scheme.
    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]                                                                     # X-velocity (D coefficient).
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]                                                                     # Y-velocity (E coefficient).
    
    # 3. Apply Boundary and Initial Conditions
    if callable(f):                                                                                                                     # If data is a function.
        for k in range(t):                                                                                                              # Loop through all time steps.
            u_ap[boun_n, k] = np.asarray(f(p[boun_n, 0], p[boun_n, 1], T[k], coef))                                                     # Boundary condition (Dirichlet).
        u_ap[:, 0] = np.asarray(f(p[:, 0], p[:, 1], T[0], coef))                                                                        # Initial condition across all nodes.
    elif isinstance(f, np.ndarray):                                                                                                     # If data is an array.
        if f.ndim == 2 and f.shape == (m, t):                                                                                           # Spatiotemporal data array.
            u_ap[boun_n, :] = f[boun_n, :]                                                                                              # Spatiotemporal boundary conditions.
            u_ap[:, 0] = f[:, 0]                                                                                                        # Initial conditions.
        elif f.ndim == 1 and f.shape[0] == m:                                                                                           # Spatial constant data array.
            for k in range(t):                                                                                                          # Loop over time.
                u_ap[boun_n, k] = f[boun_n]                                                                                             # Constant boundary conditions.
            u_ap[:, 0] = f                                                                                                              # Initial conditions.
        if isinstance(f, np.ndarray) and f.shape not in [(m, t), (m,)]:                                                                 # Check if numeric matrix shape matches dimensions.
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")                                       # Raise explicit error.
    elif isinstance(f, (int, float)):                                                                                                   # If data is a constant scalar.
        u_ap[boun_n, :] = f                                                                                                             # Constant boundary condition.
        u_ap[:, 0] = f                                                                                                                  # Constant initial condition.
    
    # 4. Neighbor search
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind:                                                                                                                      # If an Upwind stencil is requested.
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                                     # Upwind-biased neighbor selection.
        else:                                                                                                                           # All the other cases.
            vec = Neighbors.compute_neighbors(p, nvec)                                                                                  # Standard distance-based neighbors.

    # 5. Compute differentiation matrix (K)
    L         = operator[:-1]                                                                                                           # Original operator weights.
    if matrix_free:
        from scipy.sparse.linalg import LinearOperator                                                                                  # Import LinearOperator.
        K_matvec = Gammas.compute_K_matvec(p, vec, L)                                                                                   # Generate the on-the-fly matrix-vector closure.
    else:
        K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                             # Build sparse spatial differentiation matrix.
        K         = dt * K_spatial                                                                                                      # Scale by time step.
    
    # 6. Time Integration (Generalized Finite Differences)
    converged = True                                                                                                                    # Assume convergence by default.
    if matrix_free and preconditioner is not None:
        raise ParameterError("Preconditioners are not currently supported in matrix_free=True mode.")                                   # Preconditioners need explicit matrix elements.
    
    if not implicit:
        # Explicit scheme (Forward Euler)
        if matrix_free:
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                un              = u_ap[:, k-1] + dt * K_matvec(u_ap[:, k-1])                                                            # Explicit matrix-free multiplication.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
        else:
            K2 = eye(m) + K                                                                                                             # LHS Explicit Matrix.
            if device == "cuda":                                                                                                        # If GPU acceleration is requested.
                try:
                    import cupy as cp                                                                                                   # Import CuPy dynamically.
                    from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix                                                          # Import CuPy CSR matrix.
                except ImportError:                                                                                                     # If not installed.
                    raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]'.")                         # Friendly error message.
                K2 = cp_csr_matrix(K2)                                                                                                  # Transfer explicit matrix to GPU.
                u_ap = cp.asarray(u_ap)                                                                                                 # Transfer solution matrix to GPU.
                inne_n = cp.asarray(inne_n)                                                                                             # Transfer inner nodes to GPU.
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                un              = K2.dot(u_ap[:, k-1])                                                                                  # Explicit matrix-vector multiplication.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
    else:
        # Implicit scheme (Theta Method)
        if matrix_free:
            if linear_solver == "spsolve":
                raise ParameterError("Direct solver 'spsolve' is incompatible with matrix_free=True.")                                  # Direct solvers need explicit matrix.
                
            def A_matvec(x):                                                                                                            # Define the matvec closure for LHS matrix.
                y = x - lam * dt * K_matvec(x)                                                                                          # Theta parameter applied to inner nodes.
                y[boun_n] = x[boun_n]                                                                                                   # Identity applied to boundary nodes.
                return y                                                                                                                # Return the result.
            
            A = LinearOperator(shape=(m, m), matvec=A_matvec, dtype=np.float64)                                                         # type: ignore
            M = None                                                                                                                    # No preconditioner in matrix-free.
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                u_prev = u_ap[:, k-1]                                                                                                   # Previous step solution.
                RHS = u_prev + (1 - lam) * dt * K_matvec(u_prev)                                                                        # Matrix-free RHS calculation.
                RHS[boun_n] = u_ap[boun_n, k]                                                                                           # Inject exact boundary conditions.
                if linear_solver == "bicgstab":                                                                                         # Iterative solver (BiCGStab).
                    un, info = bicgstab(A, RHS, x0=u_prev, M=M)                                                                         # Solve with previous step as initial guess.
                elif linear_solver == "gmres":                                                                                          # Iterative solver (GMRES).
                    un, info = gmres(A, RHS, x0=u_prev, M=M)                                                                            # Solve with previous step as initial guess.
                if info != 0:                                                                                                           # If not strictly converged.
                    converged = False                                                                                                   # Mark convergence failure.
                    if verbose: logger.warning(f"{linear_solver.upper()} did not converge perfectly (code {info}) at time {k}.")        # Warn on convergence issues.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.
        else:
            # Boundary rows must be isolated so they don't corrupt the implicit solve for internal nodes.
            Id_inner = diags(inne_n.astype(float))                                                                                      # Diagonal mask for inner nodes.
            Id_bound = diags(boun_n.astype(float))                                                                                      # Diagonal mask for boundary nodes.
            
            A        = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                         # LHS Matrix: Theta parameter applied to inner, Identity to boundary.
            A        = A.tocsc()                                                                                                        # Convert to CSC format for efficient SuperLU factorization.
            B        = Id_inner @ (eye(m) + (1 - lam) * K)                                                                              # RHS Matrix: Zeros for boundaries, explicit part for inner.
            
            if device == "cuda":                                                                                                        # If GPU acceleration is requested.
                try:
                    import cupy as cp                                                                                                   # Import CuPy dynamically.
                    from cupyx.scipy.sparse import csc_matrix as cp_csc_matrix, csr_matrix as cp_csr_matrix                             # Import CuPy sparse matrices.
                    from cupyx.scipy.sparse.linalg import spsolve as cp_spsolve, bicgstab as cp_bicgstab, gmres as cp_gmres             # Import CuPy solvers.
                except ImportError:                                                                                                     # If not installed.
                    raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]'.")                         # Friendly error message.
                A = cp_csc_matrix(A)                                                                                                    # Transfer LHS to GPU.
                B = cp_csr_matrix(B)                                                                                                    # Transfer RHS operator to GPU.
                u_ap = cp.asarray(u_ap)                                                                                                 # Transfer solution matrix to GPU.
                boun_n = cp.asarray(boun_n)                                                                                             # Transfer boundary nodes to GPU.
                inne_n = cp.asarray(inne_n)                                                                                             # Transfer inner nodes to GPU.
            
            M        = Gammas.compute_preconditioner(A, preconditioner) if preconditioner else None                                     # type: ignore
            
            if linear_solver == "spsolve":                                                                                              # Direct pre-factorized solver.
                if device == "cuda":                                                                                                    # GPU Factorization.
                    from cupyx.scipy.sparse.linalg import factorized as cp_factorized                                                   # Import CuPy factorized.
                    solve = cp_factorized(A)                                                                                            # Pre-factorize LHS on GPU.
                else:                                                                                                                   # CPU Factorization.
                    solve = factorized(A)                                                                                               # Pre-factorize LHS for fast repeated solves.
            elif linear_solver not in ["bicgstab", "gmres"]:                                                                            # Invalid iterative solver choice.
                raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")       # Raise explicit error.
                
            for k in range(1, t):                                                                                                       # Loop over all time steps.
                RHS             = B.dot(u_ap[:, k-1])                                                                                   # Right-hand side from previous step.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary conditions.
                
                if linear_solver == "spsolve":                                                                                          # Direct pre-factorized solver.
                    un       = solve(RHS)                                                                                               # Solve global system for time level k.
                elif linear_solver == "bicgstab":                                                                                       # Iterative solver (BiCGStab).
                    if device == "cuda":                                                                                                # GPU Iterative Solver.
                        un, info = cp_bicgstab(A, RHS, x0=u_ap[:, k-1], M=M)                                                            # Solve on GPU.
                    else:                                                                                                               # CPU Iterative Solver.
                        un, info = bicgstab(A, RHS, x0=u_ap[:, k-1], M=M)                                                               # Solve with previous step as initial guess.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")                   # Warn on convergence issues.
                elif linear_solver == "gmres":                                                                                          # Iterative solver (GMRES).
                    if device == "cuda":                                                                                                # GPU Iterative Solver.
                        un, info = cp_gmres(A, RHS, x0=u_ap[:, k-1], M=M)                                                               # Solve on GPU.
                    else:                                                                                                               # CPU Iterative Solver.
                        un, info = gmres(A, RHS, x0=u_ap[:, k-1], M=M)                                                                  # Solve with previous step as initial guess.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")                      # Warn on convergence issues.
                    
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes for time level k.     
    if device == "cuda":                                                                                                                # If array is in VRAM.
        u_ap = getattr(u_ap, "get")()                                                                                                   # Retrieve array to CPU RAM.
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info("\tSolver finished successfully.")                                                                                  # Print completion message.
        
    compute_time = time.perf_counter() - start_time                                                                                     # Calculate total compute time.
    u_ap_packed = repack_solution(p_orig, u_ap)                                                                                         # Repack into original Pandas/Xarray format.
    return SolverResult(solution=u_ap_packed, neighbors=vec, converged=converged, compute_time=compute_time)                            # Return structured result.
