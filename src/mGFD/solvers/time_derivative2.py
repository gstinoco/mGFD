"""
TimeDerivative2 — Second-order transient PDEs solver

Overview:
    Numerical solver for PDEs with a second-order time derivative using a Meshless Generalized Finite
    Difference (mGFD) scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    TimeDerivative2     Main solver function for second-order transient problems.

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

from scipy.sparse import eye, diags                                                                                                     # Sparse matrix generation.
from scipy.sparse.linalg import factorized, bicgstab, gmres, spsolve                                                                    # Direct and iterative sparse linear solvers.
from typing import Callable, Optional, Tuple, List, Union, Any                                                                          # Type hinting.

from mGFD.exceptions import CloudShapeError, InputTypeError, DimensionMismatchError, OperatorFormatError, ParameterError                # Custom exceptions.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.core.adapters import extract_cloud, repack_solution                                                                           # Pandas/Xarray adapters.
import mGFD.core.gammas as Gammas                                                                                                       # Gammas calculation and sparse matrix builder.
import mGFD.core.neighbors as Neighbors                                                                                                 # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def TimeDerivative2(p: Union[np.ndarray, Any], 
                    f: Union[Callable, np.ndarray, float, int, Any], 
                    g: Union[Callable, np.ndarray, float, int, Any], 
                    t: int, 
                    coef: List[float], 
                    operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2], [0]]), 
                    upwind: bool = False, 
                    vec: Optional[np.ndarray] = None, 
                    nvec: int = 20, 
                    implicit: bool = False, 
                    lam: float = 0.5, 
                    linear_solver: str = "spsolve",
                    verbose: bool = True) -> SolverResult:
    """
    Numerical solution of partial differential equations with second-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \\frac{\\partial^2 u}{\\partial t^2} = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        f                           Callable        Function declared with the initial condition for u.
        g                           Callable        Function declared with the initial condition for du/dt.
        t                           int             Number of time steps to compute.
        coef                        List            Physical coefficients for the problem formulation.
        operator    6 x 1           ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        implicit                    bool            If True, uses an implicit (theta) time integration scheme.
        lam                         float           Theta parameter for implicit scheme (0=explicit, 1=fully implicit, 0.5=Newmark).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
        verbose                     bool            If True, prints solver progress.
    
    Output:
        SolverResult                Returns a structured result containing the approximation and the neighbor list.
    """
    
    start_time = time.perf_counter()                                                                                                    # Start solver timer.
    
    # 0. Input validation
    p_orig = p                                                                                                                          # Preserve original format.
    p = extract_cloud(p)                                                                                                                # Extract NumPy array if Pandas/Xarray.
    f = extract_cloud(f) if not callable(f) and not isinstance(f, (float, int)) else f                                                  # Extract array from Pandas/Xarray for f.
    g = extract_cloud(g) if not callable(g) and not isinstance(g, (float, int)) else g                                                  # Extract array from Pandas/Xarray for g.
    
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                  # Raise explicit error on bad input.
    if not (callable(f) or isinstance(f, (np.ndarray, float, int))):                                                                    # Validate initial condition data type.
        raise InputTypeError("Initial condition function 'f' must be a callable, ndarray, or numeric constant.")                        # Raise explicit error on bad input.
    if not (callable(g) or isinstance(g, (np.ndarray, float, int))):                                                                    # Validate initial velocity data type.
        raise InputTypeError("Initial velocity function 'g' must be a callable, ndarray, or numeric constant.")                         # Raise explicit error on bad input.
    if not isinstance(t, int) or t <= 0:                                                                                                # Validate time steps.
        raise ParameterError("Number of time steps 't' must be a positive integer.")                                                    # Raise explicit error on bad input.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients.")                                       # Raise explicit error on bad input.

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
    # Evaluate f (initial position and boundaries)
    if callable(f):                                                                                                                     # If data is a function.
        for k in range(t):                                                                                                          # Loop through all time steps.
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
        else:                                                                                                                           # Invalid shape.
            raise DimensionMismatchError(f"Data array 'f' must have shape ({m}, {t}) or ({m},).")                                       # Raise explicit error.
    elif isinstance(f, (int, float)):                                                                                                   # If data is a constant scalar.
        u_ap[boun_n, :] = f                                                                                                             # Constant boundary condition.
        u_ap[:, 0]      = f                                                                                                             # Constant initial condition.

    # Evaluate g (initial velocity)
    if callable(g):                                                                                                                     # If data is a function.
        g_eval = np.asarray(g(p[:, 0], p[:, 1], T[0], coef))                                                                            # Evaluate initial velocity function.
    elif isinstance(g, np.ndarray):                                                                                                     # If data is an array.
        if g.ndim == 1 and g.shape[0] == m:                                                                                             # Check correct spatial dimensions.
            g_eval = g                                                                                                                  # Use provided initial velocity array.
        else:                                                                                                                           # Invalid shape.
            raise DimensionMismatchError(f"Data array 'g' must have shape ({m},).")                                                     # Raise explicit error.
    else:                                                                                                                               # If data is a constant scalar.
        g_eval = np.full(m, float(g))                                                                                                   # Use scalar value.
    
    # 4. Neighbor search
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind:                                                                                                                      # If an Upwind stencil is requested.
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                                     # Upwind-biased neighbor selection.
        else:                                                                                                                           # All the other cases.
            vec = Neighbors.compute_neighbors(p, nvec)                                                                                  # Standard distance-based neighbors.

    # 5. Compute differentiation matrix (K)
    L         = operator[:-1]                                                                                                           # Original operator weights.
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                                                                                 # Build sparse spatial differentiation matrix.
    K         = (dt**2) * K_spatial                                                                                                     # Scale by time step squared.
    
    # 6. Time Integration (Generalized Finite Differences)
    converged = True                                                                                                                    # Assume convergence by default.
    
    if not implicit:                                                                                                                    # If an explicit scheme is requested.
        # Explicit scheme
        K2 = eye(m) + (1/2)*K                                                                                                           # Explicit matrix for k = 1.
        K4 = 2*eye(m) + K                                                                                                               # Explicit matrix for k = 2, ..., t.
        
        for k in range(1, t):                                                                                                           # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step logic (k=1).
                g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                      # Evaluate velocity for first step.
                un              = K2.dot(u_ap[:, k - 1]) + dt*g_val                                                                     # New time-level (k=1) using du/dt.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
            else:                                                                                                                       # Following time steps (k>=2).
                un              = K4.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]                                                               # New time-level (k>=2) using Leapfrog/Verlet.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
    else:                                                                                                                               # If an implicit scheme is requested.
        # Implicit scheme
        Id_inner = diags(inne_n.astype(float))                                                                                          # Diagonal mask for inner nodes.
        Id_bound = diags(boun_n.astype(float))                                                                                          # Diagonal mask for boundary nodes.
        
        # Step k = 1 Formulation (with du/dt condition)
        A1       = Id_inner @ (eye(m) - lam * (1/2) * K) + Id_bound                                                                     # LHS Matrix: Theta parameter applied to inner, Identity to boundary (k=1).
        B1       = Id_inner @ (eye(m) + (1 - lam) * (1/2) * K)                                                                          # RHS Matrix: Zeros for boundaries, explicit part for inner (k=1).
        
        # Step k >= 2 Formulation
        A2       = Id_inner @ (eye(m) - lam * K) + Id_bound                                                                             # LHS Matrix: Theta parameter applied to inner, Identity to boundary (k>=2).
        A2       = A2.tocsc()                                                                                                           # Convert to CSC format for efficient SuperLU factorization.
        B2       = Id_inner @ (2*eye(m) + (1 - lam) * K)                                                                                # RHS Matrix: Zeros for boundaries, explicit part for inner (k>=2).
        
        if linear_solver == "spsolve":                                                                                                  # Direct pre-factorized solver.
            solve1   = factorized(A1.tocsc())                                                                                           # Pre-factorize implicit matrix for k = 1.
            solve2   = factorized(A2)                                                                                                   # Pre-factorize implicit matrix for k = 2, ..., t.
        elif linear_solver not in ["bicgstab", "gmres"]:                                                                                # Invalid iterative solver choice.
            raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")           # Raise explicit error.

        for k in range(1, t):                                                                                                       # Loop over all time steps.
            if k == 1:                                                                                                                  # Initial time step logic (k=1).
                g_val = g(p[:, 0], p[:, 1], T[k], coef) if callable(g) else g_eval                                                      # Evaluate velocity for first step.
                RHS             = B1.dot(u_ap[:, k - 1]) + dt*g_val                                                                     # Right-hand side from previous step + du/dt.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary condition for time k=1.
                
                if linear_solver == "spsolve":                                                                                                  # Direct pre-factorized solver.
                    un       = solve1(RHS)                                                                                                      # Solve global system for time level k=1.
                elif linear_solver == "bicgstab":                                                                                               # Iterative solver (BiCGStab).
                    un, info = bicgstab(A1, RHS)                                                                                                # Solve system.
                    if info != 0:                                                                                                               # If not strictly converged.
                        converged = False                                                                                                       # Mark convergence failure.
                        if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")                           # Warn on convergence issues.
                elif linear_solver == "gmres":                                                                                                  # Iterative solver (GMRES).
                    un, info = gmres(A1, RHS)                                                                                                   # Solve system.
                    if info != 0:                                                                                                               # If not strictly converged.
                        converged = False                                                                                                       # Mark convergence failure.
                        if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")                              # Warn on convergence issues.
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.
            else:                                                                                                                       # Following time steps (k>=2).
                RHS             = B2.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]                                                               # Right-hand side from previous steps.
                RHS[boun_n]     = u_ap[boun_n, k]                                                                                       # Inject exact boundary condition for time k>=2.
                
                if linear_solver == "spsolve":                                                                                          # Direct solver (spsolve).
                    un = solve2(RHS)                                                                                                    # Solve system.
                elif linear_solver == "bicgstab":                                                                                       # Iterative solver (BiCGStab).
                    un, info = bicgstab(A2, RHS, x0=u_ap[:, k-1])                                                                       # Solve system.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}) at time {k}.")                   # Warn on convergence issues.
                elif linear_solver == "gmres":                                                                                          # Iterative solver (GMRES).
                    un, info = gmres(A2, RHS, x0=u_ap[:, k-1])                                                                          # Solve system.
                    if info != 0:                                                                                                       # If not strictly converged.
                        converged = False                                                                                               # Mark convergence failure.
                        if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}) at time {k}.")                      # Warn on convergence issues.
                    
                u_ap[inne_n, k] = un[inne_n]                                                                                            # Update interior nodes.

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info("\tSolver finished successfully.")                                                                                  # Print completion message.
        
    compute_time = time.perf_counter() - start_time                                                                                     # Calculate total compute time.
    u_ap_packed = repack_solution(p_orig, u_ap)                                                                                         # Repack into original Pandas/Xarray format.
    return SolverResult(solution=u_ap_packed, neighbors=vec, converged=converged, compute_time=compute_time)                            # Return structured result.
