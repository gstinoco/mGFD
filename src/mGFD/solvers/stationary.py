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

from scipy.sparse.linalg import spsolve, bicgstab, gmres                                                                                # Sparse linear solvers.
from typing import Callable, Optional, Tuple, List, Union, Any                                                                          # Type hinting.

from mGFD.exceptions import CloudShapeError, InputTypeError, DimensionMismatchError, OperatorFormatError, ParameterError                # Custom exceptions.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.core.adapters import extract_cloud, repack_solution                                                                           # Pandas/Xarray adapters.
import mGFD.core.gammas as Gammas                                                                                                       # Gammas calculation and sparse matrix builder.
import mGFD.core.neighbors as Neighbors                                                                                                 # Neighbor search routines.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def Stationary(p: Union[np.ndarray, Any], 
               phi: Union[Callable, np.ndarray, float, int, Any],
               f: Union[Callable, np.ndarray, float, int, Any],
               operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2], [0]]),
               upwind: bool = False,
               vec: Optional[np.ndarray] = None,
               nvec: int = 20,
               linear_solver: str = "spsolve",
               preconditioner: Optional[str] = None,
               matrix_free: bool = False,
               device: str = "cpu",
               verbose: bool = True) -> SolverResult:
    """
    Numerical solution of partial differential equations with no time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu = -f(x,y)
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        phi                         Union           Boundary condition phi(x, y) or constant/array data.
        f                           Union           Right-hand side f(x, y) or constant/array data.
        operator    6 x 1           ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
        linear_solver               str             Algebraic backend: 'spsolve', 'bicgstab', or 'gmres'.
        preconditioner              Optional[str]   Preconditioner method: 'ilu', 'jacobi', or None.
        matrix_free                 bool            If True, uses on-the-fly matrix-vector multiplication (requires iterative solver).
        device                      str             'cpu' or 'cuda'. If 'cuda', offloads the sparse solvers to the GPU (requires CuPy).
        verbose                     bool            If True, prints solver progress.
    
    Output:
        SolverResult                Returns a structured result containing the approximation and the neighbor list.
    """
    
    start_time = time.perf_counter()                                                                                                    # Start solver timer.
    
    # 0. Input validation
    p_orig = p                                                                                                                          # Preserve original format.
    p = extract_cloud(p)                                                                                                                # Extract NumPy array if Pandas/Xarray.
    phi = extract_cloud(phi) if not callable(phi) and not isinstance(phi, (float, int)) else phi                                        # Extract array from Pandas/Xarray for phi.
    f = extract_cloud(f) if not callable(f) and not isinstance(f, (float, int)) else f                                                  # Extract array from Pandas/Xarray for f.
    
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                  # Raise explicit error on bad input.
    if not (callable(phi) or isinstance(phi, (np.ndarray, float, int))):                                                                # Validate boundary condition type.
        raise InputTypeError("Boundary condition 'phi' must be a callable, ndarray, or numeric constant.")                              # Raise explicit error on bad input.
    if not (callable(f) or isinstance(f, (np.ndarray, float, int))):                                                                    # Validate RHS type.
        raise InputTypeError("Right-hand side 'f' must be a callable, ndarray, or numeric constant.")                                   # Raise explicit error on bad input.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients.")                                       # Raise explicit error on bad input.
    if device not in ["cpu", "cuda"]:                                                                                                   # Validate device choice.
        raise ParameterError("Unsupported device. Choose 'cpu' or 'cuda'.")                                                             # Raise explicit error.
    if device == "cuda" and matrix_free:                                                                                                # Validate matrix_free and cuda.
        raise ParameterError("matrix_free=True is incompatible with device='cuda'.")                                                    # CuPy requires explicit sparse matrices.

    # 1. Variable initialization
    m      = len(p[:, 0])                                                                                                               # The total number of nodes is calculated.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f"Solving Stationary problem for {m} nodes...")                                                                     # Print solver progress.
        
    u_ap   = np.zeros([m])                                                                                                              # u_ap initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                                                            # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                                                                               # Save the inner nodes.

    # 2. Extract advection velocities for Upwind scheme.
    if upwind:                                                                                                                          # If an Upwind stencil is requested.
        a = operator[0][0] if operator.ndim == 2 else operator[0]                                                                       # Value of the velocity on x.
        b = operator[1][0] if operator.ndim == 2 else operator[1]                                                                       # Value of the velocity on y.

    # 3. Apply Boundary Conditions
    if callable(phi):                                                                                                                   # If boundary condition is a function.
        u_ap[boun_n] = np.asarray(phi(p[boun_n, 0], p[boun_n, 1]))                                                                      # The boundary condition is assigned via function.
    elif isinstance(phi, np.ndarray):                                                                                                   # If boundary condition is an array.
        u_ap[boun_n] = phi[boun_n]                                                                                                      # The boundary condition is assigned via array.
    elif isinstance(phi, (int, float)):                                                                                                 # If boundary condition is a constant.
        u_ap[boun_n] = phi                                                                                                              # The boundary condition is assigned via constant.
    
    # 4. Neighbor search for all the nodes.
    if vec is None:                                                                                                                     # If no neighbor list is provided.
        if upwind:                                                                                                                      # If an Upwind stencil is requested.
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                                                                     # Neighbor search with the proper routine.
        else:                                                                                                                           # All the other cases.
            vec = Neighbors.compute_neighbors(p, nvec)                                                                                  # Neighbor search with the proper routine.

    # 5. Computation of Gamma values
    L = operator[:-1]                                                                                                                   # The values of the differential operator are assigned.
    if matrix_free:                                                                                                                     # If matrix-free computation is requested.
        from scipy.sparse.linalg import LinearOperator                                                                                  # Import LinearOperator.
        K_matvec = Gammas.compute_K_matvec(p, vec, L)                                                                                   # Generate the on-the-fly matrix-vector closure.
        K = LinearOperator(shape=(m, m), matvec=K_matvec, dtype=np.float64)                                                             # type: ignore
    else:                                                                                                                               # Standard dense matrix allocation.
        K = Gammas.compute_sparse_matrix(p, vec, L)                                                                                     # K computation with the required Gammas (Sparse).
        if device == "cuda":                                                                                                            # If GPU acceleration is requested.
            try:
                import cupy as cp                                                                                                       # Import CuPy array library.
                from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix                                                              # Import CuPy sparse matrix.
            except ImportError:                                                                                                         # If not installed.
                raise ImportError("CuPy is not installed. Please install it with 'pip install mGFD[gpu]' to use device='cuda'.")        # Friendly error message.
            K = cp_csr_matrix(K)                                                                                                        # Transfer explicit matrix K to the GPU.
    
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                                                                           # Right-hand side of the equation.
    if device == "cuda":                                                                                                                # If GPU acceleration is requested.
        R = cp.asarray(R)                                                                                                               # Transfer RHS vector R to the GPU.
    
    converged = True                                                                                                                    # Assume convergence by default.
    if matrix_free and preconditioner is not None:                                                                                      # Validate preconditioner in matrix-free mode.
        raise ParameterError("Preconditioners are not currently supported in matrix_free=True mode.")                                   # Preconditioners need explicit matrix elements.
        
    M = None if matrix_free else (Gammas.compute_preconditioner(K, preconditioner) if preconditioner else None)                         # type: ignore

    # 6. Solution of the linear system (Generalized Finite Differences)
    if linear_solver == "spsolve":                                                                                                      # Direct SciPy solver.
        if matrix_free:                                                                                                                 # If matrix-free computation is enabled.
            raise ParameterError("Direct solver 'spsolve' is incompatible with matrix_free=True. Use 'gmres' or 'bicgstab'.")           # Direct solvers need the explicit sparse matrix.
        if device == "cuda":                                                                                                            # GPU Direct Solver.
            from cupyx.scipy.sparse.linalg import spsolve as cp_spsolve                                                                 # Import CuPy spsolve.
            un = cp_spsolve(K, R)                                                                                                       # Solve via CuPy.
        else:                                                                                                                           # CPU Direct Solver.
            un = spsolve(K, R)                                                                                                          # Direct sparse solve of the linear system.
    elif linear_solver == "bicgstab":                                                                                                   # Iterative SciPy solver (BiCGStab).
        if device == "cuda":                                                                                                            # GPU Iterative Solver.
            from cupyx.scipy.sparse.linalg import bicgstab as cp_bicgstab                                                               # Import CuPy bicgstab.
            un, info = cp_bicgstab(K, R, M=M)                                                                                           # Iterative solve via CuPy.
        else:                                                                                                                           # CPU Iterative Solver.
            un, info = bicgstab(K, R, M=M)                                                                                              # Iterative sparse solve of the linear system.
        if info != 0:                                                                                                                   # If not strictly converged.
            converged = False                                                                                                           # Mark convergence failure.
            if verbose: logger.warning(f"BiCGStab did not converge perfectly (code {info}).")                                           # Warn user.
    elif linear_solver == "gmres":                                                                                                      # Iterative SciPy solver (GMRES).
        if device == "cuda":                                                                                                            # GPU Iterative Solver.
            from cupyx.scipy.sparse.linalg import gmres as cp_gmres                                                                     # Import CuPy gmres.
            un, info = cp_gmres(K, R, M=M)                                                                                              # Iterative solve via CuPy.
        else:                                                                                                                           # CPU Iterative Solver.
            un, info = gmres(K, R, M=M)                                                                                                 # Iterative sparse solve of the linear system.
        if info != 0:                                                                                                                   # If not strictly converged.
            converged = False                                                                                                           # Mark convergence failure.
            if verbose: logger.warning(f"GMRES did not converge perfectly (code {info}).")                                              # Warn user.
    else:                                                                                                                               # Unsupported solver request.
        raise ParameterError(f"Unsupported linear_solver '{linear_solver}'. Choose from 'spsolve', 'bicgstab', 'gmres'.")               # Raise explicit error.
        
    if device == "cuda":                                                                                                                # If a GPU was used.
        un = getattr(un, "get")()                                                                                                       # Pull the solution vector back to CPU RAM.

    u_ap[inne_n] = un[inne_n]                                                                                                           # Save the computed solution to the interior nodes.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info("\tSolver finished successfully.")                                                                                  # Print completion message.
        
    compute_time = time.perf_counter() - start_time                                                                                     # Calculate total compute time.
    u_ap_packed = repack_solution(p_orig, u_ap)                                                                                         # Repack into original Pandas/Xarray format.
    return SolverResult(solution=u_ap_packed, neighbors=vec, converged=converged, compute_time=compute_time)                            # Return structured result.
