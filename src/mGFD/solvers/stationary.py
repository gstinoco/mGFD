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


from typing import Callable, Optional, Tuple, List, Union, Any                                                                          # Type hinting.

from mGFD.exceptions import CloudShapeError, InputTypeError, DimensionMismatchError, OperatorFormatError, ParameterError                # Custom exceptions.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.utils.adapters import extract_cloud, repack_solution                                                                          # Pandas/Xarray adapters.



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

    # 1. Enrutamiento (Dispatch)
    if device == "cpu":
        from mGFD.solvers._backends.cpu.stationary import solve_cpu
        u_ap, vec, converged = solve_cpu(p, phi, f, operator, upwind, vec, nvec, linear_solver, preconditioner, matrix_free, verbose)
    else:
        from mGFD.solvers._backends.cuda.stationary import solve_cuda
        u_ap, vec, converged = solve_cuda(p, phi, f, operator, upwind, vec, nvec, linear_solver, preconditioner, matrix_free, verbose)
        
    compute_time = time.perf_counter() - start_time                                                                                     # Calculate total compute time.
    u_ap_packed = repack_solution(p_orig, u_ap)                                                                                         # Repack into original Pandas/Xarray format.
    return SolverResult(solution=u_ap_packed, neighbors=vec, converged=converged, compute_time=compute_time)                            # Return structured result.
