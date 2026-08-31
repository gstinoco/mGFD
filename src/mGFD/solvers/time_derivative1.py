"""
TimeDerivative1 — First-order Transient PDEs solver

Overview:
    Numerical solver for first-order transient PDEs (like the Heat equation) using a Meshless
    Generalized Finite Difference scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    TimeDerivative1     Main solver function for parabolic/first-order transient problems.

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
import time                                                                                                                             # Timing for SolverResult.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Callable, Optional, Tuple, List, Union, Any                                                                          # Type hinting.

from mGFD.temporal.cfl import estimate_cfl_dt                                                                                        # Temporal CFL and time discretization utility.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.utils.adapters import extract_cloud, repack_solution                                                                          # Pandas/Xarray adapters.
from mGFD.exceptions import CloudShapeError, InputTypeError, OperatorFormatError, ParameterError                                         # Custom exceptions.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def TimeDerivative1(p: Union[np.ndarray, Any], 
                    f: Union[Callable, np.ndarray, float, int, Any],
                    t: Optional[int] = None,
                    coef: List[float] = [],
                    operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2], [0]]),
                    upwind: bool = False,
                    vec: Optional[np.ndarray] = None,
                    nvec: int = 20,
                    implicit: bool = False,
                    lam: float = 0.5,
                    device: str = "cpu",
                    verbose: bool = True,
                    cfl: Optional[float] = 0.5,
                    dt: Optional[float] = None) -> SolverResult:
    """
    Numerical solution of partial differential equations with first-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \\frac{\\partial u}{\\partial t} = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        f                           Union           Boundary/Initial condition f(x, y, t) or data array.
        t                           Optional[int]   Number of time steps to compute (optional if cfl/dt provided).
        coef                        List            Physical coefficients for the problem formulation.
        operator    6 x 1           ndarray         Array with the weights for the operator.
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
        implicit                    bool            If True, uses an implicit (theta) time integration scheme.
        lam                         float           Theta parameter for implicit scheme (0=explicit, 1=fully implicit, 0.5=Crank-Nicolson).
        device                      str             'cpu' or 'cuda'.
        verbose                     bool            If True, prints solver progress.
        cfl                         Optional[float] Target CFL number for adaptive time-stepping (default 0.5).
        dt                          Optional[float] Explicit time step size (optional).
    
    Output:
        SolverResult                Returns a structured result containing the approximation, neighbor list, and CFL metrics.
    """
    
    start_time = time.perf_counter()                                                                                                    # Start solver timer.
    
    # 0. Input validation
    p_orig = p                                                                                                                          # Preserve original format.
    p = extract_cloud(p)                                                                                                                # Extract NumPy array if Pandas/Xarray.
    f = extract_cloud(f) if not callable(f) and not isinstance(f, (float, int)) else f                                                  # Extract array from Pandas/Xarray for f.
    
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                  # Raise explicit error on bad input.
    if not (callable(f) or isinstance(f, (np.ndarray, float, int))):                                                                    # Validate initial/boundary condition type.
        raise InputTypeError("Forcing function 'f' must be a callable, ndarray, or numeric constant.")                                 # Raise explicit error on bad input.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients.")                                       # Raise explicit error on bad input.
    if device not in ["cpu", "cuda"]:                                                                                                   # Validate device choice.
        raise ParameterError("Unsupported device. Choose 'cpu' or 'cuda'.")                                                             # Raise explicit error.

    # 1. CFL and Time-stepping calculation
    target_cfl = cfl if cfl is not None else 0.5                                                                                        # Default target CFL.
    dt_est, t_est, actual_cfl = estimate_cfl_dt(p, operator, cfl=target_cfl, order=1, vec=vec, t_end=1.0)                              # Estimate CFL baseline metrics.

    if dt is not None:
        if dt <= 0:
            raise ParameterError("Time step 'dt' must be positive.")
        dt_use = float(dt)
        t_use = max(1, int(np.ceil(1.0 / dt_use)))
        actual_cfl = actual_cfl * (dt_use / dt_est) if dt_est > 0 else target_cfl
    elif t is not None:
        if not isinstance(t, int) or t <= 0:
            raise ParameterError("Number of time steps 't' must be a positive integer.")
        t_use = t
        dt_use = 1.0 / t_use
        actual_cfl = actual_cfl * (dt_use / dt_est) if dt_est > 0 else target_cfl
    else:
        dt_use = dt_est
        t_use = t_est

    if verbose:
        logger.info(f"Adaptive Time-Stepping: t={t_use}, dt={dt_use:.6e}, CFL={actual_cfl:.4f}")

    # 2. Dispatch solver execution
    if device == "cpu":
        from mGFD.solvers._backends.cpu.time_derivative1 import solve_cpu
        u_ap, vec, converged = solve_cpu(p, f, t_use, coef, operator, upwind, vec, nvec, implicit, lam, verbose)
    else:
        from mGFD.solvers._backends.cuda.time_derivative1 import solve_cuda
        u_ap, vec, converged = solve_cuda(p, f, t_use, coef, operator, upwind, vec, nvec, implicit, lam, verbose)
        
    compute_time = time.perf_counter() - start_time                                                                                     # Calculate total compute time.
    u_ap_packed = repack_solution(p_orig, u_ap)                                                                                         # Repack into original Pandas/Xarray format.
    return SolverResult(solution=u_ap_packed, neighbors=vec, converged=converged, compute_time=compute_time,
                        dt=dt_use, cfl=actual_cfl, t_steps=t_use)                                                                      # Return structured result with CFL.
