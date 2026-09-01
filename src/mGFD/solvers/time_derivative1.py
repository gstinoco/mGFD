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

from mGFD.temporal.cfl import estimate_cfl_dt                                                                                           # Temporal CFL and time discretization utility.
from mGFD.solvers.results import SolverResult                                                                                           # Standard solver output structure.
from mGFD.utils.adapters import extract_cloud, repack_solution                                                                          # Pandas/Xarray adapters.
from mGFD.exceptions import CloudShapeError, InputTypeError, OperatorFormatError, ParameterError                                        # Custom exceptions.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def TimeDerivative1(p: Union[np.ndarray, Any], 
                    f: Union[Callable, np.ndarray, float, int, Any] = 0.0,                                                              # Assign f: Union[Callable, np.ndarray, float, int, Any].
                    t: Optional[int] = None,                                                                                            # Assign t: Optional[int].
                    coef: List[float] = [],                                                                                             # Assign coef: List[float].
                    operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2], [0]]),                                                   # Assign operator: np.ndarray.
                    upwind: bool = False,                                                                                               # Assign upwind: bool.
                    vec: Optional[np.ndarray] = None,                                                                                   # Assign vec: Optional[np.ndarray].
                    nvec: int = 20,                                                                                                     # Assign nvec: int.
                    implicit: bool = False,                                                                                             # Assign implicit: bool.
                    lam: float = 0.5,                                                                                                   # Assign lam: float.
                    device: str = "cpu",                                                                                                # Assign device: str.
                    verbose: bool = True,                                                                                               # Assign verbose: bool.
                    cfl: Optional[float] = 0.5,                                                                                         # Assign cfl: Optional[float].
                    dt: Optional[float] = None,                                                                                         # Assign dt: Optional[float].
                    ic: Optional[Union[Callable, np.ndarray, float, int, Any]] = None,                                                  # Assign ic: Optional[Union[Callable, np.ndarray, float, int, Any]].
                    bc: Optional[Union[Callable, np.ndarray, float, int, Any]] = None,                                                  # Assign bc: Optional[Union[Callable, np.ndarray, float, int, Any]].
                    source: Optional[Union[Callable, np.ndarray, float, int, Any]] = None,                                              # Assign source: Optional[Union[Callable, np.ndarray, float, int, Any]].
                    t_span: Tuple[float, float] = (0.0, 1.0)) -> SolverResult:                                                          # Assign t_span: Tuple[float, float].
    """
    Numerical solution of partial differential equations with first-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \\frac{\\partial u}{\\partial t} = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u + F_{source}(x, y, t)
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        f                           Union           Legacy/unified Boundary & Initial condition f(x, y, t) or data array.
        t                           Optional[int]   Number of time steps to compute (optional if cfl/dt provided).
        coef                        List            Physical coefficients for the problem formulation.
        operator    6 x 1           ndarray         Array with the weights for the operator [D, E, A, B, C, F_react].
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
        ic                          Optional        Initial condition u0(x, y) at t_span[0] across all nodes.
        bc                          Optional        Boundary condition phi(x, y, t) on boundary nodes over time.
        source                      Optional        Independent heat/forcing source term F(x, y, t).
        t_span      Tuple[float]                    Physical time domain tuple (t_start, t_end). Default (0.0, 1.0).
    
    Output:
        SolverResult                Returns a structured result containing the approximation, neighbor list, and CFL metrics.
    """
    
    start_time = time.perf_counter()                                                                                                    # Start solver timer.
    
    # 0. Input validation
    p_orig = p                                                                                                                          # Preserve original format.
    p = extract_cloud(p)                                                                                                                # Extract NumPy array if Pandas/Xarray.
    f = extract_cloud(f) if not callable(f) and not isinstance(f, (float, int)) else f                                                  # Extract array from Pandas/Xarray for f.
    ic = extract_cloud(ic) if ic is not None and not callable(ic) and not isinstance(ic, (float, int)) else ic                          # Extract array from Pandas/Xarray for ic.
    bc = extract_cloud(bc) if bc is not None and not callable(bc) and not isinstance(bc, (float, int)) else bc                          # Extract array from Pandas/Xarray for bc.
    source = extract_cloud(source) if source is not None and not callable(source) and not isinstance(source, (float, int)) else source  # Extract array from Pandas/Xarray for source.
    
    if not isinstance(p, np.ndarray) or p.ndim != 2 or p.shape[1] != 3:                                                                 # Validate point cloud array shape and type.
        raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with 3 columns (x, y, flag).")                                  # Raise explicit error on bad input.
    if f is None and ic is None and bc is None:                                                                                         # Validate condition presence.
        raise InputTypeError("Must provide either 'f' or both 'ic' and 'bc'.")                                                          # Raise explicit error.
    if f is not None and not (callable(f) or isinstance(f, (np.ndarray, float, int))):                                                  # Validate forcing function type.
        raise InputTypeError("Forcing function 'f' must be a callable, ndarray, or numeric constant.")                                  # Raise explicit error.
    if ic is not None and not (callable(ic) or isinstance(ic, (np.ndarray, float, int))):                                               # Validate initial condition type.
        raise InputTypeError("Initial condition function 'ic' must be a callable, ndarray, or numeric constant.")                       # Raise explicit error.
    if bc is not None and not (callable(bc) or isinstance(bc, (np.ndarray, float, int))):                                               # Validate boundary condition type.
        raise InputTypeError("Boundary condition function 'bc' must be a callable, ndarray, or numeric constant.")                      # Raise explicit error.
    if not isinstance(operator, np.ndarray) or operator.shape[0] < 5:                                                                   # Validate operator array.
        raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients.")                                       # Raise explicit error on bad input.
    if device not in ["cpu", "cuda"]:                                                                                                   # Validate device choice.
        raise ParameterError("Unsupported device. Choose 'cpu' or 'cuda'.")                                                             # Raise explicit error.

    # 1. CFL and Time-stepping calculation
    t_domain_length = t_span[1] - t_span[0]                                                                                             # Length of physical time domain.
    target_cfl      = cfl if cfl is not None else 0.5                                                                                   # Default target CFL.
    dt_est, t_est, actual_cfl = estimate_cfl_dt(p, operator, cfl=target_cfl, order=1, vec=vec, t_end=t_domain_length)                   # Estimate CFL baseline metrics across t_span.

    if dt is not None:                                                                                                                  # Check if explicit dt provided.
        if dt <= 0:                                                                                                                     # Validate dt positive.
            raise ParameterError("Time step 'dt' must be positive.")                                                                    # Raise error for non-positive dt.
        dt_use = float(dt)                                                                                                              # Convert dt to float.
        t_use  = max(1, int(np.ceil(t_domain_length / dt_use)))                                                                         # Compute required time step count.
        actual_cfl = actual_cfl * (dt_use / dt_est) if dt_est > 0 else target_cfl                                                       # Recalculate actual CFL.
    elif t is not None:                                                                                                                 # Check if explicit step count provided.
        if not isinstance(t, int) or t <= 0:                                                                                            # Validate integer step count.
            raise ParameterError("Number of time steps 't' must be a positive integer.")                                                # Raise error for invalid step count.
        t_use  = t                                                                                                                      # Use explicit step count.
        dt_use = t_domain_length / t_use                                                                                                # Calculate corresponding dt.
        actual_cfl = actual_cfl * (dt_use / dt_est) if dt_est > 0 else target_cfl                                                       # Recalculate actual CFL.
    else:                                                                                                                               # Automatic adaptive time-stepping.
        dt_use = dt_est                                                                                                                 # Use estimated dt.
        t_use  = t_est                                                                                                                  # Use estimated step count.

    if verbose:                                                                                                                         # Log adaptive time-stepping info.
        logger.info(f"Adaptive Time-Stepping: t={t_use}, dt={dt_use:.6e}, CFL={actual_cfl:.4f}, t_span={t_span}")                       # Report configuration parameters.

    # 2. Dispatch solver execution
    if device == "cpu":                                                                                                                 # Execute CPU backend solver.
        from mGFD.solvers._backends.cpu.time_derivative1 import solve_cpu                                                               # Import CPU solver.
        u_ap, vec, converged = solve_cpu(p, f, t_use, coef, operator, upwind, vec, nvec, implicit, lam, verbose=verbose, ic=ic, bc=bc, source=source, t_span=t_span) # Assign u_ap, vec, converged.
    else:                                                                                                                               # Execute CUDA backend solver.
        from mGFD.solvers._backends.cuda.time_derivative1 import solve_cuda                                                             # Import CUDA solver.
        u_ap, vec, converged = solve_cuda(p, f, t_use, coef, operator, upwind, vec, nvec, implicit, lam, verbose=verbose, ic=ic, bc=bc, source=source, t_span=t_span) # Assign u_ap, vec, converged.

    compute_time = time.perf_counter() - start_time                                                                                     # Measure total solver execution duration.

    # 3. Repackage solution format
    solution = repack_solution(p_orig, u_ap)                                                                                            # Format output matching input container type.

    return SolverResult(                                                                                                                # Return standardized solver result structure.
        solution=solution,                                                                                                              # Assign solution.
        neighbors=vec,                                                                                                                  # Assign neighbors.
        converged=converged,                                                                                                            # Assign converged.
        compute_time=compute_time,                                                                                                      # Assign compute_time.
        dt=dt_use,                                                                                                                      # Assign dt.
        cfl=actual_cfl,                                                                                                                 # Assign cfl.
        t_steps=t_use                                                                                                                   # Assign t_steps.
    )                                                                                                                                   # Execute statement.
