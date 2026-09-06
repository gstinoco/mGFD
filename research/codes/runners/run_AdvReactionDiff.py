"""
run_AdvReactionDiff — Reference batch runner for Advection-Reaction-Diffusion PDEs with Forcing Source Term

Overview:
    This script runs an Advection-Reaction-Diffusion problem with an independent forcing source term F_source(x, y, t)
    and reaction coefficient F_react u across all point clouds under Data/, using the meshless mGFD solver.

PDE Formulation:
    du/dt + a * du/dx + b * du/dy = nu * (d2u/dx2 + d2u/dy2) + F_react * u + F_source(x, y, t)

Manufactured Exact Solution:
    u_ex(x, y, t) = cos(pi * x) * cos(pi * y) * exp(-lambda_val * t)

Forcing Term:
    F_source(x, y, t) = (2 * nu * pi^2 - lambda_val - F_react) * u_ex - pi * exp(-lambda_val * t) * (a * sin(pi * x) * cos(pi * y) + b * cos(pi * x) * sin(pi * y))

Public API:
    process_cloud       Process a single point cloud for the Advection-Reaction-Diffusion equation.
    main                Batch runner entry point for the Advection-Reaction-Diffusion problem.

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
    August, 2026.
Last Modification:
    September, 2026.
"""

## Library importation.
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.path manipulation.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Numerical arrays and math.

from typing import List, Any, Optional, Tuple                                                                                           # Type hinting.

from mGFD import Cloud, Dirichlet, PDE, Solver                                                                                          # OOP Architecture classes.
from mGFD.solvers.results import SolverResult                                                                                           # Standardized SolverResult dataclass.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_transient                                                                                               # Plotting utilities for results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory.
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for runs.
from utils.batch_utils import run_batch_suite, save_metrics                                                                             # Dataset loading + batch metrics helpers.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.

## Problem parameters.
a: float                    = 0.1                                                                                                       # Advection coefficient x.
b: float                    = 0.1                                                                                                       # Advection coefficient y.
nu: float                   = 0.2                                                                                                       # Diffusion coefficient.
F_react: float              = 0.5                                                                                                       # Reaction coefficient.
lambda_val: float           = 0.5                                                                                                       # Temporal decay rate.
t_span: Tuple[float, float] = (0.0, 2.0)                                                                                                # Physical time domain [t_start, t_end].

def f_analytical(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Manufactured analytical exact solution for the Advection-Reaction-Diffusion system.

    Evaluates:
        u_ex(x, y, t) = cos(pi * x) * cos(pi * y) * exp(-lambda * t)

    Input:
        x               m           ndarray         x coordinates of the point cloud.
        y               m           ndarray         y coordinates of the point cloud.
        t_val                       float           Time value at evaluation.
        coef            5           list            Parameters [nu, a, b, lambda_val, F_react].

    Output:
        u_ex            m           ndarray         Evaluated exact state array.
    """
    return np.cos(np.pi * x) * np.cos(np.pi * y) * np.exp(-lambda_val * t_val)                                                          # Return output values.

def ic(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Initial condition for the Advection-Reaction-Diffusion problem.

    Evaluates u_0(x, y) = f_analytical(x, y, t_start, coef) at the start of simulation.

    Input:
        x               m           ndarray         x coordinates of the point cloud.
        y               m           ndarray         y coordinates of the point cloud.
        t_val                       float           Time value (t_start).
        coef            5           list            Parameters [nu, a, b, lambda_val, F_react].

    Output:
        u0              m           ndarray         Initial state vector.
    """
    return f_analytical(x, y, t_span[0], coef)                                                                                          # Return output values.

def bc(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Boundary condition for the Advection-Reaction-Diffusion problem.

    Evaluates Dirichlet boundary values phi(x, y, t) on boundary nodes.

    Input:
        x               m           ndarray         x coordinates of boundary nodes.
        y               m           ndarray         y coordinates of boundary nodes.
        t_val                       float           Time value at current time step.
        coef            5           list            Parameters [nu, a, b, lambda_val, F_react].

    Output:
        phi_val         m           ndarray         Evaluated boundary state vector.
    """
    return f_analytical(x, y, t_val, coef)                                                                                              # Return output values.

def source(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """
    Non-homogeneous forcing source term F_source(x, y, t).

    Satisfies the manufactured solution balance equation:
        F_source = (2 * nu * pi^2 - lambda - F_react) * u - (a * u_x + b * u_y)

    Input:
        x               m           ndarray         x coordinates of the point cloud.
        y               m           ndarray         y coordinates of the point cloud.
        t_val                       float           Time value at current step.
        coef            5           list            Parameters [nu, a, b, lambda_val, F_react].

    Output:
        f_src           m           ndarray         Evaluated source forcing vector.
    """
    decay    = np.exp(-lambda_val * t_val)                                                                                              # Compute temporal exponential decay.
    u_val    = np.cos(np.pi * x) * np.cos(np.pi * y) * decay                                                                            # Analytical state.
    adv_part = np.pi * decay * (a * np.sin(np.pi * x) * np.cos(np.pi * y) + b * np.cos(np.pi * x) * np.sin(np.pi * y))                  # Advective derivative combination.
    return (2 * nu * np.pi**2 - lambda_val - F_react) * u_val - adv_part                                                                # Return balanced source values.

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the Advection-Reaction-Diffusion problem with source term on a single point cloud file.

    Input:
        dataset         1           str             Name of the geometry region/dataset (e.g. 'Cuitzeo').
        scale           1           str             Scale refinement identifier (e.g. '1').
        cloud_path      1           str             Absolute path to the point cloud CSV file.
        results_path    1           str             Base path to the Results output directory.
        save            1           bool            If True, saves graphical animation outputs.
        verbose         1           bool            If True, prints execution logs.
        **kwargs                    Any             Optional configuration overrides (nvec, device, upwind, cfl, t).

    Output:
        None
    """
    # 0. Input validation
    if not isinstance(dataset, str):                                                                                                    # Validate dataset argument.
        raise TypeError("Dataset name must be a string.")                                                                               # Raise explicit error on bad input.
    if not isinstance(scale, str):                                                                                                      # Validate scale argument.
        raise TypeError("Scale must be a string.")                                                                                      # Raise explicit error on bad input.
    if not isinstance(cloud_path, str) or not os.path.exists(cloud_path):                                                               # Validate cloud path.
        raise ValueError("Cloud path must be a valid existing file path.")                                                              # Raise explicit error on bad input.

    # 1. Variable initialization
    region_id = f'{dataset}/{scale}'                                                                                                    # Region identifier.
    out_dir   = os.path.join(results_path, 'AdvReactionDiff', dataset)                                                                  # Destination path.
    os.makedirs(out_dir, exist_ok=True)                                                                                                 # Create directory.

    if verbose:                                                                                                                         # Verbosity check.
        logger.info(f'Working on region: {region_id}')                                                                                  # Log active region.

    # 2. Data Loading & Neighbor Cache
    nvec            = kwargs.get('nvec', 20)                                                                                            # Assign nvec.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Assign verbose_solvers.
    device          = kwargs.get('device', 'cpu')                                                                                       # Assign device.
    upwind          = kwargs.get('upwind', True)                                                                                        # Assign upwind.
    cfl_param       = kwargs.get('cfl', 0.5)                                                                                            # Assign cfl_param.
    t_param         = kwargs.get('t', None)                                                                                             # Assign t_param.
    config_id       = f'nvec_{nvec}_{device}_upwind_{upwind}'                                                                           # Assign config_id.

    # 2. Data Loading
    p         = load_points(cloud_path, verbose=False)                                                                                  # Assign p.
    L         = np.vstack([[-a], [-b], [2 * nu], [0], [2 * nu], [F_react]])                                                             # Assign differential operator L.

    # 3. Solver Execution via OOP Interface
    cloud     = Cloud.from_array(p)                                                                                                     # Instantiate Cloud.
    domain    = cloud.set_boundary(Dirichlet(bc))                                                                                       # Set Dirichlet boundary.
    pde       = PDE(operator=L, ic=ic, source=source, order=1, coef=[nu, a, b])                                                         # Formulate PDE.
    solver    = Solver(domain, pde, device=device, cfl=cfl_param, nvec=nvec, upwind=upwind, verbose=verbose_solvers, implicit=True)     # Create Solver.
    dt_input  = ((t_span[1] - t_span[0]) / int(t_param)) if t_param is not None else None                                               # Explicit dt or None for adaptive CFL.
    result    = solver.solve(t_span=t_span, dt=dt_input)                                                                                # Solve PDE with adaptive CFL or custom step.
    u_ap, vec = result.solution, result.neighbors                                                                                       # Unpack solution and neighbors.
    comp_time = result.compute_time                                                                                                     # Get compute time.

    # 4. Exact Solution and Metrics
    t_used    = result.t_steps if result.t_steps is not None else (int(t_param) if t_param is not None else 2)                          # Retrieve simulated time steps from solver.
    T_arr     = np.linspace(t_span[0], t_span[1], t_used)                                                                               # Reconstruct time vector.
    X         = p[:, 0, None]                                                                                                           # Node X coordinates broadcast shape (m, 1).
    Y         = p[:, 1, None]                                                                                                           # Node Y coordinates broadcast shape (m, 1).
    T         = T_arr[None, :]                                                                                                          # Time steps broadcast shape (1, t).
    u_ex      = np.cos(np.pi * X) * np.cos(np.pi * Y) * np.exp(-lambda_val * T)                                                         # Vectorized exact spatio-temporal solution.
    metrics   = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)                                            # Compute comprehensive transient error metrics.

    metrics['CFL']        = result.cfl                                                                                                  # Store CFL number.
    metrics['dt']         = result.dt                                                                                                   # Store time step size.
    metrics['Time_Steps'] = result.t_steps                                                                                              # Store total time step count.

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                                                # Log output message.

    # 5. Output persistence
    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Save metrics using the common utility.

    # 6. Graphical rendering
    if save:                                                                                                                            # Save graphical outputs if requested.
        plot_scales = kwargs.get('plot_scales', ['2'])                                                                                  # Extract target plotting scales.
        if plot_scales is None or 'all' in plot_scales or scale in [str(s) for s in plot_scales]:                                       # Filter rendering for selected scale(s).
            cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                               # Extract clean cloud name.
            appx_nom   = os.path.join(out_dir, f'Appx_{config_id}_scale_{scale}_{cloud_name}')                                          # Define numerical approximation plot path.
            plot_transient(p, u_ap, save=True, nom=appx_nom, title='Transient Approximation', verbose=verbose, t_span=t_span)           # Save transient animation of approximation.

            exact_nom = os.path.join(out_dir, f'Exact_scale_{scale}_{cloud_name}')                                                      # Define theoretical exact solution plot path.
            if not (os.path.exists(exact_nom + '.mp4') or os.path.exists(exact_nom + '.gif')):                                          # Avoid regenerating exact animation if already present.
                plot_transient(p, u_ex, save=True, nom=exact_nom, title='Theoretical Solution', verbose=verbose, t_span=t_span)         # Save exact transient animation.

def main(**kwargs: Any) -> None:
    """
    main
    Run the batch execution for Advection-Reaction-Diffusion problem across all dataset clouds.

    Input:
        **kwargs                    Any             Configuration parameters passed by the parameter sweep orchestrator.

    Output:
        None
    """
    run_batch_suite(process_cloud, DATA_ROOT, RESULTS_ROOT, **kwargs)                                                                   # Execute universal batch orchestrator.

if __name__ == '__main__':                                                                                                              # Evaluate condition.
    main()                                                                                                                              # Execute statement.
