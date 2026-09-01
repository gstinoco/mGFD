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

Date:
    August, 2026.
Last Modification:
    August, 2026.
"""

## Library importation.
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.path manipulation.
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Numerical arrays and math.
import pandas as pd                                                                                                                     # DataFrames for adapter tests.
from typing import List, Any                                                                                                            # Type hinting.

from mGFD import TimeDerivative1                                                                                                        # First-order transient solver.
from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.
from mGFD.viz.graph import plot_transient                                                                                               # Plotting utilities for results.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory.
sys.path.append(BASE_DIR)                                                                                                               # Allow importing from research/codes/utils/.

import utils.metrics as Errors                                                                                                          # Error metrics for runs.
from utils.batch_utils import load_neighbors, save_neighbors, run_batch_suite, save_metrics                                             # Dataset loading + neighbor cache helpers.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

DATA_ROOT: str    = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                     # Input dataset root directory.
RESULTS_ROOT: str = os.path.join(os.path.dirname(BASE_DIR), 'results')                                                                  # Output results root directory.

## Problem parameters.
a: float          = 0.1                                                                                                                 # Advection coefficient x.
b: float          = 0.1                                                                                                                 # Advection coefficient y.
nu: float         = 0.2                                                                                                                 # Diffusion coefficient.
F_react: float    = 0.5                                                                                                                 # Reaction coefficient.
lambda_val: float = 0.5                                                                                                                 # Temporal decay rate.
t_span            = (0.0, 2.0)                                                                                                          # Physical time domain [t_start, t_end].

def f_analytical(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """Analytical manufactured solution u_ex(x, y, t)."""
    return np.cos(np.pi * x) * np.cos(np.pi * y) * np.exp(-lambda_val * t_val)                                                          # Return output values.

def ic(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """Initial condition u_0(x, y) at t_start."""
    return f_analytical(x, y, t_span[0], coef)                                                                                          # Return output values.

def bc(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """Boundary condition phi(x, y, t) on boundary nodes."""
    return f_analytical(x, y, t_val, coef)                                                                                              # Return output values.

def source(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:
    """Forcing source term F_source(x, y, t)."""
    u_val = np.cos(np.pi * x) * np.cos(np.pi * y) * np.exp(-lambda_val * t_val)                                                         # Assign u_val.
    adv_part = np.pi * np.exp(-lambda_val * t_val) * (a * np.sin(np.pi * x) * np.cos(np.pi * y) + b * np.cos(np.pi * x) * np.sin(np.pi * y)) # Assign adv_part.
    return (2 * nu * np.pi**2 - lambda_val - F_react) * u_val - adv_part                                                                # Return output values.

def process_cloud(dataset: str, scale: str, cloud_path: str, results_path: str, save: bool, verbose: bool = True, **kwargs: Any) -> None:
    """
    process_cloud
    Run the Advection-Reaction-Diffusion problem with source term on a single point cloud file.

    Input:
        dataset                     str             Name of the geometry region/dataset (e.g. 'Cuitzeo').
        scale                       str             Scale refinement identifier (e.g. '1').
        cloud_path                  str             Absolute path to the point cloud CSV file.
        results_path                str             Base path to the Results output directory.
        save                        bool            If True, saves graphical animation outputs.
        verbose                     bool            If True, prints execution logs.
        **kwargs                    Any             Optional configuration overrides (nvec, device, upwind, cfl, t).

    Output:
        None
    """
    
    # 1. Output directory setup
    region_id = f'{dataset}/{scale}'                                                                                                    # Region identifier.
    out_dir   = os.path.join(results_path, 'AdvReactionDiff', dataset)                                                                  # Destination path.
    os.makedirs(out_dir, exist_ok=True)                                                                                                 # Create directory.
    
    if verbose:                                                                                                                         # Verbosity check.
        logger.info(f'Working on region: {region_id}')                                                                                  # Log active region.

    nvec            = kwargs.get('nvec', 20)                                                                                            # Assign nvec.
    verbose_solvers = kwargs.get('verbose_solvers', False)                                                                              # Assign verbose_solvers.
    device          = kwargs.get('device', 'cpu')                                                                                       # Assign device.
    upwind          = kwargs.get('upwind', True)                                                                                        # Assign upwind.
    input_types     = kwargs.get('input_types', ['callable'])                                                                           # Assign input_types.
    cfl_param       = kwargs.get('cfl', 0.5)                                                                                            # Assign cfl_param.
    t_param         = kwargs.get('t', None)                                                                                             # Assign t_param.
    config_id       = f'nvec_{nvec}_{device}_upwind_{upwind}'                                                                           # Assign config_id.

    p    = load_points(cloud_path, verbose=False)                                                                                       # Assign p.
    vec0 = load_neighbors(cloud_path, nvec)                                                                                             # Assign vec0.
    
    # Differential operator L = [-a, -b, 2*nu, 0, 2*nu, F_react]
    L = np.vstack([[-a], [-b], [2 * nu], [0], [2 * nu], [F_react]])                                                                     # Assign L.
    
    u_ap, vec = None, vec0                                                                                                              # Assign u_ap, vec.
    comp_time = 0.0                                                                                                                     # Assign comp_time.
    res_main  = None                                                                                                                    # Assign res_main.

    # --- A. Using Callable ---
    if 'callable' in input_types:                                                                                                       # Evaluate condition.
        res_call = TimeDerivative1(                                                                                                     # Assign res_call.
            p,                                                                                                                          # Execute statement.
            ic=ic,                                                                                                                      # Assign ic.
            bc=bc,                                                                                                                      # Assign bc.
            source=source,                                                                                                              # Assign source.
            t_span=t_span,                                                                                                              # Assign t_span.
            t=t_param,                                                                                                                  # Assign t.
            cfl=cfl_param,                                                                                                              # Assign cfl.
            coef=[a, b, nu, F_react],                                                                                                   # Assign coef.
            operator=L,                                                                                                                 # Assign operator.
            upwind=upwind,                                                                                                              # Assign upwind.
            vec=vec0,                                                                                                                   # Assign vec.
            nvec=nvec,                                                                                                                  # Assign nvec.
            implicit=True,                                                                                                              # Assign implicit.
            lam=0.5,                                                                                                                    # Assign lam.
            device=device,                                                                                                              # Assign device.
            verbose=verbose_solvers                                                                                                     # Assign verbose.
        )                                                                                                                               # Execute statement.
        u_ap, vec  = res_call.solution, res_call.neighbors                                                                              # Assign u_ap, vec.
        comp_time  = res_call.compute_time                                                                                              # Assign comp_time.
        res_main   = res_call                                                                                                           # Assign res_main.

    t_used = res_main.t_steps if res_main is not None and res_main.t_steps is not None else (t_param if t_param is not None else 200)   # Evaluate condition.

    # Compute analytical exact solution array across time domain for RMSE calculation
    T_arr = np.linspace(t_span[0], t_span[1], t_used)                                                                                   # Assign T_arr.
    X     = p[:, 0, None]                                                                                                               # Assign X.
    Y     = p[:, 1, None]                                                                                                               # Assign Y.
    T     = T_arr[None, :]                                                                                                              # Assign T.
    u_ex  = np.cos(np.pi * X) * np.cos(np.pi * Y) * np.exp(-lambda_val * T)                                                             # Assign u_ex.

    # --- B. Using Numpy Arrays ---
    if 'array' in input_types:                                                                                                          # Evaluate condition.
        res_arr = TimeDerivative1(                                                                                                      # Assign res_arr.
            p,                                                                                                                          # Execute statement.
            ic=u_ex[:, 0],                                                                                                              # Assign ic.
            bc=u_ex,                                                                                                                    # Assign bc.
            source=source,                                                                                                              # Assign source.
            t_span=t_span,                                                                                                              # Assign t_span.
            t=t_used,                                                                                                                   # Assign t.
            coef=[a, b, nu, F_react],                                                                                                   # Assign coef.
            operator=L,                                                                                                                 # Assign operator.
            upwind=upwind,                                                                                                              # Assign upwind.
            vec=vec0,                                                                                                                   # Assign vec.
            nvec=nvec,                                                                                                                  # Assign nvec.
            implicit=True,                                                                                                              # Assign implicit.
            lam=0.5,                                                                                                                    # Assign lam.
            device=device,                                                                                                              # Assign device.
            verbose=False                                                                                                               # Assign verbose.
        )                                                                                                                               # Execute statement.
        if u_ap is not None:                                                                                                            # Evaluate condition.
            assert np.allclose(u_ap, res_arr.solution), "Mismatch between Callable and Array solver outputs."                           # Execute statement.
        else:                                                                                                                           # Execute fallback branch.
            u_ap, vec  = res_arr.solution, res_arr.neighbors                                                                            # Assign u_ap, vec.
            comp_time  = res_arr.compute_time                                                                                           # Assign comp_time.
            res_main   = res_arr                                                                                                        # Assign res_main.

    # --- C. Using Pandas DataFrames ---
    if 'pandas' in input_types:                                                                                                         # Evaluate condition.
        bc_pd = pd.DataFrame(u_ex)                                                                                                      # Assign bc_pd.
        res_pd = TimeDerivative1(                                                                                                       # Assign res_pd.
            p,                                                                                                                          # Execute statement.
            ic=u_ex[:, 0],                                                                                                              # Assign ic.
            bc=bc_pd,                                                                                                                   # Assign bc.
            source=source,                                                                                                              # Assign source.
            t_span=t_span,                                                                                                              # Assign t_span.
            t=t_used,                                                                                                                   # Assign t.
            coef=[a, b, nu, F_react],                                                                                                   # Assign coef.
            operator=L,                                                                                                                 # Assign operator.
            upwind=upwind,                                                                                                              # Assign upwind.
            vec=vec0,                                                                                                                   # Assign vec.
            nvec=nvec,                                                                                                                  # Assign nvec.
            implicit=True,                                                                                                              # Assign implicit.
            lam=0.5,                                                                                                                    # Assign lam.
            device=device,                                                                                                              # Assign device.
            verbose=False                                                                                                               # Assign verbose.
        )                                                                                                                               # Execute statement.
        if u_ap is not None:                                                                                                            # Evaluate condition.
            assert np.allclose(u_ap, res_pd.solution), "Mismatch between Array/Callable and Pandas solver outputs."                     # Execute statement.
        else:                                                                                                                           # Execute fallback branch.
            u_ap, vec  = res_pd.solution, res_pd.neighbors                                                                              # Assign u_ap, vec.
            comp_time  = res_pd.compute_time                                                                                            # Assign comp_time.
            res_main   = res_pd                                                                                                         # Assign res_main.

    if u_ap is None:                                                                                                                    # Evaluate condition.
        raise ValueError("No valid input_types were specified.")                                                                        # Raise exception.

    # Calculate comprehensive transient error metrics
    metrics = Errors.Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=comp_time)                                              # Assign metrics.
    
    if res_main is not None:                                                                                                            # Evaluate condition.
        metrics['CFL']        = res_main.cfl                                                                                            # Assign metrics['CFL'].
        metrics['dt']         = res_main.dt                                                                                             # Assign metrics['dt'].
        metrics['Time_Steps'] = res_main.t_steps                                                                                        # Assign metrics['Time_Steps'].
    
    if 'array' in input_types and 'res_arr' in locals(): metrics['Time_Array']  = res_arr.compute_time                                  # Evaluate condition.
    if 'pandas' in input_types and 'res_pd' in locals(): metrics['Time_Pandas'] = res_pd.compute_time                                   # Evaluate condition.
    
    if verbose:                                                                                                                         # Evaluate condition.
        logger.info(f'\tError (Mean RMSE): {metrics["Time_Mean_RMSE"]}')                                                                # Log output message.

    if vec0 is None:                                                                                                                    # Evaluate condition.
        save_neighbors(cloud_path, nvec, vec)                                                                                           # Execute statement.

    save_metrics(out_dir, metrics, config_id=config_id, scale=scale, p=p)                                                               # Assign save_metrics(out_dir, metrics, config_id.
    if save:                                                                                                                            # Evaluate condition.
        cloud_name = os.path.basename(cloud_path).replace('.csv', '')                                                                   # Assign cloud_name.
        if scale == '3':                                                                                                                # Evaluate condition.
            if config_id.startswith('nvec_20_spsolve') or kwargs.get('plot_approximations', False):                                     # Evaluate condition.
                plot_transient(p, u_ap, save=True, nom=os.path.join(out_dir, f'Appx_{config_id}_{cloud_name}'),                         # Assign plot_transient(p, u_ap, save.
                                            title='Transient Approximation', verbose=verbose, t_span=t_span)                            # Assign title.
            exact_nom = os.path.join(out_dir, f'Exact_{cloud_name}')                                                                    # Assign exact_nom.
            if not (os.path.exists(exact_nom + '.mp4') or os.path.exists(exact_nom + '.gif')):                                          # Evaluate condition.
                plot_transient(p, u_ex, save=True, nom=exact_nom,                                                                       # Assign plot_transient(p, u_ex, save.
                                            title='Theoretical Solution', verbose=verbose, t_span=t_span)                               # Assign title.

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
