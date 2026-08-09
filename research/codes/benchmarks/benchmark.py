"""
benchmark — Runtime and memory benchmarking for mGFD reference problems

Overview:
    This script benchmarks the main solvers in mGFD on the available input point clouds under Data/.
    It measures:
    - Wall-clock execution time for each solve
    - Process memory usage before/after the solve (RSS)
    - Numerical error via Scripts/Errors.py

Reference problems:
    - Poisson (stationary)
    - Heat (first-order transient)
    - Advection–Diffusion (first-order transient)
    - Wave (second-order transient)

Usage:
    python batches/benchmark.py

Outputs:
    Creates a benchmark_results/ folder in the current working directory and writes:
    - JSON with all raw results and benchmark metadata
    - CSV with results table
    - TXT summary (aggregated per equation)

Dependencies:
    Required:
        numpy
    Optional (required to run this script fully):
        psutil   (process memory)
        pandas   (CSV table export and aggregation)
"""

## Library importation.
import sys                                                                                                                      # sys.path manipulation so this script can import project modules.
import os                                                                                                                       # Filesystem and path utilities.
import time                                                                                                                     # Wall-clock timing.
import json                                                                                                                     # JSON serialization for benchmark results.
import logging                                                                                                                  # Standard logging module.
from datetime import datetime                                                                                                   # Timestamping benchmark outputs.
from typing import Optional, List, Dict, Callable, Any, Tuple                                                                   # Type hinting.
import numpy as np                                                                                                              # Numerical arrays and math.

try:                                                                                                                            # Optional dependency: memory measurements.
    import psutil                                                                                                               # Process memory usage (RSS).
except ImportError:                                                                                                             # psutil not installed.
    psutil = None                                                                                                               # Mark dependency as unavailable.

try:                                                                                                                            # Optional dependency: CSV output + aggregation.
    import pandas as pd                                                                                                         # Tabular exports and aggregation.
except ImportError:                                                                                                             # pandas not installed.
    pd = None                                                                                                                   # Mark dependency as unavailable.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                     # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                                                       # Enable imports like "from mGFD import Stationary".

import utils.metrics as Errors                                                                                                  # Error metrics for stationary/transient runs.
from mGFD import Stationary, TimeDerivative1, TimeDerivative2                                                                   # Core solvers to benchmark.
from mGFD.io.io import load_points                                                                                              # Point cloud loading utility.
from utils.batch_utils import iter_clouds                                                                                       # Dataset loading and traversal.

logger = logging.getLogger(__name__)                                                                                            # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                   # Basic logger configuration.

def _require_optional_deps() -> None:                                                                                           # Validate optional dependencies for benchmarking/reporting.
    """
    _require_optional_deps
    Ensure optional dependencies required by the benchmark script are installed.

    Raises:
        ImportError: If psutil or pandas are not available.
    """
    if psutil is None:                                                                                                          # Memory measurement dependency check.
        raise ImportError('psutil is required for benchmark.py. Install with: pip install psutil')                              # Provide installation hint.
    if pd is None:                                                                                                              # Table export/aggregation dependency check.
        raise ImportError('pandas is required for benchmark.py. Install with: pip install pandas')                              # Provide installation hint.

def measure_performance(func: Callable, *args: Any, **kwargs: Any) -> Tuple[Any, float, float, float]:                          # Measure time and memory for a solver call.
    """
    measure_performance
    Measure runtime and memory usage for a function call.

    Input:
        func                        callable        Function to execute.
        *args                                       Positional arguments passed to func.
        **kwargs                                    Keyword arguments passed to func.

    Output:
        result                      any             Return value from func.
        execution_time              float           Wall-clock execution time (seconds).
        memory_used                 float           RSS memory delta (MB).
        peak_memory                 float           RSS memory after execution (MB).
    """
    _require_optional_deps()                                                                                                    # Ensure psutil/pandas are available for benchmarking.
    assert psutil is not None                                                                                                   # Satisfy static type checkers.

    process = psutil.Process()                                                                                                  # Handle for current process metrics.
    memory_before = process.memory_info().rss / 1024 / 1024                                                                     # Initial RSS memory (MB).
    
    start_time = time.time()                                                                                                    # Start wall-clock timer.
    result = func(*args, **kwargs)                                                                                              # Execute target function.
    end_time = time.time()                                                                                                      # Stop wall-clock timer.
    
    memory_after = process.memory_info().rss / 1024 / 1024                                                                      # Final RSS memory (MB).
    
    execution_time = end_time - start_time                                                                                      # Total execution time.
    memory_used = memory_after - memory_before                                                                                  # Delta memory during execution.
    peak_memory = memory_after                                                                                                  # Peak approximation as final RSS.
    
    return result, execution_time, memory_used, peak_memory                                                                     # Return performance measurements and result.

def benchmark_poisson_equation(p: np.ndarray) -> Dict[str, Any]:                                                                # Benchmark stationary Poisson problem on one point cloud.
    """
    benchmark_poisson_equation
    Benchmark the stationary Poisson-type problem on a given point cloud.

    The setup matches batches/run_Poisson.py.

    Input:
        p                           (m, 3) ndarray  Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    def phi(x: np.ndarray, y: np.ndarray) -> np.ndarray:                                                                        # Boundary condition.
        return 2 * np.exp(2 * x + y)

    def f(x: np.ndarray, y: np.ndarray) -> np.ndarray:                                                                          # RHS forcing term.
        return 10 * np.exp(2 * x + y)

    L = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                               # Operator (matches run_Poisson.py).
    
    (u_ap, vec), exec_time, memory_used, peak_memory = measure_performance(                                                     # Measure solver performance on this cloud.
        Stationary, p, phi, f, operator = L, verbose = False                                                                    # Call stationary solver.
    )                                                                                                                           # Unpack solution and measurements.
    
    u_ex = phi(p[:, 0], p[:, 1])                                                                                                # Compute exact theoretical solution.
    er = Errors.compute_rmse_stationary(p, vec, u_ap, u_ex)                                                                     # Per-node error metric.
    avg_error = float(np.mean(er))                                                                                              # Average numerical error for summary.
    
    return {                                                                                                                    # Return a structured metrics record for this run.
        'equation': 'Poisson',                                                                                                  # Equation label.
        'num_points': len(p),                                                                                                   # Number of nodes in the cloud.
        'boundary_condition': 'phi = 2*exp(2*x + y)',                                                                           # Boundary condition description.
        'rhs_function': 'f = 10*exp(2*x + y)',                                                                                  # RHS description.
        'operator': 'L = [0, 0, 2, 0, 2, 0]',                                                                                   # Operator coefficients string.
        'execution_time_seconds': exec_time,                                                                                    # Measured time.
        'memory_used_mb': memory_used,                                                                                          # Memory delta.
        'peak_memory_mb': peak_memory,                                                                                          # Peak approximation.
        'avg_numerical_error': avg_error                                                                                        # Mean error.
    }                                                                                                                           # Return metrics dict.

def benchmark_heat_equation(p: np.ndarray) -> Dict[str, Any]:                                                                   # Benchmark transient Heat problem on one point cloud.
    """
    benchmark_heat_equation
    Benchmark the first-order transient Heat problem on a given point cloud.

    The setup matches batches/run_Heat.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray  Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    v: float = 0.2                                                                                                              # Diffusion coefficient.
    t: int = 2000                                                                                                               # Number of time steps.
    
    def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                         # Exact solution / forcing generator.
        return np.exp(-2 * np.pi**2 * coef[0] * t_val) * np.cos(np.pi * x) * np.cos(np.pi * y)
        
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                       # Operator (matches run_Heat.py).
    
    (u_ap, vec), exec_time, memory_used, peak_memory = measure_performance(                                                     # Measure solver performance on this cloud.
        TimeDerivative1, p, f, t, [v], operator = L, implicit = True, lam = 0.5, verbose = False                                # Call transient solver (1st order).
    )                                                                                                                           # Unpack solution and measurements.
    
    T_arr = np.linspace(0, 1, t)                                                                                                # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                                                # Initialize exact solution matrix.
    for k in range(t):                                                                                                          # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T_arr[k], [v])                                                                         # Compute exact theoretical solution.
    er = Errors.compute_rmse_transient(p, vec, u_ap, u_ex)                                                                      # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                                              # Average numerical error for summary.
    
    return {                                                                                                                    # Return a structured metrics record for this run.
        'equation': 'Heat',                                                                                                     # Equation label.
        'num_points': len(p),                                                                                                   # Number of nodes in the cloud.
        'diffusion_coefficient': v,                                                                                             # Diffusion parameter.
        'time_steps': t,                                                                                                        # Time-step count.
        'function': 'f = exp(-2*pi^2*v*t)*cos(pi*x)*cos(pi*y)',                                                                 # Function description.
        'operator': f'L = [0, 0, {2 * v}, 0, {2 * v}, 0]',                                                                      # Operator coefficients string.
        'implicit': False,                                                                                                      # Integration mode flag.
        'lambda': 0.5,                                                                                                          # Stabilization parameter used by the solver.
        'execution_time_seconds': exec_time,                                                                                    # Measured time.
        'memory_used_mb': memory_used,                                                                                          # Memory delta.
        'peak_memory_mb': peak_memory,                                                                                          # Peak approximation.
        'time_per_step_ms': (exec_time * 1000) / t,                                                                             # Time per step (ms).
        'avg_numerical_error': avg_error                                                                                        # Mean error.
    }                                                                                                                           # Return metrics dict.

def benchmark_advdif_equation(p: np.ndarray) -> Dict[str, Any]:                                                                 # Benchmark transient Advection–Diffusion problem on one point cloud.
    """
    benchmark_advdif_equation
    Benchmark the first-order transient Advection–Diffusion problem on a given point cloud.

    The setup matches batches/run_AdvDif.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray  Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    v: float = 0.1                                                                                                              # Diffusion coefficient.
    a: float = 0.3                                                                                                              # Transport velocity in x.
    b: float = 0.2                                                                                                              # Transport velocity in y.
    t: int = 2000                                                                                                               # Number of time steps.
    
    def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                         # Exact solution / forcing generator.
        return (1 / (4 * t_val + 1)) * np.exp(
            - (x - coef[1] * t_val - 0.5)**2 / (coef[0] * (4 * t_val + 1)) - (y - coef[2] * t_val - 0.5)**2 / (coef[0] * (4 * t_val + 1))
        )
        
    L = np.vstack([[-a], [-b], [2 * v], [0], [2 * v], [0]])                                                                     # Operator (matches run_AdvDif.py).
    
    (u_ap, vec), exec_time, memory_used, peak_memory = measure_performance(                                                     # Measure solver performance on this cloud.
        TimeDerivative1, p, f, t, [v, a, b], operator = L, implicit = True, lam = 0.5, upwind = True, verbose = False           # Call transient solver (1st order).
    )                                                                                                                           # Unpack solution and measurements.
    
    T_arr = np.linspace(0, 1, t)                                                                                                # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                                                # Initialize exact solution matrix.
    for k in range(t):                                                                                                          # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T_arr[k], [v, a, b])                                                                   # Compute exact theoretical solution.
    er = Errors.compute_rmse_transient(p, vec, u_ap, u_ex)                                                                      # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                                              # Average numerical error for summary.
    
    return {                                                                                                                    # Return a structured metrics record for this run.
        'equation': 'Advection-Diffusion',                                                                                      # Equation label.
        'num_points': len(p),                                                                                                   # Number of nodes in the cloud.
        'diffusion_coefficient': v,                                                                                             # Diffusion parameter.
        'transport_velocity_x': a,                                                                                              # Transport velocity x.
        'transport_velocity_y': b,                                                                                              # Transport velocity y.
        'time_steps': t,                                                                                                        # Time-step count.
        'function': 'f = (1/(4*t+1))*exp(-(x-a*t-0.5)^2/(v*(4*t+1)) - (y-b*t-0.5)^2/(v*(4*t+1)))',                              # Function description.
        'operator': f'L = [{-a}, {-b}, {2 * v}, 0, {2 * v}, 0]',                                                                # Operator coefficients string.
        'implicit': True,                                                                                                       # Integration mode flag.
        'lambda': 0.5,                                                                                                          # Stabilization parameter used by the solver.
        'execution_time_seconds': exec_time,                                                                                    # Measured time.
        'memory_used_mb': memory_used,                                                                                          # Memory delta.
        'peak_memory_mb': peak_memory,                                                                                          # Peak approximation.
        'time_per_step_ms': (exec_time * 1000) / t,                                                                             # Time per step (ms).
        'avg_numerical_error': avg_error                                                                                        # Mean error.
    }                                                                                                                           # Return metrics dict.

def benchmark_advection_equation(p: np.ndarray) -> Dict[str, Any]:                                                              # Benchmark the first-order transient Pure Advection problem on a given point cloud.
    """
    benchmark_advection_equation
    Benchmark the first-order transient Pure Advection problem on a given point cloud.
    
    Input:
        p                           (m, 3) ndarray  Point cloud [x, y, flag].
        
    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    a: float = 0.3                                                                                                              # Transport velocity in x.
    b: float = 0.2                                                                                                              # Transport velocity in y.
    t: int = 2000                                                                                                               # Number of time steps.
    
    def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                         # Exact solution / forcing generator.
        return np.exp(-((x - coef[0] * t_val - 0.5)**2) / 0.05 - ((y - coef[1] * t_val - 0.5)**2) / 0.05)
        
    L = np.vstack([[-a], [-b], [0], [0], [0], [0]])                                                                             # Operator.
    
    (u_ap, vec), exec_time, memory_used, peak_memory = measure_performance(                                                     # Measure solver performance on this cloud.
        TimeDerivative1, p, f, t, [a, b], operator = L, implicit = True, lam = 1.0, upwind = True, verbose = False              # Call transient solver (1st order).
    )                                                                                                                           # Unpack solution and measurements.
    
    T_arr = np.linspace(0, 1, t)                                                                                                # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                                                # Initialize exact solution matrix.
    for k in range(t):                                                                                                          # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T_arr[k], [a, b])                                                                      # Compute exact theoretical solution.
    er = Errors.compute_rmse_transient(p, vec, u_ap, u_ex)                                                                      # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                                              # Average numerical error for summary.
    
    return {                                                                                                                    # Return a structured metrics record for this run.
        'equation': 'Advection',                                                                                                # Equation label.
        'num_points': len(p),                                                                                                   # Number of nodes in the cloud.
        'transport_velocity_x': a,                                                                                              # Transport velocity x.
        'transport_velocity_y': b,                                                                                              # Transport velocity y.
        'time_steps': t,                                                                                                        # Time-step count.
        'function': 'f = exp(-((x - a*t - 0.5)^2)/0.05 - ((y - b*t - 0.5)^2)/0.05)',                                            # Function description.
        'operator': f'L = [{-a}, {-b}, 0, 0, 0, 0]',                                                                            # Operator coefficients string.
        'implicit': True,                                                                                                       # Integration mode flag.
        'lambda': 1.0,                                                                                                          # Stabilization parameter used by the solver.
        'execution_time_seconds': exec_time,                                                                                    # Measured time.
        'memory_used_mb': memory_used,                                                                                          # Memory delta.
        'peak_memory_mb': peak_memory,                                                                                          # Peak approximation.
        'time_per_step_ms': (exec_time * 1000) / t,                                                                             # Time per step (ms).
        'avg_numerical_error': avg_error                                                                                        # Mean error.
    }                                                                                                                           # Return metrics dict.

def benchmark_wave_equation(p: np.ndarray) -> Dict[str, Any]:                                                                   # Benchmark transient Wave problem on one point cloud.
    """
    benchmark_wave_equation
    Benchmark the second-order transient Wave problem on a given point cloud.

    The setup matches batches/run_Wave.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray  Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    c: float = float(np.sqrt(1 / 2))                                                                                            # Wave propagation coefficient.
    t: int = 2000                                                                                                               # Number of time steps.
    
    def f(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                         # Initial displacement / forcing generator.
        return np.cos(np.pi * t_val) * np.sin(np.pi * (x + y))
        
    def g(x: np.ndarray, y: np.ndarray, t_val: float, coef: List[float]) -> np.ndarray:                                         # Initial velocity / forcing generator.
        return -np.pi * np.sin(np.pi * t_val) * np.sin(np.pi * (x + y))
        
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                                                 # Operator (matches run_Wave.py).
    
    (u_ap, vec), exec_time, memory_used, peak_memory = measure_performance(                                                     # Measure solver performance on this cloud.
        TimeDerivative2, p, f, g, t, [c], operator = L, implicit = True, lam = 1, verbose = False                               # Call transient solver (2nd order).
    )                                                                                                                           # Unpack solution and measurements.
    
    T_arr = np.linspace(0, 1, t)                                                                                                # Reconstruct time vector.
    u_ex = np.zeros([len(p), t])                                                                                                # Initialize exact solution matrix.
    for k in range(t):                                                                                                          # Loop over all time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T_arr[k], [c])                                                                         # Compute exact theoretical solution.
    er = Errors.compute_rmse_transient(p, vec, u_ap, u_ex)                                                                      # Per-node error metric across time.
    avg_error = np.mean(er)                                                                                                     # Average numerical error for summary.
    
    return {                                                                                                                    # Return a structured metrics record for this run.
        'equation': 'Wave',                                                                                                     # Equation label.
        'num_points': len(p),                                                                                                   # Number of nodes in the cloud.
        'wave_coefficient': c,                                                                                                  # Wave parameter.
        'time_steps': t,                                                                                                        # Time-step count.
        'function_f': 'f = cos(pi*t)*sin(pi*(x+y))',                                                                            # Function f description.
        'function_g': 'g = -pi*sin(pi*t)*sin(pi*(x+y))',                                                                        # Function g description.
        'operator': f'L = [0, 0, {2 * c**2}, 0, {2 * c**2}, 0]',                                                                # Operator coefficients string.
        'implicit': False,                                                                                                      # Integration mode flag.
        'lambda': 1,                                                                                                            # Stabilization parameter used by the solver.
        'execution_time_seconds': exec_time,                                                                                    # Measured time.
        'memory_used_mb': memory_used,                                                                                          # Memory delta.
        'peak_memory_mb': peak_memory,                                                                                          # Peak approximation.
        'time_per_step_ms': (exec_time * 1000) / t,                                                                             # Time per step (ms).
        'avg_numerical_error': avg_error                                                                                        # Mean error.
    }                                                                                                                           # Return metrics dict.

def run_comprehensive_benchmark(equations: Optional[List[str]] = None, save_results: bool = True, data_root: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    run_comprehensive_benchmark
    Run the full benchmark suite across all point clouds and selected equations.

    Input:
        equations                    list[str]|None  Which equations to benchmark. Defaults to all.
        save_results                 bool            Whether to write JSON/CSV/TXT outputs to disk.
        data_root                    str|None        Root folder containing Data/; defaults to <repo>/Data.

    Output:
        results                      list[dict]      List with one entry per (cloud, equation) run.
    """
    _require_optional_deps()                                                                                                    # Ensure psutil/pandas are available for this workflow.

    if equations is None:                                                                                                       # Default equation set when not provided.
        equations = ['Poisson', 'Heat', 'Adv', 'AdvDif', 'Wave']                                                                # Default: run all available benchmark cases.
    
    equation_functions = {                                                                                                      # Map equation labels to benchmark functions.
        'Poisson': benchmark_poisson_equation,                                                                                  # Stationary Poisson benchmark.
        'Heat': benchmark_heat_equation,                                                                                        # Transient heat benchmark.
        'Adv': benchmark_advection_equation,                                                                                    # Transient advection benchmark.
        'AdvDif': benchmark_advdif_equation,                                                                                    # Transient advection–diffusion benchmark.
        'Wave': benchmark_wave_equation                                                                                         # Transient wave benchmark.
    }                                                                                                                           # Mapping from label to benchmark function.
    
    results: List[Dict[str, Any]] = []                                                                                          # Accumulator for per-run results.
    
    logger.info("Iniciando benchmark completo...")                                                                              # Console header.
    logger.info(f"Ecuaciones: {equations}")                                                                                     # Report selected equation set.
    
    if data_root is None:                                                                                                       # Default data root when not provided.
        data_root = os.path.join(os.path.dirname(BASE_DIR), 'data')                                                                              # Default to <repo_root>/Data.
    SCALES = ('1', '2', '3', '4', '5')
    clouds = list(iter_clouds(data_root, scales=SCALES))                                                                        # Enumerate all cloud CSV files.
    total_combinations = len(clouds) * len(equations)                                                                           # Total runs expected.
    current = 0                                                                                                                 # Progress counter.
    
    benchmark_start_time = time.time()                                                                                          # Start total benchmark timer.
    
    for dataset, scale, cloud_path in clouds:                                                                                   # Iterate all discovered cloud CSVs.
        p = load_points(cloud_path, verbose=False)                                                                              # Load point cloud (m, 3).
        region_id = f"{dataset}/{scale}"                                                                                        # Region identifier.
        for equation in equations:                                                                                              # Run selected equation benchmarks for each cloud.
            current += 1                                                                                                        # Increment progress.
            logger.info(f"\nProcesando {current}/{total_combinations}: {region_id} - {equation}")                               # Report progress line.

            try:                                                                                                                # Protect the benchmark loop from single-run failures.
                benchmark_func = equation_functions[equation]                                                                   # Resolve benchmark function.
                result = benchmark_func(p)                                                                                      # Execute benchmark on this cloud.
                if result:                                                                                                      # Only record successful results.
                    result['region_id'] = region_id                                                                             # Attach region identifier to output record.
                    result['cloud_file'] = os.path.basename(cloud_path)                                                         # Attach source CSV name.
                    results.append(result)                                                                                      # Store record.
                    logger.info(f"  Error promedio: {result['avg_numerical_error']:.2e}")                                       # Print error summary.
                    logger.info(f"  Tiempo: {result['execution_time_seconds']:.3f}s")                                           # Print time summary.
                    logger.info(f"  Memoria: {result['memory_used_mb']:.1f}MB (pico: {result['peak_memory_mb']:.1f}MB)")        # Print memory summary.
                    if 'time_per_step_ms' in result:                                                                            # Per-step timing exists only for transient cases.
                        logger.info(f"  Tiempo por paso: {result['time_per_step_ms']:.3f}ms")                                   # Print per-step time when applicable.
            except Exception as e:                                                                                              # Catch solver or post-processing failures.
                logger.error(f"  Error en {equation}: {e}")                                                                     # Report failure for this run.
    
    benchmark_end_time = time.time()                                                                                            # Stop total benchmark timer.
    total_benchmark_time = benchmark_end_time - benchmark_start_time                                                            # Total wall time.
    
    if save_results and results:                                                                                                # Write output artifacts when requested and results exist.
        output_dir = os.path.join(os.path.dirname(BASE_DIR), 'results', 'benchmarks')                                           # Define output directory.
        os.makedirs(output_dir, exist_ok = True)                                                                                # Create output directory if needed.
        
        equation_names = {                                                                                                      # Display names used in summary grouping.
            'Poisson': 'Poisson',                                                                                               # Display name for Poisson.
            'Heat': 'Heat',                                                                                                     # Display name for Heat.
            'Adv': 'Advection',                                                                                                 # Display name for Adv.
            'AdvDif': 'Advection-Diffusion',                                                                                    # Display name for AdvDif.
            'Wave': 'Wave'                                                                                                      # Display name for Wave.
        }                                                                                                                       # Used to group results in summary.

        benchmark_stats = {                                                                                                     # Global benchmark metadata saved alongside raw results.
            'total_benchmark_time_seconds': total_benchmark_time,                                                               # Total elapsed time.
            'total_tests': len(results),                                                                                        # Total collected results.
            'successful_tests': len([r for r in results if r is not None]),                                                     # Count of non-null records.
            'timestamp': datetime.now().isoformat(),                                                                            # ISO timestamp.
            'data_root': data_root,                                                                                             # Data root used for discovery.
            'equations_tested': list(equations)                                                                                 # Equations included in this run.
        }                                                                                                                       # Metadata stored alongside results.
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")                                                                    # Timestamp used in output filenames.
        json_file = os.path.join(output_dir, f'benchmark_{timestamp}.json')                                                     # JSON output path.
        
        output_data = {                                                                                                         # JSON payload containing metadata + raw results.
            'benchmark_stats': benchmark_stats,                                                                                 # Global benchmark metadata.
            'results': results                                                                                                  # Raw per-run results.
        }                                                                                                                       # JSON payload.
        
        with open(json_file, 'w') as f_out:                                                                                     # Write JSON output file.
            json.dump(output_data, f_out, indent = 2)                                                                           # Write pretty-printed JSON.
        
        csv_file = os.path.join(output_dir, f'benchmark_{timestamp}.csv')                                                       # CSV output path.
        assert pd is not None                                                                                                   # Satisfy static type checkers.
        df = pd.DataFrame(results)                                                                                              # Convert to DataFrame.
        df.to_csv(csv_file, index = False)                                                                                      # Write CSV table.
        
        summary_file = os.path.join(output_dir, f'benchmark_summary_{timestamp}.txt')                                           # TXT summary output path.
        with open(summary_file, 'w') as f_out:                                                                                  # Write summary file.
            f_out.write("RESUMEN DEL BENCHMARK\n")                                                                              # Summary header.
            f_out.write("=====================\n\n")                                                                            # Separator.
            f_out.write(f"Tiempo total del benchmark: {total_benchmark_time:.2f} segundos\n")                                   # Total benchmark wall time.
            f_out.write(f"Pruebas exitosas: {len(results)} de {total_combinations}\n\n")                                        # Success count for expected runs.
            
            if results:                                                                                                         # Only aggregate statistics when results exist.
                df = pd.DataFrame(results)                                                                                      # Build DataFrame for aggregation.
                f_out.write("ESTADÍSTICAS POR ECUACIÓN:\n")                                                                     # Group header.
                for eq in equations:                                                                                            # Aggregate by equation label.
                    eq_name = equation_names.get(eq, eq)                                                                        # Display name mapping.
                    eq_results = df[df['equation'] == eq_name]                                                                  # Filter results by equation.
                    if not eq_results.empty:                                                                                    # Only write section when data is available.
                        f_out.write(f"\n{eq_name}:\n")                                                                          # Equation section header.
                        f_out.write(f"  Tiempo promedio: {eq_results['execution_time_seconds'].mean():.3f}s\n")                 # Mean time across clouds.
                        f_out.write(f"  Memoria promedio: {eq_results['memory_used_mb'].mean():.1f}MB\n")                       # Mean memory delta across clouds.
                        f_out.write(f"  Error promedio: {eq_results['avg_numerical_error'].mean():.2e}\n")                      # Mean error across clouds.
        
        logger.info("\nResultados guardados en:")                                                                               # Output location header.
        logger.info(f"  JSON: {json_file}")                                                                                     # JSON path.
        logger.info(f"  CSV: {csv_file}")                                                                                       # CSV path.
        logger.info(f"  Resumen: {summary_file}")                                                                               # TXT summary path.
    
    logger.info(f"\nBenchmark completado en {total_benchmark_time:.2f} segundos")                                               # Total time report.
    logger.info(f"Total de resultados: {len(results)}")                                                                         # Results count report.
    return results                                                                                                              # Return results list.

def main() -> None:
    """
    main
    Entry point for the benchmark script.
    """
    try:                                                                                                                        # Validate optional dependencies before running.
        _require_optional_deps()                                                                                                # Ensure psutil and pandas are available.
    except ImportError as e:                                                                                                    # Missing dependencies.
        logger.error(f"Error: {e}")                                                                                             # Print actionable error message.
        raise                                                                                                                   # Propagate error for non-zero exit.

    _ = run_comprehensive_benchmark()                                                                                           # Execute benchmark with default settings.

if __name__ == "__main__":                                                                                                      # Script entry point.
    main()
