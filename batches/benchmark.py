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
import sys                                                                                              # sys.path manipulation so this script can import project modules.
import os                                                                                               # Filesystem and path utilities.
import time                                                                                             # Wall-clock timing.
import json                                                                                             # JSON serialization for benchmark results.
from datetime import datetime                                                                           # Timestamping benchmark outputs.
import numpy as np                                                                                      # Numerical arrays and math.

try:                                                                                                    # Optional dependency: memory measurements.
    import psutil                                                                                       # Process memory usage (RSS).
except ImportError:                                                                                     # psutil not installed.
    psutil = None                                                                                       # Mark dependency as unavailable.

try:                                                                                                    # Optional dependency: CSV output + aggregation.
    import pandas as pd                                                                                 # Tabular exports and aggregation.
except ImportError:                                                                                     # pandas not installed.
    pd = None                                                                                           # Mark dependency as unavailable.

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Repository root from batches/ folder.
sys.path.append(BASE_DIR)                                                                               # Enable imports like "from mGFD import Stationary".

import Scripts.Errors as Errors                                                                         # Error metrics for stationary/transient runs.
from mGFD import Stationary, TimeDerivative1, TimeDerivative2                                           # Core solvers to benchmark.
from Scripts.IO import load_points, iter_clouds                                                         # Dataset loading and traversal.


def _require_optional_deps():                                                                           # Validate optional dependencies for benchmarking/reporting.
    """
    _require_optional_deps
    Ensure optional dependencies required by the benchmark script are installed.

    Raises:
        ImportError: If psutil or pandas are not available.
    """
    if psutil is None:                                                                                  # Memory measurement dependency check.
        raise ImportError('psutil is required for benchmark.py. Install with: pip install psutil')      # Provide installation hint.
    if pd is None:                                                                                      # Table export/aggregation dependency check.
        raise ImportError('pandas is required for benchmark.py. Install with: pip install pandas')      # Provide installation hint.

def measure_performance(func, *args, **kwargs):                                                         # Measure time and memory for a solver call.
    """
    measure_performance
    Measure runtime and memory usage for a function call.

    Input:
        func                        callable        Function to execute.
        *args                                     Positional arguments passed to func.
        **kwargs                                  Keyword arguments passed to func.

    Output:
        result                      any             Return value from func.
        execution_time              float           Wall-clock execution time (seconds).
        memory_used                 float           RSS memory delta (MB).
        peak_memory                 float           RSS memory after execution (MB).
    """
    _require_optional_deps()                                                                            # Ensure psutil/pandas are available for benchmarking.

    process = psutil.Process()                                                                          # Handle for current process metrics.
    memory_before = process.memory_info().rss / 1024 / 1024                                             # Initial RSS memory (MB).
    
    start_time = time.time()                                                                            # Start wall-clock timer.
    result = func(*args, **kwargs)                                                                      # Execute target function.
    end_time = time.time()                                                                              # Stop wall-clock timer.
    
    memory_after = process.memory_info().rss / 1024 / 1024                                              # Final RSS memory (MB).
    
    execution_time = end_time - start_time                                                              # Total execution time.
    memory_used = memory_after - memory_before                                                          # Delta memory during execution.
    peak_memory = memory_after                                                                          # Peak approximation as final RSS.
    
    return result, execution_time, memory_used, peak_memory                                             # Return performance measurements and result.

def benchmark_poisson_equation(p):                                                                      # Benchmark stationary Poisson problem on one point cloud.
    """
    benchmark_poisson_equation
    Benchmark the stationary Poisson-type problem on a given point cloud.

    The setup matches batches/run_Poisson.py.

    Input:
        p                           (m, 3) ndarray Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    
    phi = lambda x, y: 2 * np.exp(2 * x + y)                                                            # Boundary condition.
    f   = lambda x, y: 10 * np.exp(2 * x + y)                                                           # RHS forcing term.
    L   = np.vstack([[0], [0], [2], [0], [2], [0]])                                                     # Operator (matches run_Poisson.py).
    
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(                       # Measure solver performance on this cloud.
        Stationary, p, phi, f, operator = L                                                             # Call stationary solver.
    )                                                                                                   # Unpack solution and measurements.
    
    er = Errors.Cloud_Stationary(p, vec, u_ap, u_ex)                                                    # Per-node error metric.
    avg_error = float(np.mean(er))                                                                      # Average numerical error for summary.
    
    return {                                                                                            # Return a structured metrics record for this run.
        'equation': 'Poisson',                                                                          # Equation label.
        'num_points': int(len(p)),                                                                      # Number of nodes in the cloud.
        'boundary_condition': 'phi = 2*exp(2*x + y)',                                                   # Boundary condition description.
        'rhs_function': 'f = 10*exp(2*x + y)',                                                          # RHS description.
        'operator': 'L = [0, 0, 2, 0, 2, 0]',                                                           # Operator coefficients string.
        'execution_time_seconds': float(exec_time),                                                     # Measured time.
        'memory_used_mb': float(memory_used),                                                           # Memory delta.
        'peak_memory_mb': float(peak_memory),                                                           # Peak approximation.
        'avg_numerical_error': float(avg_error)                                                         # Mean error.
    }                                                                                                   # Return metrics dict.

def benchmark_heat_equation(p):                                                                         # Benchmark transient Heat problem on one point cloud.
    """
    benchmark_heat_equation
    Benchmark the first-order transient Heat problem on a given point cloud.

    The setup matches batches/run_Heat.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    
    v = 0.2                                                                                             # Diffusion coefficient.
    t = 2000                                                                                            # Number of time steps.
    f = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)
                                                                                                        # Exact solution / forcing generator.
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                               # Operator (matches run_Heat.py).
    
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(                       # Measure solver performance on this cloud.
        TimeDerivative1, p, f, t, [v], operator = L, implicit = False, lam = 0.5                        # Call transient solver (1st order).
    )                                                                                                   # Unpack solution and measurements.
    
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)                                                     # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                      # Average numerical error for summary.
    
    return {                                                                                            # Return a structured metrics record for this run.
        'equation': 'Heat',                                                                             # Equation label.
        'num_points': int(len(p)),                                                                      # Number of nodes in the cloud.
        'diffusion_coefficient': float(v),                                                              # Diffusion parameter.
        'time_steps': int(t),                                                                           # Time-step count.
        'function': 'f = exp(-2*pi^2*v*t)*cos(pi*x)*cos(pi*y)',                                         # Function description.
        'operator': f'L = [0, 0, {2 * v}, 0, {2 * v}, 0]',                                              # Operator coefficients string.
        'implicit': False,                                                                              # Integration mode flag.
        'lambda': 0.5,                                                                                  # Stabilization parameter used by the solver.
        'execution_time_seconds': float(exec_time),                                                     # Measured time.
        'memory_used_mb': float(memory_used),                                                           # Memory delta.
        'peak_memory_mb': float(peak_memory),                                                           # Peak approximation.
        'time_per_step_ms': float((exec_time * 1000) / t),                                              # Time per step (ms).
        'avg_numerical_error': float(avg_error)                                                         # Mean error.
    }                                                                                                   # Return metrics dict.

def benchmark_advdif_equation(p):                                                                       # Benchmark transient Advection–Diffusion problem on one point cloud.
    """
    benchmark_advdif_equation
    Benchmark the first-order transient Advection–Diffusion problem on a given point cloud.

    The setup matches batches/run_AdvDif.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    
    v = 0.1                                                                                             # Diffusion coefficient.
    a = 0.3                                                                                             # Transport velocity in x.
    b = 0.2                                                                                             # Transport velocity in y.
    t = 2000                                                                                            # Number of time steps.
    f = lambda x, y, t, coef: (1 / (4 * t + 1)) * np.exp(                                               # Exact solution / forcing generator.
        - (x - coef[1] * t - 0.5)**2 / (coef[0] * (4 * t + 1)) - (y - coef[2] * t - 0.5)**2 / (coef[0] * (4 * t + 1))
                                                                                                        # Exponent of the advected Gaussian.
    )                                                                                                   # Function matches run_AdvDif.py.
    L = np.vstack([[-a], [-b], [2 * v], [0], [2 * v], [0]])                                             # Operator (matches run_AdvDif.py).
    
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(                       # Measure solver performance on this cloud.
        TimeDerivative1, p, f, t, [v, a, b], operator = L, implicit = False, lam = 0.5                  # Call transient solver (1st order).
    )                                                                                                   # Unpack solution and measurements.
    
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)                                                     # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                      # Average numerical error for summary.
    
    return {                                                                                            # Return a structured metrics record for this run.
        'equation': 'Advection-Diffusion',                                                              # Equation label.
        'num_points': int(len(p)),                                                                      # Number of nodes in the cloud.
        'diffusion_coefficient': float(v),                                                              # Diffusion parameter.
        'transport_velocity_x': float(a),                                                               # Transport velocity x.
        'transport_velocity_y': float(b),                                                               # Transport velocity y.
        'time_steps': int(t),                                                                           # Time-step count.
        'function': 'f = (1/(4*t+1))*exp(-(x-a*t-0.5)^2/(v*(4*t+1)) - (y-b*t-0.5)^2/(v*(4*t+1)))',      # Function description.
        'operator': f'L = [{-a}, {-b}, {2 * v}, 0, {2 * v}, 0]',                                        # Operator coefficients string.
        'implicit': False,                                                                              # Integration mode flag.
        'lambda': 0.5,                                                                                  # Stabilization parameter used by the solver.
        'execution_time_seconds': float(exec_time),                                                     # Measured time.
        'memory_used_mb': float(memory_used),                                                           # Memory delta.
        'peak_memory_mb': float(peak_memory),                                                           # Peak approximation.
        'time_per_step_ms': float((exec_time * 1000) / t),                                              # Time per step (ms).
        'avg_numerical_error': float(avg_error)                                                         # Mean error.
    }                                                                                                   # Return metrics dict.

def benchmark_wave_equation(p):                                                                         # Benchmark transient Wave problem on one point cloud.
    """
    benchmark_wave_equation
    Benchmark the second-order transient Wave problem on a given point cloud.

    The setup matches batches/run_Wave.py (explicit time integration).

    Input:
        p                           (m, 3) ndarray Point cloud [x, y, flag].

    Output:
        metrics                     dict            Performance and accuracy metrics for this run.
    """
    
    c = float(np.sqrt(1 / 2))                                                                           # Wave propagation coefficient.
    t = 2000                                                                                            # Number of time steps.
    f = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                               # Initial displacement / forcing generator.
    g = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))                      # Initial velocity / forcing generator.
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                         # Operator (matches run_Wave.py).
    
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(                       # Measure solver performance on this cloud.
        TimeDerivative2, p, f, g, t, [c], operator = L, implicit = False, lam = 1                       # Call transient solver (2nd order).
    )                                                                                                   # Unpack solution and measurements.
    
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)                                                     # Per-node error metric across time.
    avg_error = float(np.mean(er))                                                                      # Average numerical error for summary.
    
    return {                                                                                            # Return a structured metrics record for this run.
        'equation': 'Wave',                                                                             # Equation label.
        'num_points': int(len(p)),                                                                      # Number of nodes in the cloud.
        'wave_coefficient': float(c),                                                                   # Wave parameter.
        'time_steps': int(t),                                                                           # Time-step count.
        'function_f': 'f = cos(pi*t)*sin(pi*(x+y))',                                                    # Function f description.
        'function_g': 'g = -pi*sin(pi*t)*sin(pi*(x+y))',                                                # Function g description.
        'operator': f'L = [0, 0, {2 * c**2}, 0, {2 * c**2}, 0]',                                        # Operator coefficients string.
        'implicit': False,                                                                              # Integration mode flag.
        'lambda': 1,                                                                                    # Stabilization parameter used by the solver.
        'execution_time_seconds': float(exec_time),                                                     # Measured time.
        'memory_used_mb': float(memory_used),                                                           # Memory delta.
        'peak_memory_mb': float(peak_memory),                                                           # Peak approximation.
        'time_per_step_ms': float((exec_time * 1000) / t),                                              # Time per step (ms).
        'avg_numerical_error': float(avg_error)                                                         # Mean error.
    }                                                                                                   # Return metrics dict.

def run_comprehensive_benchmark(equations=None, save_results=True, data_root=None):                     # Run the benchmark across all clouds and equations.
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
    _require_optional_deps()                                                                            # Ensure psutil/pandas are available for this workflow.

    if equations is None:                                                                               # Default equation set when not provided.
        equations = ['Poisson', 'Heat', 'AdvDif', 'Wave']                                               # Default: run all available benchmark cases.
    
    equation_functions = {                                                                              # Map equation labels to benchmark functions.
        'Poisson': benchmark_poisson_equation,                                                          # Stationary Poisson benchmark.
        'Heat': benchmark_heat_equation,                                                                # Transient heat benchmark.
        'AdvDif': benchmark_advdif_equation,                                                            # Transient advection–diffusion benchmark.
        'Wave': benchmark_wave_equation                                                                 # Transient wave benchmark.
    }                                                                                                   # Mapping from label to benchmark function.
    
    results = []                                                                                        # Accumulator for per-run results.
    
    print("Iniciando benchmark completo...")                                                            # Console header.
    print(f"Ecuaciones: {equations}")                                                                   # Report selected equation set.
    
    if data_root is None:                                                                               # Default data root when not provided.
        data_root = os.path.join(BASE_DIR, 'Data')                                                      # Default to <repo_root>/Data.
    clouds = list(iter_clouds(data_root))                                                                  # Enumerate all cloud CSV files.
    total_combinations = int(len(clouds) * len(equations))                                              # Total runs expected.
    current = 0                                                                                         # Progress counter.
    
    benchmark_start_time = time.time()                                                                  # Start total benchmark timer.
    
    for dataset, scale, variant, cloud_path in clouds:                                                  # Iterate all discovered cloud CSVs.
        p = load_points(cloud_path)                                                                     # Load point cloud (m, 3).
        region_id = f"{dataset}/{scale}/{variant}"                                                      # Human-readable region identifier.
        for equation in equations:                                                                      # Run selected equation benchmarks for each cloud.
            current += 1                                                                                # Increment progress.
            print(f"\nProcesando {current}/{total_combinations}: {region_id} - {equation}")             # Report progress line.

            try:                                                                                        # Protect the benchmark loop from single-run failures.
                benchmark_func = equation_functions[equation]                                           # Resolve benchmark function.
                result = benchmark_func(p)                                                              # Execute benchmark on this cloud.
                if result:                                                                              # Only record successful results.
                    result['region_id'] = region_id                                                     # Attach region identifier to output record.
                    result['cloud_file'] = os.path.basename(cloud_path)                                 # Attach source CSV name.
                    results.append(result)                                                              # Store record.
                    print(f"  Error promedio: {result['avg_numerical_error']:.2e}")                     # Print error summary.
                    print(f"  Tiempo: {result['execution_time_seconds']:.3f}s")                         # Print time summary.
                    print(f"  Memoria: {result['memory_used_mb']:.1f}MB (pico: {result['peak_memory_mb']:.1f}MB)")
                                                                                                        # Print memory summary.
                    if 'time_per_step_ms' in result:                                                    # Per-step timing exists only for transient cases.
                        print(f"  Tiempo por paso: {result['time_per_step_ms']:.3f}ms")                 # Print per-step time when applicable.
            except Exception as e:                                                                      # Catch solver or post-processing failures.
                print(f"  Error en {equation}: {e}")                                                    # Report failure for this run.
    
    benchmark_end_time = time.time()                                                                    # Stop total benchmark timer.
    total_benchmark_time = float(benchmark_end_time - benchmark_start_time)                             # Total wall time.
    
    if save_results and results:                                                                        # Write output artifacts when requested and results exist.
        os.makedirs('benchmark_results', exist_ok = True)                                               # Create output directory if needed.
        
        equation_names = {                                                                              # Display names used in summary grouping.
            'Poisson': 'Poisson',                                                                       # Display name for Poisson.
            'Heat': 'Heat',                                                                             # Display name for Heat.
            'AdvDif': 'Advection-Diffusion',                                                            # Display name for AdvDif.
            'Wave': 'Wave'                                                                              # Display name for Wave.
        }                                                                                               # Used to group results in summary.

        benchmark_stats = {                                                                             # Global benchmark metadata saved alongside raw results.
            'total_benchmark_time_seconds': float(total_benchmark_time),                                # Total elapsed time.
            'total_tests': int(len(results)),                                                           # Total collected results.
            'successful_tests': int(len([r for r in results if r is not None])),                        # Count of non-null records.
            'timestamp': datetime.now().isoformat(),                                                    # ISO timestamp.
            'data_root': data_root,                                                                     # Data root used for discovery.
            'equations_tested': list(equations)                                                    # Equations included in this run.
        }                                                                                               # Metadata stored alongside results.
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")                                            # Timestamp used in output filenames.
        json_file = f'benchmark_results/benchmark_{timestamp}.json'                                     # JSON output path.
        
        output_data = {                                                                                 # JSON payload containing metadata + raw results.
            'benchmark_stats': benchmark_stats,                                                         # Global benchmark metadata.
            'results': results                                                                          # Raw per-run results.
        }                                                                                               # JSON payload.
        
        with open(json_file, 'w') as f:                                                                 # Write JSON output file.
            json.dump(output_data, f, indent = 2)                                                       # Write pretty-printed JSON.
        
        csv_file = f'benchmark_results/benchmark_{timestamp}.csv'                                       # CSV output path.
        df = pd.DataFrame(results)                                                                      # Convert to DataFrame.
        df.to_csv(csv_file, index = False)                                                              # Write CSV table.
        
        summary_file = f'benchmark_results/benchmark_summary_{timestamp}.txt'                           # TXT summary output path.
        with open(summary_file, 'w') as f:                                                              # Write human-readable summary file.
            f.write("RESUMEN DEL BENCHMARK\n")                                                          # Summary header.
            f.write("=====================\n\n")                                                        # Separator.
            f.write(f"Tiempo total del benchmark: {total_benchmark_time:.2f} segundos\n")               # Total benchmark wall time.
            f.write(f"Pruebas exitosas: {len(results)} de {total_combinations}\n\n")                    # Success count for expected runs.
            
            if results:                                                                                 # Only aggregate statistics when results exist.
                df = pd.DataFrame(results)                                                              # Build DataFrame for aggregation.
                f.write("ESTADÍSTICAS POR ECUACIÓN:\n")                                                 # Group header.
                for eq in equations:                                                                    # Aggregate by equation label.
                    eq_name = equation_names.get(eq, eq)                                                # Display name mapping.
                    eq_results = df[df['equation'] == eq_name]                                          # Filter results by equation.
                    if not eq_results.empty:                                                            # Only write section when data is available.
                        f.write(f"\n{eq_name}:\n")                                                      # Equation section header.
                        f.write(f"  Tiempo promedio: {eq_results['execution_time_seconds'].mean():.3f}s\n")
                                                                                                        # Mean time across clouds.
                        f.write(f"  Memoria promedio: {eq_results['memory_used_mb'].mean():.1f}MB\n")   # Mean memory delta across clouds.
                        f.write(f"  Error promedio: {eq_results['avg_numerical_error'].mean():.2e}\n")  # Mean error across clouds.
        
        print("\nResultados guardados en:")                                                             # Output location header.
        print(f"  JSON: {json_file}")                                                                   # JSON path.
        print(f"  CSV: {csv_file}")                                                                     # CSV path.
        print(f"  Resumen: {summary_file}")                                                             # TXT summary path.
    
    print(f"\nBenchmark completado en {total_benchmark_time:.2f} segundos")                             # Total time report.
    print(f"Total de resultados: {len(results)}")                                                       # Results count report.
    return results                                                                                      # Return results list.

if __name__ == "__main__":                                                                              # Script entry point.
    try:                                                                                                # Validate optional dependencies before running.
        _require_optional_deps()                                                                        # Ensure psutil and pandas are available.
    except ImportError as e:                                                                            # Missing dependencies.
        print(f"Error: {e}")                                                                            # Print actionable error message.
        raise                                                                                           # Propagate error for non-zero exit.

    results = run_comprehensive_benchmark()                                                             # Execute benchmark with default settings.
