"""
Sweep — Parameter Sweep Orchestrator

Overview:
    This script reads sweep_config.json and executes all the runners
    iterating through combinatorial configurations of nvec and linear_solver.

Notes:
    - Each batch script is executed by importing it as a module with dynamically injected kwargs.
    - Supports concurrent CPU and CUDA pipeline parallelization when parallel=True in sweep_config.json.

Public API:
    import_module_from_file     Dynamically load and execute a batch runner module.
    run_job                     Worker execution wrapper for process pool parallelization.
    main                        Execute multi-cloud parameter sweep.

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
    September, 2026.
"""

import os                                                                                                                               # OS interfaces for file/directory paths.
import json                                                                                                                             # JSON formatting.
import logging                                                                                                                          # Standard logging module.
import itertools                                                                                                                        # Combinatorial iterator.
import sys                                                                                                                              # sys.modules access for imported modules.
import importlib.util                                                                                                                   # Dynamic module import from a file path.
import concurrent.futures                                                                                                               # Process pool executor for parallel sweep execution.
from time import time                                                                                                                   # Time tracking.
from typing import List, Tuple, Any                                                                                                     # Type hints.

logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Configure standard logger format.
logger = logging.getLogger(__name__)                                                                                                    # Create module logger.

def import_module_from_file(file_path: str, verbose: bool = True, **kwargs: Any) -> bool:
    """
    import_module_from_file
    Dynamically loads and executes a Python module from its absolute filepath.

    Input:
        file_path       1           str             Absolute path to the Python file to load.
        verbose         1           bool            If True, prints errors to console.
        **kwargs        1           Any             Additional arguments to pass to the module's main().

    Output:
        ok              1           bool            True when the import succeeds, False otherwise.
    """
    if not isinstance(file_path, str):                                                                                                  # Validate file_path.
        raise TypeError("file_path must be a string.")                                                                                  # Explicit error on bad input.
    
    if not os.path.exists(file_path):                                                                                                   # Validate file exists.
        raise ValueError(f"File not found: {file_path}")                                                                                # Explicit error on bad input.

    try:                                                                                                                                # Catch import-time errors to keep the batch runner alive.
        script_dir = os.path.dirname(os.path.abspath(file_path))                                                                        # Get directory of the script.
        
        if script_dir not in sys.path:                                                                                                  # Prevent duplicate paths.
            sys.path.insert(0, script_dir)                                                                                              # Add script's directory to sys.path so it can find its local imports.
            
        module_name = os.path.splitext(os.path.basename(file_path))[0]                                                                  # Derive a stable module name from filename.
        spec        = importlib.util.spec_from_file_location(module_name, file_path)                                                    # Build import spec from file path.
        
        if spec is None or spec.loader is None:                                                                                         # Ensure a valid spec and loader were found.
            raise ImportError(f"Cannot load module from {file_path}")                                                                   # Raise an exception if spec is missing.
        
        module                   = importlib.util.module_from_spec(spec)                                                                # Create module object from spec.
        sys.modules[module_name] = module                                                                                               # Register module to allow intra-import references.
        spec.loader.exec_module(module)                                                                                                 # Execute the module (runs top-level batch script).
        
        if hasattr(module, 'main'):                                                                                                     # Explicitly call main() if it exists.
            module.main(**kwargs)                                                                                                       # Execute statement.
            
        return True                                                                                                                     # Report success.
    
    except Exception as e:                                                                                                              # Catch any failure while importing/executing the script.
        if verbose:                                                                                                                     # Check if verbosity is enabled.
            import traceback; logger.error(f'Error importing {file_path}: {traceback.format_exc()}')                                    # Print an actionable error message.
        
        return False                                                                                                                    # Report failure.


def run_job(job_tuple: Tuple[str, str, int, Any, dict, bool]) -> bool:
    """
    run_job
    Worker execution wrapper for process pool parallelization.

    Input:
        job_tuple                   Tuple           (file_path, run_file, nvec, upwind, kwargs, verbose)

    Output:
        ok                          bool            True when the execution succeeds, False otherwise.
    """
    file_path, run_file, nvec, upwind, base_kwargs, verbose = job_tuple                                                                 # Unpack job configuration tuple.
    
    if upwind is not None:                                                                                                              # Determine if upwind applies.
        config_str = f'nvec={nvec}, dev={base_kwargs.get("device")}, up={upwind}'                                                       # Format config with upwind.
    else:                                                                                                                               # Upwind doesn't apply.
        config_str = f'nvec={nvec}, dev={base_kwargs.get("device")}'                                                                    # Format config without upwind.

    if verbose:                                                                                                                         # Check verbosity.
        logger.info(f'\nExecuting {run_file} with [{config_str}]...')                                                                   # Log current execution.

    start_time = time()                                                                                                                 # Start per-script timer.
    run_kwargs = base_kwargs.copy()                                                                                                     # Copy global kwargs to avoid side-effects.
    run_kwargs['nvec'] = nvec                                                                                                           # Inject neighbor count into kwargs.

    if upwind is not None:                                                                                                              # If upwind is required for this runner.
        run_kwargs['upwind'] = upwind                                                                                                   # Inject upwind flag into kwargs.

    ok = import_module_from_file(file_path, verbose=verbose, **run_kwargs)                                                              # Dynamically load and run script with kwargs.
    execution_time = time() - start_time                                                                                                # Compute per-script runtime.

    if ok:                                                                                                                              # Check execution success.
        if verbose:                                                                                                                     # Check verbosity.
            logger.info(f'Completed {run_file} [{config_str}] in {execution_time:.2f} seconds')                                         # Log completion.
    else:                                                                                                                               # Script execution failed.
        if verbose:                                                                                                                     # Check verbosity.
            logger.error(f'Error executing {run_file} [{config_str}]')                                                                  # Report failure.

    return ok                                                                                                                           # Return execution result.


def main(verbose: bool = True, **kwargs: Any) -> None:
    """
    main
    Entry point for the parameter sweep orchestrator.

    Input:
        verbose         1           bool            If True, prints execution details to the console.
        **kwargs        1           Any             Additional arguments passed directly to the runners.

    Output:
        None
    """
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sweep_config.json')                                         # Resolve path to sweep configuration.
    if not os.path.exists(config_path):                                                                                                 # Check if the config file exists.
        logger.error(f'Error: Sweep configuration not found at {config_path}')                                                          # Report missing config.
        return                                                                                                                          # Abort sweep.
        
    with open(config_path, 'r') as f:                                                                                                   # Open the configuration file.
        config = json.load(f)                                                                                                           # Parse JSON into a dictionary.
        
    nvec_list    = config.get('nvec', [12])                                                                                             # Extract neighbor list, default to [12].
    device_list  = config.get('device', ['cpu'])                                                                                        # Extract device list.
    upwind_list  = config.get('upwind', [True, False])                                                                                  # Extract upwind list, default to [True, False].
    scales_cfg   = config.get('scales', ['1', '2', '3'])                                                                                # Extract scales list, default to ['1', '2', '3'].
    t_cfg        = config.get('t', None)                                                                                                # Extract explicit t time steps if provided.
    cfl_cfg      = config.get('cfl', None)                                                                                              # Extract CFL target override if provided.
    input_types  = config.get('input_types', ['callable'])                                                                              # Extract input types, default to callable only.
    save_cfg     = config.get('save', False)                                                                                            # Extract save config.
    save_outputs = save_cfg[0] if isinstance(save_cfg, list) else save_cfg                                                              # Safely extract boolean.
    plot_cfg     = config.get('plot_approximations', False)                                                                             # Extract plot config.
    plot_appx    = plot_cfg[0] if isinstance(plot_cfg, list) else plot_cfg                                                              # Safely extract boolean.
    parallel_cfg = config.get('parallel', True)                                                                                         # Extract parallel config (default True).
    cpu_workers  = config.get('cpu_workers', 2)                                                                                         # Extract CPU worker count (default 2).

    if t_cfg is not None:                                                                                                               # If explicit t parameter is provided.
        kwargs['t'] = t_cfg[0] if isinstance(t_cfg, list) else t_cfg                                                                    # Inject explicit t into kwargs.
    if cfl_cfg is not None:                                                                                                             # If cfl parameter override is provided.
        kwargs['cfl'] = cfl_cfg[0] if isinstance(cfl_cfg, list) else cfl_cfg                                                            # Inject cfl into kwargs.
    kwargs['scales'] = scales_cfg                                                                                                       # Inject scales list into kwargs.
    kwargs['save']   = save_outputs                                                                                                     # Inject save flag into kwargs.
    kwargs['input_types'] = input_types                                                                                                 # Inject input types into kwargs.
    kwargs['plot_approximations'] = plot_appx                                                                                           # Inject plot approximations flag into kwargs.
    
    runners_cfg  = config.get('runners', ['run_Poisson.py', 'run_AdvReactionDiff.py', 'run_Wave.py'])                                   # Extract runners list from config.
    run_files: List[str] = runners_cfg if isinstance(runners_cfg, list) else [runners_cfg]                                              # Resolve list of runners to execute.
    
    current_dir: str = os.path.dirname(os.path.abspath(__file__))                                                                       # Directory where this script is located.
    
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('Starting mGFD Parameter Sweep...')                                                                                 # Console header.
        logger.info(f'Parallel execution: {parallel_cfg} (cpu_workers={cpu_workers})')                                                  # Log parallel configuration.
        logger.info('=' * 50)                                                                                                           # Visual separator.

    total_start_time = time()                                                                                                           # Start total timer.
    successful_runs  = 0                                                                                                                # Counter for successful runs.
    
    jobs: List[Tuple[str, str, int, Any, dict, bool]] = []                                                                              # Initialize jobs list.

    for run_file in run_files:                                                                                                          # Iterate over runners.
        file_path = os.path.join(current_dir, 'runners', run_file)                                                                      # Resolve path under runners/ folder.
        if not os.path.exists(file_path):                                                                                               # Check if the runner script exists.
            if verbose:                                                                                                                 # Check verbosity.
                logger.warning(f'Warning: File {run_file} not found')                                                                   # Warn if runner is missing.
            continue                                                                                                                    # Move to the next script.

        if run_file in ['run_AdvDif.py', 'run_Perturbation.py']:                                                                        # These scripts require upwind parameter.
            combinations = list(itertools.product(nvec_list, device_list, upwind_list))                                                 # Generate combinations including upwind.
        else:                                                                                                                           # Other scripts don't use upwind.
            combinations = list(itertools.product(nvec_list, device_list, [None]))                                                      # Generate combinations without upwind.

        for nvec, device, upwind in combinations:                                                                                       # Iterate over configurations.
            run_kwargs = kwargs.copy()                                                                                                  # Copy global kwargs.
            run_kwargs['device'] = device                                                                                               # Inject device into kwargs.
            jobs.append((file_path, run_file, nvec, upwind, run_kwargs, verbose))                                                       # Append job tuple to queue.

    total_runs = len(jobs)                                                                                                              # Total number of execution jobs.

    if parallel_cfg and total_runs > 1:                                                                                                 # Execute jobs in parallel mode.
        max_workers = cpu_workers + (1 if 'cuda' in device_list else 0)                                                                 # Calculate total process pool size.
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:                                               # Initialize process pool executor.
            futures = [executor.submit(run_job, job) for job in jobs]                                                                   # Submit all sweep jobs to executor.
            for future in concurrent.futures.as_completed(futures):                                                                     # Gather completed job futures.
                if future.result():                                                                                                     # Check job execution success.
                    successful_runs += 1                                                                                                # Increment success counter.
    else:                                                                                                                               # Fallback to sequential sweep mode.
        for job in jobs:                                                                                                                # Iterate sequentially over jobs.
            if run_job(job):                                                                                                            # Execute job sequentially.
                successful_runs += 1                                                                                                    # Increment success counter.

    total_time = time() - total_start_time                                                                                              # Compute total runtime.
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('\n' + '=' * 50)                                                                                                    # Visual separator.
        logger.info(f'Sweep completed: {successful_runs} of {total_runs} combinations successfully executed')                           # Report success count.
        logger.info(f'Total execution time: {total_time:.2f} seconds')                                                                  # Report total runtime.
        
    # Generate final summary CSV always
    try:                                                                                                                                # Execute statement.
        sys.path.append(os.path.join(current_dir, 'utils'))                                                                             # Execute statement.
        from batch_utils import generate_sweep_summary                                                                                  # Import module dependency.
        results_root = os.path.join(os.path.dirname(current_dir), 'results')                                                            # Assign results_root.
        generate_sweep_summary(results_root, verbose=verbose)                                                                           # Assign generate_sweep_summary(results_root, verbose.
    except Exception as e:                                                                                                              # Execute statement.
        if verbose:                                                                                                                     # Evaluate condition.
            logger.error(f"Failed to generate final sweep summary CSV: {e}")                                                            # Log output message.

if __name__ == '__main__':                                                                                                              # Check if script is run directly.
    main()                                                                                                                              # Execute main orchestrator.
