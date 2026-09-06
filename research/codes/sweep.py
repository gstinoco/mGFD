"""
Sweep — Parameter Sweep Orchestrator

Overview:
    This script reads sweep_config.json and executes all the runners
    iterating through combinatorial configurations of nvec, device, and upwind.

Notes:
    - Each batch script is executed by importing it as a module with dynamically injected kwargs.
    - Supports concurrent CPU and CUDA pipeline parallelization when parallel=True in sweep_config.json.

Public API:
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

## Library importation.
import os                                                                                                                               # OS interfaces for file/directory paths.
import sys                                                                                                                              # System specific parameters and path manipulation.
import json                                                                                                                             # JSON formatting.
import logging                                                                                                                          # Standard logging module.
import itertools                                                                                                                        # Combinatorial iterator.
import concurrent.futures                                                                                                               # Process pool executor for parallel sweep execution.

from time import time                                                                                                                   # High-resolution wall-clock timer.
from typing import List, Tuple, Any, Dict                                                                                               # Type hints.

_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))                                               # Resolve repository root directory.
_src_dir   = os.path.join(_repo_root, 'src')                                                                                            # Resolve core library source directory.
if _src_dir not in sys.path:                                                                                                            # Ensure src directory is on system path.
    sys.path.insert(0, _src_dir)                                                                                                        # Add src directory to system path.

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils'))                                                   # Ensure utils directory is on path.
from batch_utils import run_job, generate_sweep_summary                                                                                 # Multiprocessing job wrapper and summary generator.

logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Configure standard logger format.
logger = logging.getLogger(__name__)                                                                                                    # Create module logger.

def main(verbose: bool = True, **kwargs: Any) -> None:
    """
    Entry point for the parameter sweep orchestrator.

    Loads the sweep configuration, builds combinatorial job matrices, dispatches
    tasks across an isolated process pool, and compiles consolidated CSV summaries.

    Input:
        verbose         1           bool            If True, prints execution details to the console.
        **kwargs                    Any             Additional arguments passed directly to the runners.

    Output:
        None
    """
    # 0. Configuration loading
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sweep_config.json')                                         # Resolve path to sweep configuration.
    if not os.path.exists(config_path):                                                                                                 # Check if the config file exists.
        logger.error(f'Error: Sweep configuration not found at {config_path}')                                                          # Report missing config.
        return                                                                                                                          # Abort sweep.
        
    with open(config_path, 'r') as f:                                                                                                   # Open the configuration file.
        config = json.load(f)                                                                                                           # Parse JSON into a dictionary.
        
    # 1. Parameter parsing and kwargs preparation
    nvec_list       = kwargs.get('nvec', config.get('nvec', [12]))                                                                      # Extract neighbor list, default to [12].
    device_list     = kwargs.get('device', config.get('device', ['cpu']))                                                               # Extract device list.
    upwind_list     = kwargs.get('upwind', config.get('upwind', [True, False]))                                                         # Extract upwind list, default to [True, False].
    scales_cfg      = kwargs.get('scales', config.get('scales', ['1', '2', '3', '4']))                                                  # Extract scales list, default to 1-4.
    t_cfg           = kwargs.get('t', config.get('t', None))                                                                            # Extract explicit t time steps if provided.
    cfl_cfg         = kwargs.get('cfl', config.get('cfl', None))                                                                        # Extract CFL target override if provided.
    damping_cfg     = kwargs.get('damping', config.get('damping', None))                                                                # Extract numerical damping if provided.
    alpha_cfg       = kwargs.get('alpha', config.get('alpha', None))                                                                    # Extract HHT-alpha dissipation if provided.
    save_cfg        = kwargs.get('save', config.get('save', False))                                                                     # Extract save config.
    save_outputs    = bool(save_cfg[0]) if isinstance(save_cfg, list) and save_cfg else bool(save_cfg)                                  # Safely extract boolean.
    plot_scales_cfg = kwargs.get('plot_scales', config.get('plot_scales', ['2']))                                                       # Extract plot scales configuration.
    plot_scales     = plot_scales_cfg if isinstance(plot_scales_cfg, list) else [plot_scales_cfg]                                       # Ensure plot scales is a list.
    plot_scales     = [str(s) for s in plot_scales]                                                                                     # Cast all scale identifiers to string.
    parallel_cfg    = kwargs.get('parallel', config.get('parallel', True))                                                              # Extract parallel config (default True).
    parallel_run    = bool(parallel_cfg[0]) if isinstance(parallel_cfg, list) and parallel_cfg else bool(parallel_cfg)                  # Safely extract boolean.
    workers_cfg     = kwargs.get('cpu_workers', config.get('cpu_workers', 2))                                                           # Extract CPU worker count (default 2).
    worker_val      = workers_cfg[0] if isinstance(workers_cfg, list) and workers_cfg else workers_cfg                                  # Extract worker count value.
    cpu_workers     = int(worker_val) if worker_val is not None else 2                                                                  # Safely extract integer worker count.

    nvec_list    = nvec_list if isinstance(nvec_list, list) else [nvec_list]                                                            # Ensure neighbor counts are in a list.
    device_list  = device_list if isinstance(device_list, list) else [device_list]                                                      # Ensure devices are in a list.
    upwind_list  = upwind_list if isinstance(upwind_list, list) else [upwind_list]                                                      # Ensure upwind options are in a list.

    if t_cfg is not None:                                                                                                               # If explicit t parameter is provided.
        kwargs['t'] = t_cfg[0] if isinstance(t_cfg, list) else t_cfg                                                                    # Inject explicit t into kwargs.
    if cfl_cfg is not None:                                                                                                             # If cfl parameter override is provided.
        kwargs['cfl'] = cfl_cfg[0] if isinstance(cfl_cfg, list) else cfl_cfg                                                            # Inject cfl into kwargs.
    if damping_cfg is not None:                                                                                                         # If damping parameter override is provided.
        kwargs['damping'] = damping_cfg[0] if isinstance(damping_cfg, list) else damping_cfg                                            # Inject damping into kwargs.
    if alpha_cfg is not None:                                                                                                           # If alpha parameter override is provided.
        kwargs['alpha'] = alpha_cfg[0] if isinstance(alpha_cfg, list) else alpha_cfg                                                    # Inject alpha into kwargs.

    kwargs['scales']      = scales_cfg                                                                                                  # Inject scales list into kwargs.
    kwargs['save']        = save_outputs                                                                                                # Inject save flag into kwargs.
    kwargs['plot_scales'] = plot_scales                                                                                                 # Inject plot scales into kwargs.
    
    default_runs = ['run_Poisson.py', 'run_AdvReactionDiff.py', 'run_Wave.py']                                                          # Default suite runners.
    runners_cfg  = kwargs.get('runners', None) or config.get('runners', default_runs)                                                   # Extract runners configuration.
    run_files: List[str] = [str(r) for r in runners_cfg] if isinstance(runners_cfg, list) else [str(runners_cfg)]                       # Resolve list of runners to execute.
    
    current_dir: str = os.path.dirname(os.path.abspath(__file__))                                                                       # Directory where this script is located.
    
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('Starting mGFD Parameter Sweep...')                                                                                 # Console header.
        logger.info(f'Parallel execution: {parallel_run} (cpu_workers={cpu_workers})')                                                  # Log parallel configuration.
        logger.info('=' * 50)                                                                                                           # Visual separator.

    total_start_time = time()                                                                                                           # Start total timer.
    successful_runs  = 0                                                                                                                # Counter for successful runs.
    
    # 2. Job queue assembly
    jobs: List[Tuple[str, str, Any, Any, Dict[str, Any], bool]] = []                                                                    # Initialize jobs list.

    for run_file in run_files:                                                                                                          # Iterate over runners.
        file_path = os.path.join(current_dir, 'runners', run_file)                                                                      # Resolve path under runners/ folder.
        if not os.path.exists(file_path):                                                                                               # Check if the runner script exists.
            if verbose:                                                                                                                 # Check verbosity.
                logger.warning(f'Warning: File {run_file} not found')                                                                   # Warn if runner is missing.
            continue                                                                                                                    # Move to the next script.

        if run_file in ['run_AdvDif.py', 'run_AdvReactionDiff.py', 'run_Perturbation.py']:                                              # These scripts require upwind parameter.
            combinations = list(itertools.product(nvec_list, device_list, upwind_list))                                                 # Generate combinations including upwind.
        else:                                                                                                                           # Other scripts don't use upwind.
            combinations = list(itertools.product(nvec_list, device_list, [None]))                                                      # Generate combinations without upwind.

        for nvec, device, upwind in combinations:                                                                                       # Iterate over configurations.
            run_kwargs           = kwargs.copy()                                                                                        # Copy global kwargs.
            run_kwargs['device'] = device                                                                                               # Inject device into kwargs.
            jobs.append((file_path, run_file, nvec, upwind, run_kwargs, verbose))                                                       # Append job tuple to queue.

    total_runs = len(jobs)                                                                                                              # Total number of execution jobs.

    # 3. Process pool parallel execution
    max_workers = (cpu_workers + (1 if 'cuda' in device_list else 0)) if (parallel_run and total_runs > 1) else 1                       # Calculate total process pool size.
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, max_tasks_per_child=1) as executor:                            # Initialize isolated process pool.
        futures = [executor.submit(run_job, job) for job in jobs]                                                                       # Submit all sweep jobs to executor.
        for future in concurrent.futures.as_completed(futures):                                                                         # Gather completed job futures.
            if future.result():                                                                                                         # Check job execution success.
                successful_runs += 1                                                                                                    # Increment success counter.

    total_time = time() - total_start_time                                                                                              # Compute total runtime.
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('\n' + '=' * 50)                                                                                                    # Visual separator.
        logger.info(f'Sweep completed: {successful_runs} of {total_runs} combinations successfully executed')                           # Report success count.
        logger.info(f'Total execution time: {total_time:.2f} seconds')                                                                  # Report total runtime.
        
    # 4. Summary CSV aggregation
    try:                                                                                                                                # Wrap summary generation in try-except block.
        results_root = os.path.join(os.path.dirname(current_dir), 'results')                                                            # Resolve results directory path.
        generate_sweep_summary(results_root, verbose=verbose)                                                                           # Build consolidated summary CSV.
    except Exception as e:                                                                                                              # Catch error during summary generation.
        if verbose:                                                                                                                     # Check verbosity.
            logger.error(f"Failed to generate final sweep summary CSV: {e}")                                                            # Log error message.

if __name__ == '__main__':                                                                                                              # Check if script is run directly.
    main()                                                                                                                              # Execute main orchestrator.
