"""
Sweep — Parameter Sweep Orchestrator

Overview:
    This script reads sweep_config.json and executes all the runners
    iterating through combinatorial configurations of nvec and linear_solver.

Notes:
    - Each batch script is executed by importing it as a module with dynamically injected kwargs.
    - This is intended for parameter exploration across the mGFD suite.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
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
    August, 2026.
"""

import os                                                                                                                               # OS interfaces for file/directory paths.
import json                                                                                                                             # JSON formatting.
import logging                                                                                                                          # Standard logging module.
import itertools                                                                                                                        # Combinatorial iterator.
from time import time                                                                                                                   # Time tracking.
from typing import List, Any                                                                                                            # Type hints.
import sys                                                                                                                              # sys.modules access for imported modules.
import importlib.util                                                                                                                   # Dynamic module import from a file path.

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
            module.main(**kwargs)
            
        return True                                                                                                                     # Report success.
    
    except Exception as e:                                                                                                              # Catch any failure while importing/executing the script.
        if verbose:                                                                                                                     # Check if verbosity is enabled.
            logger.error(f'Error importing {file_path}: {str(e)}')                                                                      # Print an actionable error message.
        
        return False                                                                                                                    # Report failure.

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
    solver_list  = config.get('linear_solver', ['spsolve'])                                                                             # Extract solver list, default to ['spsolve'].
    upwind_list  = config.get('upwind', [True, False])                                                                                  # Extract upwind list, default to [True, False].
    save_outputs = config.get('save', False)                                                                                            # Extract save flag, default to False.
    kwargs['save'] = save_outputs                                                                                                       # Inject save flag into kwargs.
    
    run_files: List[str] = [                                                                                                            # List of files to execute.
        'run_Poisson.py',                                                                                                               # Stationary Poisson reference.
        'run_Heat.py',                                                                                                                  # Transient Heat reference.
        'run_Wave.py',                                                                                                                  # Transient Wave reference.
        'run_AdvDif.py',                                                                                                                # Transient Advection-Diffusion reference.
        'run_Perturbation.py'                                                                                                           # Stationary perturbation case.
    ]                                                                                                                                   # End of list.
    
    current_dir: str = os.path.dirname(os.path.abspath(__file__))                                                                       # Directory where this script is located.
    
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('Starting mGFD Parameter Sweep...')                                                                                 # Console header.
        logger.info('=' * 50)                                                                                                           # Visual separator.

    total_start_time = time()                                                                                                           # Start total timer.
    successful_runs  = 0                                                                                                                # Counter for successful runs.
    total_runs       = 0                                                                                                                # Counter for total executions.

    for run_file in run_files:                                                                                                          # Iterate over runners.
        file_path = os.path.join(current_dir, 'runners', run_file)                                                                      # Resolve path under runners/ folder.
        if not os.path.exists(file_path):                                                                                               # Check if the runner script exists.
            if verbose:                                                                                                                 # Check verbosity.
                logger.warning(f'Warning: File {run_file} not found')                                                                   # Warn if runner is missing.
            continue                                                                                                                    # Move to the next script.

        if run_file in ['run_AdvDif.py', 'run_Perturbation.py']:                                                                        # These scripts require upwind parameter.
            combinations = list(itertools.product(nvec_list, solver_list, upwind_list))                                                 # Generate combinations including upwind.
        else:                                                                                                                           # Other scripts don't use upwind.
            combinations = list(itertools.product(nvec_list, solver_list, [None]))                                                      # Generate combinations without upwind.
            
        total_runs += len(combinations)                                                                                                 # Accumulate total executions.
        
        for nvec, solver, upwind in combinations:                                                                                       # Iterate over configurations.
            if upwind is not None:                                                                                                      # Determine if upwind applies.
                config_str = f'nvec={nvec}, solver={solver}, upwind={upwind}'                                                           # Format config with upwind.
            else:                                                                                                                       # Upwind doesn't apply.
                config_str = f'nvec={nvec}, solver={solver}'                                                                            # Format config without upwind.
                
            if verbose:                                                                                                                 # Check verbosity.
                logger.info(f'\nExecuting {run_file} with [{config_str}]...')                                                           # Log current execution.
            
            start_time = time()                                                                                                         # Start per-script timer.
            run_kwargs = kwargs.copy()                                                                                                  # Copy global kwargs to avoid side-effects.
            run_kwargs['nvec']          = nvec                                                                                          # Inject neighbor count into kwargs.
            run_kwargs['linear_solver'] = solver                                                                                        # Inject solver backend into kwargs.
            if upwind is not None:                                                                                                      # If upwind is required for this runner.
                run_kwargs['upwind']    = upwind                                                                                        # Inject upwind flag into kwargs.
            
            if import_module_from_file(file_path, verbose=verbose, **run_kwargs):                                                       # Dynamically load and run script with kwargs.
                execution_time = time() - start_time                                                                                    # Compute per-script runtime.
                if verbose:                                                                                                             # Check verbosity.
                    logger.info(f'Completed {run_file} [{config_str}] in {execution_time:.2f} seconds')                                 # Log completion.
                successful_runs += 1                                                                                                    # Increment success counter.
            else:                                                                                                                       # Script execution failed.
                if verbose:                                                                                                             # Check verbosity.
                    logger.error(f'Error executing {run_file} [{config_str}]')                                                          # Report failure.

    total_time = time() - total_start_time                                                                                              # Compute total runtime.
    if verbose:                                                                                                                         # Check verbosity.
        logger.info('\n' + '=' * 50)                                                                                                    # Visual separator.
        logger.info(f'Sweep completed: {successful_runs} of {total_runs} combinations successfully executed')                           # Report success count.
        logger.info(f'Total execution time: {total_time:.2f} seconds')                                                                  # Report total runtime.

if __name__ == '__main__':                                                                                                              # Check if script is run directly.
    main()                                                                                                                              # Execute main orchestrator.
