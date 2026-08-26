"""
main — Batch runner for mGFD reference scripts

Overview:
    This script executes the reference batch scripts under runners/ in a predefined order, reporting
    per-script runtime and a global summary at the end.

Notes:
    - Each batch script is executed by importing it as a module, which triggers its top-level code.
    - This is intended for batch-style scripts that run immediately on import.

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
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.modules access for imported modules.
import logging                                                                                                                          # Standard logging module.
import importlib.util                                                                                                                   # Dynamic module import from a file path.

import json                                                                                                                             # JSON parser for configuration loading.
import argparse                                                                                                                         # Command-line arguments parser.

from time import time                                                                                                                   # Wall-clock timing for runtimes.
from typing import List, Dict, Any                                                                                                      # Type hints.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

def import_module_from_file(file_path: str, verbose: bool = True, **kwargs: Any) -> bool:
    """
    import_module_from_file
    Dynamically loads and executes a Python module from its absolute filepath. and execute it as a module.

    Input:
        file_path       1           str             Absolute path to the Python file to load.
        verbose         1           bool            If True, prints errors to console.
        **kwargs        1           Any             Additional arguments to pass to the module's main().

    Output:
        ok              1           bool            True when the import succeeds, False otherwise.
    """
    # 0. Input validation
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
        elif hasattr(module, f'test_single_cloud'):
            module.test_single_cloud()
            
        return True                                                                                                                     # Report success.
    
    except Exception as e:                                                                                                              # Catch any failure while importing/executing the script.
        if verbose:                                                                                                                     # Check if verbosity is enabled.
            logger.error(f'Error importing {file_path}: {str(e)}')                                                                      # Print an actionable error message.
        
        return False                                                                                                                    # Report failure.

def main(verbose: bool = True, **kwargs: Any) -> None:
    """
    main
    Run the predefined list of reference batch scripts under runners/ and print a summary.

    Input:
        verbose         1           bool            If True, prints execution details to the console.
        **kwargs        1           Any             Additional arguments passed directly to the runners.

    Output:
        None
    """
    run_files: List[str] = [                                                                                                            # List of files to execute in order.
        'run_Poisson.py',                                                                                                               # Stationary Poisson reference.
        'run_Heat.py',                                                                                                                  # Transient Heat reference.
        'run_Wave.py',                                                                                                                  # Transient Wave reference.
        'run_AdvDif.py',                                                                                                                # Transient Advection–Diffusion reference.
        'run_Perturbation.py'                                                                                                           # Stationary perturbation case (Adv=True).
    ]                                                                                                                                   # End of run list.

    current_dir: str = os.path.dirname(os.path.abspath(__file__))                                                                       # Directory where this script is located (repo root).

    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info('Starting mGFD scripts execution...')                                                                               # Console header for batch run.
        logger.info('=' * 50)                                                                                                           # Visual separator.

    total_start_time     = time()                                                                                                       # Start total timer.
    successful_runs: int = 0                                                                                                            # Count successful script executions.

    for run_file in run_files:                                                                                                          # Iterate scripts in predefined order.
        file_path = os.path.join(current_dir, 'runners', run_file)                                                                      # Resolve path under runners/ folder.
        
        if not os.path.exists(file_path):                                                                                               # Skip missing files without stopping the suite.
            if verbose:                                                                                                                 # Check if verbosity is enabled.
                logger.warning(f'Warning: File {run_file} not found')                                                                   # Report missing script.
            continue                                                                                                                    # Move to next script.

        if verbose:                                                                                                                     # Check if verbosity is enabled.
            logger.info(f'\nExecuting {run_file}...')                                                                                   # Report script execution start.
        start_time = time()                                                                                                             # Start per-script timer.

        if import_module_from_file(file_path, verbose=verbose, **kwargs):                                                               # Import and execute the script module.
            end_time = time()                                                                                                           # Stop per-script timer.
            execution_time = end_time - start_time                                                                                      # Compute per-script runtime.
            
            if verbose:                                                                                                                 # Check if verbosity is enabled.
                logger.info(f'Completed {run_file} in {execution_time:.2f} seconds')                                                    # Report per-script runtime.
            
            successful_runs += 1                                                                                                        # Count successful execution.
        else:                                                                                                                           # Script import/execution failed.
            if verbose:                                                                                                                 # Check if verbosity is enabled.
                logger.error(f'Error executing {run_file}')                                                                             # Report failure while keeping the suite running.

    total_time = time() - total_start_time                                                                                              # Compute total elapsed time.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info('\n' + '=' * 50)                                                                                                    # Visual separator before summary.
        logger.info(f'Execution completed: {successful_runs} of {len(run_files)} scripts successfully executed')                        # Print success count.
        logger.info(f'Total execution time: {total_time:.2f} seconds')                                                                  # Print total runtime.

def parse_and_run() -> None:
    """
    parse_and_run
    Parses CLI arguments, loads the JSON configuration, merges them, and runs main().
    """
    parser = argparse.ArgumentParser(description="mGFD Batch Runner Orchestrator")                                                      # Create CLI parser.
    parser.add_argument('--solver', type=str, default=None, choices=['spsolve', 'bicgstab', 'gmres'], help="Linear solver backend")     # Add solver flag.
    parser.add_argument('--nvec', type=int, default=None, help="Number of neighbors (e.g. 12, 16, 20)")                                 # Add nvec flag.
    parser.add_argument('--upwind', type=str, default=None, choices=['true', 'false'], help="Enable/disable upwind (for AdvDif)")       # Add upwind flag.
    parser.add_argument('--verbose-solvers', action='store_true', help="Make individual solvers print their iterations")                # Add verbosity flag.
    
    args = parser.parse_args()                                                                                                          # Parse terminal arguments.
    
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')                                               # Find config.json.
    config: Dict[str, Any] = {}                                                                                                         # Initialize configuration dict.
    
    if os.path.exists(config_path):                                                                                                     # If config file exists.
        try:
            with open(config_path, 'r') as f:                                                                                           # Read file.
                config = json.load(f)                                                                                                   # Parse JSON dict.
        except Exception as e:                                                                                                          # Catch malformed JSON.
            logger.error(f"Error loading config.json: {e}")                                                                             # Warn about broken JSON.
            
    # Overwrite config with explicitly provided CLI arguments
    if args.solver is not None:
        config['linear_solver'] = args.solver
    if args.nvec is not None:
        config['nvec'] = args.nvec
    if args.upwind is not None:
        config['upwind'] = (args.upwind.lower() == 'true')
    if args.verbose_solvers:
        config['verbose_solvers'] = True
        
    main(verbose=True, **config)                                                                                                        # Run the batch suite with the merged config.

if __name__ == '__main__':                                                                                                              # Script entry point.
    parse_and_run()                                                                                                                     # Execute parse sequence.
