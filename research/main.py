"""
main — Batch runner for mGFD reference scripts

Overview:
    This script executes the reference batch scripts under batches/ in a predefined order, reporting
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
import os                                                                                               # Filesystem and path utilities.
import sys                                                                                              # sys.modules access for imported modules.
import importlib.util                                                                                   # Dynamic module import from a file path.

from time import time                                                                                   # Wall-clock timing for runtimes.


def import_module_from_file(file_path, verbose=True):
    """
    import_module_from_file
    Dynamically loads and executes a Python module from its absolute filepath. and execute it as a module.

    Input:
        file_path                   str             Absolute path to the Python file to load.
        verbose                     bool            If True, prints errors to console.

    Output:
        ok                          bool            True when the import succeeds, False otherwise.
    """
    try:                                                                                                # Catch import-time errors to keep the batch runner alive.
        script_dir = os.path.dirname(os.path.abspath(file_path))                                        # Get directory of the script.
        if script_dir not in sys.path:                                                                  # Prevent duplicate paths.
            sys.path.insert(0, script_dir)                                                              # Add script's directory to sys.path so it can find its local imports.
            
        module_name = os.path.splitext(os.path.basename(file_path))[0]                                  # Derive a stable module name from filename.
        spec = importlib.util.spec_from_file_location(module_name, file_path)                           # Build import spec from file path.
        if spec is None or spec.loader is None:                                                         # Ensure a valid spec and loader were found.
            raise ImportError(f"Cannot load module from {file_path}")                                   # Raise an exception if spec is missing.
        module = importlib.util.module_from_spec(spec)                                                  # Create module object from spec.
        sys.modules[module_name] = module                                                               # Register module to allow intra-import references.
        spec.loader.exec_module(module)                                                                 # Execute the module (runs top-level batch script).
        return True                                                                                     # Report success.
    except Exception as e:                                                                              # Catch any failure while importing/executing the script.
        if verbose:
            print(f'Error importing {file_path}: {str(e)}')                                             # Print an actionable error message.
        return False                                                                                    # Report failure.


def main(verbose=True):
    """
    main
    Run the predefined list of reference batch scripts under batches/ and print a summary.

    Output:
        None
    """
    run_files = [                                                                                       # List of files to execute in order.
#        'run_Poisson.py',                                                                               # Stationary Poisson reference.
#        'run_Heat.py',                                                                                  # Transient Heat reference.
#        'run_Wave.py',                                                                                  # Transient Wave reference.
#        'run_AdvDif.py',                                                                                # Transient Advection–Diffusion reference.
        'run_Perturbation.py',                                                                          # Stationary perturbation case (Adv=True).
        'run_Perturbation2.py'                                                                          # Stationary perturbation case (Adv=False).
    ]                                                                                                   # End of run list.

    current_dir = os.path.dirname(os.path.abspath(__file__))                                            # Directory where this script is located (repo root).

    if verbose:
        print('Starting mGFD scripts execution...')                                                     # Console header for batch run.
        print('=' * 50)                                                                                 # Visual separator.

    total_start_time = time()                                                                           # Start total timer.
    successful_runs = 0                                                                                 # Count successful script executions.

    for run_file in run_files:                                                                          # Iterate scripts in predefined order.
        file_path = os.path.join(current_dir, 'batches', run_file)                                      # Resolve path under batches/ folder.
        if not os.path.exists(file_path):                                                               # Skip missing files without stopping the suite.
            if verbose:
                print(f'Warning: File {run_file} not found')                                            # Report missing script.
            continue                                                                                    # Move to next script.

        if verbose:
            print(f'\nExecuting {run_file}...')                                                         # Report script execution start.
        start_time = time()                                                                             # Start per-script timer.

        if import_module_from_file(file_path, verbose=verbose):                                         # Import and execute the script module.
            end_time = time()                                                                           # Stop per-script timer.
            execution_time = end_time - start_time                                                      # Compute per-script runtime.
            if verbose:
                print(f'Completed {run_file} in {execution_time:.2f} seconds')                          # Report per-script runtime.
            successful_runs += 1                                                                        # Count successful execution.
        else:                                                                                           # Script import/execution failed.
            if verbose:
                print(f'Error executing {run_file}')                                                    # Report failure while keeping the suite running.

    total_time = time() - total_start_time                                                              # Compute total elapsed time.
    if verbose:
        print('\n' + '=' * 50)                                                                          # Visual separator before summary.
        print(f'Execution completed: {successful_runs} of {len(run_files)} scripts successfully executed') # Print success count.
        print(f'Total execution time: {total_time:.2f} seconds')                                        # Print total runtime.


if __name__ == '__main__':                                                                              # Script entry point.
    main(verbose=True)                                                                                  # Run the batch suite.
