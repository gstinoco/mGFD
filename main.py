"""
Meshless Generalized Finite Differences - Main Execution Script
"""

import os
import importlib.util
import sys
from time import time

def import_module_from_file(file_path):
    try:
        module_name = os.path.splitext(os.path.basename(file_path))[0]
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return True
    except Exception as e:
        print(f'Error importing {file_path}: {str(e)}')
        return False

def main():
    # List of files to execute in order
    run_files = [
        'run_Poisson.py',
        'run_Heat.py',
        'run_Wave.py',
        'run_AdvDif.py',
        'run_Perturbation.py',
        'run_Perturbation2.py'
    ]
    
    # Current directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    print('Starting mGFD scripts execution...')
    print('=' * 50)
    
    total_start_time = time()
    successful_runs = 0
    
    for run_file in run_files:
        file_path = os.path.join(current_dir, run_file)
        if not os.path.exists(file_path):
            print(f'Warning: File {run_file} not found')
            continue
            
        print(f'\nExecuting {run_file}...')
        start_time = time()
        
        if import_module_from_file(file_path):
            end_time = time()
            execution_time = end_time - start_time
            print(f'Completed {run_file} in {execution_time:.2f} seconds')
            successful_runs += 1
        else:
            print(f'Error executing {run_file}')
            
    total_time = time() - total_start_time
    print('\n' + '=' * 50)
    print(f'Execution completed: {successful_runs} of {len(run_files)} scripts successfully executed')
    print(f'Total execution time: {total_time:.2f} seconds')

if __name__ == '__main__':
    main()