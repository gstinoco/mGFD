"""
Batch Utils — I/O and helper functions for research batches

Overview:
    Utilities designed exclusively for iterating through the experimental datasets
    (Data folder), caching/retrieving neighbor lists, executing dynamic batch runners,
    and aggregating parameter sweep benchmarks across the mGFD research laboratory.

Public API:
    project_root                Get the root directory of the research batches.
    neighbors_path              Build the filename for storing the neighbor cache.
    load_neighbors              Load cached neighbor lists from disk.
    save_neighbors              Save computed neighbor lists to disk.
    import_module_from_file     Dynamically load and execute a Python runner module.
    run_job                     Worker execution wrapper for process pool parallelization.
    iter_clouds                 Iterate through all dataset CSVs for given scales.
    run_batch_suite             Universal orchestrator for all runner batch scripts.
    save_metrics                Record and serialize execution time and numerical error metrics.
    generate_sweep_summary      Aggregate individual batch JSON metric files into a master CSV.

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
import glob                                                                                                                             # File pattern matching.
import json                                                                                                                             # JSON formatting.
import logging                                                                                                                          # Standard logging module.
import importlib.util                                                                                                                   # Dynamic module import from a file path.

import pandas as pd                                                                                                                     # Tabular data manipulation.
import numpy as np                                                                                                                      # Core numerical operations.

from time import time                                                                                                                   # High-resolution wall-clock timer.
from datetime import datetime                                                                                                           # Datetime stamp generation.
from typing import Optional, Union, Tuple, Iterator, Dict, Any, Callable                                                                # Type hints.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))                                               # Resolve repository root directory.
_src_dir   = os.path.join(_repo_root, 'src')                                                                                            # Resolve core library source directory.
if _src_dir not in sys.path:                                                                                                            # Ensure src directory is on system path.
    sys.path.insert(0, _src_dir)                                                                                                        # Add src directory to system path.


def project_root() -> str:
    """
    project_root
    Return the absolute path to the research directory.

    Input:
        None

    Output:
        path            1           str             Absolute path to the research folder.
    """
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                                  # Compute root from file location.

def neighbors_path(cloud_path: str, nvec: int, tag: Optional[str] = None) -> str:
    """
    neighbors_path
    Build the filename for storing the neighbor cache.

    Input:
        cloud_path      1           str             Path to the original point cloud CSV.
        nvec            1           int             Number of neighbors.
        tag             1           str             (Optional) Tag to insert in the filename.
        
    Output:
        path            1           str             Constructed path for the neighbor cache.
    """
    if not isinstance(cloud_path, str):                                                                                                 # Validate cloud_path argument.
        raise TypeError("cloud_path must be a string.")                                                                                 # Raise explicit error on bad input.
    
    if not isinstance(nvec, int):                                                                                                       # Validate nvec argument.
        raise TypeError("nvec must be an integer.")                                                                                     # Raise explicit error on bad input.

    base, ext = os.path.splitext(cloud_path)                                                                                            # Split filename and extension.
    
    if tag:                                                                                                                             # If a tag is provided.
        return f'{base}_{tag}_neighbors_{nvec}{ext}'                                                                                    # Return path with tag.
    
    return f'{base}_neighbors_{nvec}{ext}'                                                                                              # Return path without tag.

def load_neighbors(cloud_path: str, nvec: int, tag: Optional[str] = None) -> Optional[np.ndarray]:
    """
    load_neighbors
    Load cached neighbor lists from disk if they exist.

    Input:
        cloud_path      1           str             Path to the original point cloud CSV.
        nvec            1           int             Number of neighbors expected.
        tag             1           str             (Optional) Tag for the cache filename.
        
    Output:
        vec             m x nvec    ndarray         Array containing the neighbor indices, or None if not found.
    """
    if not isinstance(cloud_path, str):                                                                                                 # Validate cloud_path argument.
        raise TypeError("cloud_path must be a string.")                                                                                 # Raise explicit error on bad input.
    
    if not isinstance(nvec, int):                                                                                                       # Validate nvec argument.
        raise TypeError("nvec must be an integer.")                                                                                     # Raise explicit error on bad input.

    path = neighbors_path(cloud_path, nvec, tag=tag)                                                                                    # Get cache path.
    
    if not os.path.exists(path):                                                                                                        # Check if cache exists.
        return None                                                                                                                     # Return None if it doesn't.
    
    try:                                                                                                                                # Execute statement.
        vec = np.loadtxt(path, delimiter=',', dtype=np.int32)                                                                           # Load matrix from CSV.
    except Exception:                                                                                                                   # Execute statement.
        return None                                                                                                                     # Return None on parse failure.
    
    if vec.ndim == 1:                                                                                                                   # Handle 1D edge case.
        vec = vec.reshape(1, -1)                                                                                                        # Reshape to 2D.
    
    if vec.shape[1] != nvec:                                                                                                            # Validate dimensions.
        return None                                                                                                                     # Return None if mismatch.
    
    return vec                                                                                                                          # Return cached matrix.

def import_module_from_file(file_path: str, verbose: bool = True, **kwargs: Any) -> bool:
    """
    Dynamically load and execute a Python module from its absolute filepath.

    Loads the module at runtime using importlib, adds its directory to sys.path,
    registers it in sys.modules, executes its top-level code, and invokes its main()
    entry point with injected keyword arguments.

    Input:
        file_path       1           str             Absolute file path of the runner script.
        verbose         1           bool            If True, logs warnings and errors.
        **kwargs                    Any             Optional arguments forwarded directly to the module's main().

    Output:
        success         1           bool            True if execution succeeded without exception, False otherwise.
    """
    if not isinstance(file_path, str):                                                                                                  # Validate file_path argument.
        raise TypeError("file_path must be a string.")                                                                                  # Explicit error on bad input.
    if not os.path.exists(file_path):                                                                                                   # Validate file existence.
        raise ValueError(f"File not found: {file_path}")                                                                                # Explicit error on bad input.
    try:                                                                                                                                # Catch import-time errors to keep the batch runner alive.
        script_dir = os.path.dirname(os.path.abspath(file_path))                                                                        # Get directory of the script.
        if script_dir not in sys.path:                                                                                                  # Prevent duplicate paths.
            sys.path.insert(0, script_dir)                                                                                              # Add script's directory to sys.path so it can find its local imports.
        src_dir = os.path.join(os.path.dirname(project_root()), 'src')                                                                  # Resolve core library source directory.
        if src_dir not in sys.path:                                                                                                     # Ensure src directory is on system path.
            sys.path.insert(0, src_dir)                                                                                                 # Add src directory to system path.
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
            import traceback                                                                                                            # Import traceback for diagnostic dump.
            logger.error(f'Error importing {file_path}: {traceback.format_exc()}')                                                      # Print an actionable error message.
        return False                                                                                                                    # Report failure.

def run_job(job_tuple: tuple) -> bool:
    """
    Worker execution wrapper for process pool parallelization.

    Unpacks the runner configuration tuple, measures per-script wall-clock execution
    time, and safely invokes import_module_from_file in a subprocess worker.

    Input:
        job_tuple       6           tuple           Tuple containing (file_path, run_file, nvec, upwind, base_kwargs, verbose).

    Output:
        success         1           bool            True if the runner finished successfully, False on error.
    """
    file_path, run_file, nvec, upwind, base_kwargs, verbose = job_tuple                                                                 # Unpack job configuration tuple.
    if upwind is not None:                                                                                                              # Determine if upwind applies.
        config_str = f'nvec={nvec}, dev={base_kwargs.get("device")}, up={upwind}'                                                       # Format config with upwind.
    else:                                                                                                                               # Upwind doesn't apply.
        config_str = f'nvec={nvec}, dev={base_kwargs.get("device")}'                                                                    # Format config without upwind.
    if verbose:                                                                                                                         # Check verbosity.
        logger.info(f'\nExecuting {run_file} with [{config_str}]...')                                                                   # Log current execution.
    
    start_time         = time()                                                                                                         # Start per-script timer.
    run_kwargs         = base_kwargs.copy()                                                                                             # Copy global kwargs to avoid side-effects.
    run_kwargs['nvec'] = nvec                                                                                                           # Inject neighbor count into kwargs.
    
    if upwind is not None:                                                                                                              # If upwind is required for this runner.
        run_kwargs['upwind'] = upwind                                                                                                   # Inject upwind flag into kwargs.
    
    ok             = import_module_from_file(file_path, verbose=verbose, **run_kwargs)                                                  # Dynamically load and run script with kwargs.
    execution_time = time() - start_time                                                                                                # Compute per-script runtime.
    
    if ok:                                                                                                                              # Check execution success.
        if verbose:                                                                                                                     # Check verbosity.
            logger.info(f'Completed {run_file} [{config_str}] in {execution_time:.2f} seconds')                                         # Log completion.
    else:                                                                                                                               # Script execution failed.
        if verbose:                                                                                                                     # Check verbosity.
            logger.error(f'Error executing {run_file} [{config_str}]')                                                                  # Report failure.
    return ok                                                                                                                           # Return execution result.

def save_neighbors(cloud_path: str, nvec: int, vec: Union[np.ndarray, list], tag: Optional[str] = None) -> str:
    """
    save_neighbors
    Save computed neighbor lists to disk.

    Input:
        cloud_path      1           str             Path to the original point cloud CSV.
        nvec            1           int             Number of neighbors.
        vec             m x nvec    ndarray         Matrix of neighbor indices to save.
        tag             1           str             (Optional) Tag for the cache filename.
        
    Output:
        path            1           str             The path where the cache was saved.
    """
    if not isinstance(cloud_path, str):                                                                                                 # Validate cloud_path argument.
        raise TypeError("cloud_path must be a string.")                                                                                 # Raise explicit error on bad input.
    
    if not isinstance(nvec, int):                                                                                                       # Validate nvec argument.
        raise TypeError("nvec must be an integer.")                                                                                     # Raise explicit error on bad input.
    
    if not isinstance(vec, (np.ndarray, list)):                                                                                         # Validate vec argument.
        raise TypeError("vec must be an array-like object.")                                                                            # Raise explicit error on bad input.

    path = neighbors_path(cloud_path, nvec, tag=tag)                                                                                    # Get cache path.
    vec_arr = np.asarray(vec)                                                                                                           # Ensure vec is an array.
    
    if vec_arr.ndim != 2:                                                                                                               # Validate 2D structure.
        raise ValueError('vec must be a 2D array')                                                                                      # Error if not 2D.
    
    if vec_arr.shape[1] != nvec:                                                                                                        # Validate correct columns.
        raise ValueError('vec has incorrect number of columns')                                                                         # Error if mismatch.
    
    np.savetxt(path, vec_arr, delimiter=',', fmt='%d')                                                                                  # Save as integer CSV.
    
    return path                                                                                                                         # Return save location.

def iter_clouds(data_root: Optional[str] = None, scales: Union[Tuple[str, ...], str, list] = ('1', '2', '3', '4', '5')) -> Iterator[Tuple[str, str, str]]:
    """
    iter_clouds
    Iterate through all dataset CSVs for given scales.

    Input:
        data_root       1           str             Root folder for the Data/. Defaults to research/Data.
        scales          1           tuple           Tuple of scale subfolders to process.
        
    Output:
        dataset         1           str             The name of the dataset (e.g. Zempoala).
        scale           1           str             The scale/density (e.g. 1).
        cloud_path      1           str             Absolute path to the point cloud CSV.
    """
    if data_root is None:                                                                                                               # Default to research/Data if None.
        data_root = os.path.join(project_root(), 'Data')                                                                                # Build default path.
    
    if isinstance(scales, str):                                                                                                         # Convert string to list if necessary.
        scales = [scales]                                                                                                               # Wrap in list.
    
    elif not isinstance(scales, (list, tuple)):                                                                                         # Validate scales argument.
        raise TypeError("scales must be a string, list, or tuple.")                                                                     # Raise explicit error on bad input.

    for dataset in sorted(os.listdir(data_root)):                                                                                       # Iterate dataset folders.
        dataset_path = os.path.join(data_root, dataset)                                                                                 # Build full dataset path.
        
        if not os.path.isdir(dataset_path):                                                                                             # Ignore non-directories.
            continue                                                                                                                    # Execute statement.
        
        for scale in scales:                                                                                                            # Iterate requested scales.
            scale_path = os.path.join(dataset_path, scale)                                                                              # Build full scale path.
            
            if not os.path.isdir(scale_path):                                                                                           # Ignore if scale is missing.
                continue                                                                                                                # Execute statement.
            csv_paths = sorted(glob.glob(os.path.join(scale_path, '*.csv')))                                                            # Find all CSVs in scale folder.
            
            for cloud_path in csv_paths:                                                                                                # Iterate over CSVs.
                base = os.path.basename(cloud_path)                                                                                     # Get filename.
                
                if base.endswith('_cloud.csv'):                                                                                         # Filter for cloud files.
                    yield dataset, scale, cloud_path                                                                                    # Yield dataset info and path.

def run_batch_suite(
    process_func: Callable,                                                                                                             # Point cloud processing function.
    data_root:    str,                                                                                                                  # Input dataset root directory.
    results_root: str,                                                                                                                  # Output results root directory.
    scales:       Optional[Union[tuple, list, str]] = None,                                                                             # Target refinement scales.
    **kwargs:     Any                                                                                                                   # Dynamic keyword arguments.
) -> None:                                                                                                                              # Return None.
    """
    Universal orchestrator for all runner batch scripts.

    Iterates over point clouds from data_root for given scales and processes each with process_func.

    Input:
        process_func    1           Callable        The function to process each point cloud.
        data_root       1           str             Input dataset root directory.
        results_root    1           str             Output results root directory.
        scales          1           tuple/list/str  Tuple or list of scale subfolders to process (optional).
        **kwargs                    Any             Dynamic configuration forwarded from main orchestrator.
        
    Output:
        None
    """
    # 1. Scale configuration resolution
    if scales is None:                                                                                                                  # Evaluate condition.
        scales = kwargs.pop('scales', None)                                                                                             # Extract scales from kwargs.
    if scales is None:                                                                                                                  # If still None, inspect sweep config.
        config_path = os.path.join(project_root(), 'codes', 'sweep_config.json')                                                        # Resolve config path.
        if os.path.exists(config_path):                                                                                                 # Check config existence.
            try:                                                                                                                        # Try parsing config.
                with open(config_path, 'r') as f:                                                                                       # Open configuration file.
                    cfg    = json.load(f)                                                                                               # Load JSON.
                    scales = cfg.get('scales', ('1', '2', '3', '4'))                                                                    # Extract scales list.
            except Exception:                                                                                                           # Handle JSON read error.
                scales = ('1', '2', '3', '4')                                                                                           # Fallback to scales 1-4.
        else:                                                                                                                           # Config file not present.
            scales = ('1', '2', '3', '4')                                                                                               # Default to scales 1-4.

    if isinstance(scales, (int, str)):                                                                                                  # Single scale scalar.
        scales = (str(scales),)                                                                                                         # Wrap in tuple.
    else:                                                                                                                               # Multiple scales collection.
        scales = tuple(str(s) for s in scales)                                                                                          # Cast all elements to str.

    verbose = kwargs.pop('verbose', True)                                                                                               # Verbosity flag.
    save    = kwargs.pop('save', True)                                                                                                  # Save outputs flag.
    found   = 0                                                                                                                         # Match counter.
    
    # 2. Dataset iteration and execution
    if verbose:                                                                                                                         # Evaluate verbosity condition.
        logger.info(f'Processing point clouds from {data_root} (scales={scales}).')                                                     # Log active processing batch.
        
    for dataset, scale, cloud_path in iter_clouds(data_root, scales):                                                                   # Iterate through point clouds.
        found += 1                                                                                                                      # Increment match counter.
        process_func(dataset, scale, cloud_path, results_root, save, verbose=verbose, **kwargs)                                         # Execute cloud processing function.
        
    # 3. Completion check
    if found == 0:                                                                                                                      # If no matching clouds were located.
        if verbose:                                                                                                                     # Check verbosity.
            logger.warning(f'No point clouds found in {data_root} for scales {scales}.')                                                # Warn user.

def save_metrics(out_dir: str, metrics: Dict[str, float], config_id: Optional[str] = None, scale: Optional[str] = None, p: Optional[np.ndarray] = None) -> None:
    """
    save_metrics
    Dumps the error metrics dictionary into a JSON file in the output directory.
    If config_id is provided, appends the metrics to the JSON file under that key.

    Input:
        out_dir         1           str             Directory to save the metrics in.
        metrics         1           dict            Dictionary with calculated errors/times.
        config_id       1           str             Optional key for parameter sweeps.
        scale           1           str             Optional scale identifier.
        p               m x 3       ndarray         Optional point cloud to extract statistics.
        
    Output:
        None
    """
    exec_time = metrics.pop('Compute_Time_Secs', None)                                                                                  # Extract primary solver execution time if present.
    if exec_time is not None:                                                                                                           # If primary compute time present.
        metrics['Time_Secs'] = exec_time                                                                                                # Set standardized execution time metric.
        
    metrics.pop('Max_Abs_Error', None)                                                                                                  # Standardize metrics.
    metrics.pop('Mean_Abs_Error', None)                                                                                                 # Standardize metrics.
    if 'Time_Mean_Abs_Error' in metrics:                                                                                                # Standardize metrics.
        metrics.pop('Time_Mean_Abs_Error', None)                                                                                        # Standardize metrics.
        metrics.pop('Time_Max_Abs_Error', None)                                                                                         # Standardize metrics.
    os.makedirs(out_dir, exist_ok=True)                                                                                                 # Ensure output directory exists.
    metrics_path = os.path.join(out_dir, 'Metrics.json')                                                                                # Define the absolute path for the metrics file.
    if config_id:                                                                                                                       # Check if a configuration ID was provided.
        data = {}                                                                                                                       # Initialize empty dictionary for sweep data.
        if os.path.exists(metrics_path):                                                                                                # Check if the JSON file already exists.
            with open(metrics_path, 'r') as f:                                                                                          # Open the existing file in read mode.
                try:                                                                                                                    # Attempt to load the JSON content.
                    data = json.load(f)                                                                                                 # Load data into the dictionary.
                except json.JSONDecodeError:                                                                                            # Catch decoding errors (e.g., empty file).
                    pass                                                                                                                # Ignore error and start fresh.
        if scale:                                                                                                                       # If scale is provided, nest it.
            scale_key = f'Scale_{scale}'                                                                                                # Format scale key.
            if scale_key not in data:                                                                                                   # Create scale if missing.
                data[scale_key] = {}                                                                                                    # Initialize empty dictionary.
            if p is not None and "Cloud_Statistics" not in data[scale_key]:                                                             # Add statistics.
                data[scale_key]["Cloud_Statistics"] = {                                                                                 # Statistics block.
                    "Total_Nodes": len(p),                                                                                              # Total count.
                    "Interior_Nodes": int(np.sum(p[:, 2] == 0)),                                                                        # Interior count.
                    "Boundary_Nodes": int(np.sum(p[:, 2] != 0)),                                                                        # Boundary count.
                    "X_Min": float(np.min(p[:, 0])),                                                                                    # Min X.
                    "X_Max": float(np.max(p[:, 0])),                                                                                    # Max X.
                    "Y_Min": float(np.min(p[:, 1])),                                                                                    # Min Y.
                    "Y_Max": float(np.max(p[:, 1]))                                                                                     # Max Y.
                }                                                                                                                       # End block.
            if "Metrics" not in data[scale_key]:                                                                                        # Initialize metrics block.
                data[scale_key]["Metrics"] = {}                                                                                         # Initialize empty dict.
            data[scale_key]["Metrics"][config_id] = metrics                                                                             # Add metrics.
        else:                                                                                                                           # Legacy mode.
            data[config_id] = metrics                                                                                                   # Add metrics.
        with open(metrics_path, 'w') as f:                                                                                              # Open the file in write mode.
            json.dump(data, f, indent=4)                                                                                                # Dump the updated data to the file.
    else:                                                                                                                               # Legacy mode (no config_id).
        with open(metrics_path, 'w') as f:                                                                                              # Open the file in write mode.
            json.dump(metrics, f, indent=4)                                                                                             # Dump the single metrics dictionary.

def generate_sweep_summary(results_root: str, verbose: bool = True) -> None:
    """
    Crawls all Metrics.json files in the results directory and generates a master CSV summary.

    This allows for rapid comparison of GPU vs CPU times and solver combinations.

    Input:
        results_root    1           str             Output results root directory.
        verbose         1           bool            If True, prints progress.
        
    Output:
        None
    """
    # 1. Results root verification
    all_records = []                                                                                                                    # Initialize list for all records.
    
    if not os.path.exists(results_root):                                                                                                # Check if results root exists.
        if verbose:                                                                                                                     # Check verbosity.
            logger.warning(f"Results root {results_root} does not exist. Cannot generate summary.")                                     # Log warning.
        return                                                                                                                          # Exit function.
        
    # 2. Metric crawling across equations, regions, and scales
    for equation in os.listdir(results_root):                                                                                           # Iterate through equations.
        eq_path = os.path.join(results_root, equation)                                                                                  # Build equation path.
        if not os.path.isdir(eq_path) or equation == 'benchmarks':                                                                      # Check if it is a valid directory.
            continue                                                                                                                    # Skip if not.
            
        for dataset in os.listdir(eq_path):                                                                                             # Iterate through datasets.
            dataset_path = os.path.join(eq_path, dataset)                                                                               # Build dataset path.
            if not os.path.isdir(dataset_path):                                                                                         # Check if it is a directory.
                continue                                                                                                                # Skip if not.
                
            metrics_file = os.path.join(dataset_path, 'Metrics.json')                                                                   # Build metrics file path.
            if not os.path.exists(metrics_file):                                                                                        # Check if metrics file exists.
                continue                                                                                                                # Skip if not.
                
            try:                                                                                                                        # Try to open and parse JSON.
                with open(metrics_file, 'r') as f:                                                                                      # Open file in read mode.
                    data = json.load(f)                                                                                                 # Parse JSON data.
            except Exception as e:                                                                                                      # Catch any parsing errors.
                if verbose:                                                                                                             # Check verbosity.
                    logger.error(f"Failed to read {metrics_file}: {e}")                                                                 # Log error.
                continue                                                                                                                # Skip to next file.
                
            for scale_key, scale_data in data.items():                                                                                  # Iterate through scales.
                if not scale_key.startswith('Scale_'):                                                                                  # Validate scale key.
                    continue                                                                                                            # Skip if invalid.
                scale = scale_key.replace('Scale_', '')                                                                                 # Extract scale number.
                
                metrics_dict = scale_data.get('Metrics', {})                                                                            # Get metrics dictionary.
                for config_id, metrics in metrics_dict.items():                                                                         # Iterate through configurations.
                    record = {                                                                                                          # Initialize record dictionary.
                        'Equation':  equation,                                                                                          # Store equation name.
                        'Dataset':   dataset,                                                                                           # Store dataset name.
                        'Scale':     scale,                                                                                             # Store scale string.
                        'Config_ID': config_id                                                                                          # Store configuration ID.
                    }                                                                                                                   # End record initialization.
                    
                    for k, v in metrics.items():                                                                                        # Iterate through metrics.
                        record[k] = v                                                                                                   # Append metric to record.
                        
                    all_records.append(record)                                                                                          # Append record to list.
                    
    # 3. Export to consolidated CSV
    if not all_records:                                                                                                                 # Check if records list is empty.
        if verbose:                                                                                                                     # Check verbosity.
            logger.warning("No metrics found to summarize.")                                                                            # Log warning.
        return                                                                                                                          # Exit function.
        
    df        = pd.DataFrame(all_records)                                                                                               # Convert records to DataFrame.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")                                                                                # Generate timestamp string.
    out_file  = os.path.join(results_root, f"sweep_summary_{timestamp}.csv")                                                            # Build output file path.
    df.to_csv(out_file, index=False)                                                                                                    # Export DataFrame to CSV.
    
    if verbose:                                                                                                                         # Check verbosity.
        logger.info(f"\nSuccessfully generated master sweep summary CSV at:")                                                           # Log success message.
        logger.info(f" -> {out_file}")                                                                                                  # Log output path.
        logger.info(f"Total aggregated configurations: {len(df)}")                                                                      # Log total records.
