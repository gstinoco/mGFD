"""
Batch Utils — I/O and helper functions for research batches

Overview:
    Utilities designed exclusively for iterating through the experimental datasets
    (Data folder) and caching/retrieving neighbor lists. These functions are intended
    for research benchmarks and are not part of the public mGFD API.

Public API:
    project_root        Get the root directory of the research batches.
    neighbors_path      Build the filename for the neighbor cache.
    load_neighbors      Load cached neighbor lists from disk.
    save_neighbors      Save computed neighbor lists to disk.
    iter_clouds         Iterate through all dataset CSVs for given scales.

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

Date:
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
import os                                                                                                                               # OS interfaces for file/directory paths.
import glob                                                                                                                             # File pattern matching.
import numpy as np                                                                                                                      # Core numerical operations.
import json                                                                                                                             # JSON formatting.

from typing import Optional, Union, Tuple, Iterator, Dict, Any, Callable                                                                # Type hints.

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

import json                                                                                                                             # Import module dependency.
import logging                                                                                                                          # Import module dependency.
from typing import Callable, Any, Dict                                                                                                  # Import module dependency.

logger = logging.getLogger(__name__)                                                                                                    # Assign logger.

def run_batch_suite(
    process_func: Callable,                                                                                                             # Execute statement.
    data_root: str,                                                                                                                     # Execute statement.
    results_root: str,                                                                                                                  # Execute statement.
    scales: Optional[Union[tuple, list, str]] = None,                                                                                   # Assign scales: Optional[Union[tuple, list, str]].
    **kwargs: Any
) -> None:                                                                                                                              # Execute statement.
    """
    run_batch_suite
    Universal orchestrator for all runner batch scripts.
    Iterates all point clouds in the datasets and invokes the solver logic.

    Input:
        process_func    1           Callable        The function to process each point cloud.
        data_root       1           str             Input dataset root directory.
        results_root    1           str             Output results root directory.
        scales          1           tuple/list/str  Tuple or list of scale subfolders to process (optional).
        **kwargs        Any         dict            Dynamic configuration from main orchestrator.
        
    Output:
        None
    """
    if scales is None:                                                                                                                  # Evaluate condition.
        scales = kwargs.pop('scales', None)                                                                                             # Assign scales.
    if scales is None:                                                                                                                  # Evaluate condition.
        config_path = os.path.join(project_root(), 'codes', 'sweep_config.json')                                                        # Assign config_path.
        if os.path.exists(config_path):                                                                                                 # Evaluate condition.
            try:                                                                                                                        # Execute statement.
                with open(config_path, 'r') as f:                                                                                       # Execute statement.
                    cfg = json.load(f)                                                                                                  # Assign cfg.
                    scales = cfg.get('scales', ('1', '2', '3'))                                                                         # Assign scales.
            except Exception:                                                                                                           # Execute statement.
                scales = ('1', '2', '3')                                                                                                # Assign scales.
        else:                                                                                                                           # Execute fallback branch.
            scales = ('1', '2', '3')                                                                                                    # Assign scales.

    if isinstance(scales, (int, str)):                                                                                                  # Evaluate condition.
        scales = (str(scales),)                                                                                                         # Assign scales.
    else:                                                                                                                               # Execute fallback branch.
        scales = tuple(str(s) for s in scales)                                                                                          # Iterate over collection.

    verbose = kwargs.pop('verbose', True)                                                                                               # Assign verbose.
    save = kwargs.pop('save', True)                                                                                                     # Assign save.
    found = 0                                                                                                                           # Assign found.
    
    if verbose:                                                                                                                         # Evaluate condition.
        logger.info(f'Processing point clouds from {data_root} (scales={scales}).')                                                     # Assign logger.info(f'Processing point clouds from {data_root} (scales.
        
    for dataset, scale, cloud_path in iter_clouds(data_root, scales):                                                                   # Iterate over collection.
        found += 1                                                                                                                      # Assign found +.
        process_func(dataset, scale, cloud_path, results_root, save, verbose=verbose, **kwargs)                                         # Assign process_func(dataset, scale, cloud_path, results_root, save, verbose.
        
    if found == 0:                                                                                                                      # Evaluate condition.
        if verbose:                                                                                                                     # Evaluate condition.
            logger.warning(f'No point clouds found in {data_root} for scales {scales}.')                                                # Iterate over collection.

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
        if 'Time_Callable' not in metrics and 'Time_Array' not in metrics and 'Time_Pandas' not in metrics:                             # If no type-specific key was created.
            metrics['Time_Callable'] = exec_time                                                                                        # Fallback to Callable label for backward compatibility.
    elif 'Time_Callable' in metrics:                                                                                                    # If Callable time is present.
        metrics['Time_Secs'] = metrics['Time_Callable']                                                                                 # Mirror to unified execution time metric.
    elif 'Time_Array' in metrics:                                                                                                       # If Array time is present.
        metrics['Time_Secs'] = metrics['Time_Array']                                                                                    # Mirror to unified execution time metric.
    elif 'Time_Pandas' in metrics:                                                                                                      # If Pandas time is present.
        metrics['Time_Secs'] = metrics['Time_Pandas']                                                                                   # Mirror to unified execution time metric.
        
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

import pandas as pd                                                                                                                     # Import module dependency.
from datetime import datetime                                                                                                           # Import module dependency.

def generate_sweep_summary(results_root: str, verbose: bool = True) -> None:
    """
    generate_sweep_summary
    Crawls all Metrics.json files in the results directory and generates a master CSV summary.
    This allows for rapid comparison of GPU vs CPU times and solver combinations.

    Input:
        results_root    1           str             Output results root directory.
        verbose         1           bool            If True, prints progress.
        
    Output:
        None
    """
    all_records = []                                                                                                                    # Initialize list for all records.
    
    if not os.path.exists(results_root):                                                                                                # Check if results root exists.
        if verbose:                                                                                                                     # Check verbosity.
            logger.warning(f"Results root {results_root} does not exist. Cannot generate summary.")                                     # Log warning.
        return                                                                                                                          # Exit function.
        
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
                        'Equation': equation,                                                                                           # Store equation name.
                        'Dataset': dataset,                                                                                             # Store dataset name.
                        'Scale': scale,                                                                                                 # Store scale string.
                        'Config_ID': config_id                                                                                          # Store configuration ID.
                    }                                                                                                                   # End record initialization.
                    
                    for k, v in metrics.items():                                                                                        # Iterate through metrics.
                        record[k] = v                                                                                                   # Append metric to record.
                        
                    all_records.append(record)                                                                                          # Append record to list.
                    
    if not all_records:                                                                                                                 # Check if records list is empty.
        if verbose:                                                                                                                     # Check verbosity.
            logger.warning("No metrics found to summarize.")                                                                            # Log warning.
        return                                                                                                                          # Exit function.
        
    df = pd.DataFrame(all_records)                                                                                                      # Convert records to DataFrame.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")                                                                                # Generate timestamp string.
    out_file = os.path.join(results_root, f"sweep_summary_{timestamp}.csv")                                                             # Build output file path.
    df.to_csv(out_file, index=False)                                                                                                    # Export DataFrame to CSV.
    
    if verbose:                                                                                                                         # Check verbosity.
        logger.info(f"\nSuccessfully generated master sweep summary CSV at:")                                                           # Log success message.
        logger.info(f" -> {out_file}")                                                                                                  # Log output path.
        logger.info(f"Total aggregated configurations: {len(df)}")                                                                      # Log total records.
