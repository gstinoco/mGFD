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

from typing import Optional, Union, Tuple, Iterator                                                                                     # Type hints.

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
    
    try:
        vec = np.loadtxt(path, delimiter=',', dtype=np.int32)                                                                           # Load matrix from CSV.
    except Exception:
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
            continue
        
        for scale in scales:                                                                                                            # Iterate requested scales.
            scale_path = os.path.join(dataset_path, scale)                                                                              # Build full scale path.
            
            if not os.path.isdir(scale_path):                                                                                           # Ignore if scale is missing.
                continue
            csv_paths = sorted(glob.glob(os.path.join(scale_path, '*.csv')))                                                            # Find all CSVs in scale folder.
            
            for cloud_path in csv_paths:                                                                                                # Iterate over CSVs.
                base = os.path.basename(cloud_path)                                                                                     # Get filename.
                
                if base.endswith('_cloud.csv'):                                                                                         # Filter for cloud files.
                    yield dataset, scale, cloud_path                                                                                    # Yield dataset info and path.
