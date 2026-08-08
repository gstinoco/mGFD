"""
IO — Data and neighbor-list I/O helpers

Overview:
    Small utilities to:
    - Locate the project root
    - Load a point cloud from CSV into the (m, 3) convention used by mGFD (x, y, flag)
    - Build/load/save cached neighbor lists vec as CSV
    - Iterate over datasets in the repository Data/ folder structure

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    project_root        Resolve the repository root directory from this file location.
    neighbors_path      Compute the canonical neighbor CSV filename for a given cloud and nvec.
    load_neighbors      Load a cached vec from disk (returns None if unavailable/invalid).
    save_neighbors      Save vec to disk after basic shape validation.
    load_points         Load a point cloud from CSV into the (m, 3) format used by the solvers.
    iter_clouds         Yield (dataset, scale, cloud_path) for available input CSVs under Data/.

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
import os                                                                                               # Filesystem path utilities.
import glob                                                                                             # Globbing for dataset traversal.
import numpy as np                                                                                      # CSV loading and array utilities.

def project_root():
    """
    project_root
    Return the absolute path to the repository root directory.
    
    This is computed by taking the parent directory of the Scripts/ folder (i.e. two levels above
    this file).
    
    Output:
        root                        str             Absolute path to the project root.
    """
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                  # Go two directories up: Scripts/ -> repo root.


def neighbors_path(cloud_path, nvec, tag=None):
    """
    neighbors_path
    Build the filename used to cache neighbor lists for a cloud.
    
    The generated filename is based on cloud_path, inserting a suffix:
        <base>[_<tag>]_neighbors_<nvec><ext>
    
    Input:
        cloud_path                   str            Path to the original cloud CSV file.
        nvec                         int            Number of neighbors per node (vec width).
        tag                          str|None       Optional extra tag to disambiguate variants.
    
    Output:
        path                         str            Path where the neighbor CSV is expected/saved.
    """
    base, ext = os.path.splitext(cloud_path)                                                            # Split input path into base name and extension.
    if tag:                                                                                             # Optional disambiguation tag (e.g. 'adv', 'cloud_exterior').
        return f'{base}_{tag}_neighbors_{int(nvec)}{ext}'                                               # Insert tag and neighbor count before extension.
    
    return f'{base}_neighbors_{int(nvec)}{ext}'                                                         # Insert neighbor count before extension.


def load_neighbors(cloud_path, nvec, tag=None):
    """
    load_neighbors
    Load a cached neighbor list vec from disk.
    
    This function prevents recalculating the neighborhood and distances. If the canonical cache file
    exists, it will be loaded directly, improving performance significantly on repeated runs.
    
    Input:
        cloud_path                   str            Path to the original cloud CSV file.
        nvec                         int            Expected number of neighbors per node.
        tag                          str|None       Optional extra tag to disambiguate variants.
    
    Output:
        vec                          ndarray|None   Loaded neighbor array (m, nvec) if valid, else None.
    
    Notes:
        Returns None if the file does not exist, cannot be parsed, or the loaded array does not
        match the expected number of columns.
    """
    path = neighbors_path(cloud_path, nvec, tag = tag)                                                  # Compute expected cache path.
    if not os.path.exists(path):                                                                        # No cache file found.
        return None                                                                                     # Signal cache miss.
    try:                                                                                                # Try to parse CSV as int32 neighbor indices.
        vec = np.loadtxt(path, delimiter = ',', dtype = np.int32)                                       # Load neighbor matrix (may be 1D if single row).
    except Exception:                                                                                   # Parse or I/O error.
        return None                                                                                     # Treat invalid cache as missing.
    if vec.ndim == 1:                                                                                   # Ensure 2D shape even for single-row files.
        vec = vec.reshape(1, -1)                                                                        # Promote to (1, nvec).
    if vec.shape[1] != int(nvec):                                                                       # Validate expected neighbor count.
        return None                                                                                     # Reject cache with wrong width.
    
    return vec                                                                                          # Return loaded neighbor list.


def save_neighbors(cloud_path, nvec, vec, tag=None):
    """
    save_neighbors
    Save a neighbor list vec to disk as a CSV file.
    
    Input:
        cloud_path                   str            Path to the original cloud CSV file.
        nvec                         int            Number of neighbors per node (vec width).
        vec                          array-like     Neighbor indices (m, nvec).
        tag                          str|None       Optional extra tag to disambiguate variants.
    
    Output:
        path                         str            Path to the saved neighbor CSV file.
    
    Raises:
        ValueError: If vec is not 2D or does not have exactly nvec columns.
    """
    path = neighbors_path(cloud_path, nvec, tag = tag)                                                  # Compute output cache path.
    vec  = np.asarray(vec)                                                                              # Normalize input to ndarray.
    if vec.ndim != 2:                                                                                   # vec must be a matrix (m, nvec).
        raise ValueError('vec must be a 2D array')                                                      # Reject invalid shape.
    if vec.shape[1] != int(nvec):                                                                       # Validate expected neighbor count.
        raise ValueError('vec has incorrect number of columns')                                         # Reject mismatched width.
    np.savetxt(path, vec, delimiter = ',', fmt = '%d')                                                  # Save as comma-separated integers.
    
    return path                                                                                         # Return the saved file path.


def load_points(file_path):
    """
    load_points
    Load a point cloud from a CSV file into the (m, 3) format used by mGFD.
    
    Supported formats:
        - CSV with header columns including x, y, and optionally flag
        - CSV with header columns x, y, and optionally classification (mapped to flag)
        - Headerless CSV with 2 columns (x, y) or 3 columns (x, y, flag)
    
    Input:
        file_path                    str            Path to a CSV file.
    
    Output:
        p                            ndarray         Array with columns [x, y, flag] (float, float, int).
    
    Notes:
        classification mapping (case-insensitive):
            'boundary' -> flag=1
            'hole', 'interior_boundary', 'inner_boundary' -> flag=2
            otherwise -> flag=0
    """
    try:                                                                                                # First attempt: header-aware CSV parsing.
        data = np.genfromtxt(file_path, delimiter = ',', names = True, dtype = None, encoding = 'utf-8')# Load with column names when available.
        if data.dtype.names and 'x' in data.dtype.names and 'y' in data.dtype.names:                    # Require at least x and y columns.
            x = np.asarray(data['x'], dtype = float)                                                    # Extract x coordinates.
            y = np.asarray(data['y'], dtype = float)                                                    # Extract y coordinates.
            if 'flag' in data.dtype.names:                                                              # Direct flag column available.
                flag = np.asarray(data['flag'], dtype = int)                                            # Use provided integer flag.
            elif 'classification' in data.dtype.names:                                                  # Text-based classification column available.
                cls  = np.asarray(data['classification']).astype(str)                                   # Normalize to string array.
                cls  = np.char.lower(cls)                                                               # Lowercase for case-insensitive matching.
                flag = np.zeros(cls.shape, dtype = int)                                                 # Default to interior (0).
                flag[np.isin(cls, ['boundary'])] = 1                                                    # Map boundary class to 1.
                flag[np.isin(cls, ['hole', 'interior_boundary', 'inner_boundary'])] = 2                 # Map interior boundary/hole to 2.
            else:                                                                                       # No flag information present.
                flag = np.zeros(x.shape, dtype = int)                                                   # Default to interior (0).
            return np.column_stack([x, y, flag])                                                        # Return (m,3) point array.
    except Exception:                                                                                   # Any parsing error: fall back to headerless reader.
        pass                                                                                            # Continue to headerless parsing.

    p = np.genfromtxt(file_path, delimiter = ',', skip_header = 0)                                      # Fallback: load numeric CSV without named columns.
    if p.ndim == 1:                                                                                     # Single row case becomes 1D.
        p = np.vstack([p])                                                                              # Promote to (1, ncols).
    if p.shape[1] == 2:                                                                                 # If only x,y provided, append flag=0.
        p = np.column_stack([p, np.zeros(p.shape[0])])                                                  # Add default flag column.
    
    return p                                                                                            # Return parsed point array.


def iter_clouds(data_root=None, scales=('0.50x', '1.00x', '1.50x', '2.00x')):
    """
    iter_clouds
    Iterate over cloud CSV files under the repository Data/ folder structure.
    
    This generator yields tuples:
        (dataset, scale, cloud_path)
    
    Input:
        data_root                    str|None       Root directory to search. Defaults to <project_root>/Data.
        scales                       tuple|list|str Scale folders to traverse (e.g., '1.00x', '2.00x').
    
    Output:
        dataset                      str            Dataset folder name.
        scale                        str            Scale folder name.
        cloud_path                   str            Path to the cloud CSV file.
    
    Notes:
        Neighbor cache files (those containing 'neighbors' in the filename) are skipped.
    """
    if data_root is None:                                                                               # Default to repository Data/ directory.
        data_root = os.path.join(project_root(), 'Data')                                                # Resolve <repo_root>/Data.
    if isinstance(scales, str):                                                                         # Allow passing a single scale string.
        scales = [scales]                                                                               # Normalize to list for iteration.

    for dataset in sorted(os.listdir(data_root)):                                                       # Iterate dataset folders under Data/.
        dataset_path = os.path.join(data_root, dataset)                                                 # Full path to dataset directory.
        if not os.path.isdir(dataset_path):                                                             # Guard against non-directory entries.
            continue                                                                                    # Skip invalid dataset entry.
        for scale in scales:                                                                            # Iterate requested scale folders (e.g., 2x/3x).
            scale_path = os.path.join(dataset_path, scale)                                              # Full path to scale directory.
            if not os.path.isdir(scale_path):                                                           # Skip missing scales.
                continue                                                                                # Skip to next scale.
            csv_paths = sorted(glob.glob(os.path.join(scale_path, '*.csv')))                            # Collect CSV files in scale folder.
            for cloud_path in csv_paths:                                                                # Iterate CSV candidates.
                base = os.path.basename(cloud_path)                                                     # File name for filtering.
                if base.endswith('_cloud.csv'):
                    yield dataset, scale, cloud_path                                                    # Yield tuple for downstream processing.
                else:
                    continue                                                                            # Skip to next file.
