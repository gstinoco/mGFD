"""
IO — Data and neighbor-list I/O helpers

Overview:
    Small utilities to:
    - Load a point cloud from CSV into the (m, 3) convention used by mGFD (x, y, flag)

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.

Public API:
    load_points         Load a point cloud from CSV into the (m, 3) format used by the solvers.

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
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.

def load_points(file_path: str, verbose: bool = False) -> np.ndarray:
    """
    load_points
    Load a point cloud from a CSV file into the (m, 3) format used by mGFD.
    
    Supported formats:
        - CSV with header columns including x, y, and optionally flag
        - CSV with header columns x, y, and optionally classification (mapped to flag)
        - Headerless CSV with 2 columns (x, y) or 3 columns (x, y, flag)
    
    Input:
        file_path                   str             Path to a CSV file.
        verbose                     bool            If True, prints confirmation of loaded points.
    
    Output:
        p           m x 3           ndarray         Array with columns [x, y, flag] (float, float, int).
    """
    # 1. Header parsing
    try:                                                                                                                                # Attempt loading with headers.
        data = np.genfromtxt(file_path, delimiter = ',', names = True, dtype = None, encoding = 'utf-8')                                # Load CSV data with headers.
        
        if data.dtype.names and 'x' in data.dtype.names and 'y' in data.dtype.names:                                                    # Check for required coordinate headers.
            x = np.asarray(data['x'], dtype = float)                                                                                    # Extract X coordinates.
            y = np.asarray(data['y'], dtype = float)                                                                                    # Extract Y coordinates.
            
            if 'flag' in data.dtype.names:                                                                                              # Check for explicit boundary flags.
                flag = np.asarray(data['flag'], dtype = int)                                                                            # Extract boundary flags.
            elif 'classification' in data.dtype.names:                                                                                  # Check for string classifications.
                cls  = np.asarray(data['classification']).astype(str)                                                                   # Extract classifications as strings.
                cls  = np.char.lower(cls)                                                                                               # Convert to lowercase for comparison.
                flag = np.zeros(cls.shape, dtype = int)                                                                                 # Initialize integer flags array.
                flag[np.isin(cls, ['boundary'])] = 1                                                                                    # Map outer boundaries to flag 1.
                flag[np.isin(cls, ['hole', 'interior_boundary', 'inner_boundary'])] = 2                                                 # Map internal boundaries to flag 2.
            else:                                                                                                                       # Missing flag or classification.
                flag = np.zeros(x.shape, dtype = int)                                                                                   # Default to interior points.
                
            if verbose:                                                                                                                 # Check if verbosity is enabled.
                print(f"Loaded {x.shape[0]} points from {file_path}")                                                                   # Print load confirmation.
                
            return np.column_stack([x, y, flag])                                                                                        # Return loaded point cloud.
            
    except Exception:                                                                                                                   # Fallback if headers fail.
        pass                                                                                                                            # Do nothing.

    # 2. Headerless parsing fallback
    p = np.genfromtxt(file_path, delimiter = ',', skip_header = 0)                                                                      # Load CSV data as raw matrix.
    
    if p.ndim == 1:                                                                                                                     # Check if loaded data is 1D.
        p = np.vstack([p])                                                                                                              # Reshape into a 2D array.
        
    if p.shape[1] == 2:                                                                                                                 # Check if only 2 columns exist.
        p = np.column_stack([p, np.zeros(p.shape[0])])                                                                                  # Append default flags column.
    
    # 3. Load confirmation and return
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        print(f"Loaded {p.shape[0]} points from {file_path}")                                                                           # Print load confirmation.
    
    return p                                                                                                                            # Return loaded point cloud.
