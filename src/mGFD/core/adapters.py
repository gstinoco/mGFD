"""
Adapters — Data structure format integrators

Overview:
    Data structure adapters for integrating pandas and xarray transparently into the mGFD ecosystem.
    Provides logic to seamlessly wrap and unwrap native industry formats into C-contiguous NumPy arrays.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    extract_cloud       Extracts NumPy arrays from pandas DataFrames or xarray DataArrays.
    repack_solution     Packs NumPy solutions back into the original data format.

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
    August, 2026.
Last Modification:
    August, 2026.
"""

## Library importation.
import sys                                                                                                                              # System-specific parameters and functions.
import numpy as np                                                                                                                      # Core numerical operations.

def extract_cloud(p: np.ndarray) -> np.ndarray:
    """
    Safely extracts NumPy arrays from pandas DataFrames or xarray DataArrays using a soft-dependency (duck typing) approach.
    
    This adapter acts as the ingestion layer for the mGFD suite. It ensures that the core mathematical 
    components (Sparse Matrices, Krylov Solvers, Numba JIT loops) can operate at maximum C-contiguous 
    performance using raw NumPy arrays, without requiring heavy libraries like pandas or xarray to be 
    added as hard dependencies in the project environment.
    
    Input:
        p           m x 3           ndarray         The point cloud input, potentially a DataFrame or DataArray.
        
    Output:
        p           m x 3           ndarray         The NumPy representation of the point cloud.
    """
    
    # 1. Check for Pandas DataFrame
    if "pandas" in sys.modules:                                                                                                         # If pandas is currently loaded in memory.
        import pandas as pd                                                                                                             # Import pandas locally to avoid global dependency.
        if isinstance(p, (pd.DataFrame, pd.Series)):                                                                                    # If input is a pandas DataFrame or Series.
            return np.asarray(p.values)                                                                                                 # Extract and return the underlying numpy array.

    # 2. Check for Xarray DataArray
    if "xarray" in sys.modules:                                                                                                         # If xarray is currently loaded in memory.
        import xarray as xr                                                                                                             # Import xarray locally to avoid global dependency.
        if isinstance(p, xr.DataArray):                                                                                                 # If input is an xarray DataArray.
            return np.asarray(p.values)                                                                                                 # Extract and return the underlying numpy array.

    return np.asarray(p)                                                                                                                # Return default numpy array fallback.

def repack_solution(p_orig: np.ndarray, u_ap: np.ndarray):
    """
    Packs the computed NumPy solution back into the original data format (DataFrame or DataArray) if applicable.
    
    Once the numerical PDE approximation is finished, this function intercepts the raw NumPy solution 
    and appends it as new columns or variables into a copy of the original dataset structure.
    For stationary solvers, it attaches a single `u_ap` column. For transient solvers, it dynamically 
    generates a new column/variable for each computed time step (e.g., `u_ap_t0`, `u_ap_t1`, ...).
    
    Input:
        p_orig      m x 3           ndarray         The original point cloud format, potentially a DataFrame or DataArray.
        u_ap        m x t           ndarray         The computed solution array generated by the solver.
        
    Output:
        packed                      Any             The original packed object (or raw ndarray) with the solution appended.
    """
    
    # 1. Pack into Pandas DataFrame
    if "pandas" in sys.modules:                                                                                                         # If pandas is currently loaded in memory.
        import pandas as pd                                                                                                             # Import pandas locally.
        if isinstance(p_orig, pd.DataFrame):                                                                                            # If the original input was a DataFrame.
            df = p_orig.copy()                                                                                                          # Generate a mutable copy of the DataFrame.
            if u_ap.ndim == 1:                                                                                                          # If the solution is 1D (Stationary problem).
                df['u_ap'] = u_ap                                                                                                       # type: ignore # Assign solution to a new column.
            else:                                                                                                                       # If the solution is 2D (Transient problem).
                for k in range(u_ap.shape[1]):                                                                                          # Iterate through each computed time step.
                    df[f'u_ap_t{k}'] = u_ap[:, k]                                                                                       # type: ignore # Assign time step to new column.
            return df                                                                                                                   # Return appended DataFrame.
            
    # 2. Pack into Xarray DataArray
    if "xarray" in sys.modules:                                                                                                         # If xarray is currently loaded in memory.
        import xarray as xr                                                                                                             # Import xarray locally.
        if isinstance(p_orig, xr.DataArray):                                                                                            # If the original input was a DataArray.
            xa = p_orig.copy()                                                                                                          # Generate a mutable copy of the DataArray.
            ds = xr.Dataset({"cloud": xa})                                                                                              # Initialize an xarray Dataset wrapper.
            if u_ap.ndim == 1:                                                                                                          # If the solution is 1D (Stationary problem).
                ds["u_ap"] = (xa.dims[0], u_ap)                                                                                         # Append solution mapping to spatial dimension.
            else:                                                                                                                       # If the solution is 2D (Transient problem).
                ds["u_ap"] = ((xa.dims[0], "time"), u_ap)                                                                               # Append solution mapping to space-time dimensions.
            return ds                                                                                                                   # Return appended Dataset.
            
    return u_ap                                                                                                                         # Return raw array if no specialized format matched.
