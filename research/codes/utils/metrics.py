"""
Metrics — Area-weighted error metrics on point clouds

Overview:
    Utilities to compute an area-weighted RMSE between an approximate solution (u_ap) and a reference
    solution (u_ex) on an unstructured 2D point cloud.

    The local area weight for each node i is estimated as the area of a polygon constructed from the
    coordinates of its immediate neighbors vec[i, :]. The polygon area is computed with the shoelace
    formula (PolyArea).

Notes:
    The polygon area calculation assumes that the neighbor vertices are ordered around the node
    (counterclockwise or clockwise). If vec does not provide an ordered ring, the computed area can
    be inaccurate (including self-intersections).

Public API:
    compute_rmse_transient      Area-weighted RMSE per time step.
    compute_rmse_stationary     Area-weighted RMSE for a single snapshot.
    Compute_Metrics_Stationary  Comprehensive error metrics for a single snapshot.
    Compute_Metrics_Transient   Comprehensive error metrics for transient solutions.

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
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Optional, Dict, Union                                                                                                # Type hints.

from mGFD.utils.core_utils import poly_area                                                                                             # Polygon area calculator.

def compute_rmse_transient(p: Union[np.ndarray, list], vec: Union[np.ndarray, list], u_ap: Union[np.ndarray, list], u_ex: Union[np.ndarray, list]) -> np.ndarray:
    """
    compute_rmse_transient
    Compute the area-weighted RMSE on an unstructured point cloud for a transient solution.
    
    For each node i, a local weight area[i] is computed as the polygon area formed by its neighbor
    coordinates p[vec[i, :], :]. For each time step k, the error is computed as:
        er[k] = sqrt(mean( (u_ap[:, k] - u_ex[:, k])^2 * area[:] ))
    
    Input:
        p               m x 2       ndarray         Node coordinates [x, y].
        vec             m x nvec    ndarray         Neighbor indices per node (unused slots padded with -1).
        u_ap            m x t       ndarray         Approximate solution values on nodes over time.
        u_ex            m x t       ndarray         Reference/exact solution values on nodes over time.
    
    Output:
        er              t           ndarray         Area-weighted RMSE at each time step.
    
    Notes:
        The neighbor list vec is assumed to define a valid polygonal ring around each node. If vec
        is unordered, poly_area may compute an inaccurate area and the weighting becomes unreliable.
    """
    # 0. Input validation
    if not isinstance(p, (np.ndarray, list)):                                                                                           # Validate p.
        raise TypeError("p must be an array-like object.")                                                                              # Raise explicit error on bad input.
    
    if not isinstance(vec, (np.ndarray, list)):                                                                                         # Validate vec.
        raise TypeError("vec must be an array-like object.")                                                                            # Raise explicit error on bad input.
    
    if not isinstance(u_ap, (np.ndarray, list)):                                                                                        # Validate u_ap.
        raise TypeError("u_ap must be an array-like object.")                                                                           # Raise explicit error on bad input.
    
    if not isinstance(u_ex, (np.ndarray, list)):                                                                                        # Validate u_ex.
        raise TypeError("u_ex must be an array-like object.")                                                                           # Raise explicit error on bad input.

    # 1. Variable initialization
    p_arr    = np.asarray(p)                                                                                                            # Normalize coordinates to ndarray.
    vec_arr  = np.asarray(vec)                                                                                                          # Normalize neighbors to ndarray.
    u_ap_arr = np.asarray(u_ap)                                                                                                         # Normalize approximate solution to ndarray.
    u_ex_arr = np.asarray(u_ex)                                                                                                         # Normalize reference solution to ndarray.
    m, t     = p_arr.shape[0], u_ap_arr.shape[1]                                                                                        # Number of nodes and time steps.
    er       = np.zeros(t)                                                                                                              # Per-time-step RMSE accumulator.
    area     = np.zeros(m)                                                                                                              # Per-node area weights.

    # 2. Area computation for each node
    for i in range(m):                                                                                                                  # Loop through all nodes.
        nvec_i  = int(np.sum(vec_arr[i, :] != -1))                                                                                      # Count valid neighbors (skip padding -1).
        nindex  = vec_arr[i, :nvec_i].astype(int)                                                                                       # Neighbor indices for node i.
        polix   = p_arr[nindex, 0]                                                                                                      # x-coordinates of polygon vertices.
        poliy   = p_arr[nindex, 1]                                                                                                      # y-coordinates of polygon vertices.
        area[i] = poly_area(polix, poliy)                                                                                               # Local area weight for node i.

    # 3. Error computation
    for k in range(t):                                                                                                                  # For each time step.
        diff2 = np.square(u_ap_arr[:, k] - u_ex_arr[:, k])                                                                              # Pointwise squared difference.
        err   = diff2 * area                                                                                                            # Area-weighted squared error.
        er[k] = np.sqrt(np.mean(err))                                                                                                   # Area-weighted RMSE for this time step.
    
    return er                                                                                                                           # Return per-time-step RMSE.

def compute_rmse_stationary(p: Union[np.ndarray, list], vec: Union[np.ndarray, list], u_ap: Union[np.ndarray, list], u_ex: Union[np.ndarray, list]) -> float:
    """
    compute_rmse_stationary
    Compute the area-weighted RMSE on an unstructured point cloud for a stationary (single snapshot) solution.
    
    The local area weights are computed in the same way as in compute_rmse_transient().
    
    Input:
        p               m x 2       ndarray         Node coordinates [x, y].
        vec             m x nvec    ndarray         Neighbor indices per node (unused slots padded with -1).
        u_ap            m           ndarray         Approximate solution values on nodes.
        u_ex            m           ndarray         Reference/exact solution values on nodes.
    
    Output:
        er              1           float           Area-weighted RMSE of the snapshot.
    
    Notes:
        The neighbor list vec is assumed to define a valid polygonal ring around each node. If vec
        is unordered, poly_area may compute an inaccurate area and the weighting becomes unreliable.
    """
    # 0. Input validation
    if not isinstance(p, (np.ndarray, list)):                                                                                           # Validate p.
        raise TypeError("p must be an array-like object.")                                                                              # Raise explicit error on bad input.
    
    if not isinstance(vec, (np.ndarray, list)):                                                                                         # Validate vec.
        raise TypeError("vec must be an array-like object.")                                                                            # Raise explicit error on bad input.
    
    if not isinstance(u_ap, (np.ndarray, list)):                                                                                        # Validate u_ap.
        raise TypeError("u_ap must be an array-like object.")                                                                           # Raise explicit error on bad input.
    
    if not isinstance(u_ex, (np.ndarray, list)):                                                                                        # Validate u_ex.
        raise TypeError("u_ex must be an array-like object.")                                                                           # Raise explicit error on bad input.

    # 1. Variable initialization
    p_arr    = np.asarray(p)                                                                                                            # Normalize coordinates to ndarray.
    vec_arr  = np.asarray(vec)                                                                                                          # Normalize neighbors to ndarray.
    u_ap_arr = np.asarray(u_ap)                                                                                                         # Normalize approximate solution to ndarray.
    u_ex_arr = np.asarray(u_ex)                                                                                                         # Normalize reference solution to ndarray.
    m        = p_arr.shape[0]                                                                                                           # Number of nodes.
    area     = np.zeros(m)                                                                                                              # Per-node area weights.

    # 2. Area computation for each node
    for i in range(m):                                                                                                                  # Loop through all nodes.
        nvec_i  = int(np.sum(vec_arr[i, :] != -1))                                                                                      # Count valid neighbors (skip padding -1).
        nindex  = vec_arr[i, :nvec_i].astype(int)                                                                                       # Neighbor indices for node i.
        polix   = p_arr[nindex, 0]                                                                                                      # x-coordinates of polygon vertices.
        poliy   = p_arr[nindex, 1]                                                                                                      # y-coordinates of polygon vertices.
        area[i] = poly_area(polix, poliy)                                                                                               # Local area weight for node i.

    # 3. Error computation
    diff2 = np.square(u_ap_arr[:] - u_ex_arr[:])                                                                                        # Pointwise squared difference.
    err   = diff2 * area                                                                                                                # Area-weighted squared error.
    er    = np.sqrt(np.mean(err))                                                                                                       # Area-weighted RMSE.
    
    return float(er)                                                                                                                    # Return scalar RMSE.

def Compute_Metrics_Stationary(p: Union[np.ndarray, list], vec: Union[np.ndarray, list], u_ap: Union[np.ndarray, list], u_ex: Union[np.ndarray, list], compute_time: Optional[float] = None) -> Dict[str, float]:
    """
    Compute_Metrics_Stationary
    Compute comprehensive error metrics for a stationary solution.
    
    This function computes the area-weighted RMSE as well as the maximum and mean
    absolute errors. It optionally includes the computation time.
    
    Input:
        p               m x 2       ndarray         Node coordinates [x, y].
        vec             m x nvec    ndarray         Neighbor indices per node (unused slots padded with -1).
        u_ap            m           ndarray         Approximate solution values on nodes.
        u_ex            m           ndarray         Reference/exact solution values on nodes.
        compute_time    1           float           (Optional) Time spent computing the solution in seconds.
        
    Output:
        metrics         1           dict            Dictionary containing 'RMSE', 'Max_Abs_Error', 'Mean_Abs_Error', and optionally 'Compute_Time_Secs'.
    """
    # 0. Input validation
    if not isinstance(p, (np.ndarray, list)):                                                                                           # Validate p.
        raise TypeError("p must be an array-like object.")                                                                              # Raise explicit error on bad input.
    
    if not isinstance(vec, (np.ndarray, list)):                                                                                         # Validate vec.
        raise TypeError("vec must be an array-like object.")                                                                            # Raise explicit error on bad input.
    
    if not isinstance(u_ap, (np.ndarray, list)):                                                                                        # Validate u_ap.
        raise TypeError("u_ap must be an array-like object.")                                                                           # Raise explicit error on bad input.
    
    if not isinstance(u_ex, (np.ndarray, list)):                                                                                        # Validate u_ex.
        raise TypeError("u_ex must be an array-like object.")                                                                           # Raise explicit error on bad input.

    # 1. Computation
    u_ap_arr = np.asarray(u_ap)                                                                                                         # Normalize approximate solution to ndarray.
    u_ex_arr = np.asarray(u_ex)                                                                                                         # Normalize reference solution to ndarray.

    rmse     = compute_rmse_stationary(p, vec, u_ap_arr, u_ex_arr)                                                                      # Compute the scalar RMSE for the stationary result.
    abs_diff = np.abs(u_ap_arr[:] - u_ex_arr[:])                                                                                        # Pointwise absolute difference.
    
    metrics: Dict[str, float] = {                                                                                                       # Dictionary to store the metrics.
        "RMSE":           rmse,                                                                                                         # Store Area-weighted RMSE.
        "Max_Abs_Error":  float(np.max(abs_diff)),                                                                                      # Store maximum absolute error.
        "Mean_Abs_Error": float(np.mean(abs_diff))                                                                                      # Store mean absolute error.
    }                                                                                                                                   # Execute statement.
    
    if compute_time is not None:                                                                                                        # If computation time is provided.
        metrics["Compute_Time_Secs"] = compute_time                                                                                     # Store computation time in seconds.
    
    return metrics                                                                                                                      # Return the populated metrics dictionary.

def Compute_Metrics_Transient(p: Union[np.ndarray, list], vec: Union[np.ndarray, list], u_ap: Union[np.ndarray, list], u_ex: Union[np.ndarray, list], compute_time: Optional[float] = None) -> Dict[str, float]:
    """
    Compute_Metrics_Transient
    Compute comprehensive error metrics for a transient solution.
    
    This function computes the area-weighted RMSE for the entire time domain, as well
    as the maximum absolute error across all nodes and time steps, and the RMSE of the final step.
    It optionally includes the computation time.
    
    Input:
        p               m x 2       ndarray         Node coordinates [x, y].
        vec             m x nvec    ndarray         Neighbor indices per node (unused slots padded with -1).
        u_ap            m x t       ndarray         Approximate solution values on nodes over time.
        u_ex            m x t       ndarray         Reference/exact solution values on nodes over time.
        compute_time    1           float           (Optional) Time spent computing the solution in seconds.
        
    Output:
        metrics         1           dict            Dictionary containing 'Time_Mean_RMSE', 'Max_Abs_Error', and optionally 'Compute_Time_Secs'.
    """
    # 0. Input validation
    if not isinstance(p, (np.ndarray, list)):                                                                                           # Validate p.
        raise TypeError("p must be an array-like object.")                                                                              # Raise explicit error on bad input.
    
    if not isinstance(vec, (np.ndarray, list)):                                                                                         # Validate vec.
        raise TypeError("vec must be an array-like object.")                                                                            # Raise explicit error on bad input.
    
    if not isinstance(u_ap, (np.ndarray, list)):                                                                                        # Validate u_ap.
        raise TypeError("u_ap must be an array-like object.")                                                                           # Raise explicit error on bad input.
    
    if not isinstance(u_ex, (np.ndarray, list)):                                                                                        # Validate u_ex.
        raise TypeError("u_ex must be an array-like object.")                                                                           # Raise explicit error on bad input.

    # 1. Computation
    u_ap_arr = np.asarray(u_ap)                                                                                                         # Normalize approximate solution to ndarray.
    u_ex_arr = np.asarray(u_ex)                                                                                                         # Normalize reference solution to ndarray.

    rmse_array = compute_rmse_transient(p, vec, u_ap_arr, u_ex_arr)                                                                     # Compute array of RMSE values across time steps.
    abs_diff   = np.abs(u_ap_arr - u_ex_arr)                                                                                            # Pointwise absolute difference over time and space.
    
    metrics: Dict[str, float] = {                                                                                                       # Dictionary to store the metrics.
        "Time_Mean_RMSE": float(np.mean(rmse_array)),                                                                                   # Store mean RMSE over all time steps.
        "Max_Abs_Error":  float(np.max(abs_diff))                                                                                       # Store absolute error peak.
    }                                                                                                                                   # Execute statement.
    if compute_time is not None:                                                                                                        # If computation time is provided.
        metrics["Compute_Time_Secs"] = compute_time                                                                                     # Store computation time in seconds.
    
    return metrics                                                                                                                      # Return the populated metrics dictionary.
