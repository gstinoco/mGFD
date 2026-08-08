"""
Errors — Area-weighted error metrics on point clouds

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
    PolyArea             Polygon area by the shoelace formula.
    compute_rmse_transient      Area-weighted RMSE per time step.
    compute_rmse_stationary     Area-weighted RMSE for a single snapshot.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
        Aula CIMNE-Morelia. México
        SIIIA-MATH: Soluciones de Ingeniería. México

    Date:
        May, 2024.
    Last Modification:
        August, 2026.
"""
## Library importation.
import numpy as np                                                                                      # Core numerical operations.

from mGFD.core.utils import poly_area                                                                     # Polygon area calculator.

def compute_rmse_transient(p, vec, u_ap, u_ex):
    '''
    compute_rmse_transient
    Compute the area-weighted RMSE on an unstructured point cloud for a transient solution.
    
    For each node i, a local weight area[i] is computed as the polygon area formed by its neighbor
    coordinates p[vec[i, :], :]. For each time step k, the error is computed as:
        er[k] = sqrt(mean( (u_ap[:, k] - u_ex[:, k])^2 * area[:] ))
    
    Input:
        p           m x 2           array-like      Node coordinates [x, y].
        vec         m x nvec        array-like[int] Neighbor indices per node (unused slots padded with -1).
        u_ap        m x t           array-like      Approximate solution values on nodes over time.
        u_ex        m x t           array-like      Reference/exact solution values on nodes over time.
    
    Output:
        er          t               ndarray         Area-weighted RMSE at each time step.
    
    Notes:
        The neighbor list vec is assumed to define a valid polygonal ring around each node. If vec
        is unordered, PolyArea may compute an inaccurate area and the weighting becomes unreliable.
    '''
    ## Variable initialization.
    p     = np.asarray(p)                                                                               # Normalize coordinates to ndarray.
    vec   = np.asarray(vec)                                                                             # Normalize neighbors to ndarray.
    u_ap  = np.asarray(u_ap)                                                                            # Normalize approximate solution to ndarray.
    u_ex  = np.asarray(u_ex)                                                                            # Normalize reference solution to ndarray.
    m, t  = p.shape[0], u_ap.shape[1]                                                                   # Number of nodes and time steps.
    er    = np.zeros(t)                                                                                 # Per-time-step RMSE accumulator.
    area  = np.zeros(m)                                                                                 # Per-node area weights.

    ## Area computation for each node.
    for i in np.arange(m):
        nvec_i   = int(np.sum(vec[i, :] != -1))                                                         # Count valid neighbors (skip padding -1).
        nindex   = vec[i, :nvec_i].astype(int)                                                          # Neighbor indices for node i.
        polix    = p[nindex, 0]                                                                         # x-coordinates of polygon vertices.
        poliy    = p[nindex, 1]                                                                         # y-coordinates of polygon vertices.
        area[i]  = poly_area(polix, poliy)                                                              # Local area weight for node i.

    ## Error computation.
    for k in np.arange(t):                                                                              # For each time step.
        diff2 = np.square(u_ap[:, k] - u_ex[:, k])                                                      # Pointwise squared difference.
        err   = diff2 * area                                                                            # Area-weighted squared error.
        er[k] = np.sqrt(np.mean(err))                                                                   # Area-weighted RMSE for this time step.
    
    return er                                                                                           # Return per-time-step RMSE.

def compute_rmse_stationary(p, vec, u_ap, u_ex):
    '''
    compute_rmse_stationary
    Compute the area-weighted RMSE on an unstructured point cloud for a stationary (single snapshot) solution.
    
    The local area weights are computed in the same way as Cloud_Transient().
    
    Input:
        p           m x 2           array-like      Node coordinates [x, y].
        vec         m x nvec        array-like[int] Neighbor indices per node (unused slots padded with -1).
        u_ap        m               array-like      Approximate solution values on nodes.
        u_ex        m               array-like      Reference/exact solution values on nodes.
    
    Output:
        er                          float           Area-weighted RMSE of the snapshot.
    
    Notes:
        The neighbor list vec is assumed to define a valid polygonal ring around each node. If vec
        is unordered, PolyArea may compute an inaccurate area and the weighting becomes unreliable.
    '''
    ## Variable initialization.
    p     = np.asarray(p)                                                                               # Normalize coordinates to ndarray.
    vec   = np.asarray(vec)                                                                             # Normalize neighbors to ndarray.
    u_ap  = np.asarray(u_ap)                                                                            # Normalize approximate solution to ndarray.
    u_ex  = np.asarray(u_ex)                                                                            # Normalize reference solution to ndarray.
    m     = p.shape[0]                                                                                  # Number of nodes.
    area  = np.zeros(m)                                                                                 # Per-node area weights.

    ## Area computation for each node.
    for i in np.arange(m):
        nvec_i   = int(np.sum(vec[i, :] != -1))                                                         # Count valid neighbors (skip padding -1).
        nindex   = vec[i, :nvec_i].astype(int)                                                          # Neighbor indices for node i.
        polix    = p[nindex, 0]                                                                         # x-coordinates of polygon vertices.
        poliy    = p[nindex, 1]                                                                         # y-coordinates of polygon vertices.
        area[i]  = poly_area(polix, poliy)                                                              # Local area weight for node i.

    ## Error computation.
    diff2 = np.square(u_ap[:] - u_ex[:])                                                                # Pointwise squared difference.
    err   = diff2 * area                                                                                # Area-weighted squared error.
    er    = np.sqrt(np.mean(err))                                                                       # Area-weighted RMSE.
    
    return float(er)                                                                                    # Return scalar RMSE.

def Compute_Metrics_Stationary(p, vec, u_ap, u_ex, compute_time=None):
    '''
    Compute_Metrics_Stationary
    Compute comprehensive error metrics for a stationary solution.
    
    This function computes the area-weighted RMSE as well as the maximum and mean
    absolute errors. It optionally includes the computation time.
    
    Input:
        p                           ndarray         Node coordinates [x, y].
        vec                         ndarray[int]    Neighbor indices per node (unused slots padded with -1).
        u_ap                        ndarray         Approximate solution values on nodes.
        u_ex                        ndarray         Reference/exact solution values on nodes.
        compute_time                float           (Optional) Time spent computing the solution in seconds.
        
    Output:
        metrics                     dict            Dictionary containing 'RMSE', 'Max_Abs_Error', 'Mean_Abs_Error', and optionally 'Compute_Time_Secs'.
    '''
    rmse = compute_rmse_stationary(p, vec, u_ap, u_ex)                                                  # Compute the scalar RMSE for the stationary result.
    abs_diff = np.abs(u_ap[:] - u_ex[:])                                                                # Pointwise absolute difference.
    metrics = {                                                                                         # Dictionary to store the metrics.
        "RMSE": float(rmse),                                                                            # Store Area-weighted RMSE.
        "Max_Abs_Error": float(np.max(abs_diff)),                                                       # Store maximum absolute error.
        "Mean_Abs_Error": float(np.mean(abs_diff))                                                      # Store mean absolute error.
    }
    if compute_time is not None:                                                                        # If computation time is provided.
        metrics["Compute_Time_Secs"] = float(compute_time)                                              # Store computation time in seconds.
    
    return metrics                                                                                      # Return the populated metrics dictionary.

def Compute_Metrics_Transient(p, vec, u_ap, u_ex, compute_time=None):
    '''
    Compute_Metrics_Transient
    Compute comprehensive error metrics for a transient solution.
    
    This function computes the area-weighted RMSE for the entire time domain, as well
    as the maximum absolute error across all nodes and time steps, and the RMSE of the final step.
    It optionally includes the computation time.
    
    Input:
        p                           ndarray         Node coordinates [x, y].
        vec                         ndarray[int]    Neighbor indices per node (unused slots padded with -1).
        u_ap                        ndarray         Approximate solution values on nodes over time.
        u_ex                        ndarray         Reference/exact solution values on nodes over time.
        compute_time                float           (Optional) Time spent computing the solution in seconds.
        
    Output:
        metrics                     dict            Dictionary containing 'Time_Mean_RMSE', 'Max_Abs_Error', 'Final_Step_RMSE', and optionally 'Compute_Time_Secs'.
    '''
    rmse_array = compute_rmse_transient(p, vec, u_ap, u_ex)                                             # Compute array of RMSE values across time steps.
    abs_diff = np.abs(u_ap - u_ex)                                                                      # Pointwise absolute difference over time and space.
    metrics = {                                                                                         # Dictionary to store the metrics.
        "Time_Mean_RMSE": float(np.mean(rmse_array)),                                                   # Store mean RMSE over all time steps.
        "Max_Abs_Error": float(np.max(abs_diff)),                                                       # Store absolute error peak.
        "Final_Step_RMSE": float(rmse_array[-1]) if len(rmse_array) > 0 else 0.0                        # Store RMSE at the very last step.
    }
    if compute_time is not None:                                                                        # If computation time is provided.
        metrics["Compute_Time_Secs"] = float(compute_time)                                              # Store computation time in seconds.
    
    return metrics                                                                                      # Return the populated metrics dictionary.