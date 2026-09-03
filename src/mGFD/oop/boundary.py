"""
mGFD.oop.boundary — Boundary Condition Abstractions for mGFD OOP Interface

Overview:
    This module provides object-oriented abstractions for boundary conditions in mGFD,
    including Dirichlet and Neumann boundary conditions, supporting callables, scalars,
    NumPy arrays, and Pandas Series/DataFrames.

Public API:
    BoundaryCondition   Base class for boundary conditions.
    Dirichlet           Dirichlet boundary condition class (u = val).
    Neumann             Neumann boundary condition class (du/dn = val).

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

Date:
    September, 2026.
"""

## Library importation.
import numpy as np                                                                                                      # Core numerical operations.
from scipy.spatial import cKDTree                                                                                       # Fast spatial coordinate lookup.
from typing import Union, Callable, Any, Optional                                                                       # Type hinting interfaces.
from mGFD.exceptions import InputTypeError                                                                              # Custom input exception.

class BoundaryCondition:                                                                                                # Base boundary condition class.
    """
    Base class for boundary condition representations in mGFD.
    """
    def __init__(self, value: Union[float, Callable, np.ndarray, Any] = 0.0) -> None:                                   # Initialize base boundary condition.
        """
        __init__
        Initialize boundary condition with value, callable, array, or Pandas series/dataframe.
        
        Input:
            value       Union[float, Callable, ndarray, Any]    Boundary value specification. Default is 0.0.
        """
        # Validate value type
        valid_scalar = isinstance(value, (int, float, np.number))                                                       # Check scalar.
        valid_callable = callable(value)                                                                                # Check callable.
        valid_array = hasattr(value, 'to_numpy') or hasattr(value, 'toarray') or isinstance(value, (np.ndarray, list, tuple))# Check array/pandas.
        
        if not (valid_scalar or valid_callable or valid_array):                                                         # If invalid type.
            raise InputTypeError("Boundary condition 'value' must be a callable, ndarray, or numeric constant.")        # Raise exception.

        self.value  = value                                                                                             # Store raw value or callable.
        self.points = None                                                                                              # Optional attached points cloud array (m, 3).
        self.t_span = (0.0, 1.0)                                                                                        # Optional physical time span domain.

    def __call__(self, x: np.ndarray, y: np.ndarray, t: Any = 0.0, coef: Optional[Any] = None) -> np.ndarray:           # Evaluates boundary condition at (x, y, t).
        """
        __call__
        Evaluates boundary condition at coordinates (x, y) at physical time t.
        
        Input:
            x           ndarray     X coordinates array.
            y           ndarray     Y coordinates array.
            t           Any         Physical time value or time array. Default is 0.0.
            coef        Any         Optional extra coefficient parameter.
            
        Output:
            bc_vals     ndarray     Vector of boundary condition values matching x shape.
        """
        if callable(self.value):                                                                                        # If value is a function callable.
            try:                                                                                                        # Attempt evaluation with (x, y, t, coef).
                res = self.value(x, y, t, coef) if coef is not None else self.value(x, y, t)                            # Evaluate callable function.
            except TypeError:                                                                                           # Fallback if callable takes (x, y).
                res = self.value(x, y)                                                                                  # Evaluate (x, y) signature.
            return np.asarray(res, dtype=np.float64)                                                                    # Return as NumPy array.
        elif isinstance(self.value, (int, float, np.number)):                                                           # If value is numeric scalar.
            return np.full_like(x, float(self.value), dtype=np.float64)                                                 # Return array filled with scalar.
        else:                                                                                                           # Array or Pandas DataFrame/Series handling.
            val_arr = self.value.to_numpy() if hasattr(self.value, 'to_numpy') else np.asarray(self.value, dtype=np.float64)# Convert to numpy array.
            x_arr   = np.asarray(x).ravel()                                                                             # Flatten x coordinates.
            y_arr   = np.asarray(y).ravel()                                                                             # Flatten y coordinates.

            # Determine spatial indexing
            if val_arr.shape[0] == len(x_arr):                                                                          # Exact spatial match.
                spatial_idx = slice(None)                                                                               # Use all rows.
            elif self.points is not None:                                                                               # KDTree spatial lookup.
                tree = cKDTree(self.points[:, :2])                                                                      # Build tree on full cloud points.
                _, spatial_idx = tree.query(np.column_stack([x_arr, y_arr]))                                            # Query nearest point indices.
            else:                                                                                                       # Fallback slice.
                spatial_idx = slice(0, len(x_arr))                                                                      # Truncate to length of x.

            # Determine temporal indexing
            if val_arr.ndim == 1:                                                                                       # 1D array.
                return val_arr[spatial_idx]                                                                             # Return 1D slice.
            elif val_arr.ndim == 2:                                                                                     # 2D spatiotemporal array.
                t0, t1 = self.t_span if self.t_span is not None else (0.0, 1.0)                                         # Extract t_span.
                t_len  = t1 - t0 if t1 > t0 else 1.0                                                                    # Physical time domain length.
                t_np   = np.asarray(t, dtype=np.float64)                                                                # Convert t to numpy array.
                frac   = np.clip((t_np - t0) / t_len, 0.0, 1.0)                                                         # Fractional time.
                k_idx  = np.minimum(np.round(frac * (val_arr.shape[1] - 1)).astype(int), val_arr.shape[1] - 1)          # Discrete time step indices.
                
                if k_idx.ndim == 0:                                                                                     # Scalar time index.
                    return val_arr[spatial_idx, int(k_idx)]                                                             # Return 1D vector.
                elif k_idx.ndim == 1:                                                                                   # 1D time index vector.
                    return val_arr[spatial_idx, k_idx]                                                                  # Return 2D matrix slice.
                else:                                                                                                   # 2D time index matrix (vectorized space-time).
                    k_flat = k_idx.ravel()                                                                              # Flatten time indices.
                    if isinstance(spatial_idx, slice):                                                                  # If all spatial rows.
                        return val_arr[:, k_flat]                                                                       # Return 2D matrix slice.
                    else:                                                                                               # If selected spatial rows.
                        return val_arr[spatial_idx][:, k_flat]                                                          # Return 2D matrix slice.
            else:                                                                                                       # Higher dimension fallback.
                return val_arr.flatten()[:len(x_arr)]                                                                   # Flatten and truncate.


class Dirichlet(BoundaryCondition):                                                                                     # Dirichlet boundary condition class.
    """
    Dirichlet Boundary Condition (u = val on boundary).
    """
    def __repr__(self) -> str:                                                                                          # String representation.
        return f"Dirichlet(value={self.value})"                                                                         # Return string format.


class Neumann(BoundaryCondition):                                                                                       # Neumann boundary condition class.
    """
    Neumann Boundary Condition (du/dn = val on boundary).
    """
    def __repr__(self) -> str:                                                                                          # String representation.
        return f"Neumann(value={self.value})"                                                                           # Return string format.
