"""
Preconditioners — CUDA preconditioner utilities
"""

import numpy as np
from typing import Optional, Union, Any

try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix
    from cupyx.scipy.sparse.linalg import spilu as cp_spilu
    from cupyx.scipy.sparse.linalg import LinearOperator as cp_LinearOperator
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None
    cp_csr_matrix = None
    cp_spilu = None
    cp_LinearOperator = None

def compute_preconditioner(K: Any, method: str) -> Any:
    """
    compute_preconditioner
    Generates a Krylov preconditioner (LinearOperator) from a sparse matrix K to accelerate GMRES or BiCGStab on GPU.

    Input:
        K           m x m           cp_csr_matrix       The original sparse stiffness matrix.
        method                      str                 The preconditioning method: 'ilu' or 'jacobi'.

    Output:
        M           m x m           cp_LinearOperator   The preconditioner object (or None if method is None).
    """
    if not CUPY_AVAILABLE:
        raise ImportError("CuPy is not available.")

    if method is None or method.lower() == "none":                                                                                      # If no preconditioner is requested.
        return None                                                                                                                     # Return None.

    method = method.lower()                                                                                                             # Normalize string.
    m = K.shape[0]                                                                                                                      # Matrix dimension.

    if method == "ilu":                                                                                                                 # Incomplete LU Factorization.
        try:
            ilu_decomp = cp_spilu(K.tocsc())                                                                                            # type: ignore
            
            def M_x(x):                                                                                                                 # Define the solve operator.
                """
                M_x
                
                ILU preconditioner application on GPU. Solves the system M * y = x using the ILU decomposition.
                
                Input:
                    x           m               ndarray         Input vector to be preconditioned.
                    
                Output:
                    y           m               ndarray         Result of applying the ILU preconditioner to x.
                """
                return ilu_decomp.solve(x)                                                                                              # Return the solved value.
            
            return cp_LinearOperator((m, m), M_x)                                                                                       # type: ignore
        except RuntimeError:                                                                                                            # If ILU factorization fails.
            return None                                                                                                                 # Fallback to no preconditioner.

    elif method == "jacobi":                                                                                                            # Jacobi / Diagonal Scaling.
        diag = K.diagonal()                                                                                                             # Extract the main diagonal.
        diag[diag == 0] = 1e-12                                                                                                         # Prevent division by zero.
        inv_diag = 1.0 / diag                                                                                                           # Compute the inverse diagonal.
        
        def M_x(x):                                                                                                                     # Define the scaling operator.
            """
            M_x
            
            Jacobi preconditioner application on GPU. Scales the input vector x by the inverse of the matrix diagonal.
            
            Input:
                x           m               ndarray         Input vector to be preconditioned.
                
            Output:
                y           m               ndarray         Result of applying the Jacobi preconditioner to x.
            """
            return inv_diag * x                                                                                                         # Return the scaled value.
        
        return cp_LinearOperator((m, m), M_x)                                                                                           # type: ignore

    else:
        raise ValueError(f"Unknown preconditioning method: {method}. Choose 'ilu' or 'jacobi'.")                                        # Raise error on bad input.
