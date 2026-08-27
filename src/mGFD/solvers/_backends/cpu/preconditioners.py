"""
Preconditioners — CPU preconditioner utilities
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spilu, LinearOperator
from typing import Optional, Union, Any

def compute_preconditioner(K: Union[csr_matrix, Any], method: str) -> Optional[Union[LinearOperator, Any]]:
    """
    compute_preconditioner
    Generates a Krylov preconditioner (LinearOperator) from a sparse matrix K to accelerate GMRES or BiCGStab.

    Input:
        K           m x m           csr_matrix          The original sparse stiffness matrix.
        method                      str                 The preconditioning method: 'ilu' or 'jacobi'.

    Output:
        M           m x m           LinearOperator      The preconditioner object (or None if method is None).
    """
    if method is None or method.lower() == "none":                                                                                      # If no preconditioner is requested.
        return None                                                                                                                     # Return None.

    method = method.lower()                                                                                                             # Normalize string.
    m = K.shape[0]                                                                                                                      # Matrix dimension.

    if method == "ilu":                                                                                                                 # Incomplete LU Factorization.
        try:
            ilu_decomp = spilu(K.tocsc())                                                                                               # Compute ILU decomposition (requires CSC).
            
            def M_x(x):                                                                                                                 # Define the solve operator.
                """
                M_x
                
                ILU preconditioner application. Solves the system M * y = x using the ILU decomposition.
                
                Input:
                    x           m               ndarray         Input vector to be preconditioned.
                    
                Output:
                    y           m               ndarray         Result of applying the ILU preconditioner to x.
                """
                return ilu_decomp.solve(x)                                                                                              # Return the solved value.
            
            return LinearOperator((m, m), M_x)                                                                                          # type: ignore
        except RuntimeError:                                                                                                            # If ILU factorization fails (e.g., exactly singular).
            return None                                                                                                                 # Fallback to no preconditioner.

    elif method == "jacobi":                                                                                                            # Jacobi / Diagonal Scaling.
        diag = K.diagonal()                                                                                                             # Extract the main diagonal.
        diag[diag == 0] = 1e-12                                                                                                         # Prevent division by zero.
        inv_diag = 1.0 / diag                                                                                                           # Compute the inverse diagonal.
        
        def M_x(x):                                                                                                                     # Define the scaling operator.
            """
            M_x
            
            Jacobi preconditioner application. Scales the input vector x by the inverse of the matrix diagonal.
            
            Input:
                x           m               ndarray         Input vector to be preconditioned.
                
            Output:
                y           m               ndarray         Result of applying the Jacobi preconditioner to x.
            """
            return inv_diag * x                                                                                                         # Return the scaled value.
        
        return LinearOperator((m, m), M_x)                                                                                              # type: ignore

    else:
        raise ValueError(f"Unknown preconditioning method: {method}. Choose 'ilu' or 'jacobi'.")                                        # Raise error on bad input.
