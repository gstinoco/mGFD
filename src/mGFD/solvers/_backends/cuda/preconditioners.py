"""
Preconditioners — CUDA preconditioner utilities

Overview:
    Utilities to compute Krylov preconditioners on the GPU.

Public API:
    compute_preconditioner      Generates a CuPy LinearOperator preconditioner.

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
    August, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.

from typing import Optional, Union, Any                                                                                                 # Type hinting.

try:                                                                                                                                    # Attempt CuPy import.
    import cupy as cp                                                                                                                   # CuPy GPU array operations.
    from cupyx.scipy.sparse.linalg import spilu as cp_spilu                                                                             # CuPy sparse linear solvers.
    from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix                                                                          # CuPy sparse matrices.
    from cupyx.scipy.sparse.linalg import LinearOperator as cp_LinearOperator                                                           # CuPy sparse linear solvers.
    CUPY_AVAILABLE = True                                                                                                               # CuPy availability flag.
except ImportError:                                                                                                                     # Catch missing CuPy.
    CUPY_AVAILABLE = False                                                                                                              # CuPy availability flag.
    cp = None                                                                                                                           # Fallback dummy.
    cp_spilu = None                                                                                                                     # Fallback dummy.
    cp_csr_matrix = None                                                                                                                # Fallback dummy.
    cp_LinearOperator = None                                                                                                            # Fallback dummy.

def compute_preconditioner(K: Any, method: str) -> Any:                                                                                 # Function definition.
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
        try:                                                                                                                            # Attempt CuPy import.
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
