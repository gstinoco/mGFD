"""
Preconditioners — CPU preconditioner utilities

Overview:
    Utilities to compute Krylov preconditioners on the CPU.

Public API:
    compute_preconditioner      Generates a SciPy LinearOperator preconditioner.

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

from scipy.sparse import csr_matrix                                                                                                     # SciPy sparse matrices.
from typing import Optional, Union, Any                                                                                                 # Type hinting.
from scipy.sparse.linalg import spilu, LinearOperator                                                                                   # SciPy sparse linear solvers.

def compute_preconditioner(K: Union[csr_matrix, Any], method: str) -> Optional[Union[LinearOperator, Any]]:                             # Function definition.
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
        try:                                                                                                                            # Attempt CuPy import.
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
