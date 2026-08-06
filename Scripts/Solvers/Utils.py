"""
Utils — Utility functions for mGFD solvers

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
import numpy as np
from Scripts import Neighbors
import sys

__all__ = ["_import_scipy_stationary", "_import_scipy_implicit", "_node_masks", "_prepare_neighbors", "_next_nvec"]

def _import_scipy_stationary():
    '''
    _import_scipy_stationary
    Lazy import of the SciPy sparse utilities used by Stationary().
    
    Notes:
        This keeps SciPy as an optional dependency for the module. If SciPy is not installed,
        this helper raises ImportError and the caller may fall back to non-SciPy paths.
    
    Output:
        LinearOperator            SciPy LinearOperator class.
        bicgstab                  SciPy BiCGStab iterative solver.
        spsolve                   SciPy direct sparse solver.
        coo_matrix                SciPy COO sparse matrix constructor.
    '''
    try:                                                                                                # Import SciPy lazily to keep SciPy optional for the module.
        from scipy.sparse.linalg import LinearOperator, bicgstab, spsolve                               # Stationary: matrix-free operator + iterative/direct solvers.
        from scipy.sparse import coo_matrix                                                             # Stationary fallback: sparse matrix assembly.
    except ImportError as e:                                                                            # SciPy is not installed in the environment.
        raise ImportError('SciPy is required for this operation') from e                                # Provide a clear dependency error.
    return LinearOperator, bicgstab, spsolve, coo_matrix                                                # Return imported symbols to the caller.


def _import_scipy_implicit():
    '''
    _import_scipy_implicit
    Lazy import of the SciPy sparse utilities used by implicit schemes in TimeDerivative1/2.
    
    Notes:
        This keeps explicit schemes usable without SciPy. If SciPy is not installed, this helper
        raises ImportError and the implicit branches should surface a clear dependency message.
    
    Output:
        coo_matrix                SciPy COO sparse matrix constructor.
        eye                       SciPy sparse identity constructor.
        bicgstab                  SciPy BiCGStab iterative solver.
        spsolve                   SciPy direct sparse solver.
        splu                      SciPy sparse LU factorization.
    '''
    try:                                                                                                # Import SciPy lazily to keep explicit schemes usable without SciPy.
        from scipy.sparse import coo_matrix, eye                                                        # Sparse matrix utilities (assembly + identity).
        from scipy.sparse.linalg import bicgstab, spsolve, splu                                         # Sparse solvers (iterative/direct/LU).
    except ImportError as e:                                                                            # SciPy is not installed in the environment.
        raise ImportError('SciPy is required for implicit schemes') from e                              # Provide a clear dependency error.
    return coo_matrix, eye, bicgstab, spsolve, splu                                                     # Return imported symbols to the caller.


def _node_masks(p):
    '''
    _node_masks
    Build boolean masks for boundary and interior nodes.
    
    Input:
        p           m x 3           ndarray         Array with node coordinates and a node flag in column 2.
                                                    flag = 0 for interior, flag = 1/2 for boundary.
    
    Output:
        boun_n      m               ndarray         Boolean mask: True on boundary nodes (flag 1 or 2).
        inne_n      m               ndarray         Boolean mask: True on interior nodes (flag 0).
    '''
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                                            # Boundary nodes mask.
    inne_n = p[:, 2] == 0                                                                               # Interior nodes mask.
    return boun_n, inne_n                                                                               # Return both masks.


def _prepare_neighbors(p, vec, nvec, Adv, operator):
    '''
    _prepare_neighbors
    Validate/load/build the neighbor array used by the stencil.
    
    Input:
        p           m x 3           ndarray         Array with node coordinates and type flag.
        vec                         ndarray|None    Optional precomputed neighbors (m x nvec, padded with -1).
        nvec                        int             Requested maximum number of neighbors.
        Adv                         bool            If True, use upwind neighbor selection (CloudAdv).
        operator                    ndarray         Operator coefficients (used to extract velocities when Adv=True).
                                                    It assumes operator[0] and operator[1] correspond to D and E.
    
    Output:
        vec         m x nvec        ndarray         Neighbor indices per node (int32).
        nvec                        int             Effective nvec (taken from vec if vec provided).
    '''
    m    = int(p.shape[0])                                                                              # Number of nodes.
    nvec = int(nvec)                                                                                    # Ensure integer nvec.

    if vec is None:                                                                                     # Compute neighbors when not provided.
        if Adv:                                                                                         # Upwind neighbor selection (advection problems).
            a   = float(operator[0][0])                                                                 # Velocity in x (for upwind).
            b   = float(operator[1][0])                                                                 # Velocity in y (for upwind).
            vec = Neighbors.CloudAdv(p, a, b, nvec)                                                     # Compute neighbors with upwind preference.
        else:                                                                                           # Standard neighbor selection (no upwind bias).
            vec = Neighbors.Cloud(p, nvec)                                                              # Compute standard neighbors.
    else:                                                                                               # Validate and normalize a provided vec.
        vec = np.asarray(vec)                                                                           # Normalize input to ndarray.
        if vec.ndim == 1:                                                                               # Allow passing a single row as 1D.
            vec = vec.reshape(1, -1)                                                                    # Support single-row vec passed as 1D.
        if vec.shape[0] != m:                                                                           # Ensure one neighbor row per node.
            raise ValueError('vec has incorrect number of rows')                                        # Shape mismatch.
        nvec = int(vec.shape[1])                                                                        # Trust nvec from provided vec.

    vec = vec.astype(np.int32, copy = False)                                                            # Ensure fast indexing dtype.
    return vec, nvec                                                                                    # Return neighbors and effective nvec.


def _next_nvec(nvec):
    '''
    _next_nvec
    Choose the next neighbor count when expanding stencils.
    
    Input:
        nvec                        int             Current nvec.
    
    Output:
        nvec_next                   int|None        Next nvec in {12,16,20,30} or None if no expansion left.
    '''
    nvec = int(nvec)                                                                                    # Ensure integer nvec.
    if nvec < 12:                                                                                       # Promote to the first expansion size.
        return 12                                                                                       # First expansion.
    if nvec < 16:                                                                                       # Promote to the second expansion size.
        return 16                                                                                       # Second expansion.
    if nvec < 20:                                                                                       # Promote to the third expansion size.
        return 20                                                                                       # Third expansion.
    if nvec < 30:                                                                                       # Promote to the largest supported expansion size.
        return 30                                                                                       # Max expansion used in batches.
    return None                                                                                         # No further expansion possible.

