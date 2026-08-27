"""
Test Stationary — Unit tests for the Stationary solver routines

Overview:
    This file contains the unit tests for validating the Stationary solver, ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

Public API:
    test_stationary_input_validation    Validates input argument checking for Stationary.
    test_stationary_poisson             Validates Stationary solver with a Poisson equation.

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
import pytest                                                                                                                           # Unit testing framework.
import numpy as np                                                                                                                      # Core numerical operations.
from mGFD import Stationary                                                                                                             # Stationary solver.
from mGFD.exceptions import CloudShapeError, InputTypeError, OperatorFormatError                                                        # Custom exceptions.

def generate_square_cloud(n: int = 11) -> np.ndarray:
    """
    generate_square_cloud
    Helper function to generate a normalized square point cloud for testing.
    
    Input:
        n                           int             Number of nodes per side of the square.
        
    Output:
        p           m x 3           ndarray         Generated point cloud [x, y, flag].
    """
    # 1. Geometry generation
    x    = np.linspace(0, 1, n)                                                                                                         # X-axis coordinates.
    y    = np.linspace(0, 1, n)                                                                                                         # Y-axis coordinates.
    X, Y = np.meshgrid(x, y)                                                                                                            # Mesh grid generation.
    X    = X.flatten()                                                                                                                  # Flatten X array.
    Y    = Y.flatten()                                                                                                                  # Flatten Y array.
    flag = np.zeros(len(X))                                                                                                             # Node flag initialization.
    
    # 2. Boundary detection
    boun       = (X == 0) | (Y == 0) | (X == 1) | (Y == 1)                                                                              # Logical mask for boundaries.
    flag[boun] = 1                                                                                                                      # Flag boundary nodes as 1.
    
    return np.column_stack([X, Y, flag])                                                                                                # Return assembled point cloud.

def test_stationary_input_validation() -> None:
    """
    test_stationary_input_validation
    Validates input argument checking for Stationary.
    
    Systematically feeds invalid data structures, malformed point clouds, and incorrect
    operator lengths to the solver to ensure that our custom exception suite intercepts
    every erroneous state before it reaches numerical execution.
    """
    # 1. Validation of Cloud Shape Error
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                                              # Expect CloudShapeError.
        Stationary(p=[1, 2, 3], phi=lambda x,y: 0, f=lambda x,y: 0)  # type: ignore                                                     # Execute with invalid type.

    # 2. Validation of Input Type Errors
    p = generate_square_cloud()                                                                                                         # Generate point cloud.
    with pytest.raises(InputTypeError, match="Boundary condition 'phi' must be a callable, ndarray, or numeric constant."):             # Expect InputTypeError.
        Stationary(p=p, phi="invalid_string", f=lambda x,y: 0)  # type: ignore                                                          # Execute with invalid type.

    with pytest.raises(InputTypeError, match="Right-hand side 'f' must be a callable, ndarray, or numeric constant."):                  # Expect InputTypeError.
        Stationary(p=p, phi=lambda x,y: 0, f="invalid_string")  # type: ignore                                                          # Execute with invalid type.

    # 3. Validation of Operator Format Error
    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):                       # Expect OperatorFormatError.
        Stationary(p=p, phi=lambda x,y: 0, f=lambda x,y: 0, operator=np.array([1, 2]))  # type: ignore                                  # Execute with invalid type.

def test_stationary_poisson() -> None:
    """
    test_stationary_poisson
    Validates Stationary solver with a Poisson equation.
    
    Sets up a classic Poisson equation over a unit square. Validates that the approximated
    solution computed through spsolve, bicgstab, and gmres converges securely towards
    the exact analytical solution within the expected order of magnitude.
    """
    # 1. Test initialization
    p   = generate_square_cloud(21)                                                                                                     # Generate point cloud.
    phi = lambda x, y: 2 * np.exp(2 * x + y)                                                                                            # Boundary condition.
    f   = lambda x, y: 10 * np.exp(2 * x + y)                                                                                           # Forcing term.
    L   = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                                     # Laplacian operator.
    u_ex = phi(p[:, 0], p[:, 1])                                                                                                        # Exact solution for error check.
    
    # 2. Execution and Assertion - SPSOLVE
    u_ap, vec = Stationary(p, phi, f, operator=L, linear_solver="spsolve", verbose=False)                                               # Execute solver.
    error_spsolve = np.sqrt(np.mean((u_ap - u_ex)**2))                                                                                  # Compute RMS error.
    assert error_spsolve < 5e-2                                                                                                         # Verify error bound.

    # 3. Execution and Assertion - BICGSTAB
    u_ap_bicgstab, _ = Stationary(p, phi, f, operator=L, linear_solver="bicgstab", verbose=False)                                       # Execute solver.
    error_bicgstab = np.sqrt(np.mean((u_ap_bicgstab - u_ex)**2))                                                                        # Compute RMS error.
    assert error_bicgstab < 5e-2                                                                                                        # Verify error bound.

    # 4. Execution and Assertion - GMRES
    u_ap_gmres, _ = Stationary(p, phi, f, operator=L, linear_solver="gmres", verbose=False)                                             # Execute solver.
    error_gmres = np.sqrt(np.mean((u_ap_gmres - u_ex)**2))                                                                              # Compute RMS error.
    assert error_gmres < 5e-2                                                                                                           # Verify error bound.
