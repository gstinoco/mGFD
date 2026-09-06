"""
Test Stationary — Unit tests for the Stationary PDE OOP solver routines

Overview:
    This file contains the unit tests for validating the OOP Stationary solver (PDE order 0), ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

Public API:
    generate_square_cloud                           Generates a square point cloud for testing.
    test_stationary_input_validation                Validates input validation and exception raising in Stationary solver.
    test_stationary_poisson                         Validates numerical convergence and accuracy for Poisson problem.

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
    September, 2026.
Last Modification:
    September, 2026.
"""

## Library importation.
import pytest                                                                                                                           # Unit testing framework.
import numpy as np                                                                                                                      # Core numerical operations.

from mGFD import Cloud, Dirichlet, Domain, PDE, Solver                                                                                  # OOP Architecture classes.
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
    Validates input argument checking for Cloud and Dirichlet.
    """
    # 1. Validation of Cloud Shape Error
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                                              # Expect CloudShapeError.
        Cloud.from_array([1, 2, 3])                                                                                                     # type: ignore.

    # 2. Validation of Input Type Errors
    p = generate_square_cloud()                                                                                                         # Generate point cloud.
    cloud = Cloud.from_array(p)                                                                                                         # Create Cloud.
    with pytest.raises(InputTypeError, match="Boundary condition 'value' must be a callable, ndarray, or numeric constant."):           # Expect InputTypeError.
        cloud.set_boundary(Dirichlet("invalid_string"))                                                                                 # type: ignore

    # 3. Validation of Operator Format Error
    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):                       # Expect OperatorFormatError.
        PDE(source=lambda x, y: 0, operator=np.array([1, 2]), order=0)                                                                  # type: ignore

def test_stationary_poisson() -> None:
    """
    test_stationary_poisson
    Validates OOP Stationary solver with a Poisson equation.
    """
    # 1. Test initialization
    p    = generate_square_cloud(21)                                                                                                    # Generate point cloud.
    phi  = lambda x, y: 2 * np.exp(2 * x + y)                                                                                           # Boundary condition.
    f    = lambda x, y: 10 * np.exp(2 * x + y)                                                                                          # Forcing term.
    u_ex = phi(p[:, 0], p[:, 1])                                                                                                        # Exact solution for error check.
    
    # 2. OOP Execution and Assertion
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(phi))                                                                                         # Set Dirichlet boundary.
    pde    = PDE(operator=[0, 0, 2, 0, 2, 0], source=f, order=0)                                                                        # Formulate Stationary PDE.
    solver = Solver(domain, pde, verbose=False)                                                                                         # Instantiate Solver.
    result = solver.solve()                                                                                                             # Execute solver.

    u_ap      = result.solution                                                                                                         # Extract solution array.
    error_rms = np.sqrt(np.mean((u_ap - u_ex)**2))                                                                                      # Compute RMS error.
    assert error_rms < 5e-2                                                                                                             # Verify error bound.
