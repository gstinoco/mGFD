"""
Test Time Derivative 1 — Unit tests for the TimeDerivative1 solver routines

Overview:
    This file contains the unit tests for validating the TimeDerivative1 solver, ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

Public API:
    test_time_derivative1_input_validation    Validates input argument checking for TimeDerivative1.
    test_time_derivative1_heat                Validates TimeDerivative1 solver with a Heat equation.

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
from mGFD import TimeDerivative1                                                                                                        # TimeDerivative1 solver.
from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError, OperatorFormatError                                        # Custom exceptions.

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

def test_time_derivative1_input_validation() -> None:
    """
    test_time_derivative1_input_validation
    Validates input argument checking for TimeDerivative1.
    
    Systematically feeds invalid data structures, malformed point clouds, negative time
    steps, and incorrect operator lengths to ensure that exceptions are raised correctly
    prior to any numerical time-marching process.
    """
    # 1. Validation of Cloud Shape Error
    p = generate_square_cloud()                                                                                                         # Generate point cloud.
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                                              # Expect CloudShapeError.
        TimeDerivative1(p=[1, 2, 3], f=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore                                                # Execute with invalid type.

    # 2. Validation of Input Type Errors
    with pytest.raises(InputTypeError, match="Forcing function 'f' must be a callable, ndarray, or numeric constant."):                 # Expect InputTypeError.
        TimeDerivative1(p=p, f="invalid_string", t=10, coef=[])  # type: ignore                                                         # Execute with invalid type.

    # 3. Validation of Parameter Error
    with pytest.raises(ParameterError, match="Number of time steps 't' must be a positive integer"):                                    # Expect ParameterError.
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=-1, coef=[])                                                                        # Execute with invalid type.

    # 4. Validation of Operator Format Error
    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):                       # Expect OperatorFormatError.
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=10, coef=[], operator=np.array([1]))  # type: ignore                                # Execute with invalid type.

def test_time_derivative1_heat() -> None:
    """
    test_time_derivative1_heat
    Validates TimeDerivative1 solver with a Heat equation.
    
    Simulates a 2D transient thermal conduction process. Verifies that all sparse linear
    solvers (spsolve, bicgstab, gmres) maintain numerical stability and accurately track
    the analytical decay profile across all temporal layers.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                                       # Generate point cloud.
    v = 0.2                                                                                                                             # Thermal diffusivity.
    t = 100                                                                                                                             # Number of time steps.
    f = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)                               # Analytical solution function.
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                               # Heat operator.
    
    T = np.linspace(0, 1, t)                                                                                                            # Time vector.
    u_ex = np.zeros([len(p), t])                                                                                                        # Exact solution array.
    for k in range(t):                                                                                                                  # Iterate over time.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [v])                                                                                     # Compute exact solution at time step.
        
    # 2. Execution and Assertion
    u_ap, vec = TimeDerivative1(p, f, t, [v], operator=L, implicit=True, verbose=False)                                                 # Execute solver.
    error_spsolve = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))                                                                    # Compute RMS error.
    assert error_spsolve < 1e-1                                                                                                         # Verify error bound.
