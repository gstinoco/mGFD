"""
Test Time Derivative 2 — Unit tests for the TimeDerivative2 solver routines

Overview:
    This file contains the unit tests for validating the TimeDerivative2 solver, ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

Public API:
    test_time_derivative2_input_validation    Validates input argument checking for TimeDerivative2.
    test_time_derivative2_wave                Validates TimeDerivative2 solver with a Wave equation.

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
from mGFD import TimeDerivative2                                                                                                        # TimeDerivative2 solver.
from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError                                                             # Custom exceptions.

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

def test_time_derivative2_input_validation() -> None:
    """
    test_time_derivative2_input_validation
    Validates input argument checking for TimeDerivative2.
    
    Systematically feeds invalid data structures, malformed point clouds, zero time steps,
    and incorrect functional definitions for position and velocity to ensure the solver
    intercepts them through the custom exception framework.
    """
    # 1. Validation of Cloud Shape Error
    p = generate_square_cloud()                                                                                                         # Generate point cloud.
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                                              # Expect CloudShapeError.
        TimeDerivative2(p=[1, 2], f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore                              # Execute with invalid type.

    # 2. Validation of Input Type Errors
    with pytest.raises(InputTypeError, match="Initial condition function 'f' must be a callable, ndarray, or numeric constant."):       # Expect InputTypeError.
        TimeDerivative2(p=p, f="invalid_string", g=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore                                    # Execute with invalid type.

    with pytest.raises(InputTypeError, match="Initial velocity function 'g' must be a callable, ndarray, or numeric constant."):        # Expect InputTypeError.
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g="invalid_string", t=10, coef=[])  # type: ignore                                    # Execute with invalid type.

    # 3. Validation of Parameter Error
    with pytest.raises(ParameterError, match="Number of time steps 't' must be a positive integer"):                                    # Expect ParameterError.
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=0, coef=[])                                                    # Execute with invalid type.

def test_time_derivative2_wave() -> None:
    """
    test_time_derivative2_wave
    Validates TimeDerivative2 solver with a Wave equation.
    
    Evaluates a vibrating 2D membrane using a standard hyperbolic wave model. Ensures
    that the solver effectively incorporates the initial position and velocity arrays
    and propagates the signal forward without artificial numerical damping.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                                       # Generate point cloud.
    c = float(np.sqrt(1 / 2))                                                                                                           # Wave speed.
    t = 100                                                                                                                             # Number of time steps.
    f = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                                                               # Initial condition function.
    g = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))                                                      # Initial velocity function.
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                                                         # Wave operator.
    
    T = np.linspace(0, 1, t)                                                                                                            # Time vector.
    u_ex = np.zeros([len(p), t])                                                                                                        # Exact solution array.
    for k in range(t):                                                                                                                  # Iterate over time.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [c])                                                                                     # Compute exact solution at time step.

    # 2. Execution and Assertion
    u_ap, vec = TimeDerivative2(p, f, g, t, [c], operator=L, implicit=True, verbose=False)                                             # Execute solver.
    error_spsolve = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))                                                                    # Compute RMS error.
    assert error_spsolve < 2e-1                                                                                                         # Verify error bound.
