"""
Test Time Derivative 1 — Unit tests for the HeatEquation / 1st-order transient PDE OOP solver routines

Overview:
    This file contains the unit tests for validating the OOP 1st-order transient solver (HeatEquation), ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

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
import pytest                                                                                                           # Unit testing framework.
import numpy as np                                                                                                      # Core numerical operations.
from mGFD import Cloud, Dirichlet, HeatEquation, Solver                                                                 # OOP Architecture classes.
from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError, OperatorFormatError                        # Custom exceptions.

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
    x    = np.linspace(0, 1, n)                                                                                         # X-axis coordinates.
    y    = np.linspace(0, 1, n)                                                                                         # Y-axis coordinates.
    X, Y = np.meshgrid(x, y)                                                                                            # Mesh grid generation.
    X    = X.flatten()                                                                                                  # Flatten X array.
    Y    = Y.flatten()                                                                                                  # Flatten Y array.
    flag = np.zeros(len(X))                                                                                             # Node flag initialization.
    
    # 2. Boundary detection
    boun       = (X == 0) | (Y == 0) | (X == 1) | (Y == 1)                                                              # Logical mask for boundaries.
    flag[boun] = 1                                                                                                      # Flag boundary nodes as 1.
    
    return np.column_stack([X, Y, flag])                                                                                # Return assembled point cloud.

def test_time_derivative1_input_validation() -> None:
    """
    test_time_derivative1_input_validation
    Validates input argument checking for Cloud and HeatEquation.
    """
    # 1. Validation of Cloud Shape Error
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                              # Expect CloudShapeError.
        Cloud.from_array([1, 2, 3])                                                                                     # type: ignore                                                                                      # Execute with invalid type.

    # 2. Validation of Input Type Errors
    p = generate_square_cloud()                                                                                         # Generate point cloud.
    cloud = Cloud.from_array(p)                                                                                         # Instantiate Cloud.
    with pytest.raises(InputTypeError, match="Boundary condition 'value' must be a callable, ndarray, or numeric constant."):# Expect InputTypeError.
        cloud.set_boundary(Dirichlet("invalid_string"))                                                                 # type: ignore                                                                # Execute with invalid type.

    # 3. Validation of Operator Format Error
    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):       # Expect OperatorFormatError.
        HeatEquation(k=0.2, operator=np.array([1]))                                                                     # type: ignore                                                                     # Execute with invalid type.

def test_time_derivative1_heat() -> None:
    """
    test_time_derivative1_heat
    Validates OOP TimeDerivative1 solver with a Heat equation.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                       # Generate point cloud.
    v = 0.2                                                                                                             # Thermal diffusivity.
    t = 100                                                                                                             # Number of time steps.
    f = lambda x, y, t=0, coef=None: np.exp(-2 * np.pi**2 * v * t) * np.cos(np.pi * x) * np.cos(np.pi * y)              # Analytical solution function.
    ic = lambda x, y: np.cos(np.pi * x) * np.cos(np.pi * y)                                                             # Initial state.

    T = np.linspace(0, 1, t)                                                                                            # Time vector.
    u_ex = np.zeros([len(p), t])                                                                                        # Exact solution array.
    for k in range(t):                                                                                                  # Iterate over time.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [v])                                                                     # Compute exact solution at time step.
        
    # 2. Execution and Assertion
    cloud  = Cloud.from_array(p)                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(f))                                                                           # Set Dirichlet boundary condition.
    pde    = HeatEquation(k=v, ic=ic)                                                                                   # Formulate Heat PDE.
    solver = Solver(domain, pde, verbose=False)                                                                         # Instantiate Solver.
    result = solver.solve(dt=1.0/t)                                                                                     # Execute solver.

    u_ap = result.solution                                                                                              # Extract solution matrix.
    error_rms = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))                                                        # Compute RMS error.
    assert error_rms < 1e-1                                                                                             # Verify error bound.
