"""
Test Time Derivative 2 — Unit tests for the WaveEquation / 2nd-order transient PDE OOP solver routines

Overview:
    This file contains the unit tests for validating the OOP 2nd-order transient solver (WaveEquation), ensuring that the 
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
from mGFD import Cloud, Dirichlet, WaveEquation, Solver                                                                 # OOP Architecture classes.
from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError                                             # Custom exceptions.

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

def test_time_derivative2_input_validation() -> None:
    """
    test_time_derivative2_input_validation
    Validates input argument checking for Cloud and WaveEquation.
    """
    # 1. Validation of Cloud Shape Error
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):                              # Expect CloudShapeError.
        Cloud.from_array([1, 2])                                                                                        # type: ignore                                                                                        # Execute with invalid type.

    # 2. Validation of Input Type Errors
    p = generate_square_cloud()                                                                                         # Generate point cloud.
    cloud = Cloud.from_array(p)                                                                                         # Instantiate Cloud.
    with pytest.raises(InputTypeError, match="Boundary condition 'value' must be a callable, ndarray, or numeric constant."):# Expect InputTypeError.
        cloud.set_boundary(Dirichlet("invalid_string"))                                                                 # type: ignore                                                                # Execute with invalid type.

def test_time_derivative2_wave() -> None:
    """
    test_time_derivative2_wave
    Validates OOP TimeDerivative2 solver with a Wave equation.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                       # Generate point cloud.
    c = float(np.sqrt(1 / 2))                                                                                           # Wave speed.
    t = 100                                                                                                             # Number of time steps.
    f = lambda x, y, t=0, coef=None: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                                        # Boundary & analytical solution function.
    ic = lambda x, y: np.sin(np.pi * (x + y))                                                                           # Initial position u(x, y, 0).
    g = 0.0                                                                                                             # Initial velocity u_t(x, y, 0).
    
    T = np.linspace(0, 1, t)                                                                                            # Time vector.
    u_ex = np.zeros([len(p), t])                                                                                        # Exact solution array.
    for k in range(t):                                                                                                  # Iterate over time.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [c])                                                                     # Compute exact solution at time step.

    # 2. Execution and Assertion
    cloud  = Cloud.from_array(p)                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(f))                                                                           # Set Dirichlet boundary.
    pde    = WaveEquation(c=c, ic=ic, g=g)                                                                              # Formulate Wave PDE.
    solver = Solver(domain, pde, verbose=False)                                                                         # Instantiate Solver.
    result = solver.solve(dt=1.0/t)                                                                                     # Execute solver.

    u_ap = result.solution                                                                                              # Extract solution array.
    error_rms = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))                                                        # Compute RMS error.
    assert error_rms < 2e-1                                                                                             # Verify error bound.

def test_time_derivative2_hht_alpha_and_damping() -> None:
    """
    test_time_derivative2_hht_alpha_and_damping
    Validates that HHT-alpha dissipation and velocity damping execute correctly and keep solution bounded.
    """
    p = generate_square_cloud(11)                                                                                       # Generate 11x11 point cloud.
    c = 0.5                                                                                                             # Wave speed.
    ic = lambda x, y: np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))                                                       # Initial position Gaussian pulse.
    
    cloud  = Cloud.from_array(p)                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                         # Set zero Dirichlet boundary.
    pde    = WaveEquation(c=c, damping=0.01, alpha=-0.05, ic=ic, g=0.0)                                                 # Formulate Wave PDE with HHT-alpha & damping.
    solver = Solver(domain, pde, verbose=False)                                                                         # Instantiate Solver.
    res    = solver.solve(t_span=(0.0, 1.0))                                                                            # Execute solver.

    sol = res.solution                                                                                                  # Extract solution array.
    assert np.all(np.isfinite(sol))                                                                                     # Verify finite values (no NaN/Inf).
    assert np.max(np.abs(sol)) < 2.0                                                                                    # Verify solution remains bounded.

def test_time_derivative2_irregular_cloud_stability() -> None:
    """
    test_time_derivative2_irregular_cloud_stability
    Validates that WaveEquation remains bounded and stable on an irregular natural point cloud over long physical time steps.
    """
    import csv, os                                                                                                      # Standard OS and CSV interfaces.
    import pandas as pd                                                                                                 # DataFrame manipulation.
    from mGFD.cloud_generator.core.generator import generate_cloud_natural                                              # Natural point cloud generator.

    theta = np.linspace(0, 2 * np.pi, 60)                                                                               # Angular sampling array.
    r_star = 0.5 + 0.15 * np.sin(5 * theta)                                                                             # Star radius contour profile.
    star_contour = [(0.5 + r_star[i] * np.cos(theta[i]), 0.5 + r_star[i] * np.sin(theta[i])) for i in range(len(theta))]# Compute star boundary vertices.

    contour_file = 'temp_contour_test_wave.csv'                                                                         # Temporary contour file path.
    cloud_file   = 'temp_cloud_test_wave.csv'                                                                           # Temporary cloud file path.

    with open(contour_file, 'w', newline='') as f:                                                                      # Open temporary contour file.
        writer = csv.writer(f)                                                                                          # Create CSV writer object.
        writer.writerow(['x', 'y'])                                                                                     # Write CSV header row.
        for pt in star_contour: writer.writerow(pt)                                                                     # Write star contour vertices.

    cloud = Cloud.generate_natural(contour_file, cloud_file, cloud_size=0.03, save=False, show=False)                   # Generate natural point cloud via Cloud factory.

    c  = 0.8                                                                                                            # Wave propagation speed.
    ic = lambda x, y: np.exp(-100 * ((x - 0.5)**2 + (y - 0.5)**2))                                                      # Initial position pulse.

    domain = cloud.set_boundary(Dirichlet(0.0))                                                                         # Zero boundary walls.
    pde    = WaveEquation(c=c, ic=ic, g=0.0)                                                                            # Formulate Wave PDE.
    solver = Solver(domain, pde, nvec=16, verbose=False)                                                                # Instantiate Solver.
    res    = solver.solve(t_span=(0.0, 2.0))                                                                            # Execute solver.
    sol    = res.solution                                                                                               # Extract solution array.

    if os.path.exists(contour_file): os.remove(contour_file)                                                            # Clean up temporary contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                # Clean up temporary cloud file.

    assert np.all(np.isfinite(sol))                                                                                     # Verify finite values (no NaN/Inf).
    assert np.max(np.abs(sol)) < 1.5                                                                                    # Verify wave solution remains bounded < 1.5 on irregular cloud.
