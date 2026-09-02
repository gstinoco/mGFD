"""
Test Data Arrays — Unit tests for the data array handling in solvers

Overview:
    This file contains the unit tests for validating that the mGFD solvers can handle both callable functions
    and pre-computed data arrays consistently, yielding mathematically equivalent solutions.

Public API:
    test_stationary_arrays              Validates array handling in the Stationary solver.
    test_time_derivative1_arrays        Validates array handling in the TimeDerivative1 solver.
    test_time_derivative2_arrays        Validates array handling in the TimeDerivative2 solver.

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
from mGFD import Stationary, TimeDerivative1, TimeDerivative2                                                                           # mGFD solvers.

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

def test_stationary_arrays() -> None:
    """
    test_stationary_arrays
    Validates array handling in the Stationary solver.
    
    Verifies that providing pre-computed numpy arrays (or constants) for boundary
    conditions and forcing terms yields mathematically identical results as providing
    callable Python functions, within a strict tolerance of 1e-12.
    """
    # 1. Test initialization
    p = generate_square_cloud(15)                                                                                                       # Generate point cloud.
    
    phi_func = lambda x, y: x**2 + y**2                                                                                                 # Callable boundary condition.
    f_func   = lambda x, y: -4.0                                                                                                        # Callable forcing term.
    
    phi_arr  = p[:, 0]**2 + p[:, 1]**2                                                                                                  # Array boundary condition.
    f_const  = -4.0                                                                                                                     # Constant forcing term.
    
    op = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                                      # Laplacian operator.
    
    # 2. Execution
    u_ap_func, _ = Stationary(p, phi_func, f_func, operator=op, nvec=12, verbose=False)                                                 # Run with Callable.
    u_ap_arr, _  = Stationary(p, phi_arr, f_const, operator=op, nvec=12, verbose=False)                                                 # Run with Arrays and Constants.
    
    # 3. Assertions
    error = np.max(np.abs(u_ap_func - u_ap_arr))                                                                                        # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.

def test_time_derivative1_arrays() -> None:
    """
    test_time_derivative1_arrays
    Validates array handling in the TimeDerivative1 solver.
    
    Tests the substitution of a spatiotemporal callable function with a pre-computed
    dense array covering all spatial nodes and time steps, guaranteeing the internal
    interpolation mechanism handles indexing properly without loss of precision.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                                       # Generate point cloud.
    t_steps = 10                                                                                                                        # Number of time steps.
    v = 0.2                                                                                                                             # Thermal diffusivity.
    T = np.linspace(0, 1, t_steps)                                                                                                      # Time vector.
    
    f_func = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)                          # Callable forcing function.
    
    f_arr = np.zeros((len(p), t_steps))                                                                                                 # Spatiotemporal array.
    for k in range(t_steps):                                                                                                            # Iterate over time.
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k], [v])                                                                               # Populate array values.
        
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                               # Heat operator.
    
    # 2. Execution
    u_ap_func, _ = TimeDerivative1(p, f_func, t_steps, [v], operator=L, implicit=True, verbose=False)                                   # Run with Callable.
    u_ap_arr, _  = TimeDerivative1(p, f_arr, t_steps, [v], operator=L, implicit=True, verbose=False)                                    # Run with Array.
    
    # 3. Assertions
    error = np.max(np.abs(u_ap_func - u_ap_arr))                                                                                        # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.

def test_time_derivative2_arrays() -> None:
    """
    test_time_derivative2_arrays
    Validates array handling in the TimeDerivative2 solver.
    
    Checks that the hyperbolic solver correctly processes both a spatiotemporal array
    for initial position constraints and a purely spatial array for initial velocities,
    matching the functional implementations within near-zero tolerance bounds.
    """
    # 1. Test initialization
    p = generate_square_cloud(11)                                                                                                       # Generate point cloud.
    t_steps = 10                                                                                                                        # Number of time steps.
    c = float(np.sqrt(1 / 2))                                                                                                           # Wave speed.
    T = np.linspace(0, 1, t_steps)                                                                                                      # Time vector.
    
    f_func = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                                                          # Callable initial condition.
    g_func = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))                                                 # Callable initial velocity.
    
    f_arr = np.zeros((len(p), t_steps))                                                                                                 # Spatiotemporal array for f.
    g_arr = np.zeros(len(p))                                                                                                            # Spatial array for g.
    
    for k in range(t_steps):                                                                                                            # Iterate over time.
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k], [c])                                                                               # Populate array values for f.
    
    g_arr = g_func(p[:, 0], p[:, 1], T[0], [c])                                                                                         # Velocity 'g' is evaluated at t=0 (k=0) in solver.
    
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                                                         # Wave operator.
    
    # 2. Execution
    u_ap_func, _ = TimeDerivative2(p, f_func, g_func, t_steps, [c], operator=L, implicit=True, verbose=False)                           # Run with Callable.
    u_ap_arr, _  = TimeDerivative2(p, f_arr, g_arr, t_steps, [c], operator=L, implicit=True, verbose=False)                             # Run with Array.
    
    # 3. Assertions
    error = np.max(np.abs(u_ap_func - u_ap_arr))                                                                                        # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.
