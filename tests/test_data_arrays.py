"""
Test Data Arrays — Unit tests for the data array handling in OOP solvers

Overview:
    This file contains the unit tests for validating that the mGFD OOP solvers can handle both callable functions
    and pre-computed data arrays consistently, yielding mathematically equivalent solutions.

Public API:
    generate_square_cloud                   Generates a square point cloud for testing.
    test_stationary_arrays                  Validates data array handling in OOP Stationary solver.
    test_time_derivative1_arrays            Validates data array handling in OOP 1st-order transient solver.
    test_time_derivative2_arrays            Validates data array handling in OOP 2nd-order transient solver.

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
from mGFD import Cloud, Dirichlet, PDE, Solver                                                                                          # OOP Architecture classes.

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
    Validates array handling in the OOP Stationary solver.
    """
    # 1. Test initialization
    p        = generate_square_cloud(15)                                                                                                # Generate point cloud.
    
    phi_func = lambda x, y: x**2 + y**2                                                                                                 # Callable boundary condition.
    f_func   = lambda x, y: -4.0                                                                                                        # Callable forcing term.
    
    phi_arr  = p[:, 0]**2 + p[:, 1]**2                                                                                                  # Array boundary condition.
    f_const  = -4.0                                                                                                                     # Constant forcing term.
    
    # 2. Execution
    cloud    = Cloud.from_array(p)                                                                                                      # Instantiate Cloud.
    
    dom_func = cloud.set_boundary(Dirichlet(phi_func))                                                                                  # Callable boundary domain.
    pde_func = PDE(operator=[0, 0, 2, 0, 2, 0], source=f_func, order=0)                                                                 # Callable PDE.
    res_func = Solver(dom_func, pde_func, nvec=12, verbose=False).solve()                                                               # Solve with Callable.

    dom_arr  = cloud.set_boundary(Dirichlet(phi_arr))                                                                                   # Array boundary domain.
    pde_arr  = PDE(operator=[0, 0, 2, 0, 2, 0], source=f_const, order=0)                                                                # Constant PDE.
    res_arr  = Solver(dom_arr, pde_arr, nvec=12, verbose=False).solve()                                                                 # Solve with Array/Constant.
    
    # 3. Assertions
    error    = np.max(np.abs(res_func.solution - res_arr.solution))                                                                     # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.

def test_time_derivative1_arrays() -> None:
    """
    test_time_derivative1_arrays
    Validates array handling in the OOP 1st-order transient PDE solver.
    """
    # 1. Test initialization
    p       = generate_square_cloud(11)                                                                                                 # Generate point cloud.
    t_steps = 10                                                                                                                        # Number of time steps.
    v       = 0.2                                                                                                                       # Thermal diffusivity.
    T       = np.linspace(0, 1, t_steps)                                                                                                # Time vector.
    
    f_func = lambda x, y, t=0, coef=None: np.exp(-2 * np.pi**2 * v * t) * np.cos(np.pi * x) * np.cos(np.pi * y)                         # Callable boundary function.
    ic     = lambda x, y: np.cos(np.pi * x) * np.cos(np.pi * y)                                                                         # Initial state.
    
    f_arr  = np.zeros((len(p), t_steps))                                                                                                # Spatiotemporal array.
    for k in range(t_steps):                                                                                                            # Iterate over time.
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k])                                                                                    # Populate array values.
        
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    dt_val = 1.0 / (t_steps - 1)                                                                                                        # Step size matching T grid.
    
    # 2. Execution
    dom_func = cloud.set_boundary(Dirichlet(f_func))                                                                                    # Callable domain.
    pde_func = PDE(operator=[0, 0, 0.4, 0, 0.4, 0], ic=ic, order=1)                                                                     # Heat PDE.
    res_func = Solver(dom_func, pde_func, verbose=False).solve(t_span=(0.0, 1.0), dt=dt_val)                                            # Solve with Callable.

    dom_arr  = cloud.set_boundary(Dirichlet(f_arr))                                                                                     # Array domain.
    pde_arr  = PDE(operator=[0, 0, 0.4, 0, 0.4, 0], ic=ic, order=1)                                                                     # Heat PDE.
    res_arr  = Solver(dom_arr, pde_arr, verbose=False).solve(t_span=(0.0, 1.0), dt=dt_val)                                              # Solve with Array.
    
    # 3. Assertions
    error = np.max(np.abs(res_func.solution - res_arr.solution))                                                                        # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.

def test_time_derivative2_arrays() -> None:
    """
    test_time_derivative2_arrays
    Validates array handling in the OOP 2nd-order transient PDE solver.
    """
    # 1. Test initialization
    p       = generate_square_cloud(11)                                                                                                 # Generate point cloud.
    t_steps = 10                                                                                                                        # Number of time steps.
    c       = float(np.sqrt(1 / 2))                                                                                                     # Wave speed.
    T       = np.linspace(0, 1, t_steps)                                                                                                # Time vector.
    
    f_func = lambda x, y, t=0, coef=None: np.cos(np.pi * t) * np.sin(np.pi * (x + y))                                                   # Callable boundary function.
    ic     = lambda x, y: np.sin(np.pi * (x + y))                                                                                       # Initial position function.
    g_func = lambda x, y, t=0, coef=None: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))                                          # Callable velocity function.
    
    f_arr = np.zeros((len(p), t_steps))                                                                                                 # Spatiotemporal array for f.
    for k in range(t_steps):                                                                                                            # Iterate over time.
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k])                                                                                    # Populate array values.
    
    g_arr  = g_func(p[:, 0], p[:, 1], T[0])                                                                                             # Velocity 'g' at t=0.
    
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    dt_val = 1.0 / (t_steps - 1)                                                                                                        # Step size matching T grid.
    
    # 2. Execution
    dom_func = cloud.set_boundary(Dirichlet(f_func))                                                                                    # Callable domain.
    pde_func = PDE(operator=[0, 0, 1.0, 0, 1.0, 0], ic=ic, g=g_func, order=2)                                                           # Wave PDE.
    res_func = Solver(dom_func, pde_func, verbose=False).solve(t_span=(0.0, 1.0), dt=dt_val)                                            # Solve with Callable.

    dom_arr  = cloud.set_boundary(Dirichlet(f_arr))                                                                                     # Array domain.
    pde_arr  = PDE(operator=[0, 0, 1.0, 0, 1.0, 0], ic=ic, g=g_arr, order=2)                                                            # Wave PDE.
    res_arr  = Solver(dom_arr, pde_arr, verbose=False).solve(t_span=(0.0, 1.0), dt=dt_val)                                              # Solve with Array.
    
    # 3. Assertions
    error = np.max(np.abs(res_func.solution - res_arr.solution))                                                                        # Compute max difference.
    assert error < 1e-12                                                                                                                # Verify mathematical equivalence.
