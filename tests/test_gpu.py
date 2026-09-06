"""
Test GPU — Unit tests for the GPU-accelerated routines of mGFD (OOP Interface)

Overview:
    This file contains the unit tests for validating the GPU-accelerated algorithms of mGFD using the OOP Architecture,
    ensuring that numerical and computational behaviors remain stable on CUDA hardware.

Public API:
    generate_square_cloud                           Generates a square point cloud for testing.
    test_gpu_stationary                             Validates GPU acceleration in OOP Stationary solver.
    test_gpu_time_derivative1                       Validates GPU acceleration in OOP 1st-order transient solver.
    test_gpu_time_derivative2                       Validates GPU acceleration in OOP 2nd-order transient solver.
    test_gpu_without_cupy_raises_importerror        Validates ImportError when CuPy is unavailable.

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

try:
    import cupy as cp                                                                                                                   # CuPy library for GPU.
    HAS_CUPY = True                                                                                                                     # CuPy availability flag.
except ImportError:
    HAS_CUPY = False                                                                                                                    # CuPy availability flag.

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

@pytest.mark.skipif(not HAS_CUPY, reason="CuPy not installed")
def test_gpu_stationary() -> None:
    """
    test_gpu_stationary
    Validates the GPU Stationary solver execution and shape correctness in OOP mode.
    """
    # 1. Test initialization
    p      = generate_square_cloud(11)                                                                                                  # Generate point cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                                         # Set boundary condition.
    pde    = PDE(operator=[0, 0, 2, 0, 2, 0], source=1.0, order=0)                                                                      # Formulate Poisson PDE.
    
    # 2. Execution
    solver = Solver(domain, pde, device="cuda", verbose=False)                                                                          # Instantiate Solver on GPU.
    result = solver.solve()                                                                                                             # Execute solver.
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape[0] == 121                                                                                              # Solution size must match nodes.

@pytest.mark.skipif(not HAS_CUPY, reason="CuPy not installed")
def test_gpu_time_derivative1() -> None:
    """
    test_gpu_time_derivative1
    Validates the GPU HeatEquation solver execution and shape correctness in OOP mode.
    """
    # 1. Test initialization
    p      = generate_square_cloud(11)                                                                                                  # Generate point cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                                         # Set boundary condition.
    pde    = PDE(operator=[0, 0, 0.4, 0, 0.4, 0], order=1)                                                                              # Formulate Heat PDE.
    
    # 2. Execution
    solver = Solver(domain, pde, device="cuda", verbose=False)                                                                          # Instantiate Solver on GPU.
    result = solver.solve(t_span=(0.0, 1.0), dt=1.0)                                                                                    # Execute solver with 2 time steps (dt=1.0).
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape == (121, 2)                                                                                            # Solution size must match (nodes, time_steps).

@pytest.mark.skipif(not HAS_CUPY, reason="CuPy not installed")
def test_gpu_time_derivative2() -> None:
    """
    test_gpu_time_derivative2
    Validates the GPU WaveEquation solver execution and shape correctness in OOP mode.
    """
    # 1. Test initialization
    p      = generate_square_cloud(11)                                                                                                  # Generate point cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                                         # Set boundary condition.
    pde    = PDE(operator=[0, 0, 0.5, 0, 0.5, 0], g=0.0, order=2)                                                                       # Formulate Wave PDE.
    
    # 2. Execution
    solver = Solver(domain, pde, device="cuda", verbose=False)                                                                          # Instantiate Solver on GPU.
    result = solver.solve(t_span=(0.0, 1.0), dt=1.0)                                                                                    # Execute solver with 2 time steps (dt=1.0).
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape == (121, 2)                                                                                            # Solution size must match (nodes, time_steps).

def test_gpu_without_cupy_raises_importerror(monkeypatch) -> None:
    """
    test_gpu_without_cupy_raises_importerror
    Validates fallback behavior when CuPy is missing but CUDA is requested.
    """
    # 1. Test initialization
    import sys                                                                                                                          # System module.
    monkeypatch.setitem(sys.modules, "cupy", None)                                                                                      # Mock CuPy unavailability.
    
    p      = generate_square_cloud(3)                                                                                                   # Generate small point cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                                         # Set boundary condition.
    pde    = PDE(operator=[0, 0, 2, 0, 2, 0], source=1.0, order=0)                                                                      # Formulate Poisson PDE.
    solver = Solver(domain, pde, device="cuda", verbose=False)                                                                          # Instantiate Solver on GPU.
    
    # 2. Execution and Assertions
    with pytest.raises(ImportError, match="CuPy is not installed"):                                                                     # Assert exception raised.
        solver.solve()                                                                                                                  # Execute solver with CUDA request.
