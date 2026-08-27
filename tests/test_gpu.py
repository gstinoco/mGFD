"""
Test GPU — Unit tests for the GPU-accelerated routines of mGFD

Overview:
    This file contains the unit tests for validating the GPU-accelerated algorithms of mGFD, ensuring that the 
    numerical and computational behaviors remain stable and compliant with the theoretical bounds of the method.

Public API:
    test_gpu_stationary                 Validates GPU Stationary solver behavior.
    test_gpu_time_derivative1           Validates GPU TimeDerivative1 solver behavior.
    test_gpu_time_derivative2           Validates GPU TimeDerivative2 solver behavior.
    test_gpu_without_cupy_raises_importerror Validates fallback behavior when CuPy is missing.
    test_gpu_matrix_free_conflict       Validates parameter validation for matrix-free with CUDA.

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
from mGFD.exceptions import ParameterError                                                                                              # Exceptions.

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
    Validates the GPU Stationary solver execution and shape correctness.
    
    This test ensures that when the device='cuda' flag is provided, the solver
    successfully performs the computation using CuPy routines without failure,
    returning a solution vector of the exact length matching the input cloud.
    """
    # 1. Test initialization
    p        = generate_square_cloud(11)                                                                                                # Generate point cloud.
    phi      = 0.0                                                                                                                      # Boundary condition.
    f        = 1.0                                                                                                                      # Forcing term.
    operator = np.array([0, 0, 1, 0, 1, 0], dtype=np.float64)                                                                           # Laplacian operator.
    
    # 2. Execution
    result   = Stationary(p, phi, f, operator=operator, linear_solver="spsolve", device="cuda", verbose=False)                          # Run Stationary solver on GPU.
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape[0] == 121                                                                                              # Solution size must match nodes.

@pytest.mark.skipif(not HAS_CUPY, reason="CuPy not installed")
def test_gpu_time_derivative1() -> None:
    """
    test_gpu_time_derivative1
    Validates the GPU TimeDerivative1 solver execution and shape correctness.
    
    This test runs a simplified parabolic equation simulation over a short duration
    to assert that the underlying matrix factorizations and sparse multiplications
    are appropriately mapped to GPU memory without crashing.
    """
    # 1. Test initialization
    p        = generate_square_cloud(11)                                                                                                # Generate point cloud.
    f_cond   = lambda x, y, t, coef: 0.0                                                                                                # Boundary condition function.
    operator = np.array([0, 0, 1, 0, 1, 0], dtype=np.float64)                                                                           # Laplacian operator.
    
    # 2. Execution
    result   = TimeDerivative1(p, f_cond, 2, [], operator=operator, implicit=True, linear_solver="bicgstab", preconditioner="jacobi", device="cuda", verbose=False) # Run TimeDerivative1 solver on GPU.
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape == (121, 2)                                                                                            # Solution size must match (nodes, time_steps).

@pytest.mark.skipif(not HAS_CUPY, reason="CuPy not installed")
def test_gpu_time_derivative2() -> None:
    """
    test_gpu_time_derivative2
    Validates the GPU TimeDerivative2 solver execution and shape correctness.
    
    This test checks the hyperbolic PDE solver on the GPU, verifying that both initial
    conditions and velocity conditions are transferred and computed correctly inside
    the time-stepping loop on device memory.
    """
    # 1. Test initialization
    p        = generate_square_cloud(11)                                                                                                # Generate point cloud.
    f_cond   = lambda x, y, t, coef: 0.0                                                                                                # Initial condition function.
    g_cond   = lambda x, y, t, coef: 0.0                                                                                                # Initial velocity function.
    operator = np.array([0, 0, 1, 0, 1, 0], dtype=np.float64)                                                                           # Laplacian operator.
    
    # 2. Execution
    result   = TimeDerivative2(p, f_cond, g_cond, 2, [], operator=operator, implicit=True, linear_solver="gmres", preconditioner="jacobi", device="cuda", verbose=False) # Run TimeDerivative2 solver on GPU.
    
    # 3. Assertions
    assert result.converged                                                                                                             # Solver must converge.
    assert result.solution.shape == (121, 2)                                                                                            # Solution size must match (nodes, time_steps).

def test_gpu_without_cupy_raises_importerror(monkeypatch) -> None:
    """
    test_gpu_without_cupy_raises_importerror
    Validates fallback behavior when CuPy is missing but CUDA is requested.
    
    By artificially removing CuPy from sys.modules, this test ensures the system
    safely traps the import error and throws a descriptive exception rather than
    failing mysteriously deep in the solver logic.
    """
    # 1. Test initialization
    import sys                                                                                                                          # System module.
    monkeypatch.setitem(sys.modules, "cupy", None)                                                                                      # Mock CuPy unavailability.
    
    p        = generate_square_cloud(3)                                                                                                 # Generate small point cloud.
    phi      = 0.0                                                                                                                      # Boundary condition.
    f        = 1.0                                                                                                                      # Forcing term.
    operator = np.array([0, 0, 1, 0, 1, 0], dtype=np.float64)                                                                           # Laplacian operator.
    
    # 2. Execution and Assertions
    with pytest.raises(ImportError, match="CuPy is not installed"):                                                                     # Assert exception raised.
        Stationary(p, phi, f, operator=operator, device="cuda", verbose=False)                                                          # Run solver with CUDA request.

def test_gpu_matrix_free_conflict() -> None:
    """
    test_gpu_matrix_free_conflict
    Validates parameter validation for matrix-free with CUDA.
    
    Currently, the matrix-free iteration strategy is optimized for Numba JIT on the CPU.
    This test verifies that asking for both 'matrix_free=True' and 'device=cuda'
    immediately triggers a ParameterError.
    """
    # 1. Test initialization
    p        = generate_square_cloud(3)                                                                                                 # Generate small point cloud.
    phi      = 0.0                                                                                                                      # Boundary condition.
    f        = 1.0                                                                                                                      # Forcing term.
    operator = np.array([0, 0, 1, 0, 1, 0], dtype=np.float64)                                                                           # Laplacian operator.
    
    # 2. Execution and Assertions
    with pytest.raises(ParameterError, match="incompatible with device='cuda'"):                                                        # Assert exception raised.
        Stationary(p, phi, f, operator=operator, device="cuda", matrix_free=True, verbose=False)                                        # Run solver with incompatible flags.
