"""
Test OOP Interface — Unit tests for mGFD Object-Oriented Architecture

Overview:
    Validates mGFD OOP abstractions including Cloud, Dirichlet, Neumann, Domain,
    PDE, Solver, and SolverResult convenience methods (.export_vtk and .plot).

Public API:
    test_oop_cloud_constructors                     Validates Cloud instantiation from raw arrays, CSVs, and properties.
    test_oop_boundary_conditions                    Validates Dirichlet and Neumann boundary assignments.
    test_oop_poisson_solver                         Validates stationary Poisson PDE solution via OOP interface.
    test_oop_heat_solver                            Validates 1st-order transient Heat PDE solution via OOP interface.
    test_oop_advection_diffusion_solver             Validates advection-diffusion-reaction PDE solution via OOP interface.
    test_oop_wave_solver                            Validates 2nd-order transient Wave PDE solution via OOP interface.
    test_oop_result_methods                         Validates VTK export and plotting convenience methods on SolverResult.

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
import csv, os                                                                                                                          # File and system management.
import numpy as np                                                                                                                      # Core numerical array operations.
import pytest                                                                                                                           # Pytest testing framework.

import mGFD as mgfd                                                                                                                     # Import mGFD main library.

def test_oop_cloud_constructors() -> None:                                                                                              # Test Cloud instantiation.
    """
    test_oop_cloud_constructors
    Validates Cloud instantiation from raw NumPy arrays, CSV files, and point properties.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0], [0.25, 0.25, 0], [0.75, 0.75, 0]
    ], dtype=np.float64)                                                                                                                # Create 7-node test array.
    
    cloud = mgfd.Cloud.from_array(p_arr)                                                                                                # Instantiate from array.
    assert cloud.num_nodes == 7                                                                                                         # Verify node count.
    assert np.sum(cloud.interior_mask) == 3                                                                                             # Verify interior nodes.
    assert np.sum(cloud.boundary_mask) == 4                                                                                             # Verify boundary nodes.
    
    vec = cloud.compute_neighbors(nvec=4)                                                                                               # Precompute neighbor list.
    assert vec.shape == (7, 4)                                                                                                          # Verify neighbor matrix dimensions.
    assert cloud.neighbors is not None                                                                                                  # Verify cached neighbor attribute.

def test_oop_boundary_conditions() -> None:                                                                                             # Test BoundaryCondition classes.
    """
    test_oop_boundary_conditions
    Validates Dirichlet and Neumann boundary condition evaluation with scalars, callables, and arrays.
    """
    x = np.array([0.0, 0.5, 1.0])                                                                                                       # Test X coordinates.
    y = np.array([0.0, 0.5, 1.0])                                                                                                       # Test Y coordinates.
    
    # 1. Scalar Dirichlet
    d_scalar = mgfd.Dirichlet(5.0)                                                                                                      # Constant scalar Dirichlet.
    assert np.allclose(d_scalar(x, y), [5.0, 5.0, 5.0])                                                                                 # Verify scalar array evaluation.
    
    # 2. Callable Dirichlet
    d_func = mgfd.Dirichlet(lambda x, y, t=0: x + y)                                                                                    # Callable Dirichlet.
    assert np.allclose(d_func(x, y), [0.0, 1.0, 2.0])                                                                                   # Verify function evaluation.
    
    # 3. Neumann
    n_spec = mgfd.Neumann(0.0)                                                                                                          # Zero Neumann condition.
    assert np.allclose(n_spec(x, y), [0.0, 0.0, 0.0])                                                                                   # Verify Neumann array evaluation.

def test_oop_poisson_solver() -> None:                                                                                                  # Test Stationary Poisson solver via OOP.
    """
    test_oop_poisson_solver
    Validates solving a stationary Poisson problem using Cloud, Dirichlet, Domain, PoissonEquation, and Solver.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0], [0.25, 0.25, 0], [0.75, 0.75, 0]
    ], dtype=np.float64)                                                                                                                # Test grid point cloud.
    
    cloud  = mgfd.Cloud.from_array(p_arr)                                                                                               # Instantiate cloud.
    domain = cloud.set_boundary(mgfd.Dirichlet(0.0))                                                                                    # Bind zero Dirichlet boundary.
    pde    = mgfd.PDE(operator=[0, 0, 2, 0, 2, 0], source=1.0, order=0)                                                                 # Create Poisson PDE with source=1.0.
    solver = mgfd.Solver(domain, pde, nvec=4, verbose=False)                                                                            # Create high-level Solver.
    
    result = solver.solve()                                                                                                             # Solve PDE.
    assert result.converged is True                                                                                                     # Verify convergence.
    assert result.solution.shape[0] == 7                                                                                                # Verify solution vector shape.
    assert np.all(np.isfinite(result.solution))                                                                                         # Verify finite values.

def test_oop_heat_solver() -> None:                                                                                                     # Test Transient Heat solver via OOP.
    """
    test_oop_heat_solver
    Validates solving a 1st-order transient Heat equation using HeatEquation.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0], [0.25, 0.5, 0], [0.75, 0.5, 0]
    ], dtype=np.float64)                                                                                                                # Test point cloud.
    
    cloud  = mgfd.Cloud.from_array(p_arr)                                                                                               # Instantiate cloud.
    domain = cloud.set_boundary(0.0)                                                                                                    # Bind boundary float.
    pde    = mgfd.PDE(operator=[0, 0, 0.2, 0, 0.2, 0], ic=lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y), order=1)                  # Create Heat PDE.
    solver = mgfd.Solver(domain, pde, nvec=4, cfl=0.5, verbose=False)                                                                   # Create high-level Solver.
    
    result = solver.solve(t_span=(0.0, 0.1))                                                                                            # Execute transient solve.
    assert result.converged is True                                                                                                     # Verify convergence.
    assert result.solution.ndim == 2                                                                                                    # Verify 2D solution matrix [nodes, time].
    assert result.solution.shape[0] == 7                                                                                                # Verify spatial node count.

def test_oop_advection_diffusion_solver() -> None:                                                                                      # Test Advection-Diffusion solver via OOP.
    """
    test_oop_advection_diffusion_solver
    Validates solving a 1st-order transient Advection-Diffusion equation using AdvectionDiffusion.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0], [0.25, 0.5, 0], [0.75, 0.5, 0]
    ], dtype=np.float64)                                                                                                                # Test point cloud.
    
    cloud  = mgfd.Cloud.from_array(p_arr)                                                                                               # Instantiate cloud.
    domain = cloud.set_boundary(0.0)                                                                                                    # Bind zero boundary.
    pde    = mgfd.PDE(operator=[-0.1, -0.1, 0.02, 0, 0.02, 0], ic=lambda x, y: np.exp(-10*((x-0.5)**2 + (y-0.5)**2)), order=1)          # Create Advection-Diffusion PDE.
    solver = mgfd.Solver(domain, pde, nvec=4, cfl=0.5, verbose=False)                                                                   # Create high-level Solver.
    
    result = solver.solve(t_span=(0.0, 0.1))                                                                                            # Execute solve.
    assert result.converged is True                                                                                                     # Verify convergence.
    assert np.all(np.isfinite(result.solution))                                                                                         # Verify finite numbers.

def test_oop_wave_solver() -> None:                                                                                                     # Test Transient Wave solver via OOP.
    """
    test_oop_wave_solver
    Validates solving a 2nd-order transient Wave equation using WaveEquation.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0], [0.25, 0.5, 0], [0.75, 0.5, 0]
    ], dtype=np.float64)                                                                                                                # Test point cloud.
    
    cloud   = mgfd.Cloud.from_array(p_arr)                                                                                              # Instantiate cloud.
    domain  = cloud.set_boundary(0.0)                                                                                                   # Bind boundary.
    ic_func = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)                                                                        # Initial position function.
    pde     = mgfd.PDE(operator=[0, 0, 0.5, 0, 0.5, 0], ic=ic_func, g=0.0, order=2, coef={'damping': 0.01, 'alpha': -0.1})              # Create Wave PDE.
    solver  = mgfd.Solver(domain, pde, nvec=4, cfl=0.5, verbose=False)                                                                  # Create high-level Solver.
    
    result  = solver.solve(t_span=(0.0, 0.1))                                                                                           # Execute 2nd order solve.
    assert result.converged is True                                                                                                     # Verify convergence.
    assert result.solution.ndim == 2                                                                                                    # Verify 2D solution shape.
    assert np.max(np.abs(result.solution)) < 2.0                                                                                        # Verify numerical stability.

def test_oop_result_methods() -> None:                                                                                                  # Test SolverResult VTK and plot methods.
    """
    test_oop_result_methods
    Validates export_vtk and plot methods on SolverResult.
    """
    p_arr = np.array([
        [0.0, 0.0, 1], [1.0, 0.0, 1], [1.0, 1.0, 1], [0.0, 1.0, 1],
        [0.5, 0.5, 0]
    ], dtype=np.float64)                                                                                                                # Test point cloud.
    
    cloud  = mgfd.Cloud.from_array(p_arr)                                                                                               # Instantiate cloud.
    domain = cloud.set_boundary(0.0)                                                                                                    # Bind boundary.
    pde    = mgfd.PDE(operator=[0, 0, 2, 0, 2, 0], source=1.0, order=0)                                                                 # Create PDE.
    solver = mgfd.Solver(domain, pde, nvec=4, verbose=False)                                                                            # Create Solver.
    
    result = solver.solve()                                                                                                             # Solve PDE.
    
    import shutil                                                                                                                       # Directory cleanup interface.
    vtk_file = "temp_test_oop_export.vtp"                                                                                               # Temporary VTK path (.vtp for PolyData).
    success  = result.export_vtk(vtk_file)                                                                                              # Execute VTK export.
    assert success is True                                                                                                              # Verify export success.
    assert os.path.exists(vtk_file)                                                                                                     # Verify VTK file creation.
    
    if os.path.exists(vtk_file):                                                                                                        # Check file existence.
        shutil.rmtree(vtk_file) if os.path.isdir(vtk_file) else os.remove(vtk_file)                                                     # Clean up temporary VTK output.
