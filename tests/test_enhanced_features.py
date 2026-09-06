"""
Test Enhanced Features — Unit tests for mGFD enhanced solver features (OOP Interface)

Overview:
    This test module validates the generalized solver capabilities introduced in mGFD,
    including flexible physical time intervals (t_span), independent initial (ic) and
    boundary (bc) conditions, source forcing terms, and reaction operator coefficients.

Public API:
    test_reaction_coefficient_gammas    Validates reaction coefficient modification on spatial matrix diagonal.
    test_t_span_custom_domain           Validates time domain scaling over custom physical intervals.
    test_independent_ic_and_bc          Validates independent evaluation of initial and boundary conditions.
    test_independent_source_term        Validates independent non-homogeneous source term injection.

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
import pytest                                                                                                                           # Test framework.
import numpy as np                                                                                                                      # Numerical array operations.

from mGFD import Cloud, Dirichlet, Domain, PDE, Solver                                                                                  # High-level OOP API classes.
from mGFD.spatial.gammas import compute_sparse_matrix                                                                                   # Spatial stencil matrix builder.

def test_reaction_coefficient_gammas() -> None:
    """
    test_reaction_coefficient_gammas
    Verify that 6th operator coefficient F_react modifies diagonal per Eq. 190: diag[i] = -sum(YY) + F_react.
    """
    # 1. Point cloud and neighbor definitions
    p = np.array([
        [0.0, 0.0, 0],                                                                                                                  # Interior node.
        [1.0, 0.0, 1],                                                                                                                  # Boundary node 1.
        [-1.0, 0.0, 1],                                                                                                                 # Boundary node 2.
        [0.0, 1.0, 1],                                                                                                                  # Boundary node 3.
        [0.0, -1.0, 1]                                                                                                                  # Boundary node 4.
    ])
    vec = np.array([
        [1, 2, 3, 4],                                                                                                                   # Neighbor connectivity for node 0.
        [-1, -1, -1, -1],                                                                                                               # Null neighbors for boundary 1.
        [-1, -1, -1, -1],                                                                                                               # Null neighbors for boundary 2.
        [-1, -1, -1, -1],                                                                                                               # Null neighbors for boundary 3.
        [-1, -1, -1, -1]                                                                                                                # Null neighbors for boundary 4.
    ])
    
    # 2. Operator assembly and matrix generation
    op5 = np.vstack([[0], [0], [1], [0], [1]])                                                                                          # 5-element Laplacian operator.
    K5  = compute_sparse_matrix(p, vec, op5).toarray()                                                                                  # Base spatial matrix.
    
    op6 = np.vstack([[0], [0], [1], [0], [1], [3.5]])                                                                                   # 6-element operator with reaction.
    K6  = compute_sparse_matrix(p, vec, op6).toarray()                                                                                  # Matrix with reaction coefficient.
    
    # 3. Verification
    assert np.isclose(K6[0, 0], K5[0, 0] + 3.5)                                                                                         # Verify diagonal shift equals F_react.

def test_t_span_custom_domain() -> None:
    """
    test_t_span_custom_domain
    Verify that t_span solves across physical time domain [0, 100].
    """
    # 1. Formulation
    p      = np.array([[0.0, 0.0, 0], [1.0, 0.0, 1], [-1.0, 0.0, 1], [0.0, 1.0, 1], [0.0, -1.0, 1]])                                    # 5-node test cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(0.0))                                                                                         # Set zero Dirichlet boundary.
    pde    = PDE(operator=[0, 0, 2, 0, 2, 0], ic=0.0, order=1)                                                                          # Formulate 1st order PDE.
    solver = Solver(domain, pde, verbose=False)                                                                                         # Create Solver.
    
    # 2. Execution and Assertion
    res    = solver.solve(t_span=(0.0, 100.0), dt=10.0)                                                                                 # Solve across custom t_span.
    assert res.solution.shape == (5, 11)                                                                                                # Verify solution matrix dimensions.
    assert res.dt is not None                                                                                                           # Verify non-null time step.
    assert np.isclose(res.dt, 100.0 / 10)                                                                                               # Verify time step equals 10.0.

def test_independent_ic_and_bc() -> None:
    """
    test_independent_ic_and_bc
    Verify independent initial condition (ic=20.0) and boundary condition (bc=100.0).
    """
    # 1. Formulation
    p      = np.array([[0.0, 0.0, 0], [1.0, 0.0, 1], [-1.0, 0.0, 1], [0.0, 1.0, 1], [0.0, -1.0, 1]])                                    # 5-node test cloud.
    cloud  = Cloud.from_array(p)                                                                                                        # Instantiate Cloud.
    domain = cloud.set_boundary(Dirichlet(100.0))                                                                                       # Set boundary condition bc=100.0.
    pde    = PDE(operator=[0, 0, 2, 0, 2, 0], ic=20.0, order=1)                                                                         # Set initial condition ic=20.0.
    solver = Solver(domain, pde, verbose=False)                                                                                         # Create Solver.
    
    # 2. Execution and Assertion
    res    = solver.solve(t_span=(0.0, 1.0), dt=0.25)                                                                                   # Solve transient system.
    sol    = res.solution                                                                                                               # Extract solution matrix.
    
    assert np.allclose(sol[:, 0], 20.0)                                                                                                 # Verify ic at t=0 across all nodes.
    assert np.allclose(sol[1:, 1:], 100.0)                                                                                              # Verify bc for t > 0 on boundary nodes.

def test_independent_source_term() -> None:
    """
    test_independent_source_term
    Verify independent source term F_source(x, y, t) = 5.0 increases interior values.
    """
    # 1. Formulation
    p               = np.array([[0.0, 0.0, 0], [1.0, 0.0, 1], [-1.0, 0.0, 1], [0.0, 1.0, 1], [0.0, -1.0, 1]])                           # 5-node test cloud.
    cloud           = Cloud.from_array(p)                                                                                               # Instantiate Cloud.
    domain          = cloud.set_boundary(Dirichlet(0.0))                                                                                # Set zero Dirichlet boundary.
    
    # 2. Solves with and without source
    pde_no_source   = PDE(operator=[0, 0, 2, 0, 2, 0], ic=0.0, source=0.0, order=1)                                                     # PDE without source.
    res_no_source   = Solver(domain, pde_no_source, verbose=False).solve(t_span=(0.0, 1.0), dt=0.25)                                    # Baseline solve.

    pde_with_source = PDE(operator=[0, 0, 2, 0, 2, 0], ic=0.0, source=5.0, order=1)                                                     # PDE with source F=5.0.
    res_with_source = Solver(domain, pde_with_source, verbose=False).solve(t_span=(0.0, 1.0), dt=0.25)                                  # Solve with source term.
    
    # 3. Verification
    assert res_with_source.solution[0, -1] > res_no_source.solution[0, -1]                                                              # Interior node value must increase.
