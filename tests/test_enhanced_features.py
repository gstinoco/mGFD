"""
Unit tests for mGFD enhanced solver features:
- Custom time domain t_span
- Independent initial (ic) and boundary (bc) conditions
- Independent source / forcing terms (source)
- Reaction terms F_react (6th operator coefficient per Eq. 190 of Tinoco24_r1 paper)
"""

import pytest
import numpy as np

from mGFD.solvers.time_derivative1 import TimeDerivative1
from mGFD.solvers.time_derivative2 import TimeDerivative2
from mGFD.spatial.gammas import compute_sparse_matrix

def test_reaction_coefficient_gammas():
    """Verify that 6th operator coefficient F_react modifies diagonal per Eq. 190: diag[i] = -sum(YY) + F_react."""
    p = np.array([
        [0.0, 0.0, 0], # Interior node
        [1.0, 0.0, 1],
        [-1.0, 0.0, 1],
        [0.0, 1.0, 1],
        [0.0, -1.0, 1]
    ])
    vec = np.array([
        [1, 2, 3, 4],
        [-1, -1, -1, -1],
        [-1, -1, -1, -1],
        [-1, -1, -1, -1],
        [-1, -1, -1, -1]
    ])
    
    # 5-element operator: Laplacian (A=1, C=1) -> [D=0, E=0, A=1, B=0, C=1]
    op5 = np.vstack([[0], [0], [1], [0], [1]])
    K5 = compute_sparse_matrix(p, vec, op5).toarray()
    
    # 6-element operator: Laplacian + Reaction F_react = 3.5 -> [D=0, E=0, A=1, B=0, C=1, F_react=3.5]
    op6 = np.vstack([[0], [0], [1], [0], [1], [3.5]])
    K6 = compute_sparse_matrix(p, vec, op6).toarray()
    
    # Check that diagonal on interior node i=0 is increased exactly by F_react=3.5
    assert np.isclose(K6[0, 0], K5[0, 0] + 3.5)

def test_t_span_custom_domain():
    """Verify that t_span solves across physical time domain [0, 100]."""
    p = np.array([
        [0.0, 0.0, 0],
        [1.0, 0.0, 1],
        [-1.0, 0.0, 1],
        [0.0, 1.0, 1],
        [0.0, -1.0, 1]
    ])
    res = TimeDerivative1(p, f=0.0, t=10, implicit=True, t_span=(0.0, 100.0), verbose=False)
    assert res.solution.shape == (5, 10)
    assert res.dt is not None
    assert np.isclose(res.dt, 100.0 / 10)

def test_independent_ic_and_bc():
    """Verify independent initial condition (ic=20.0) and boundary condition (bc=100.0)."""
    p = np.array([
        [0.0, 0.0, 0],
        [1.0, 0.0, 1],
        [-1.0, 0.0, 1],
        [0.0, 1.0, 1],
        [0.0, -1.0, 1]
    ])
    res = TimeDerivative1(p, ic=20.0, bc=100.0, t=5, implicit=True, verbose=False)
    sol = res.solution
    
    # Initial condition at t=0 across all nodes
    assert np.allclose(sol[:, 0], 20.0)
    # Boundary condition for t > 0 on boundary nodes
    assert np.allclose(sol[1:, 1:], 100.0)

def test_independent_source_term():
    """Verify independent source term F_source(x, y, t) = 5.0 increases interior values."""
    p = np.array([
        [0.0, 0.0, 0],
        [1.0, 0.0, 1],
        [-1.0, 0.0, 1],
        [0.0, 1.0, 1],
        [0.0, -1.0, 1]
    ])
    res_no_source = TimeDerivative1(p, ic=0.0, bc=0.0, source=0.0, t=5, implicit=True, verbose=False)
    res_with_source = TimeDerivative1(p, ic=0.0, bc=0.0, source=5.0, t=5, implicit=True, verbose=False)
    
    # Interior node (row 0) should be > 0 with positive heat generation source term
    assert res_with_source.solution[0, -1] > res_no_source.solution[0, -1]
