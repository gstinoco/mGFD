"""
Test Neighbors — Unit tests for the meshless neighbor search routines

Overview:
    This file contains the unit tests for validating the geometric and adaptive neighbor
    search algorithms of mGFD, ensuring that the stencil distributions remain stable and
    compliant with the theoretical bounds of the method.

Public API:
    test_find_distances                 Validates the isotropic search radius calculation.
    test_compute_neighbors_adaptive     Validates the Round-Robin selection and geometric stability.
    test_compute_upwind_neighbors       Validates the biased convective directionality search.

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
import numpy as np                                                                                                      # Core numerical operations.
import pytest                                                                                                           # Unit testing framework.

from mGFD.spatial.neighbors import compute_neighbors, compute_upwind_neighbors, find_distances                          # Neighbor search core modules.

def generate_square_cloud(n: int = 11) -> np.ndarray:
    """
    generate_square_cloud
    Helper function to generate a normalized square point cloud for geometric testing.
    
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

def test_find_distances() -> None:
    """
    test_find_distances
    Validates that the isotropic search radius calculation behaves consistently and scales properly.
    
    The algorithm returns 1.5 * max(min_distances). For a perfectly uniform square grid of spacing
    0.1, the output should exactly equal 0.15.
    """
    # 1. Test initialization
    p    = generate_square_cloud(11)                                                                                    # Generate a perfect 0.1 spaced grid.
    
    # 2. Distance evaluation
    dist = find_distances(p)                                                                                            # Evaluate critical distance.
    
    # 3. Assertions
    assert np.isclose(dist, 0.15, atol=1e-3)                                                                            # Verify optimal distance constraint.

def test_compute_neighbors_adaptive() -> None:
    """
    test_compute_neighbors_adaptive
    Validates the Round-Robin selection and geometric stability thresholds of the adaptive solver.
    
    Interior and border nodes on a well-conditioned geometric configuration should dynamically stop
    searching at exactly 9 points regardless of the provided nvec maximum limit.
    """
    # 1. Test initialization
    p    = generate_square_cloud(11)                                                                                    # Generate standard geometry.
    nvec = 20                                                                                                           # Extreme maximum limit.
    
    # 2. Adaptive neighbor search
    vec  = compute_neighbors(p, nvec)                                                                                   # Request neighbor computation.
    
    # 3. Structural validation
    assert vec.shape == (len(p), nvec)                                                                                  # Ensure matrix footprint conforms to padded expectations.
    
    # 4. Behavioral geometric validation
    nvec_active  = np.sum(vec != -1, axis=1)                                                                            # Compute active node count per stencil.
    interior_idx = np.argmin((p[:, 0] - 0.5)**2 + (p[:, 1] - 0.5)**2)                                                   # Locate an absolute interior node.
    
    assert nvec_active[interior_idx] == 9                                                                               # Interior node must forcefully stabilize at 9 neighbors.
    
    boundary_idx = np.argmin((p[:, 0] - 1.0)**2 + (p[:, 1] - 0.5)**2)                                                   # Locate a flat boundary node.
    
    assert nvec_active[boundary_idx] == 9                                                                               # Normal boundary node must also stabilize at 9 neighbors.

def test_compute_upwind_neighbors() -> None:
    """
    test_compute_upwind_neighbors
    Validates the convective directionality logic to ensure advective biased stencils.
    
    Given an extreme velocity vector, the algorithm must prioritize neighbors located in the upwind
    trajectory to properly secure the numerical stability of hyperbolic solvers.
    """
    # 1. Test initialization
    p    = generate_square_cloud(11)                                                                                    # Generate standard geometry.
    nvec = 12                                                                                                           # Request convective stencil limit.
    a    = 1.0                                                                                                          # Upwind velocity X.
    b    = 0.0                                                                                                          # Upwind velocity Y.
    
    # 2. Convective neighbor search
    vec  = compute_upwind_neighbors(p, a, b, nvec)                                                                      # Execute biased neighbor search.
    
    # 3. Structural validation
    assert vec.shape == (len(p), nvec)                                                                                  # Verify array footprint constraint.
    
    # 4. Directional verification
    idx             = np.argmin((p[:, 0] - 0.5)**2 + (p[:, 1] - 0.5)**2)                                                # Locate target node.
    neighbors       = vec[idx]                                                                                          # Retrieve corresponding neighbors.
    valid_neighbors = neighbors[neighbors != -1]                                                                        # Filter empty padding slots.
    nx              = p[valid_neighbors, 0]                                                                             # Extract X coordinates.
    
    upwind_count    = np.sum(nx <= 0.5 + 1e-10)                                                                         # Count nodes strictly located upwind.
    
    assert upwind_count > len(valid_neighbors) / 2                                                                      # Statistically biased majority must exist upwind.
