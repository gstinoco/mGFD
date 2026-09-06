"""
Test Cloud Generator — Unit tests for point generation routines

Overview:
    This file contains the unit tests for validating the point generation algorithms of mGFD,
    ensuring that the boundary and interior points are distributed according to the expected rules.

Public API:
    test_boundary_generation            Validates the boundary generation algorithm.

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
import numpy as np                                                                                                                      # Core numerical operations.
from mGFD.cloud_generator.core.point_generation.boundary import generate_boundary_points                                                # Point generation module.

def test_boundary_generation() -> None:
    """
    test_boundary_generation
    Validates that the boundary generation calculation creates the expected number of nodes.
    
    A square domain of 1x1 with a spacing of 0.5 must yield at least 8 boundary nodes.
    """
    # 1. Test initialization
    contour    = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]                                                                       # Square contour definition.
    cloud_size = 0.5                                                                                                                    # Target spacing.
    
    # 2. Execution
    p          = generate_boundary_points(contour, cloud_size)                                                                          # Generate point cloud.
    
    # 3. Assertions
    assert isinstance(p, np.ndarray)                                                                                                    # Result must be an array.
    assert p.shape[1] == 2                                                                                                              # Boundary points have 2 dimensions (x, y).
    assert p.shape[0] >= 4                                                                                                              # At least the 4 corners must be present.
