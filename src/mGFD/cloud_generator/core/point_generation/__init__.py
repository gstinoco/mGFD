"""
CloudGenerator.core.point_generation — Point Generation Module

Overview:
    This package implements the core algorithms for generating point clouds using
    both Regular (Grid-based) and Natural (Poisson Disk Sampling) distributions.
    It handles the complex logic of boundary and interior point generation for
    multi-region geometries.

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
    March, 2026.
Last Modification:
    August, 2026.
"""

from .geometry import create_fast_polygon_checker                                                                                       # Geometry fast polygon checker.
from .boundary import generate_boundary_points                                                                                          # Boundary generation routine.
from .relaxation import lloyd_relaxation                                                                                                # Lloyd relaxation algorithm.
from .regular_generation import (generate_interior_points, generate_region_cloud_with_uniform_density,                                  # Regular interior generators.
                                 generate_region_cloud_with_holes, generate_region_task, generate_interior_regions_clouds)              # Regular domain generators.
from .poisson_generation import (poisson_disk_sampling, generate_interior_points_poisson, generate_region_cloud_poisson,                # Poisson disk generators.
                                 generate_region_cloud_with_holes_poisson, generate_region_task_poisson,                                # Poisson cloud generators.
                                 generate_interior_regions_clouds_poisson)                                                              # Poisson multi-region generators.

__all__ = [                                                                                                                             # Public export symbols list.
    'create_fast_polygon_checker', 'generate_boundary_points', 'lloyd_relaxation',                                                      # Core geometric and relaxation symbols.
    'generate_interior_points', 'generate_region_cloud_with_uniform_density',                                                           # Uniform grid generation symbols.
    'generate_region_cloud_with_holes', 'generate_region_task', 'generate_interior_regions_clouds',                                     # Complex domain grid symbols.
    'poisson_disk_sampling', 'generate_interior_points_poisson', 'generate_region_cloud_poisson',                                       # Poisson sampling symbols.
    'generate_region_cloud_with_holes_poisson', 'generate_region_task_poisson',                                                         # Poisson domain symbols.
    'generate_interior_regions_clouds_poisson'                                                                                          # Poisson interior region symbols.
]                                                                                                                                       # End of public symbols.
