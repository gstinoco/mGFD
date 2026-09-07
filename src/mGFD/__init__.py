"""
mGFD — Meshless Generalized Finite Differences (OOP Architecture)

Overview:
    Commercial-grade meshless solver suite for 2D stationary and transient PDEs using generalized finite differences (GFD).
    Exposes an intuitive Object-Oriented Architecture (Cloud, Dirichlet, Neumann, Domain, PDE physics, Solver, SolverResult).

Public API:
    Cloud                   Point cloud geometry representation.
    Dirichlet               Dirichlet boundary condition class (u = val).
    Neumann                 Neumann boundary condition class (du/dn = val).
    Domain                  Pairs Cloud and BoundaryCondition.
    PDE                     Unified generalized physics definition class.
    Solver                  High-level solver orchestrator.
    SolverResult            Standardized solution result object.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI. México.
        Coordination of Scientific Research, CIC-UMSNH. México.
        Aula CIMNE-Morelia. México.
        SIIIA-MATH: Soluciones de Ingeniería. México.

Date:
    September, 2026.
"""

## Library importation.
from mGFD.solvers.results import SolverResult                                                                                           # Standardized solver result container.
from mGFD.temporal.cfl import estimate_cfl_dt                                                                                           # Temporal CFL condition calculation.
from mGFD.spatial.gammas import compute_sparse_matrix                                                                                   # GFD star coefficients matrix computation.
from mGFD.spatial.neighbors import compute_neighbors, compute_mesh_spacing                                                              # Spatial neighborhood utilities.
from mGFD.cloud_generator import generate_cloud_natural, generate_cloud_regular, reduce_points_by_region                                # Point cloud generator interfaces.
from mGFD.oop import (                                                                                                                  # High-level OOP abstractions.
    Cloud, Dirichlet, Neumann, Domain, PDE,                                                                                             # Core geometry and boundary classes.
    Solver                                                                                                                              # High-level solver orchestrator.
)                                                                                                                                       # End of OOP imports.

## Semantic version.
__version__ = "0.12.0"                                                                                                          # Library semantic version.

__all__ = [                                                                                                                             # List of public symbols exported.
    '__version__',                                                                                                                      # Library semantic version.
    'Cloud', 'Dirichlet', 'Neumann', 'Domain', 'PDE',                                                                                   # Core domain and physics classes.
    'Solver', 'SolverResult',                                                                                                           # Solvers and result containers.
    'compute_neighbors', 'compute_mesh_spacing', 'compute_sparse_matrix',                                                               # Spatial discretization functions.
    'estimate_cfl_dt',                                                                                                                  # Temporal discretization functions.
    'generate_cloud_natural', 'generate_cloud_regular', 'reduce_points_by_region'                                                       # Cloud generation functions.
]                                                                                                                                       # End of exported symbols.

