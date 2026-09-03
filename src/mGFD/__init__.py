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
    PDE                     Base PDE class.
    PoissonEquation         Stationary Poisson PDE.
    HeatEquation            1st-order transient Heat PDE.
    AdvectionDiffusion      1st-order transient Advection-Diffusion PDE.
    WaveEquation            2nd-order transient Wave PDE.
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
from mGFD.solvers.results import SolverResult

from mGFD.spatial.neighbors import compute_neighbors, compute_mesh_spacing
from mGFD.spatial.gammas import compute_sparse_matrix

from mGFD.temporal.cfl import estimate_cfl_dt

from mGFD.cloud_generator import generate_cloud_natural, generate_cloud_regular, reduce_points_by_region

from mGFD.oop import (
    Cloud, Dirichlet, Neumann, Domain, PDE,
    PoissonEquation, HeatEquation, WaveEquation, AdvectionDiffusion,
    Solver
)

__all__ = [
    'Cloud', 'Dirichlet', 'Neumann', 'Domain', 'PDE',
    'PoissonEquation', 'HeatEquation', 'WaveEquation', 'AdvectionDiffusion',
    'Solver', 'SolverResult',
    'compute_neighbors', 'compute_mesh_spacing', 'compute_sparse_matrix',
    'estimate_cfl_dt',
    'generate_cloud_natural', 'generate_cloud_regular', 'reduce_points_by_region'
]
