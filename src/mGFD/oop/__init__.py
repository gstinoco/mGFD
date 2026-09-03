"""
mGFD.oop — Object-Oriented Architecture Package

Overview:
    This package provides clean, object-oriented abstractions for building and solving
    partial differential equations with mGFD.

Public API:
    Cloud                   Point cloud geometry representation.
    BoundaryCondition       Base class for boundary conditions.
    Dirichlet               Dirichlet boundary condition class (u = val).
    Neumann                 Neumann boundary condition class (du/dn = val).
    Domain                  Pairs Cloud and BoundaryCondition.
    PDE                     Base PDE class.
    PoissonEquation         Stationary Poisson PDE.
    HeatEquation            1st-order transient Heat PDE.
    AdvectionDiffusion      1st-order transient Advection-Diffusion PDE.
    WaveEquation            2nd-order transient Wave PDE.
    Solver                  High-level solver orchestrator.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

Date:
    September, 2026.
"""

## Library importation.
from mGFD.oop.boundary import BoundaryCondition, Dirichlet, Neumann
from mGFD.oop.cloud import Cloud
from mGFD.oop.domain import Domain
from mGFD.oop.pde import PDE, PoissonEquation, HeatEquation, AdvectionDiffusion, WaveEquation
from mGFD.oop.solver import Solver

__all__ = [
    'Cloud', 'BoundaryCondition', 'Dirichlet', 'Neumann', 'Domain',
    'PDE', 'PoissonEquation', 'HeatEquation', 'AdvectionDiffusion', 'WaveEquation',
    'Solver'
]
