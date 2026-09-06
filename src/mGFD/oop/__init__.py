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
    PDE                     Unified generalized physics definition class.
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
from mGFD.oop.boundary import BoundaryCondition, Dirichlet, Neumann                                                                     # Import boundary condition classes.
from mGFD.oop.cloud import Cloud                                                                                                        # Import point cloud abstraction.
from mGFD.oop.domain import Domain                                                                                                      # Import domain abstraction.
from mGFD.oop.pde import PDE                                                                                                            # Import PDE abstraction.
from mGFD.oop.solver import Solver                                                                                                      # Import high-level solver abstraction.

__all__ = [                                                                                                                             # List of public symbols exported.
    'Cloud', 'BoundaryCondition', 'Dirichlet', 'Neumann', 'Domain',                                                                     # Core domain classes.
    'PDE',                                                                                                                              # Generalized PDE class.
    'Solver'                                                                                                                            # Numerical solver class.
]                                                                                                                                       # End of exported symbols list.

