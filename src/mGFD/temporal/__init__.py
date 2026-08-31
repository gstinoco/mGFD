"""
mGFD.temporal — Temporal Discretization & Time-Stepping Utilities

Overview:
    This module handles temporal discretization, adaptive time-stepping, and Courant-Friedrichs-Lewy (CFL) 
    calculations for transient PDE solvers in mGFD.

Public API:
    estimate_cfl_dt     Estimate maximum stable time step dt, step count t, and effective CFL number.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

Date:
    August, 2026.
"""

## Library importation.
from mGFD.temporal.cfl import estimate_cfl_dt

__all__ = [
    'estimate_cfl_dt'
]
