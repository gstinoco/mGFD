"""
mGFD.solvers — High-level Solver APIs

Overview:
    This module provides the high-level solver APIs for the mGFD method,
    including stationary solvers, and first and second order time derivative solvers.

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
"""

## Library importation.
from mGFD.solvers.stationary import Stationary
from mGFD.solvers.time_derivative1 import TimeDerivative1
from mGFD.solvers.time_derivative2 import TimeDerivative2
from mGFD.solvers.results import SolverResult

__all__ = [
    'Stationary', 'TimeDerivative1', 'TimeDerivative2', 'SolverResult'
]
