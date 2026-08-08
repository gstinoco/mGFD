"""
mGFD — Meshless Generalized Finite Differences

Overview:
    Meshless solvers for 2D stationary and transient PDEs using generalized finite differences (GFD).
    Spatial derivatives are approximated using local reconstructions over a node neighborhood (nvec).

Public API:
    Stationary          Stationary problems with Dirichlet boundary conditions.
    TimeDerivative1     First-order-in-time problems (heat / advection-diffusion family).
    TimeDerivative2     Second-order-in-time problems (wave family).

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.
            If vec is not provided, it is computed from p using Neighbors.compute_neighbors / Neighbors.compute_upwind_neighbors.

Operator conventions:
    operator                array-like
            A 6-coefficient vector [D, E, A, B, C, F] (shape (6,) or (6, 1)).
            The spatial stencil is built with L = operator[:5] = [D, E, A, B, C], interpreted as:
                D*u_x + E*u_y + A*u_xx + B*u_xy + C*u_yy
            The reaction term F*u is reserved in the last coefficient, but it is not applied by the
            current implementation. For the Laplacian, use [0, 0, 2, 0, 2, 0].
            When upwind=True, neighbor selection is upwind-biased using velocities (D, E).

Notes:
    NumPy is required. SciPy is optional, but implicit schemes require SciPy for sparse linear algebra.
    Transient solvers use a normalized time grid T = linspace(0, 1, t).
    When instability is detected, the solver may retry with expanded neighborhoods (8→12→16→20→30).

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
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
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
from mGFD.Solvers.Stationary import Stationary
from mGFD.Solvers.TimeDerivative1 import TimeDerivative1
from mGFD.Solvers.TimeDerivative2 import TimeDerivative2

__all__ = ['Stationary', 'TimeDerivative1', 'TimeDerivative2']
