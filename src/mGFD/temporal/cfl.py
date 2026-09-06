"""
CFL — Temporal CFL estimation and adaptive time-stepping for mGFD

Overview:
    Routines for calculating node spatial spacing, maximum stable time steps (dt),
    and Courant-Friedrichs-Lewy (CFL) numbers for transient PDEs in mGFD.

Public API:
    estimate_cfl_dt     Estimate maximum stable time step dt, step count t, and effective CFL number.

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
from typing import Optional, Tuple                                                                                                      # Type hinting.

from mGFD.spatial.neighbors import compute_mesh_spacing                                                                                 # Spatial mesh spacing estimation.

def estimate_cfl_dt(p: np.ndarray,                                                                                                      # Function signature part 1.
                    operator: np.ndarray,                                                                                               # Execute statement.
                    cfl: float = 0.5,                                                                                                   # Assign cfl: float.
                    order: int = 1,                                                                                                     # Assign order: int.
                    vec: Optional[np.ndarray] = None,                                                                                   # Assign vec: Optional[np.ndarray].
                    t_end: float = 1.0) -> Tuple[float, int, float]:                                                                    # Assign t_end: float.
    """
    estimate_cfl_dt
    Estimate maximum stable time step dt, required time steps t, and effective CFL number.

    Input:
        p           m x 3           ndarray         Point cloud array [x, y, flag].
        operator    6 x 1           ndarray         Differential operator weights [D, E, A, B, C, F].
        cfl                         float           Target CFL Courant safety factor (default 0.5).
        order                       int             PDE time derivative order (1 for parabolic/advective, 2 for wave).
        vec         m x nvec        ndarray         Cached neighbor matrix (optional).
        t_end                       float           Final time horizon (default 1.0).

    Output:
        dt                          float           Estimated stable time step size.
        t                           int             Recommended number of time steps.
        actual_cfl                  float           Resulting CFL number based on discrete step counts.
    """
    h_min, h_avg = compute_mesh_spacing(p, vec)                                                                                         # Compute minimum and average node spacing.
    h_char       = max(h_min, 0.25 * h_avg) if h_min > 0.0 else h_avg                                                                   # Characteristic mesh spacing (conservative lower bound).
    
    D = abs(float(operator[0][0] if operator.ndim == 2 else operator[0]))                                                               # Advection x (D).
    E = abs(float(operator[1][0] if operator.ndim == 2 else operator[1]))                                                               # Advection y (E).
    A = abs(float(operator[2][0] if operator.ndim == 2 else operator[2]))                                                               # Diffusion/Wave xx (A).
    C = abs(float(operator[4][0] if operator.ndim == 2 else operator[4]))                                                               # Diffusion/Wave yy (C).
    F = float(operator[5][0] if operator.ndim == 2 else operator[5]) if (operator.size >= 6) else 0.0                                   # Reaction coefficient (F).

    if order == 1:                                                                                                                      # First-order transient PDE (Heat/AdvDif).
        V_adv      = float(np.hypot(D, E))                                                                                              # Advective speed magnitude.
        nu         = max(A, C) / 2.0                                                                                                    # Diffusion coefficient.
        speed      = V_adv + np.sqrt(nu)                                                                                                # Characteristic propagation speed scale.
        if speed <= 0.0:                                                                                                                # Safeguard for zero operator.
            speed  = 1.0                                                                                                                # Default reference speed.
            
        dt_max     = cfl * (h_avg / speed)                                                                                              # Compute characteristic linear time step limit.
    else:                                                                                                                               # Second-order transient PDE (Wave).
        c_wave     = np.sqrt(max(A, C) / 2.0) if max(A, C) > 0.0 else 1.0                                                               # Effective wave speed.
        dt_max     = cfl * (h_char / (c_wave * np.sqrt(2)))                                                                             # Courant-Friedrichs-Lewy wave time step limit.

    t          = max(2, int(np.ceil(t_end / dt_max)))                                                                                   # Ensure at least 2 discrete time steps (initial + final).
    dt         = t_end / t                                                                                                              # Discrete time step size.
    actual_cfl = cfl * (dt / dt_max) if dt_max > 0.0 else cfl                                                                           # Actual effective CFL number.

    return dt, t, actual_cfl                                                                                                            # Return step parameters.
