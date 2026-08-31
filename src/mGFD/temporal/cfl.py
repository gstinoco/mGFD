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

Date:
    August, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.
from typing import Optional, Tuple                                                                                                      # Type hinting.

from mGFD.spatial.neighbors import compute_mesh_spacing                                                                                # Spatial mesh spacing estimation.

def estimate_cfl_dt(p: np.ndarray,
                    operator: np.ndarray,
                    cfl: float = 0.5,
                    order: int = 1,
                    vec: Optional[np.ndarray] = None,
                    t_end: float = 1.0) -> Tuple[float, int, float]:
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
    _, h_avg = compute_mesh_spacing(p, vec)                                                                                             # Compute characteristic average node spacing.
    
    # Extract operator components
    D = abs(float(operator[0][0] if operator.ndim == 2 else operator[0]))                                                               # Advection x (D).
    E = abs(float(operator[1][0] if operator.ndim == 2 else operator[1]))                                                               # Advection y (E).
    A = abs(float(operator[2][0] if operator.ndim == 2 else operator[2]))                                                               # Diffusion/Wave xx (A).
    C = abs(float(operator[4][0] if operator.ndim == 2 else operator[4]))                                                               # Diffusion/Wave yy (C).

    if order == 1:                                                                                                                      # First-order transient PDE (Heat/AdvDif).
        V_adv = float(np.hypot(D, E))                                                                                                   # Advective speed magnitude.
        nu    = max(A, C) / 2.0                                                                                                        # Diffusion coefficient.
        
        speed = V_adv + np.sqrt(nu)                                                                                                     # Characteristic propagation speed.
        if speed <= 0.0:
            speed = 1.0                                                                                                                 # Safeguard for zero operator.
            
        dt_max = cfl * (h_avg / speed)                                                                                                  # Compute safe characteristic time step limit.
    else:                                                                                                                               # Second-order transient PDE (Wave).
        c_wave = np.sqrt(max(A, C) / 2.0) if max(A, C) > 0 else 1.0                                                                     # Effective wave speed.
        dt_max = cfl * (h_avg / (c_wave * np.sqrt(2)))                                                                                  # Courant-Friedrichs-Lewy wave time step limit.

    t = max(1, int(np.ceil(t_end / dt_max)))                                                                                            # Number of discrete time steps.
    dt = t_end / t                                                                                                                      # Discrete time step size.
    actual_cfl = cfl * (dt / dt_max) if dt_max > 0 else cfl                                                                             # Actual effective CFL number.

    return dt, t, actual_cfl                                                                                                            # Return step parameters.
