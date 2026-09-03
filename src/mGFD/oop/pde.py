"""
mGFD.oop.pde — Partial Differential Equation Formulations for mGFD OOP Interface

Overview:
    This module provides object-oriented physics abstractions for PDEs in mGFD,
    including PoissonEquation, HeatEquation, AdvectionDiffusion, and WaveEquation.

Public API:
    PDE                 Base class for PDE physics formulations.
    PoissonEquation     Stationary Poisson Equation (Order 0).
    HeatEquation        First-Order Transient Heat Equation (Order 1).
    AdvectionDiffusion  First-Order Transient Advection-Diffusion Equation (Order 1).
    WaveEquation        Second-Order Transient Wave Equation (Order 2).

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
import numpy as np                                                                                                      # Core numerical operations.
from typing import Union, Callable, Any, Optional                                                                       # Type hinting interfaces.
from mGFD.exceptions import OperatorFormatError                                                                         # Custom operator exception.

class PDE:                                                                                                              # Base PDE class.
    """
    Base class for Partial Differential Equation formulations in mGFD.
    """
    def __init__(self, operator: Optional[np.ndarray] = None, source: Optional[Union[Callable, float, Any]] = None,
                 ic: Optional[Union[Callable, float, Any]] = None, g: Optional[Union[Callable, float, Any]] = None,
                 coef: Optional[Any] = None) -> None:                                                                   # Initialize PDE base.
        """
        __init__
        Initialize base PDE with differential operator, source, initial conditions, and coefficients.
        
        Input:
            operator    ndarray                 6-element operator vector [D, E, A, B, C, F].
            source      Callable | float        Forcing/source term function F(x, y, t, coef).
            ic          Callable | float        Initial position/state u(x, y, 0).
            g           Callable | float        Initial velocity u_t(x, y, 0) (for 2nd order PDEs).
            coef        Any                     Parameter coefficient tuple or list.
        """
        if operator is not None:                                                                                        # If operator supplied.
            op_arr = np.asarray(operator)                                                                               # Convert to numpy array.
            if op_arr.size < 5:                                                                                         # Validate 5+ coefficients.
                raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients")                # Raise exception.

        self.operator = operator                                                                                        # Differential operator vector.
        self.source   = source                                                                                          # Source term function F.
        self.ic       = ic                                                                                              # Initial condition function u_0.
        self.g        = g                                                                                               # Initial velocity function u_t_0.
        self.coef     = coef                                                                                            # Additional physical coefficients.
        self.order    = 1                                                                                               # Time derivative order (0=stationary, 1=1st order, 2=2nd order).


class PoissonEquation(PDE):                                                                                             # Poisson PDE class.
    """
    Stationary Poisson Equation: L u = F(x, y)
    Defaults to Laplacian operator \\nabla^2 u = u_{xx} + u_{yy} ([0, 0, 2, 0, 2, 0]).
    """
    def __init__(self, source: Optional[Union[Callable, float, Any]] = None, coef: Optional[Any] = None, 
                 operator: Optional[np.ndarray] = None) -> None:                                                        # Constructor.
        if operator is None:                                                                                            # If custom operator not supplied.
            operator = np.vstack([[0], [0], [2], [0], [2], [0]])                                                        # Default Laplacian operator [0, 0, 2, 0, 2, 0].
        super().__init__(operator=operator, source=source, coef=coef)                                                   # Initialize base PDE.
        self.order = 0                                                                                                  # Stationary PDE (0th time order).

    def __repr__(self) -> str:                                                                                          # String representation.
        op_str = self.operator.ravel()[:5] if self.operator is not None else None                                       # Safely extract operator.
        return f"PoissonEquation(operator={op_str})"                                                                    # Return summary string.


class HeatEquation(PDE):                                                                                                # Transient Heat PDE class.
    """
    First-Order Transient Heat / Diffusion Equation:
        u_t = k * (u_{xx} + u_{yy}) + F(x, y, t)
    """
    def __init__(self, k: float = 1.0, source: Optional[Union[Callable, float, Any]] = None,
                 ic: Optional[Union[Callable, float, Any]] = None,
                 operator: Optional[np.ndarray] = None) -> None:                                                        # Constructor.
        if operator is None:                                                                                            # If operator not supplied.
            operator = np.vstack([[0], [0], [2 * k], [0], [2 * k], [0]])                                                # Thermal diffusion operator [0, 0, 2k, 0, 2k, 0].
        super().__init__(operator=operator, source=source, ic=ic, coef=[k])                                             # Initialize base PDE.
        self.k = k                                                                                                      # Store thermal diffusivity k.
        self.order = 1                                                                                                  # First-order in time.

    def __repr__(self) -> str:                                                                                          # String representation.
        return f"HeatEquation(k={self.k})"                                                                              # Return summary string.


class AdvectionDiffusion(PDE):                                                                                          # Advection-Diffusion PDE class.
    """
    First-Order Transient Advection-Diffusion Equation:
        u_t = v * (u_{xx} + u_{yy}) - (v_x * u_x + v_y * u_y) + F(x, y, t)
    """
    def __init__(self, v: float = 0.01, v_x: float = 0.0, v_y: float = 0.0,
                 source: Optional[Union[Callable, float, Any]] = None,
                 ic: Optional[Union[Callable, float, Any]] = None,
                 operator: Optional[np.ndarray] = None) -> None:                                                        # Constructor.
        if operator is None:                                                                                            # If operator not supplied.
            operator = np.vstack([[-v_x], [-v_y], [2 * v], [0], [2 * v], [0]])                                          # Advection-diffusion operator [-v_x, -v_y, 2v, 0, 2v, 0].
        super().__init__(operator=operator, source=source, ic=ic, coef=[v, v_x, v_y])                                   # Initialize base PDE.
        self.v   = v                                                                                                    # Store diffusion coefficient v.
        self.v_x = v_x                                                                                                  # Store x-velocity component.
        self.v_y = v_y                                                                                                  # Store y-velocity component.
        self.order = 1                                                                                                  # First-order in time.

    def __repr__(self) -> str:                                                                                          # String representation.
        return f"AdvectionDiffusion(v={self.v}, v_x={self.v_x}, v_y={self.v_y})"                                        # Return summary string.


class WaveEquation(PDE):                                                                                                # Second-Order Wave PDE class.
    r"""
    Second-Order Transient Wave Equation:
        u_{tt} + \eta * u_t = c^2 * (u_{xx} + u_{yy}) + F(x, y, t)
    Includes HHT-\alpha numerical dissipation and physical velocity damping \eta (damping).
    """
    def __init__(self, c: float = 1.0, damping: float = 0.0, alpha: float = -0.1,
                 source: Optional[Union[Callable, float, Any]] = None,
                 ic: Optional[Union[Callable, float, Any]] = None,
                 g: Optional[Union[Callable, float, Any]] = 0.0,
                 operator: Optional[np.ndarray] = None) -> None:                                                        # Constructor.
        g_val = g if g is not None else 0.0                                                                             # Default to 0.0 if None.
        if operator is None:                                                                                            # If operator not supplied.
            operator = np.vstack([[0], [0], [2 * (c**2)], [0], [2 * (c**2)], [0]])                                      # Wave propagation operator [0, 0, 2c^2, 0, 2c^2, 0].
        super().__init__(operator=operator, source=source, ic=ic, g=g_val, coef=[c])                                    # Initialize base PDE.
        self.c       = c                                                                                                # Wave propagation speed.
        self.damping = damping                                                                                          # Velocity damping factor \eta.
        self.alpha   = alpha                                                                                            # HHT-\alpha dissipation parameter.
        self.order   = 2                                                                                                # Second-order in time.

    def __repr__(self) -> str:                                                                                          # String representation.
        return f"WaveEquation(c={self.c}, damping={self.damping}, alpha={self.alpha})"                                  # Return summary string.
