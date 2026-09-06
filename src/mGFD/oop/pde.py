"""
mGFD.oop.pde — Partial Differential Equation Formulations for mGFD OOP Interface

Overview:
    This module provides the pure generalized object-oriented physics abstraction (`PDE`) for mGFD,
    allowing the user to define arbitrary stationary or transient physics via explicit spatial operators.

Public API:
    PDE                 Unified generalized physics definition class.

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
    September, 2026.
Last Modification:
    September, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Core numerical operations.
from typing import Union, Callable, Any, Optional, List, Tuple                                                                          # Type hinting interfaces.
from mGFD.exceptions import OperatorFormatError                                                                                         # Custom operator exception.

class PDE:                                                                                                                              # Base PDE class.
    """
    Unified generalized physics definition class for Partial Differential Equation formulations in mGFD.
    
    Encapsulates the differential spatial operator (a 6-element array for 2D Taylor reconstruction), 
    initial conditions, source terms, additional physics coefficients, and time-derivative order.
    By generalizing the physics definition, it natively supports arbitrary non-standard PDE models 
    including advection-diffusion-reaction equations, pure stationary Poisson systems, or 
    second-order hyperbolic waves without requiring rigid subclassing.
    """
    def __init__(self, operator: Optional[Union[np.ndarray, List[float], Tuple[float, ...]]] = None,                                    # Parameter operator vector.
                 source: Optional[Union[Callable, float, Any]] = None,                                                                  # Source/forcing function.
                 ic: Optional[Union[Callable, float, Any]] = None,                                                                      # Initial condition function.
                 g: Optional[Union[Callable, float, Any]] = None,                                                                       # Initial velocity function.
                 coef: Optional[Any] = None, order: int = 1) -> None:                                                                   # Initialize PDE base.
        """
        __init__
        Initialize base PDE with differential operator, source, initial conditions, and coefficients.
        
        Input:
            operator    Optional[Union[np.ndarray, List[float], Tuple[float, ...]]]         6-element operator vector [D, E, A, B, C, F].
            source      Optional[Union[Callable, float, Any]]                               Forcing/source term function F(x, y, t, coef).
            ic          Optional[Union[Callable, float, Any]]                               Initial position/state u(x, y, 0).
            g           Optional[Union[Callable, float, Any]]                               Initial velocity u_t(x, y, 0) (for 2nd order PDEs).
            coef        Optional[Any]                                                       Parameter coefficient tuple or list.
            order       int                                                                 Time derivative order (0=stationary, 1=transient 1st, 2=transient 2nd).
        """
        # 1. Operator validation
        if operator is not None:                                                                                                        # If operator supplied.
            op_arr = np.asarray(operator)                                                                                               # Convert to numpy array.
            if op_arr.size < 5:                                                                                                         # Validate 5+ coefficients.
                raise OperatorFormatError("Operator must be a numpy array with at least 5 coefficients")                                # Raise exception.
            self.operator = op_arr                                                                                                      # Differential operator vector.
        else:                                                                                                                           # Fallback for undefined operator.
            self.operator = None                                                                                                        # Initialize operator attribute as None.

        # 2. Physics properties
        self.source                   = source                                                                                          # Source term function F.
        self.ic                       = ic                                                                                              # Initial condition function u_0.
        self.g                        = g                                                                                               # Initial velocity function u_t_0.
        self.coef                     = coef                                                                                            # Additional physical coefficients.
        self.order                    = order                                                                                           # Time derivative order (0=stationary, 1=1st order, 2=2nd order).
        self.damping: Optional[float] = None                                                                                            # Optional structural damping for Wave equation.
        self.alpha: Optional[float]   = None                                                                                            # Optional HHT-alpha coefficient for Wave equation.
