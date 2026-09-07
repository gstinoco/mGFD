"""
mGFD.oop.solver — Solver Orchestrator for mGFD OOP Interface

Overview:
    This module provides the high-level Solver class that binds a Domain (Cloud + Boundary)
    and a PDE (Physics), dispatches execution to CPU or CUDA backends, and returns a SolverResult.

Public API:
    Solver      High-level PDE solver class with fluid .solve() interface.

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
import time                                                                                                                             # Performance timing utilities.
import numpy as np                                                                                                                      # Core numerical operations.
from typing import Tuple, Optional, Any, Union, List                                                                                    # Type hinting interfaces.

from mGFD.oop.domain import Domain                                                                                                      # Domain abstraction class.
from mGFD.oop.pde import PDE                                                                                                            # PDE physics classes.
from mGFD.solvers.results import SolverResult                                                                                           # Standardized SolverResult object.
from mGFD.temporal.cfl import estimate_cfl_dt                                                                                           # CFL estimation module.

class Solver:                                                                                                                           # Solver class definition.
    """
    mGFD High-Level PDE Solver.
    
    Orchestrates execution of physical PDE problems over meshless domains using CPU or CUDA backends.
    Binds the geometry (Cloud), boundary conditions (Dirichlet/Neumann), and physics (PDE) into a single 
    executable unit that routes to the appropriate transient or stationary solver backend.
    """
    def __init__(self, domain: Domain, pde: PDE, device: str = "cpu", cfl: float = 0.5,                                                 # Constructor part 1.
                 implicit: bool = True, lam: float = 0.5, nvec: int = 15, upwind: Optional[bool] = None,                                # Constructor part 2.
                 verbose: bool = True, symmetric: bool = True) -> None:                                                                 # Constructor.
        """
        __init__
        Initialize Solver with Domain, PDE, hardware device, CFL parameter, and stencil options.
        
        Input:
            domain      Domain                  Physical domain (Cloud + BoundaryCondition).
            pde         PDE                     PDE physics formulation definition.
            device      str                     Execution hardware backend ('cpu' or 'cuda'). Default is 'cpu'.
            cfl         float                   Courant-Friedrichs-Lewy safety factor. Default is 0.5.
            implicit    bool                    Use implicit time stepping scheme. Default is True.
            lam         float                   Time integration parameter (Theta/Beta). Default is 0.5 (Crank-Nicolson / avg acc).
            nvec        int                     Neighborhood stencil size. Default is 15.
            upwind      Optional[bool]          Use directional upwind-biased stencils. Default is None (auto-selected).
            verbose     bool                    Print solver execution progress to console. Default is True.
        """
        self.domain   = domain                                                                                                          # Store physical domain object.
        self.pde      = pde                                                                                                             # Store PDE physics object.
        self.device   = device.lower()                                                                                                  # Store hardware device string.
        self.cfl      = cfl                                                                                                             # Store CFL safety factor.
        self.implicit = implicit                                                                                                        # Store implicit integration flag.
        self.lam      = lam                                                                                                             # Store implicit time scheme parameter (theta/beta).
        self.nvec     = nvec                                                                                                            # Store neighbor stencil size.
        self.upwind    = upwind                                                                                                         # Store upwind stencil selection.
        self.verbose   = verbose                                                                                                        # Store verbose logging flag.
        self.symmetric = symmetric                                                                                                      # Store conservative symmetrization flag.

    def solve(self, t_span: Tuple[float, float] = (0.0, 1.0), dt: Optional[float] = None) -> SolverResult:                              # Main solve method.
        """
        solve
        Executes the numerical PDE solver over the specified physical time domain t_span.
        
        Input:
            t_span      Tuple[float, float]     Physical time range (t_start, t_end). Default is (0.0, 1.0).
            dt          Optional[float]         Explicit time step override (if None, estimated automatically via CFL).
            
        Output:
            result      SolverResult            Standardized mGFD solution object.
        """
        # 1. Initialization and validation
        start_time  = time.perf_counter()                                                                                               # Start execution timer.
        p           = self.domain.points                                                                                                # Extract points array from domain cloud.
        bc_callable = self.domain.boundary                                                                                              # Extract boundary condition evaluation wrapper.
        vec         = self.domain.cloud.neighbors                                                                                       # Extract cached neighbor list.
        op          = self.pde.operator                                                                                                 # Extract differential operator array.
        if op is None:                                                                                                                  # Validate operator is set for type checker.
            raise ValueError("PDE operator is undefined. Ensure the PDE was initialized with a valid operator.")                        # Raise ValueError exception.

        # 2. Stationary PDEs (Order 0: Poisson / Perturbation)
        if self.pde.order == 0:                                                                                                         # Check if stationary PDE.
            f_src       = self.pde.source if self.pde.source is not None else 0.0                                                       # Extract source term or default float.
            upwind_flag = False if self.upwind is None else self.upwind                                                                 # Upwind flag selection.

            if self.device == "cpu":                                                                                                    # Evaluate CPU backend.
                from mGFD.solvers._backends.cpu.stationary import solve_cpu as solve_cpu_stat                                           # Import CPU backend solver.
                u_ap, vec, converged = solve_cpu_stat(p, bc_callable, f_src, op,                                                        # Call CPU solver part 1.
                                                       upwind_flag, vec, self.nvec, verbose=self.verbose)                               # Execute CPU solver.
            else:                                                                                                                       # Evaluate CUDA backend.
                from mGFD.solvers._backends.cuda.stationary import solve_cuda as solve_cuda_stat                                        # Import CUDA backend solver.
                u_ap, vec, converged = solve_cuda_stat(p, bc_callable, f_src, op,                                                       # Call CUDA solver part 1.
                                                        upwind_flag, vec, self.nvec, verbose=self.verbose)                              # Execute CUDA solver.

            comp_time = time.perf_counter() - start_time                                                                                # Compute total duration.
            return SolverResult(solution=u_ap, neighbors=vec, converged=converged, compute_time=comp_time, p=p)                         # Return structured result.

        # 3. First-Order Transient PDEs (Order 1: Heat / Advection-Diffusion)
        elif self.pde.order == 1:                                                                                                       # Check if 1st-order transient PDE.
            upwind_flag               = False if self.upwind is None else self.upwind                                                   # Auto-select isotropic stencils for 1st order.
            t_len                     = t_span[1] - t_span[0]                                                                           # Time domain length.
            dt_est, t_est, actual_cfl = estimate_cfl_dt(p, op, cfl=self.cfl, order=1, vec=vec, t_end=t_len)                             # Estimate CFL metrics.

            self.domain.boundary.t_span = t_span                                                                                        # Store physical time domain in boundary.
            if dt is not None:                                                                                                          # If explicit dt provided.
                dt_use = dt                                                                                                             # Use provided dt.
                t_use  = max(2, round(t_len / dt_use) + 1)                                                                              # Compute discrete time step points.
            else:                                                                                                                       # Automatic CFL estimation.
                dt_use = dt_est                                                                                                         # Use estimated dt.
                t_use  = t_est                                                                                                          # Use estimated steps.

            lam                   = self.lam                                                                                            # Theta parameter for implicit scheme.
            coef_val: List[float] = [float(x) for x in self.pde.coef] if isinstance(self.pde.coef, (list, tuple)) else []               # Type-safe coefficient list.

            if self.device == "cpu":                                                                                                    # Evaluate CPU backend.
                from mGFD.solvers._backends.cpu.time_derivative1 import solve_cpu as solve_cpu_td1                                      # Import CPU backend solver.
                u_ap, vec, converged = solve_cpu_td1(p, None, t_use, coef_val, op, upwind_flag, vec, self.nvec,                         # Call CPU solver.
                                                     self.implicit, lam, verbose=self.verbose, ic=self.pde.ic, bc=bc_callable,          # Pass parameters.
                                                     source=self.pde.source, t_span=t_span)                                             # Pass t_span.
            else:                                                                                                                       # Evaluate CUDA backend.
                from mGFD.solvers._backends.cuda.time_derivative1 import solve_cuda as solve_cuda_td1                                   # Import CUDA backend solver.
                u_ap, vec, converged = solve_cuda_td1(p, None, t_use, coef_val, op, upwind_flag, vec, self.nvec,                        # Call CUDA solver.
                                                      self.implicit, lam, verbose=self.verbose, ic=self.pde.ic, bc=bc_callable,         # Pass parameters.
                                                      source=self.pde.source, t_span=t_span)                                            # Pass t_span.

            comp_time = time.perf_counter() - start_time                                                                                # Compute total duration.
            return SolverResult(solution=u_ap, neighbors=vec, converged=converged, compute_time=comp_time,                              # Construct result.
                                dt=dt_use, cfl=actual_cfl, t_steps=t_use, p=p)                                                          # Return structured result.

        # 4. Second-Order Transient PDEs (Order 2: Wave Equation)
        elif self.pde.order == 2:                                                                                                       # Check if 2nd-order transient PDE.
            upwind_flag = False if self.upwind is None else self.upwind                                                                 # Default to isotropic balanced stencils for wave.
            
            damping_val = getattr(self.pde, 'damping', None)                                                                            # Try to get from subclass.
            if damping_val is None and isinstance(self.pde.coef, dict): damping_val = self.pde.coef.get('damping', 0.0)                 # Fallback to coef dict.
            elif damping_val is None: damping_val = 0.0                                                                                 # Default.
            
            alpha_val = getattr(self.pde, 'alpha', None)                                                                                # Try to get from subclass.
            if alpha_val is None and isinstance(self.pde.coef, dict): alpha_val = self.pde.coef.get('alpha', -0.1)                      # Fallback to coef dict.
            elif alpha_val is None: alpha_val = -0.1                                                                                    # Default.
            
            t_len                     = t_span[1] - t_span[0]                                                                           # Time domain length.
            dt_est, t_est, actual_cfl = estimate_cfl_dt(p, op, cfl=self.cfl, order=2, vec=vec, t_end=t_len)                             # Estimate CFL metrics for 2nd order.

            self.domain.boundary.t_span = t_span                                                                                        # Store physical time domain in boundary.
            if dt is not None:                                                                                                          # If explicit dt provided.
                dt_use = dt                                                                                                             # Convert dt to float.
                t_use  = max(2, round(t_len / dt_use) + 1)                                                                              # Compute discrete time step points.
            else:                                                                                                                       # Automatic CFL estimation.
                dt_use = dt_est                                                                                                         # Use estimated dt.
                t_use  = t_est                                                                                                          # Use estimated steps.

            g_val = getattr(self.pde, 'g', 0.0)                                                                                         # Initial velocity u_t(0).
            if g_val is None: g_val = 0.0                                                                                               # Default to 0.0.
            lam                   = self.lam                                                                                            # Newmark beta parameter for implicit scheme.
            coef_val              = [float(x) for x in self.pde.coef] if isinstance(self.pde.coef, (list, tuple)) else []               # Type-safe coefficient list.
            symmetric_val         = getattr(self.pde, 'symmetric', getattr(self, 'symmetric', True))                                    # Determine conservative symmetrization.

            if self.device == "cpu":                                                                                                    # Evaluate CPU backend.
                from mGFD.solvers._backends.cpu.time_derivative2 import solve_cpu as solve_cpu_td2                                      # Import CPU backend solver.
                u_ap, vec, converged = solve_cpu_td2(p, None, g_val, t_use, coef_val, op, upwind_flag, vec, self.nvec,                  # Call CPU solver with correct positional args.
                                                     self.implicit, lam, verbose=self.verbose, ic=self.pde.ic, bc=bc_callable,          # Pass implicit, lam, verbose, ic, bc.
                                                     source=self.pde.source, damping=damping_val, alpha=alpha_val, t_span=t_span,       # Pass source, damping, alpha, t_span.
                                                     symmetric=symmetric_val)                                                           # Pass conservative symmetrization flag.
            else:                                                                                                                       # Evaluate CUDA backend.
                from mGFD.solvers._backends.cuda.time_derivative2 import solve_cuda as solve_cuda_td2                                   # Import CUDA backend solver.
                u_ap, vec, converged = solve_cuda_td2(p, None, g_val, t_use, coef_val, op, upwind_flag, vec, self.nvec,                 # Call CUDA solver with correct positional args.
                                                      self.implicit, lam, verbose=self.verbose, ic=self.pde.ic, bc=bc_callable,         # Pass implicit, lam, verbose, ic, bc.
                                                      source=self.pde.source, damping=damping_val, alpha=alpha_val, t_span=t_span,      # Pass source, damping, alpha, t_span.
                                                      symmetric=symmetric_val)                                                          # Pass conservative symmetrization flag.

            comp_time = time.perf_counter() - start_time                                                                                # Compute total duration.
            return SolverResult(solution=u_ap, neighbors=vec, converged=converged, compute_time=comp_time,                              # Construct result.
                                dt=dt_use, cfl=actual_cfl, t_steps=t_use, p=p)                                                          # Return structured result.

        else:                                                                                                                           # Unknown PDE order fallback.
            raise ValueError(f"Unsupported PDE order: {self.pde.order}")                                                                # Raise ValueError exception.

    def __repr__(self) -> str:                                                                                                          # String representation.
        return f"Solver(domain={self.domain}, pde={self.pde}, device='{self.device}', cfl={self.cfl})"                                  # Return summary string.
