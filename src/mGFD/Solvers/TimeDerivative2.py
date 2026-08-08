"""
TimeDerivative2 — Second-order transient PDEs solver

Overview:
    Numerical solver for PDEs with a second-order time derivative using a Meshless Generalized Finite
    Difference (mGFD) scheme on a 2D cloud of points.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.

Public API:
    TimeDerivative2     Main solver function for second-order transient problems.

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
import numpy as np

from scipy.sparse import eye, diags
from scipy.sparse.linalg import factorized

from typing import Callable, Optional, Tuple, List
import mGFD.Gammas as Gammas
import mGFD.Neighbors as Neighbors

def TimeDerivative2(p: np.ndarray, f: Callable, g: Callable, t: int, coef: List[float], operator: np.ndarray = np.vstack([[0], [0], [2], [0], [2]]), implicit: bool = False, lam: float = 0.5, upwind: bool = False, vec: Optional[np.ndarray] = None, nvec: int = 12) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Numerical solution of partial differential equations with second-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \\frac{\\partial^2 u}{\\partial t^2} = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        f                           function        Function declared with the initial condition for u.
        g                           function        Function declared with the initial condition for du/dt.
        t                           int             Number of time steps to compute.
        coef                        list            Physical coefficients for the problem formulation.
        operator                    ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        implicit                    bool            If True, uses an implicit (theta) time integration scheme.
        lam                         float           Theta parameter for implicit scheme (0=explicit, 1=fully implicit, 0.5=Newmark).
        upwind                      bool            If an Upwind stencil is requested.
        vec         m x nvec        ndarray         Cached neighbor list (optional).
        nvec                        int             Maximum number of neighbors for each node.
    
    Output:
        u_ap        m x t           ndarray         Array with the approximation computed by the routine.
        u_ex        m x t           ndarray         Array with the theoretical solution.
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    
    # 1. Variable initialization
    m      = p.shape[0]                                                             # Total number of nodes.
    T      = np.linspace(0, 1, t)                                                   # Time discretization array.
    dt     = T[1] - T[0]                                                            # Time step size.
    u_ap   = np.zeros([m, t])                                                       # Numerical approximation matrix.
    u_ex   = np.zeros([m, t])                                                       # Exact solution matrix.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                        # Boolean mask for boundary nodes.
    inne_n = p[:, 2] == 0                                                           # Boolean mask for interior nodes.

    # 2. Extract advection velocities for Upwind scheme.
    if upwind:
        a = -operator[0][0] if operator.ndim == 2 else -operator[0]
        # X-velocity (D coefficient).
        b = -operator[1][0] if operator.ndim == 2 else -operator[1]
        # Y-velocity (E coefficient).

    # 3. Apply Boundary and Initial Conditions
    for k in np.arange(t):
        u_ap[boun_n, k] = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                 # Boundary condition (Dirichlet).

    u_ap[:, 0] = f(p[:, 0], p[:, 1], T[0], coef)                                    # Initial condition u(t=0).
    
    # 4. Neighbor search
    if vec is None:
        if upwind:
            vec = Neighbors.compute_upwind_neighbors(p, a, b, nvec)                 # Upwind-biased neighbor selection.
        else:
            vec = Neighbors.compute_neighbors(p, nvec)                              # Standard distance-based neighbors.

    # 5. Compute differentiation matrix (K)
    L = operator[:-1]                                                               # Original operator weights.
    K_spatial = Gammas.compute_sparse_matrix(p, vec, L)                             # Build sparse spatial differentiation matrix.
    K = (dt**2) * K_spatial                                                         # Scale by time step squared.
    
    # 6. Time Integration (Generalized Finite Differences)
    if not implicit:
        # Explicit scheme
        K2 = eye(m) + (1/2)*K                                                       # Explicit matrix for k = 1.
        K4 = 2*eye(m) + K                                                           # Explicit matrix for k = 2, ..., t.
        
        for k in np.arange(1, t):
            if k == 1:
                un = K2.dot(u_ap[:, k - 1]) + dt*g(p[:, 0], p[:, 1], T[k], coef)    # New time-level (k=1) using du/dt.
                u_ap[inne_n, k] = un[inne_n]
            else:
                un = K4.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]                        # New time-level (k>=2) using Leapfrog/Verlet.
                u_ap[inne_n, k] = un[inne_n]
    else:
        # Implicit scheme
        Id_inner = diags(inne_n.astype(float))                                      # Diagonal mask for inner nodes.
        Id_bound = diags(boun_n.astype(float))                                      # Diagonal mask for boundary nodes.
        
        # Step k = 1 Formulation (with du/dt condition)
        A1 = Id_inner @ (eye(m) - lam * (1/2) * K) + Id_bound
        B1 = Id_inner @ (eye(m) + (1 - lam) * (1/2) * K)
        
        # Step k >= 2 Formulation
        A2 = Id_inner @ (eye(m) - lam * K) + Id_bound                               # LHS Matrix: Theta parameter applied to inner, Identity to boundary.
        A2 = A2.tocsc()                                                             # Convert to CSC format for efficient SuperLU factorization.
        B2 = Id_inner @ (2*eye(m) + (1 - lam) * K)
        
        solve1 = factorized(A1.tocsc())                                             # Pre-factorize implicit matrix for k = 1.
        solve2 = factorized(A2)                                                     # Pre-factorize implicit matrix for k = 2, ..., t.

        for k in np.arange(1, t):
            if k == 1:
                RHS = B1.dot(u_ap[:, k - 1]) + dt*g(p[:, 0], p[:, 1], T[k], coef)
                RHS[boun_n] = u_ap[boun_n, k]                                       # Inject exact boundary condition for time k=1.
                un = solve1(RHS)
                u_ap[inne_n, k] = un[inne_n]
            else:
                RHS = B2.dot(u_ap[:, k - 1]) - u_ap[:, k - 2]
                RHS[boun_n] = u_ap[boun_n, k]                                       # Inject exact boundary condition for time k>=2.
                un = solve2(RHS)
                u_ap[inne_n, k] = un[inne_n]

    # 7. Compute exact theoretical solution
    for k in np.arange(t):
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], coef)

    return u_ap, u_ex, vec