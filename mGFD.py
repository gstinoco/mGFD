"""
Meshless Generalized Finite Differences

All the codes presented below were developed by:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

With the funding of:
    National Council of Science and Technology, CONACyT (Consejo Nacional de Ciencia y Tecnología, CONACyT). México.
    Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
    Aula CIMNE-Morelia. México

Date:
    May, 2024.

Last Modification:
    May 2024
"""

import numpy as np
import Scripts.Gammas as Gammas
import Scripts.Neighbors as Neighbors

def Stationary(p, phi, f, operator = np.vstack([[0], [0], [2], [0], [2]]), triangulation = False, tt = None, Adv = False):
    '''
    Numerical solution of partial differential equations with no time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu = -f(x,y)
    
    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        phi                         function        Function declared with the boundary condition.
        f                           function        Function declared with the right-hand side of the equation.
        operator                    ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        triangulation               bool            Select whether or not there is a triangulation available.
                                                        True: Triangulation available.
                                                        False: No triangulation available (Default).
        tt          m x 3           ndarray         Array with the triangulation indexes.
    
    Output:
        u_ap        m               ndarray         Array with the approximation computed by the routine.
        u_ex        m               ndarray         Array with the theoretical solution.
        vec         m x o           ndarray         Array with the correspondence of the o neighbors of each node.  
    '''

    # Variable initialization
    m      = len(p[:, 0])                                                           # The total number of nodes is calculated.
    nvec   = 8                                                                      # Maximum number of neighbors for each node.
    u_ap   = np.zeros([m])                                                          # u_ap initialization with zeros.
    u_ex   = np.zeros([m])                                                          # u_ex initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                        # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                           # Save the inner nodes.

    # Values for the velocities form the operator
    if Adv == True:                                                                 # If an Upwind stencil is requested.
        a = operator[0][0]                                                          # Value of the velocity on x.
        b = operator[1][0]                                                          # Value of the velocity on y.

    # Boundary conditions
    u_ap[boun_n] = phi(p[boun_n, 0], p[boun_n, 1])                                  # The boundary condition is assigned.
    
    ## Neighbor search for all the nodes.
    if triangulation == True:                                                       # If there are triangles available.
        vec = Neighbors.Triangulation(p, tt, nvec)                                  # Neighbor search with the proper routine.
    elif Adv == True:                                                               # If there are no triangles available and upwind required.
        vec = Neighbors.CloudAdv(p, a, b, nvec)                                     # Neighbor search with the proper routine.
    else:                                                                           # All the other cases.
        vec = Neighbors.Cloud(p, nvec)                                              # Neighbor search with the proper routine.

    # Computation of Gamma values
    L = operator[:-1]                                                               # The values of the differential operator are assigned.
    K = Gammas.Cloud(p, vec, L)                                                     # K computation with the required Gammas.
    R = Gammas.RHS(p, boun_n, inne_n, phi, f)                                       # Right-hand side of the equation.
    
    # A Generalized Finite Differences Method
    K = np.linalg.pinv(K)
    un  = K@R
    u_ap[inne_n] = un[inne_n]
    
    # Theoretical Solution
    u_ex = phi(p[:,0], p[:,1])                                                      # The theoretical solution is computed.

    return u_ap, u_ex, vec


def TimeDerivative1(p, f, t, coef, operator = np.vstack([[0], [0], [2], [0], [2]]), triangulation = False, tt = [], implicit = False, lam = 0.5, Adv = False):
    """
    Numerical solution of partial differential equations with first-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \frac{\partial u}{\partial t} = Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the flag for boundary or inner node.
        f                           function        Function declared with the boundary condition.
        t                           int             Number of time steps to be considered.
        coef                        ndarray         Coefficients of the specific problem.
        operator                    ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        triangulation               bool            Select whether or not there is a triangulation available.
                                                        True: Triangulation available.
                                                        False: No triangulation available (Default)
        tt          m x 3           ndarray         Array with the triangulation indexes.
        implicit                    bool            Select whether or not use an implicit scheme.
                                                        True: Implicit scheme used.
                                                        False: Explicit scheme used (Default).
        lam                         float           Lambda parameter for the implicit scheme.
                                                        Must be between 0 and 1 (Default: 0.5).
    
    Output:
        u_ap        m x t           ndarray         Array with the approximation computed by the routine.
        u_ex        m x t           ndarray         Array with the theoretical solution.
        vec         m x o           ndarray         Array with the correspondence of the o neighbors of each node.
    """

    # Variable initialization
    m    = len(p[:,0])                                                              # The total number of nodes is calculated.
    nvec = 8                                                                        # Maximum number of neighbors for each node.
    T    = np.linspace(0,1,t)                                                       # Time discretization.
    dt   = T[1] - T[0]                                                              # dt computation.
    u_ap = np.zeros([m,t])                                                          # u_ap initialization with zeros.
    u_ex = np.zeros([m,t])                                                          # u_ex initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                        # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                           # Save the inner nodes.

    # Values for the velocities form the operator
    if Adv == True:                                                                 # If an Upwind stencil is requested.
        a = L[0]                                                                    # Value of the velocity on x.
        b = L[1]                                                                    # Value of the velocity on y.
    
    # Boundary conditions.
    for k in np.arange(t):                                                          # For each time step.
        u_ap[boun_n, k] = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                 # The boundary condition is assigned.
  
    # Initial condition
    u_ap[:, 0] = f(p[:, 0], p[:, 1], T[0], coef)                                    # The initial condition is assigned.
    
    ## Neighbor search for all the nodes.
    if triangulation == True:                                                       # If there are triangles available.
        vec = Neighbors.Triangulation(p, tt, nvec)                                  # Neighbor search with the proper routine.
    elif Adv == True:                                                               # If there are no triangles available and upwind required.
        vec = Neighbors.CloudAdv(p, nvec)                                           # Neighbor search with the proper routine.
    else:                                                                           # All the other cases.
        vec = Neighbors.Cloud(p, nvec)                                              # Neighbor search with the proper routine.

    # Gamma computation.
    L = dt*operator[:-1]                                                            # The values of the differential operator are assigned.
    K = Gammas.Cloud(p, vec, L)                                                     # K computation with the required Gammas.
    
    # Generalized Finite Differences Method
    if implicit == False:                                                           # For the explicit scheme.
        K2 = np.identity(m) + K                                                     # Explicit formulation of K.
    else:                                                                           # For the implicit scheme.
        K2 = np.linalg.pinv(np.identity(m) - (1-lam)*K)@(np.identity(m) + lam*K)    # Implicit formulation of K.

    for k in np.arange(1,t):                                                        # For each of the time steps.
        un = K2@u_ap[:,k-1]                                                         # The new time-level is computed.
        u_ap[inne_n,k] = un[inne_n]                                                 # Save the computed solution.
        
    # Theoretical Solution
    for k in np.arange(t):                                                          # For all the time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], coef)                                # The theoretical solution is computed.

    return u_ap, u_ex, vec


def TimeDerivative2(p, f, g, t, coef, operator = np.vstack([[0], [0], [2], [0], [2]]), triangulation = False, tt = None, implicit = False, lam = 0.5, Adv = False):
    '''
    Numerical solution of partial differential equations with second-order time derivatives using a Meshless Generalized Finite Difference Scheme.
    
    The problem to solve is:
    
    \frac{\partial^2 u}{\partial t^2} = Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu
    
    Input:
        p           m x 2           ndarray         Array with the coordinates of the nodes.
        f                           function        Function declared with the boundary condition.
        g                           function        Function declared with the boundary condition.
        t                           int             Number of time steps to be considered.
        coef                        ndarray         Coefficients of the specific problem.
        operator                    ndarray         Array with the weights for the operator.
                                                        ([D, E, A, B, C, F]).
                                                        ([0, 0, 2, 0, 2, 0] is the default).
        triangulation               bool            Select whether or not there is a triangulation available.
                                                        True: Triangulation available.
                                                        False: No triangulation available (Default).
        tt          m x 3           ndarray         Array with the triangulation indexes.
        implicit                    bool            Select whether or not use an implicit scheme.
                                                        True: Implicit scheme used.
                                                        False: Explicit scheme used (Default).
        lam                         float           Lambda parameter for the implicit scheme.
                                                        Must be between 0 and 1 (Default: 0.5).
    
    Output:
        u_ap        m x t           ndarray         Array with the approximation computed by the routine.
        u_ex        m x t           ndarray         Array with the theoretical solution.
        vec         m x o           ndarray         Array with the correspondence of the o neighbors of each node.
    '''
    
    ## Variable initialization.
    m      = len(p[:, 0])                                                           # The total number of nodes is calculated.
    nvec   = 8                                                                      # Maximum number of neighbors for each node.
    T      = np.linspace(0, 1, t)                                                   # Time discretization.
    dt     = T[1] - T[0]                                                            # dt computation.
    u_ap   = np.zeros([m, t])                                                       # u_ap initialization with zeros.
    u_ex   = np.zeros([m, t])                                                       # u_ex initialization with zeros.
    boun_n = (p[:, 2] == 1) | (p[:, 2] == 2)                                        # Save the boundary nodes.
    inne_n = p[:, 2] == 0                                                           # Save the inner nodes.

    # Values for the velocities form the operator
    if Adv == True:                                                                 # If an Upwind stencil is requested.
        a = L[0]                                                                    # Value of the velocity on x.
        b = L[1]                                                                    # Value of the velocity on y.

    ## Boundary conditions.
    for k in np.arange(t):                                                          # For each time step.
        u_ap[boun_n, k] = f(p[boun_n, 0], p[boun_n, 1], T[k], coef)                 # The boundary condition is assigned.

    ## Initial condition.
    u_ap[:, 0] = f(p[:, 0], p[:, 1], T[0], coef)                                    # The initial condition is assigned.
    
    ## Neighbor search for all the nodes.
    if triangulation == True:                                                       # If there are triangles available.
        vec = Neighbors.Triangulation(p, tt, nvec)                                  # Neighbor search with the proper routine.
    elif Adv == True:                                                               # If there are no triangles available and upwind required.
        vec = Neighbors.CloudAdv(p, nvec)                                           # Neighbor search with the proper routine.
    else:                                                                           # All the other cases.
        vec = Neighbors.Cloud(p, nvec)                                              # Neighbor search with the proper routine.

    ## Gamma computation.
    L = (dt**2)*operator[:-1]                                                       # The values of the differential operator are assigned.
    K = Gammas.Cloud(p, vec, L)                                                     # K computation with the required Gammas.

    if implicit == False:                                                           # For the explicit scheme.
        K1 = np.identity(m)                                                         # Implicit formulation of K for k = 1.
        K2 = np.identity(m) + (1/2)*K                                               # Implicit formulation of K for k = 1.
        K3 = np.identity(m)                                                         # Implicit formulation of K for k = 2, ..., t.
        K4 = 2*np.identity(m) + K                                                   # Implicit formulation of K for k = 2, ..., t.
    else:                                                                           # For the implicit scheme.
        K1 = np.linalg.pinv(np.identity(m) - (1 - lam)*(1/2)*K)                     # Implicit formulation of K for k = 1.
        K2 = np.identity(m) + lam*(1/2)*K                                           # Implicit formulation of K for k = 1.
        K3 = np.linalg.pinv(np.identity(m) - (1 - lam)*K)                           # Implicit formulation of K for k = 2, ..., t.
        K4 = 2*np.identity(m) + lam*K                                               # Implicit formulation of K for k = 2, ..., t.

    ## Generalized Finite Differences Method
    for k in np.arange(1, t):                                                       # For al time levels.
        if k == 1:                                                                  # For the first time level.
            un = K1@(K2@u_ap[:, k - 1] + dt*g(p[:, 0], p[:, 1], T[k], coef))        # The new time-level is computed.
            u_ap[inne_n, k] = un[inne_n]                                            # Save the computed solution.
        else:                                                                       # For all the other time levels.
            un = K3@(K4@u_ap[:, k - 1] - u_ap[:, k - 2])                            # The new time-level is computed.
            u_ap[inne_n, k] = un[inne_n]                                            # Save the computed solution.                

    ## Theoretical Solution
    for k in np.arange(t):                                                          # For all the time steps.
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], coef)                                # The theoretical solution is computed.

    return u_ap, u_ex, vec