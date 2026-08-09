import numpy as np
import pytest
from mGFD import Stationary

def generate_square_cloud(n=11):
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    flag = np.zeros_like(X, dtype=int)
    flag[(X == 0) | (X == 1) | (Y == 0) | (Y == 1)] = 1
    p = np.column_stack((X, Y, flag))
    return p

def test_stationary_input_validation():
    with pytest.raises(ValueError, match="Point cloud 'p' must be a 2D numpy array"):
        Stationary(p=[1, 2, 3], phi=lambda x,y: 0, f=lambda x,y: 0)

    p = generate_square_cloud()
    with pytest.raises(TypeError, match="Boundary condition 'phi' must be a callable function"):
        Stationary(p=p, phi=0, f=lambda x,y: 0)

    with pytest.raises(TypeError, match="Right-hand side 'f' must be a callable function"):
        Stationary(p=p, phi=lambda x,y: 0, f=0)

    with pytest.raises(ValueError, match="Operator must be a numpy array with at least 5 coefficients"):
        Stationary(p=p, phi=lambda x,y: 0, f=lambda x,y: 0, operator=np.array([1, 2]))

def test_stationary_poisson():
    p = generate_square_cloud(21)
    
    # Exact solution and derivatives for Laplace/Poisson
    phi = lambda x, y: 2 * np.exp(2 * x + y)
    f   = lambda x, y: 10 * np.exp(2 * x + y)
    L   = np.vstack([[0], [0], [2], [0], [2], [0]])
    
    u_ap, vec = Stationary(p, phi, f, operator=L, verbose=False)
    
    u_ex = phi(p[:, 0], p[:, 1])
    error = np.sqrt(np.mean((u_ap - u_ex)**2))
    
    # Error should be reasonably small
    assert error < 5e-2
