import numpy as np
import pytest
from mGFD import TimeDerivative1

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

def test_time_derivative1_input_validation():
    p = generate_square_cloud()
    
    with pytest.raises(ValueError, match="Point cloud 'p' must be a 2D numpy array"):
        TimeDerivative1(p=[1, 2, 3], f=lambda x,y,t,c: 0, t=10, coef=[])

    with pytest.raises(TypeError, match="Forcing function 'f' must be a callable"):
        TimeDerivative1(p=p, f=0, t=10, coef=[])

    with pytest.raises(ValueError, match="Number of time steps 't' must be a positive integer"):
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=-1, coef=[])

    with pytest.raises(ValueError, match="Operator must be a numpy array with at least 5 coefficients"):
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=10, coef=[], operator=np.array([1]))

def test_time_derivative1_heat():
    p = generate_square_cloud(11)
    
    v = 0.2
    t = 100
    f = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])
    
    u_ap, vec = TimeDerivative1(p, f, t, [v], operator=L, implicit=True, verbose=False)
    
    T = np.linspace(0, 1, t)
    u_ex = np.zeros([len(p), t])
    for k in range(t):
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [v])
        
    error = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))
    assert error < 1e-1
