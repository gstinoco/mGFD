import numpy as np
import pytest
from mGFD import TimeDerivative2

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

def test_time_derivative2_input_validation():
    p = generate_square_cloud()
    
    with pytest.raises(ValueError, match="Point cloud 'p' must be a 2D numpy array"):
        TimeDerivative2(p=[1, 2], f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=10, coef=[])

    with pytest.raises(TypeError, match="Initial condition function 'f' must be a callable"):
        TimeDerivative2(p=p, f=0, g=lambda x,y,t,c: 0, t=10, coef=[])

    with pytest.raises(TypeError, match="Initial condition function 'g' must be a callable"):
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g=0, t=10, coef=[])

    with pytest.raises(ValueError, match="Number of time steps 't' must be a positive integer"):
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=0, coef=[])

def test_time_derivative2_wave():
    p = generate_square_cloud(11)
    
    c = float(np.sqrt(1 / 2))
    t = 100
    f = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))
    g = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])
    
    u_ap, vec = TimeDerivative2(p, f, g, t, [c], operator=L, implicit=True, verbose=False)
    
    T = np.linspace(0, 1, t)
    u_ex = np.zeros([len(p), t])
    for k in range(t):
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [c])
        
    error = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))
    assert error < 2e-1
