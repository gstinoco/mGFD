import numpy as np
import pytest
from mGFD import TimeDerivative2
from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError
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
    
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):
        TimeDerivative2(p=[1, 2], f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore

    with pytest.raises(InputTypeError, match="Initial condition function 'f' must be a callable, ndarray, or numeric constant."):
        TimeDerivative2(p=p, f="invalid_string", g=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore

    with pytest.raises(InputTypeError, match="Initial velocity function 'g' must be a callable, ndarray, or numeric constant."):
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g="invalid_string", t=10, coef=[])  # type: ignore

    with pytest.raises(ParameterError, match="Number of time steps 't' must be a positive integer"):
        TimeDerivative2(p=p, f=lambda x,y,t,c: 0, g=lambda x,y,t,c: 0, t=0, coef=[])

def test_time_derivative2_wave():
    p = generate_square_cloud(11)
    
    c = float(np.sqrt(1 / 2))
    t = 100
    f = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))
    g = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])
    
    T = np.linspace(0, 1, t)
    u_ex = np.zeros([len(p), t])
    for k in range(t):
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [c])

    # Test spsolve
    u_ap, vec = TimeDerivative2(p, f, g, t, [c], operator=L, implicit=True, linear_solver="spsolve", verbose=False)
    error_spsolve = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))
    assert error_spsolve < 2e-1

    # Test bicgstab
    u_ap_bicgstab, _ = TimeDerivative2(p, f, g, t, [c], operator=L, implicit=True, linear_solver="bicgstab", verbose=False)
    error_bicgstab = np.sqrt(np.mean((u_ap_bicgstab[:, -1] - u_ex[:, -1])**2))
    assert error_bicgstab < 2e-1

    # Test gmres
    u_ap_gmres, _ = TimeDerivative2(p, f, g, t, [c], operator=L, implicit=True, linear_solver="gmres", verbose=False)
    error_gmres = np.sqrt(np.mean((u_ap_gmres[:, -1] - u_ex[:, -1])**2))
    assert error_gmres < 2e-1
