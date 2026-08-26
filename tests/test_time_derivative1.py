import numpy as np
import pytest
from mGFD import TimeDerivative1

from mGFD.exceptions import CloudShapeError, InputTypeError, ParameterError, OperatorFormatError
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
    
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):
        TimeDerivative1(p=[1, 2, 3], f=lambda x,y,t,c: 0, t=10, coef=[])  # type: ignore

    with pytest.raises(InputTypeError, match="Forcing function 'f' must be a callable, ndarray, or numeric constant."):
        TimeDerivative1(p=p, f="invalid_string", t=10, coef=[])  # type: ignore

    with pytest.raises(ParameterError, match="Number of time steps 't' must be a positive integer"):
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=-1, coef=[])

    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):
        TimeDerivative1(p=p, f=lambda x,y,t,c: 0, t=10, coef=[], operator=np.array([1]))  # type: ignore

def test_time_derivative1_heat():
    p = generate_square_cloud(11)
    
    v = 0.2
    t = 100
    f = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])
    
    T = np.linspace(0, 1, t)
    u_ex = np.zeros([len(p), t])
    for k in range(t):
        u_ex[:, k] = f(p[:, 0], p[:, 1], T[k], [v])
        
    # Test spsolve
    u_ap, vec = TimeDerivative1(p, f, t, [v], operator=L, implicit=True, linear_solver="spsolve", verbose=False)
    error_spsolve = np.sqrt(np.mean((u_ap[:, -1] - u_ex[:, -1])**2))
    assert error_spsolve < 1e-1

    # Test bicgstab
    u_ap_bicgstab, _ = TimeDerivative1(p, f, t, [v], operator=L, implicit=True, linear_solver="bicgstab", verbose=False)
    error_bicgstab = np.sqrt(np.mean((u_ap_bicgstab[:, -1] - u_ex[:, -1])**2))
    assert error_bicgstab < 1e-1

    # Test gmres
    u_ap_gmres, _ = TimeDerivative1(p, f, t, [v], operator=L, implicit=True, linear_solver="gmres", verbose=False)
    error_gmres = np.sqrt(np.mean((u_ap_gmres[:, -1] - u_ex[:, -1])**2))
    assert error_gmres < 1e-1
