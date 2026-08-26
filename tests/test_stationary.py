import numpy as np
import pytest
from mGFD import Stationary
from mGFD.exceptions import CloudShapeError, InputTypeError, OperatorFormatError

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
    with pytest.raises(CloudShapeError, match="Point cloud 'p' must be a 2D numpy array"):
        Stationary(p=[1, 2, 3], phi=lambda x,y: 0, f=lambda x,y: 0)  # type: ignore

    p = generate_square_cloud()
    with pytest.raises(InputTypeError, match="Boundary condition 'phi' must be a callable, ndarray, or numeric constant."):
        Stationary(p=p, phi="invalid_string", f=lambda x,y: 0)  # type: ignore

    with pytest.raises(InputTypeError, match="Right-hand side 'f' must be a callable, ndarray, or numeric constant."):
        Stationary(p=p, phi=lambda x,y: 0, f="invalid_string")  # type: ignore

    with pytest.raises(OperatorFormatError, match="Operator must be a numpy array with at least 5 coefficients"):
        Stationary(p=p, phi=lambda x,y: 0, f=lambda x,y: 0, operator=np.array([1, 2]))  # type: ignore

def test_stationary_poisson():
    p = generate_square_cloud(21)
    
    # Exact solution and derivatives for Laplace/Poisson
    phi = lambda x, y: 2 * np.exp(2 * x + y)
    f   = lambda x, y: 10 * np.exp(2 * x + y)
    L   = np.vstack([[0], [0], [2], [0], [2], [0]])
    
    # Test spsolve
    u_ap, vec = Stationary(p, phi, f, operator=L, linear_solver="spsolve", verbose=False)
    u_ex = phi(p[:, 0], p[:, 1])
    error_spsolve = np.sqrt(np.mean((u_ap - u_ex)**2))
    assert error_spsolve < 5e-2

    # Test bicgstab
    u_ap_bicgstab, _ = Stationary(p, phi, f, operator=L, linear_solver="bicgstab", verbose=False)
    error_bicgstab = np.sqrt(np.mean((u_ap_bicgstab - u_ex)**2))
    assert error_bicgstab < 5e-2

    # Test gmres
    u_ap_gmres, _ = Stationary(p, phi, f, operator=L, linear_solver="gmres", verbose=False)
    error_gmres = np.sqrt(np.mean((u_ap_gmres - u_ex)**2))
    assert error_gmres < 5e-2
