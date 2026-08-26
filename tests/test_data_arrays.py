import numpy as np
import pytest
from mGFD import Stationary, TimeDerivative1, TimeDerivative2

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

def test_stationary_arrays():
    p = generate_square_cloud(15)
    
    # 1. Callable functions
    phi_func = lambda x, y: x**2 + y**2
    f_func = lambda x, y: -4.0
    
    # 2. Array / Constant data
    phi_arr = p[:, 0]**2 + p[:, 1]**2
    f_const = -4.0
    
    op = np.vstack([[0], [0], [2], [0], [2], [0]]) # Laplacian
    
    # Run with Callable
    u_ap_func, _ = Stationary(p, phi_func, f_func, operator=op, nvec=12, verbose=False)
    
    # Run with Arrays and Constants
    u_ap_arr, _ = Stationary(p, phi_arr, f_const, operator=op, nvec=12, verbose=False)
    
    # Verify both methods yield the same results
    error = np.max(np.abs(u_ap_func - u_ap_arr))
    assert error < 1e-12

def test_time_derivative1_arrays():
    p = generate_square_cloud(11)
    t_steps = 10
    v = 0.2
    T = np.linspace(0, 1, t_steps)
    
    # Callable function for heat equation
    f_func = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)
    
    # Create Spatiotemporal array (m, t)
    f_arr = np.zeros((len(p), t_steps))
    for k in range(t_steps):
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k], [v])
        
    L = np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])
    
    # Run with Callable
    u_ap_func, _ = TimeDerivative1(p, f_func, t_steps, [v], operator=L, implicit=True, verbose=False)
    
    # Run with Array
    u_ap_arr, _ = TimeDerivative1(p, f_arr, t_steps, [v], operator=L, implicit=True, verbose=False)
    
    error = np.max(np.abs(u_ap_func - u_ap_arr))
    assert error < 1e-12

def test_time_derivative2_arrays():
    p = generate_square_cloud(11)
    t_steps = 10
    c = float(np.sqrt(1 / 2))
    T = np.linspace(0, 1, t_steps)
    
    # Callable functions for wave equation
    f_func = lambda x, y, t, coef: np.cos(np.pi * t) * np.sin(np.pi * (x + y))
    g_func = lambda x, y, t, coef: -np.pi * np.sin(np.pi * t) * np.sin(np.pi * (x + y))
    
    # Arrays
    f_arr = np.zeros((len(p), t_steps))
    g_arr = np.zeros(len(p))
    
    for k in range(t_steps):
        f_arr[:, k] = f_func(p[:, 0], p[:, 1], T[k], [c])
    
    # g is evaluated at k=1 in solver (velocity)
    g_arr = g_func(p[:, 0], p[:, 1], T[1], [c])
    
    L = np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])
    
    # Run with Callable
    u_ap_func, _ = TimeDerivative2(p, f_func, g_func, t_steps, [c], operator=L, implicit=True, verbose=False)
    
    # Run with Array
    u_ap_arr, _ = TimeDerivative2(p, f_arr, g_arr, t_steps, [c], operator=L, implicit=True, verbose=False)
    
    error = np.max(np.abs(u_ap_func - u_ap_arr))
    assert error < 1e-12
