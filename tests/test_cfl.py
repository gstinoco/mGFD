"""
test_cfl.py — Unit tests for Adaptive Time-Stepping and CFL estimation in mGFD
"""

import pytest
import numpy as np

from mGFD.spatial.neighbors import compute_mesh_spacing
from mGFD.temporal.cfl import estimate_cfl_dt
from mGFD.solvers import TimeDerivative1, TimeDerivative2

@pytest.fixture
def sample_cloud():
    """Build a simple 2D regular grid point cloud [x, y, flag]."""
    x = np.linspace(0, 1, 10)
    y = np.linspace(0, 1, 10)
    xx, yy = np.meshgrid(x, y)
    m = xx.size
    p = np.zeros((m, 3))
    p[:, 0] = xx.ravel()
    p[:, 1] = yy.ravel()
    
    # Boundary flag: 1 if on boundary, 0 interior
    is_boundary = (p[:, 0] == 0) | (p[:, 0] == 1) | (p[:, 1] == 0) | (p[:, 1] == 1)
    p[:, 2] = np.where(is_boundary, 1, 0)
    return p

def test_compute_mesh_spacing(sample_cloud):
    """Test estimation of minimum and average spatial node spacing."""
    h_min, h_avg = compute_mesh_spacing(sample_cloud)
    assert h_min > 0.0
    assert h_avg >= h_min
    # For a 10x10 regular grid on [0,1]^2, spacing dx = dy = 1/9 ~ 0.1111
    assert pytest.approx(h_min, abs=1e-3) == 1.0 / 9.0

def test_estimate_cfl_dt(sample_cloud):
    """Test estimate_cfl_dt for parabolic and wave operators."""
    operator = np.vstack([[0], [0], [2], [0], [2], [0]])  # Heat operator: A=2, C=2
    
    # Parabolic (Order 1)
    dt1, t1, cfl1 = estimate_cfl_dt(sample_cloud, operator, cfl=0.5, order=1, t_end=1.0)
    assert dt1 > 0.0
    assert t1 > 0
    assert cfl1 > 0.0
    
    # Hyperbolic (Order 2)
    dt2, t2, cfl2 = estimate_cfl_dt(sample_cloud, operator, cfl=0.5, order=2, t_end=1.0)
    assert dt2 > 0.0
    assert t2 > 0
    assert cfl2 > 0.0

def test_adaptive_time_derivative1(sample_cloud):
    """Test TimeDerivative1 with automatic t resolution via cfl parameter."""
    def f_init(x, y, t_val, coef):
        return np.exp(-2 * coef[0] * np.pi**2 * t_val) * np.sin(np.pi * x) * np.sin(np.pi * y)
        
    v = 0.2
    # Solve without explicit t parameter
    res = TimeDerivative1(sample_cloud, f_init, coef=[v], cfl=0.5, verbose=False)
    
    assert res.converged
    assert res.dt is not None and res.dt > 0.0
    assert res.cfl is not None and res.cfl > 0.0
    assert res.t_steps is not None and res.t_steps > 0
    
    # Verify backward compatible tuple unpacking
    u_ap, vec = res
    assert u_ap.shape[0] == sample_cloud.shape[0]
    assert u_ap.shape[1] == res.t_steps

def test_adaptive_time_derivative2(sample_cloud):
    """Test TimeDerivative2 with automatic t resolution via cfl parameter."""
    def f_init(x, y, t_val, coef):
        return np.cos(np.pi * coef[0] * np.sqrt(2) * t_val) * np.sin(np.pi * x) * np.sin(np.pi * y)
    def g_init(x, y, t_val, coef):
        return np.zeros_like(x)
        
    c = 0.5
    # Solve without explicit t parameter
    res = TimeDerivative2(sample_cloud, f_init, g_init, coef=[c], cfl=0.5, verbose=False)
    
    assert res.converged
    assert res.dt is not None and res.dt > 0.0
    assert res.cfl is not None and res.cfl > 0.0
    assert res.t_steps is not None and res.t_steps > 0
    
    # Verify backward compatible tuple unpacking
    u_ap, vec = res
    assert u_ap.shape[0] == sample_cloud.shape[0]
    assert u_ap.shape[1] == res.t_steps
