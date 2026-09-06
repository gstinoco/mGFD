"""
Test CFL — Unit tests for Adaptive Time-Stepping and CFL estimation in mGFD (OOP Interface)

Overview:
    This test module validates the Courant-Friedrichs-Lewy (CFL) time step estimation routines,
    spatial mesh spacing calculations, and autonomous time-stepping resolution in the OOP Solver.

Public API:
    sample_cloud                        Pytest fixture generating a 2D regular grid point cloud.
    test_compute_mesh_spacing           Validates characteristic spatial node spacing estimation.
    test_estimate_cfl_dt                Validates CFL step estimation for 1st and 2nd order operators.
    test_adaptive_time_derivative1      Validates autonomous time-step resolution in 1st order Solver.
    test_adaptive_time_derivative2      Validates autonomous time-step resolution in 2nd order Solver.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
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
    September, 2026.
Last Modification:
    September, 2026.
"""

## Library importation.
import pytest                                                                                                                           # Test framework.
import numpy as np                                                                                                                      # Core numerical array operations.

from mGFD.spatial.neighbors import compute_mesh_spacing                                                                                 # Mesh spacing routine.
from mGFD.temporal.cfl import estimate_cfl_dt                                                                                           # CFL estimation function.
from mGFD import Cloud, Dirichlet, PDE, Solver                                                                                          # High-level OOP API classes.

@pytest.fixture
def sample_cloud() -> np.ndarray:
    """
    sample_cloud
    Build a simple 2D regular grid point cloud [x, y, flag].

    Output:
        p           m x 3           ndarray         Constructed regular grid point cloud.
    """
    x           = np.linspace(0, 1, 10)                                                                                                 # Coordinate sampling along x.
    y           = np.linspace(0, 1, 10)                                                                                                 # Coordinate sampling along y.
    xx, yy      = np.meshgrid(x, y)                                                                                                     # 2D Cartesian coordinate grid.
    m           = xx.size                                                                                                               # Total node count.
    p           = np.zeros((m, 3))                                                                                                      # Allocate point array.
    p[:, 0]     = xx.ravel()                                                                                                            # Assign x coordinates.
    p[:, 1]     = yy.ravel()                                                                                                            # Assign y coordinates.
    
    is_boundary = (p[:, 0] == 0) | (p[:, 0] == 1) | (p[:, 1] == 0) | (p[:, 1] == 1)                                                     # Flag boundary coordinates.
    p[:, 2]     = np.where(is_boundary, 1, 0)                                                                                           # Assign boundary flags.
    return p                                                                                                                            # Return point cloud array.

def test_compute_mesh_spacing(sample_cloud: np.ndarray) -> None:
    """
    test_compute_mesh_spacing
    Test estimation of minimum and average spatial node spacing.

    Input:
        sample_cloud    m x 3       ndarray         Input point cloud fixture.
    """
    # 1. Execute spacing computation
    h_min, h_avg = compute_mesh_spacing(sample_cloud)                                                                                   # Compute spatial metrics.
    
    # 2. Verify bounds and values
    assert h_min > 0.0                                                                                                                  # Verify positive min spacing.
    assert h_avg >= h_min                                                                                                               # Verify mean exceeds or equals min.
    assert pytest.approx(h_min, abs=1e-3) == 1.0 / 9.0                                                                                  # Verify regular grid spacing dx = 1/9.

def test_estimate_cfl_dt(sample_cloud: np.ndarray) -> None:
    """
    test_estimate_cfl_dt
    Test estimate_cfl_dt for parabolic and wave operators.

    Input:
        sample_cloud    m x 3       ndarray         Input point cloud fixture.
    """
    # 1. Setup operator
    operator = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                                # Heat operator: A=2, C=2.
    
    # 2. Parabolic test (Order 1)
    dt1, t1, cfl1 = estimate_cfl_dt(sample_cloud, operator, cfl=0.5, order=1, t_end=1.0)                                                # Estimate 1st order CFL.
    assert dt1 > 0.0                                                                                                                    # Verify positive time step.
    assert t1 > 0                                                                                                                       # Verify positive step count.
    assert cfl1 > 0.0                                                                                                                   # Verify positive actual CFL.
    
    # 3. Hyperbolic test (Order 2)
    dt2, t2, cfl2 = estimate_cfl_dt(sample_cloud, operator, cfl=0.5, order=2, t_end=1.0)                                                # Estimate 2nd order CFL.
    assert dt2 > 0.0                                                                                                                    # Verify positive time step.
    assert t2 > 0                                                                                                                       # Verify positive step count.
    assert cfl2 > 0.0                                                                                                                   # Verify positive actual CFL.

def test_adaptive_time_derivative1(sample_cloud: np.ndarray) -> None:
    """
    test_adaptive_time_derivative1
    Test Solver with PDE (order 1) with automatic t resolution via cfl parameter.

    Input:
        sample_cloud    m x 3       ndarray         Input point cloud fixture.
    """
    # 1. Formulation
    def f_init(x: np.ndarray, y: np.ndarray, t_val: float = 0.0, coef: list | None = None) -> np.ndarray:                               # Analytical initial function.
        return np.exp(-2 * 0.2 * np.pi**2 * t_val) * np.sin(np.pi * x) * np.sin(np.pi * y)                                              # Spatial-temporal state.
        
    v      = 0.2                                                                                                                        # Diffusion coefficient.
    cloud  = Cloud.from_array(sample_cloud)                                                                                             # Instantiate Cloud object.
    domain = cloud.set_boundary(Dirichlet(f_init))                                                                                      # Set Dirichlet boundary condition.
    pde    = PDE(operator=[0, 0, 2*v, 0, 2*v, 0], order=1)                                                                              # Formulate 1st order Heat PDE.
    solver = Solver(domain, pde, cfl=0.5, verbose=False)                                                                                # Instantiate Solver.
    
    # 2. Execution and Verification
    res    = solver.solve()                                                                                                             # Solve with adaptive CFL.
    assert res.converged                                                                                                                # Verify convergence status.
    assert res.dt is not None and res.dt > 0.0                                                                                          # Verify computed dt.
    assert res.cfl is not None and res.cfl > 0.0                                                                                        # Verify computed CFL.
    assert res.t_steps is not None and res.t_steps > 0                                                                                  # Verify computed step count.
    
    u_ap, vec = res                                                                                                                     # Unpack result tuple.
    assert u_ap.shape[0] == sample_cloud.shape[0]                                                                                       # Verify node dimension.
    assert u_ap.shape[1] == res.t_steps                                                                                                 # Verify time steps dimension.

def test_adaptive_time_derivative2(sample_cloud: np.ndarray) -> None:
    """
    test_adaptive_time_derivative2
    Test Solver with PDE (order 2) with automatic t resolution via cfl parameter.

    Input:
        sample_cloud    m x 3       ndarray         Input point cloud fixture.
    """
    # 1. Formulation
    def f_init(x: np.ndarray, y: np.ndarray, t_val: float = 0.0, coef: list | None = None) -> np.ndarray:                               # Analytical initial displacement.
        return np.cos(np.pi * 0.5 * np.sqrt(2) * t_val) * np.sin(np.pi * x) * np.sin(np.pi * y)                                         # Spatio-temporal wave state.
        
    c      = 0.5                                                                                                                        # Wave propagation speed.
    cloud  = Cloud.from_array(sample_cloud)                                                                                             # Instantiate Cloud object.
    domain = cloud.set_boundary(Dirichlet(f_init))                                                                                      # Set Dirichlet boundary condition.
    pde    = PDE(operator=[0, 0, 2*c**2, 0, 2*c**2, 0], g=0.0, order=2)                                                                 # Formulate 2nd order Wave PDE.
    solver = Solver(domain, pde, cfl=0.5, verbose=False)                                                                                # Instantiate Solver.
    
    # 2. Execution and Verification
    res    = solver.solve()                                                                                                             # Solve with adaptive CFL.
    assert res.converged                                                                                                                # Verify convergence status.
    assert res.dt is not None and res.dt > 0.0                                                                                          # Verify computed dt.
    assert res.cfl is not None and res.cfl > 0.0                                                                                        # Verify computed CFL.
    assert res.t_steps is not None and res.t_steps > 0                                                                                  # Verify computed step count.
    
    u_ap, vec = res                                                                                                                     # Unpack result tuple.
    assert u_ap.shape[0] == sample_cloud.shape[0]                                                                                       # Verify node dimension.
    assert u_ap.shape[1] == res.t_steps                                                                                                 # Verify time steps dimension.
