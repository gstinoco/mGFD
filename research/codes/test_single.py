"""
test_single — Single cloud solver validation

Overview:
    Quick diagnostic script to validate that all PDE solvers (Stationary, Heat,
    Advection, Advection-Diffusion, Wave) run without crashing on a single sample cloud.
    It utilizes the benchmarking functions to get unified metrics.

Usage:
    python research/runners/test_single.py

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
    May, 2024.
Last Modification:
    August, 2026.
"""

## Library importation.
import os                                                                                                                               # Filesystem and path utilities.
import sys                                                                                                                              # sys.path manipulation.
import logging                                                                                                                          # Standard logging module.

from typing import Callable, List, Tuple, Dict, Any                                                                                     # Type hinting.

from mGFD.io.io import load_points                                                                                                      # Point cloud loading utility.

BASE_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                                             # Research root directory (for local utils).

if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)                                                                                                           # Allow importing from research/codes/.

import benchmarks.benchmark as bmk                                                                                                      # Benchmarking functions to test.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.
logging.basicConfig(level=logging.INFO, format='%(message)s')                                                                           # Basic logger configuration.

def test_single_cloud() -> None:
    """
    test_single_cloud
    Quickly test all solvers on a single point cloud.
    Input:
        None

    Output:
        None
    """
    # 1. Load test cloud
    cloud_path: str = os.path.join(BASE_DIR, 'data', 'Titicaca', '1', 'Titicaca_cloud.csv')                                             # Path to sample cloud.
    if not os.path.exists(cloud_path):                                                                                                  # Validate file existence.
        logger.error(f"Test cloud not found at: {cloud_path}")                                                                          # Report error.
        return                                                                                                                          # Abort test.

    logger.info(f"Starting quick test on a single cloud: {cloud_path}")                                                                 # Report startup.
    p = load_points(cloud_path, verbose=False)                                                                                          # Load point cloud silently.
    logger.info(f"Cloud loaded with {len(p)} points.")                                                                                  # Report loaded point count.
    
    # 2. Define solvers to test
    tests: List[Tuple[str, Callable]] = [                                                                                               # List of benchmarks to execute.
        ('Poisson',             bmk.benchmark_poisson_equation),                                                                        # Poisson test.
        ('Heat (t=2000)',       bmk.benchmark_heat_equation),                                                                           # Heat test.
        ('Advection',           bmk.benchmark_advection_equation),                                                                      # Pure advection test.
        ('Advection-Diffusion', bmk.benchmark_advdif_equation),                                                                         # Advection-diffusion test.
        ('Wave',                bmk.benchmark_wave_equation)                                                                            # Wave test.
    ]                                                                                                                                   # End of tests list.
    
    # 3. Execute tests
    for name, func in tests:                                                                                                            # Iterate all defined tests.
        logger.info(f"\n[{name}] Executing...")                                                                                         # Report current test.
        try:                                                                                                                            # Wrap execution to catch solver crashes.
            res: Dict[str, Any] = func(p)                                                                                               # Execute benchmark function on cloud.
            logger.info(f"  [SUCCESS] Time: {res['execution_time_seconds']:.3f}s | Avg RMSE: {res['avg_numerical_error']:.3e}")         # Report metrics.
        except Exception as e:                                                                                                          # Catch any exception.
            logger.error(f"  [ERROR] Test {name} failed: {e}")                                                                          # Report failure reason.

if __name__ == '__main__':                                                                                                              # Script entry point.
    test_single_cloud()                                                                                                                 # Run the single cloud test suite.
