# 📈 mGFD Changelog

All notable changes to the **mGFD** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.11.0]
*The High-Performance & GPU Acceleration Release. (Not yet published on PyPI)*

### ✨ Added
- **⚡ GPU Acceleration (CuPy)**: Fully integrated `device="cuda"` support for all PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`). This allows offloading algebraic bottlenecks (sparse matrices and iterative solvers) directly to NVIDIA GPUs via `cupyx.scipy.sparse`.
- **🚀 Krylov Preconditioning**: Added optional preconditioning parameter (`preconditioner="ilu"` or `"jacobi"`) to all solvers for vastly improved `GMRES` and `BiCGStab` convergence.
- **🏎️ Matrix-Free Computing**: Implemented `matrix_free=True` mode for all solvers, utilizing Numba JIT-compiled on-the-fly matrix-vector multiplications via `scipy.sparse.linalg.LinearOperator`, drastically reducing memory usage for iterative solvers.
- **🛠️ Core Functionality**: Added `compute_preconditioner` inside `gammas.py` generating `scipy.sparse.linalg.LinearOperator` wrappers for ILU factorization and diagonal scaling (with dynamic CuPy compatibility).
- **🧠 Adaptive Neighborhoods (Dynamic `nvec`)**: Modified KDTree search algorithms in `neighbors.py` (`_find_neighbors_balanced_jit`) to automatically build density-aware, condition-aware dynamic stencils using `np.linalg.cond` internally. The `nvec` solver parameter is now treated as an upper bound rather than a strict count.

### 🔧 Changed
- **🔄 Benchmark Orchestration**: Completely revamped `sweep.py` and `batch_utils.py` to automatically crawl JSON metrics and generate a consolidated `sweep_summary.csv` master file. This eliminates manual data aggregation for CPU vs GPU performance analysis.

### 🐛 Fixed
- **💥 Numba Segmentation Fault**: Fixed a critical memory fragmentation issue (segfault) in `neighbors.py` that occurred on Linux CI runners. Numba JIT loops now use strictly preallocated scratch arrays (`dx_arr`, `dy_arr`) to guarantee ABI stability across all operating systems.
- **🐢 CUDA PCIe Bottleneck**: Resolved a massive performance degradation in `TimeDerivative1` and `TimeDerivative2` by explicitly transferring boundary (`boun_n`) and interior (`inne_n`) index arrays to VRAM (`cp.asarray`). This eliminates implicit host-device PCIe synchronization during the time integration loop, achieving the true theoretical speedup for explicit and factorized solvers on the GPU.

---

## [0.10.0]
*The Architecture & Standardization Release.*

### ✨ Added
- **📊 Pandas/Xarray Adapters**: `src/mGFD/core/adapters.py` to seamlessly integrate DataFrames and DataArrays using a lightweight soft-dependency approach.
- **📦 Structured Results**: `src/mGFD/solvers/results.py` module defining the new `SolverResult` dataclass to standardize the output of all PDE solvers.
- **🛡️ Custom Exceptions**: `src/mGFD/exceptions.py` introducing a robust custom exception hierarchy (`mGFDError`, `CloudShapeError`, `InputTypeError`, `DimensionMismatchError`, `OperatorFormatError`, `ParameterError`).

### 🔧 Changed
- **🔄 Batch Processing**: Restored and upgraded `research/codes/` batch processing scripts to natively utilize the v0.10.0 `SolverResult` and flexible algebraic backends.
- **🧹 Cleanup**: Removed the redundant `run_Perturbation2.py` batch script.
- **🔄 API Standardization**: Standardized the output of all PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) to return the new `SolverResult` dataclass instead of returning raw tuples.
- **🔙 Backward Compatibility**: Implemented `__iter__` on `SolverResult` to allow 100% backward-compatible tuple unpacking (`u_ap, vec = solver(...)`), avoiding breaking changes.
- **🧮 Flexible Linear Solvers**: Integrated flexible algebraic linear solvers (`spsolve`, `bicgstab`, `gmres`) in all PDE solvers, avoiding memory bottlenecks for large systems.
- **📈 Warm-Starting**: Implemented $x_0$ warm-starting (solution from step $k-1$) in transient solvers to exponentially accelerate Krylov subspace convergence.
- **🛡️ Error Handling**: Refactored solvers to raise custom exceptions instead of generic `ValueError` or `TypeError`.
- **🧪 Test Suite**: Updated test suite (`tests/`) to validate against the new custom exception hierarchy.
- **🏎️ Numba Optimizations**: Integrated Numba (`@nb.njit`) to optimize core iterative loops in `gammas.py`, achieving massive performance gains for meshless operator constructions.
- **🏎️ KDTree Optimizations**: Refactored KDTree neighbor queries in `neighbors.py` to use SciPy's vectorized bulk queries coupled with Numba JIT functions.
- **🏎️ Ray-Casting**: Optimized `cloud_generator` bottlenecks via Numba JIT compilation by replacing `matplotlib.path` logic with a fast ray-casting algorithm.
- **🏎️ Spatial Queries**: Upgraded point classification in `classification.py` to use fully vectorized GEOS routines and Shapely 2.0 `spatial_index.query` arrays.

### 🐛 Fixed
- **📝 Typing**: Corrected type hinting signatures in `TimeDerivative1` and `TimeDerivative2` to properly reflect the `SolverResult` return type, resolving static analyzer (mypy/pylance) warnings.

---

## [0.9.1]
*The PyPI & CI Release.*

### ✨ Added
- **📦 PyPI Release**: Official publication of version `0.9.1` on PyPI.
- **🤖 Continuous Integration**: CI workflow in GitHub Actions with tests for Python 3.9, 3.10, 3.11, 3.12, and 3.13.
- **📝 Typing**: Massive refactoring to comply with strict Type Hints (Mypy).
- **📝 Compatibility**: Mathematical signature correction and type compatibility (NumPy/SciPy/PyVista/OpenCV).
- **📦 TestPyPI**: Pre-release and verification on TestPyPI.

---

## [0.9.0]
*The Modularity & Validation Release.*

### ✨ Added
- **📝 Changelog**: `CHANGELOG.md` to track ongoing project modifications.
- **🧩 Modularity**: `src/mGFD/cloud_generator/core/point_generation/` package to modularize the cloud generation logic (geometry, boundary, relaxation, regular, poisson).
- **🧪 Single Testing**: `research/batches/test_single.py` for rapid validation of all PDE solvers on a single point cloud.
- **🧪 Pytest Suite**: `tests/` directory with a robust `pytest` suite for validating mathematical solvers and cloud generation logic.
- **📊 Logging**: Standard `logging` infrastructure across all modules, replacing raw `print()` statements for better external consumption.
- **🛡️ Input Validation**: Comprehensive `fail-fast` input validation in all solvers to catch dimensionality mismatches early.
- **🤖 GitHub Actions**: CI workflow (`.github/workflows/python-tests.yml`) for automated testing across Python 3.9-3.12.

### 🔧 Changed
- **🧩 Refactoring**: Refactored `point_generation.py` monolith into multiple modular files.
- **🏷️ Renaming**: Renamed `grid_generation` functions and files to `regular_generation` to better align with the "meshless" philosophy.
- **📖 Documentation**: Updated `mgfd-cloud` CLI help text and `README.md` to reflect the new `regular` terminology and provide better examples.
- **🏛️ Architecture**: Architecturally decoupled the theoretical exact solution (`u_ex`) from the solvers. They now return `Tuple[np.ndarray, np.ndarray]` (`u_ap, vec`).
- **🔄 Scripts**: Adapted all benchmarking scripts in `research/batches/` to compute `u_ex` locally instead of relying on the solvers.
- **📦 Configuration**: Modernized project configuration entirely into a PEP 621 compliant `pyproject.toml` with updated PyPI classifiers and metadata.
- **📝 Typing**: Improved Type Hinting (`mypy`) coverage throughout all core modules and utilities.

### 🐛 Fixed
- **🧹 Linting**: Addressed multiple syntax, indentation, and f-string errors flagged by `flake8`.
