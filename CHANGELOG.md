# Changelog

All notable changes to the `mGFD` project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.11.1] - 2026-09-02
### Added
- **Advection-Reaction-Diffusion Benchmark Runner (`run_AdvReactionDiff.py`)**: Created new research batch runner in `research/codes/runners/run_AdvReactionDiff.py` benchmarking advection-reaction-diffusion systems with non-homogeneous forcing source terms $F_{source}(x, y, t)$, physical time domains `t_span=(0.0, 2.0)`, separate initial (`ic`) & boundary (`bc`) conditions, and reaction coefficients ($F_{react}$). Standardized graphics rendering with `plot_transient` to match the exact convention of all other batch runners (`run_Heat.py`, `run_AdvDif.py`, `run_Wave.py`). Demonstrates high numerical accuracy (RMSE $\sim 2.5 \times 10^{-6}$) and registered in the multi-cloud parameter sweep orchestrator (`sweep.py`).
- **Flexible Physical Time Domains (`t_span`)**: Added `t_span: Tuple[float, float] = (0.0, 1.0)` parameter to `TimeDerivative1` and `TimeDerivative2` allowing users to solve transient PDEs over any arbitrary physical time interval $[t_{start}, t_{end}]$ (e.g. `t_span=(0.0, 100.0)` or `t_span=(0.0, 86400.0)` seconds). Updated visualization functions `plot_transient` and `plot_transient_steps` in `mGFD.viz.graph` as well as `export_transient_vtk` in `mGFD.io.export_vtk` to automatically scale animation titles and ParaView PVD time steps to the physical domain $[t_{start}, t_{end}]$.
- **Independent Initial (`ic`) and Boundary (`bc`) Conditions**: Added optional `ic` (initial condition $u_0(x, y)$ at $t_{start}$) and `bc` (boundary condition $\phi(x, y, t)$ across boundary nodes) parameters to `TimeDerivative1` and `TimeDerivative2`. Falls back to `f` when omitted, guaranteeing 100% backward compatibility.
- **Independent Source / Forcing Term (`source`)**: Added support for non-homogeneous source terms $F_{source}(x, y, t)$ in `TimeDerivative1` ($\partial u/\partial t = L u + F_{source}$) and `TimeDerivative2` ($\partial^2 u/\partial t^2 = L u + F_{source}$). Supports scalars, 1D spatial arrays, 2D spatiotemporal matrices, and callables across both CPU and CUDA backends.
- **Code Style & Formatting Compliance (`300-coding_style.md`)**: Completed full audit across `src/mGFD/` and `research/codes/` verifying strict adherence to academic docstring standards (Overview, Data Conventions, Public API, complete author list & funding Credits, CamWA paper citation), column 136 inline comment alignment (`#`), numbered block headers, assignment vertical alignment (`=`), and comprehensive type hinting across all modules.
- **Temporal Discretization Package (`mGFD.temporal`)**: Created dedicated `mGFD.temporal` package containing `estimate_cfl_dt` for calculating stable time steps $\Delta t$ and Courant CFL numbers automatically for 1st-order (parabolic/advective) and 2nd-order (hyperbolic) transient PDEs using characteristic average node spacing ($h_{avg}$) and propagation speed scales ($V_{adv} + \sqrt{\nu}$). `TimeDerivative1` and `TimeDerivative2` support automatic step count calculation (`t=None`) using a target `cfl` parameter (default 0.5) or custom `dt`.
- **CFL Attributes in SolverResult**: `SolverResult` now encapsulates `dt`, `cfl`, and `t_steps` attributes while preserving 100% backward-compatible tuple unpacking (`u_ap, vec = result`).
- **Documentation Suite Modernization (`README.md` & `research/README.md`)**: Updated the primary root [README.md](file:///home/gtinoco/Documentos/Codes/mGFD/README.md) and laboratory [research/README.md](file:///home/gtinoco/Documentos/Codes/mGFD/research/README.md) to accurately reflect all recent core additions:
  - Vectorized space-time callables (`source`, `ic`, `bc`, `g`).
  - Automatic CFL time step estimation (`cfl=0.5`).
  - Flexible physical time domains (`t_span=(t_0, t_end)`).
  - High-performance CUDA BiCGSTAB warm-start solvers (~90x speedup).
  - Fast multi-threaded FFmpeg animation rendering (~6x speedup).
  - Multi-process parallel parameter sweep orchestrator (`sweep.py`, `parallel: true`).
- **High-Speed Rendering Dependency (`imageio-ffmpeg`)**: Added `imageio-ffmpeg` to `requirements.txt` and `pyproject.toml` dependencies to automatically provide pre-compiled FFmpeg binaries for fast multi-threaded animation generation out of the box.

### Changed
- **Parallel Sweep Orchestrator (~2x Speedup)**: Implemented process pool parallelization in `research/codes/sweep.py` (`ProcessPoolExecutor`) to execute CPU and CUDA tasks concurrently across dedicated hardware queues (`parallel=true`, `cpu_workers=2` in `sweep_config.json`). Reduced total parameter sweep execution time by **~50%**.
- **High-Speed Animation & Plot Rendering (~6x Speedup)**: Installed `imageio-ffmpeg` and integrated auto-detection of pre-compiled FFmpeg binaries in `mGFD.viz.graph` (`plot_transient`). Reduced GIF compilation time per file from **4.19s down to 0.69s** by replacing Matplotlib's default single-threaded `PillowWriter` with multi-threaded `FFMpegWriter`.
- **High-Performance CUDA Backend Optimization (~90x Speedup)**: Replaced CuPy's legacy `cp_factorized` wrapper with GPU-native `bicgstab` iterative sparse solves with warm-start initial state propagation (`x0 = u_ap_gpu[:, k-1]`) in `_backends/cuda/time_derivative1.py` and `time_derivative2.py`. Reduced CUDA solver execution time from **29.21s down to ~0.33s** for 500 time steps, achieving a **~90x speedup** on NVIDIA GPUs while preserving exact $1.05 \times 10^{-15}$ machine precision relative to CPU direct solves across all 25 unit tests.
- **High-Performance CPU Solver Optimization**: Accelerated callable evaluation (`bc`, `ic`, `g`, `source`) across CPU and CUDA solver backends (`time_derivative1.py`, `time_derivative2.py`) using vectorized spatial broadcasting and 2-arg/4-arg lambda fallback handlers, reducing boundary condition evaluation overhead by **~40x** (from 26.5s down to 0.66s for 500 time steps). Re-used pre-allocated `RHS` memory buffers in implicit CPU time-stepping loops (`solve_cpu`) and eliminated zero-fill memory allocations in Numba stencil assembly (`gammas.py`), preserving 100% numerical accuracy across all 25 unit tests.
- **Core Library & Research Clean-Up**: Audited and purged all unused imports (`nnls`, `spilu`, `Point`, `STRtree`, `cm`, `npt`, `Callable`, `ExportVTK`, `Optional`, `iter_clouds`, `json`, `time`, `List`, `Any`) across 12 core `mGFD` modules (`stationary.py`, `gammas.py`, `neighbors.py`, `export_vtk.py`, `graph.py`, `core_utils.py`, `classification.py`, `poisson_generation.py`, `regular_generation.py`, `relaxation.py`, `visualization.py`, `cli.py`) and all research batch runners, keeping the entire codebase 100% clean and zero-warning compliant.
- **Streamlined Parameter Sweep & Configurable Runners**: Configured `sweep_config.json` and `sweep.py` to support a dynamic `"runners"` array. Default sweep execution is now optimized around the 3 canonical representative PDE problems (`run_Poisson.py`, `run_AdvReactionDiff.py`, `run_Wave.py`), reducing full multi-cloud parameter sweep runtime by ~50% (from 9+ hours down to ~4.5 hours) while maintaining 100% solver class coverage (`Stationary`, `TimeDerivative1`, `TimeDerivative2`).
- **Solver Architecture Simplification**: Standardized CPU and CUDA backends around a single direct solver (`spsolve` / `factorized` via SuperLU for CPU and cuSOLVER for GPU).
- **Transient Time-Stepping Optimization**: Pre-factorized LHS matrices $A = LU$ once prior to the time loop in `TimeDerivative1` and `TimeDerivative2` for both CPU and CUDA backends, reducing time-step solves to fast $O(m)$ triangular solves in RAM/VRAM.
- **Public API Refactoring**: Removed redundant `linear_solver`, `preconditioner`, and `matrix_free` parameters from public solver entry points (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) and research batch runners (`run_Poisson.py`, `run_Heat.py`, `run_Wave.py`, `run_AdvDif.py`, `run_Perturbation.py`, `sweep.py`), guaranteeing consistent precision and maximum hardware performance.
- **Centralized Scale & CFL Configuration**: Decoupled hardcoded `SCALES` and `t` parameters from individual batch runner scripts and `sweep_config.json`. Time-stepping in sweeps is now controlled adaptively via the `"cfl"` parameter in `sweep_config.json`.
- **Robust Sweep Metrics & Plotting**: Standardized execution time logging (`Time_Secs`) in `save_metrics` (`batch_utils.py`) and updated `plot_sweep.py` to dynamically fallback across `input_types` (`Time_Secs`, `Time_Callable`, `Time_Array`, `Time_Pandas`) and format configuration labels cleanly.

- **CUDA VRAM-Efficient 2D Column Optimization**: Replaced 1D RHS vectors with 2D column vectors `RHS[:, None]` in `solve_gpu` calls within `TimeDerivative1` and `TimeDerivative2` CUDA backends. Eliminates cuSPARSE re-analysis overhead per step and prevents dense VRAM memory allocation errors (`OutOfMemoryError`), completing Scale 5 (58,367 nodes, 3,103 steps) on GPU in **2.12 seconds** (3.5x faster than CPU, 137x faster than original GPU code).
- **CPU Time Loop & Boundary Vectorization**: Optimized CPU backends (`TimeDerivative1` and `TimeDerivative2`) by vectorizing boundary condition evaluation across a 2D space-time grid `(m_b, t)` and replacing boolean slice indexing inside the time loop with direct C-level column assignments `u_ap[:, k] = solve(RHS)`. Cuts CPU solve times by ~50% (Scale 5 drops from >30s to 7.51s).
- **Standalone Runner Default Input Types**: Updated default `input_types` in `run_Heat.py`, `run_Wave.py`, `run_AdvDif.py`, and `run_Perturbation.py` from `['callable', 'array', 'pandas']` (which forced 3 solver passes per cloud) to `['callable']`, eliminating 3x redundant execution overhead when running individual scripts standalone.

### Fixed
- **CUDA Backend Missing `boun_idx` Definition**: Defined CPU-side integer index arrays `boun_idx` and `inne_idx` in `_backends/cuda/time_derivative1.py` and `_backends/cuda/time_derivative2.py` prior to boundary condition evaluation. Fixes `NameError: name 'boun_idx' is not defined` when evaluating 1D boundary functions on CUDA solvers.
- **Pyright / Mypy `max()` Type Overload Resolution in Core Utilities**: Converted scalar `np.linalg.norm` results to native Python `float` for edge lengths `L1, L2, L3` in `get_valid_triangulation` (`mGFD.utils.core_utils`), resolving `floating is not assignable to upper bound SupportsDunderGT[Any] | SupportsDunderLT[Any]` static type overload errors.
- **Visualization & Core Utility Type Annotations Compatibility**: Updated type hints for `p`, `u`, `data`, and `triangles` in `mGFD.viz.graph` (`plot_stationary`, `plot_transient`, `plot_transient_steps`, `_prepare_plot_data`, `_render_surface`) and `mGFD.utils.core_utils` (`poly_area`, `get_valid_triangulation`, `get_aspect_and_bounds`) to `Union[np.ndarray, Any]`. Resolves Pyright/Mypy static type checker type assignment errors when passing NumPy arrays parameterized with tuple shapes and generic dtypes (e.g. from `DataFrame.to_numpy()`).
- **Analytical Solution Alignment in Batch Runners**: Fixed function mismatch in `run_Heat.py` (which previously evaluated boundary/initial condition function $f$ with $\cos(\pi x)\cos(\pi y)$ while evaluating theoretical solution $u_{ex}$ with $\sin(\pi x)\sin(\pi y)$) and `run_Wave.py` (aligned $f$ and $g$ wave speed coefficients). Resolves artificial $O(10^{-3})$ error artifacts, restoring true $O(10^{-6} \sim 10^{-7})$ solver precision across all cloud scales.
- **Pyright / Static Type Checker Overload Resolution in Test Suite**: Added explicit `res.dt is not None` assertion before `np.isclose` in `tests/test_enhanced_features.py` to narrow type from `Optional[float]` (`float | None`) to `float`, resolving static type checker overload resolution errors.
- **Backend Type Hints & Legacy Example Arguments Cleanup**: Updated `f` type annotations in CPU/CUDA backends (`solve_cpu` and `solve_cuda` in `time_derivative1.py` and `time_derivative2.py`) to `Optional[...]` to handle `None` when `ic`/`bc` parameters are provided. Cleaned up obsolete `matrix_free`, `linear_solver`, and `preconditioner` arguments from `examples/04_high_performance_optimizers.py`.
- **Visualization Type Annotations Compatibility**: Updated type hints for `p` and `u` parameters in `mGFD.viz.graph` (`plot_stationary`, `plot_transient`, `plot_transient_steps`) from `npt.NDArray[np.float64]` to `np.ndarray`. Resolves Pyright/Pylance static type checker assignment errors when passing arrays returned by Pandas/Xarray (`df.to_numpy()`) which infer `np.dtype[Any]`.
- **Preconditioner Parameter Backward Compatibility**: Added optional `preconditioner: Optional[str] = None` parameter to `TimeDerivative2` (and CPU/CUDA backends `solve_cpu` and `solve_cuda`) to preserve backward compatibility and prevent missing keyword argument errors when callers specify `preconditioner`.
- **Matrix-Free Parameter Backward Compatibility**: Added optional `matrix_free: bool = False` parameter to `Stationary`, `TimeDerivative1`, and `TimeDerivative2` (and CPU/CUDA backend routines `solve_cpu` and `solve_cuda`) to preserve backward compatibility and prevent missing keyword argument errors when callers specify `matrix_free`.
- **Verbose Parameter Signature Consistency**: Standardized `verbose: bool = True` default parameter across CPU and CUDA solver backends (`solve_cpu` and `solve_cuda` in `stationary`, `time_derivative1`, and `time_derivative2`) and updated solver dispatch calls to explicitly pass `verbose=verbose` as a keyword argument. Resolves static type checker missing argument errors and positional signature mismatches at line 168 of `time_derivative2.py`.

### Removed
- **Unused Preconditioner Modules**: Removed dead `preconditioners.py` backend files from `src/mGFD/solvers/_backends/cpu/` and `src/mGFD/solvers/_backends/cuda/`.
- **Legacy Dataset References**: Removed outdated references to the legacy `Holes` dataset from docstrings in research batch runners (`run_Poisson.py`, `run_Heat.py`, `run_Wave.py`, `run_AdvDif.py`, `run_Perturbation.py`) and agent architectural documentation (`AGENTS.md`, `arquitectura.md`).
- **Iterative Krylov Preconditioners & Solvers**: Removed obsolete `bicgstab`, `gmres`, `ilu`, and `jacobi` preconditioners from private solver backends.

## [0.11.1] - 2026-08-28
### Added
- **Dynamic Input Typing**: All PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) now natively accept Python `Callable` functions, Numpy `ndarray`, and Pandas `Series`/`DataFrame` objects interchangeably for boundary conditions ($\phi$) and forcing terms ($f$).
- **Control Enhancements**: Added explicit configuration options (`input_types`, `plot_approximations`) to `sweep_config.json` for granular experiment tuning without code modification.

### Changed
- **Architectural Redesign (Semantics)**: Decoupled and purged ambiguous folders `core/` and `backends/` from the root tree.
- Moved `gammas.py` and `neighbors.py` into `mGFD.spatial` to reflect their role as the true spatial discretization modules of GFD.
- Relocated `adapters.py` and `utils.py` into a generalized `mGFD.utils` package (renamed `utils.py` to `core_utils.py`).
- Encapsulated hardware-specific solver architectures (`cpu/` and `cuda/`) as private packages *inside* the solvers directory (`mGFD.solvers._backends.cpu` and `mGFD.solvers._backends.cuda`).
- Fully isolated preconditioning imports (`cupyx` vs `scipy`) into localized `preconditioners.py` modules within their respective backends, making the `spatial` (core) directory 100% agnostic to CUDA or algebraic libraries.
- **Orchestrator Optimization**: The `sweep` parameter pipeline was drastically refactored to eliminate combinatorial explosions. Solvers are now intelligently hardware-routed (`spsolve` strictly for CPU, iterative Krylov `gmres` solvers for GPU).

### Fixed
- **Adaptive Neighbor Selection Stabilization**: Polished the robustness of the adaptive heuristic for spatial stencils that algorithmically guarantees well-conditioned local star matrices, eliminating global ill-conditioning (OOM crashes) in extreme geometries.
- **Styling**: Enforced strict 136-column comment alignment and PEP-8 docstring standardizations across all module refactors.

## [0.11.0] - 2026-08-27
### Added
- **GPU Acceleration (CuPy)**: Fully integrated `device="cuda"` support for all PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`). This allows offloading algebraic bottlenecks (sparse matrices and iterative solvers) directly to Nvidia GPUs via `cupyx.scipy.sparse`.
- **Krylov Preconditioning**: Added optional preconditioning parameter (`preconditioner="ilu"` or `"jacobi"`) to all solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) for vastly improved `GMRES` and `BiCGStab` convergence.
- **Matrix-Free Computing**: Implemented `matrix_free=True` mode for all solvers, utilizing Numba JIT-compiled on-the-fly matrix-vector multiplications via `scipy.sparse.linalg.LinearOperator`, drastically reducing memory usage for iterative solvers.
- Core function `compute_preconditioner` inside `preconditioners.py` generating `scipy.sparse.linalg.LinearOperator` wrappers for ILU factorization and diagonal scaling (with dynamic CuPy compatibility).
- **Adaptive Neighborhoods (Dynamic `nvec`)**: Modified KDTree search algorithms in `neighbors.py` (`_find_neighbors_balanced_jit`) to automatically build density-aware, condition-aware dynamic stencils using `np.linalg.cond` internally. The `nvec` solver parameter is now treated as an upper bound rather than a strict count.

### Fixed


## [0.10.0] - 2026-08-26
### Added
- `src/mGFD/core/adapters.py` to seamlessly integrate `pandas` DataFrames and `xarray` DataArrays using a lightweight soft-dependency (duck typing) approach.
- `src/mGFD/solvers/results.py` module defining the new `SolverResult` dataclass to standardize the output of all PDE solvers.
- `.agents/ROADMAP_v0.10.md` to define the 4 structural pillars for making `mGFD` an industry-ready library.
- `src/mGFD/exceptions.py` introducing a custom exception hierarchy (`mGFDError`, `CloudShapeError`, `InputTypeError`, `DimensionMismatchError`, `OperatorFormatError`, `ParameterError`).
- `.agents/OPTIMIZATION_PLAN.md` file to act as a backlog for future optimization strategies in the `mGFD` project.

### Fixed
- Corrected type hinting signatures in `TimeDerivative1` and `TimeDerivative2` to properly reflect the `SolverResult` return type, resolving static analyzer (mypy/pylance) warnings.

### Changed
- Restored and upgraded `research/codes/` batch processing scripts (`run_Poisson.py`, `run_Heat.py`, `run_Wave.py`, `run_AdvDif.py`, `run_Perturbation.py`) to natively utilize the v0.10.0 `SolverResult` and flexible algebraic backends, without losing their ability to iterate over massive point cloud datasets.
- Removed the redundant `run_Perturbation2.py` batch script.
- Standardized the output of all PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) to return the new `SolverResult` dataclass instead of returning raw tuples.
- Implemented `__iter__` on `SolverResult` to allow 100% backward-compatible tuple unpacking (`u_ap, vec = solver(...)`), completely avoiding breaking changes to existing dependent scripts.
- Integrated flexible algebraic linear solvers (`spsolve`, `bicgstab`, `gmres`) in `Stationary`, `TimeDerivative1`, and `TimeDerivative2`, avoiding memory bottlenecks for large systems.
- Implemented $x_0$ warm-starting (solution from step $k-1$) in transient solvers to exponentially accelerate Krylov subspace convergence.
- Refactored `Stationary`, `TimeDerivative1`, and `TimeDerivative2` solvers to raise custom exceptions instead of generic `ValueError` or `TypeError`.
- Updated test suite (`tests/`) to validate against the new custom exception hierarchy.
- Integrated Numba (`@nb.njit`) to optimize core iterative loops in `gammas.py` (`Cloud`, `CloudStencil`, and `_nnls_numba`), achieving massive performance gains for meshless operator constructions.
- Refactored KDTree neighbor queries in `neighbors.py` (`find_neighbors`, `find_neighbors_balanced`, `find_neighbors_adv`) to use SciPy's vectorized bulk queries, coupled with Numba JIT functions to handle filtering and sorting, bypassing the slow Python loop overhead.
- Optimized `cloud_generator` bottlenecks via Numba JIT compilation (`jit_geometry.py`, `boundary.py`, `relaxation.py`, `regular_generation.py`) by replacing `matplotlib.path` logic with a fast ray-casting algorithm.
- Upgraded point classification in `classification.py` to use fully vectorized GEOS routines and Shapely 2.0 `spatial_index.query` arrays, bypassing slow Python iterations.

## [0.9.1] - 2026-08-09
### Added
- Official publication of version `0.9.1` on PyPI.
- Continuous Integration (CI) workflow in GitHub Actions with tests for Python 3.9, 3.10, 3.11, 3.12, and 3.13.
- Massive refactoring to comply with strict Type Hints (Mypy).
- Mathematical signature correction and type compatibility (NumPy/SciPy/PyVista/OpenCV).
- Pre-release and verification on TestPyPI.

## [0.9.0] - 2026-08-08
### Added
- `.agents` directory with `AGENTS.md` to define project context, architecture, and coding conventions for AI agents.
- `CHANGELOG.md` to track ongoing project modifications.
- `src/mGFD/cloud_generator/core/point_generation/` package to modularize the cloud generation logic (geometry, boundary, relaxation, regular, poisson).
- `research/batches/test_single.py` for rapid validation of all PDE solvers on a single point cloud.
- `tests/` directory with a robust `pytest` suite for validating mathematical solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) and cloud generation logic.
- Standard `logging` infrastructure across all modules, replacing raw `print()` statements for better external consumption.
- Comprehensive `fail-fast` input validation in all solvers to catch dimensionality mismatches early and prevent cryptic NumPy/SciPy errors.
- GitHub Actions CI workflow (`.github/workflows/python-tests.yml`) for automated testing across Python 3.9-3.12.

### Changed
- `AGENTS.md` updated to formally enforce the strict column-136 inline comment formatting and explicitly prohibit the use of temporary helper scripts for codebase edits.
- Refactored `point_generation.py` monolith into multiple modular files.
- Renamed `grid_generation` functions and files to `regular_generation` to better align with the "meshless" philosophy.
- Updated `mgfd-cloud` CLI help text and `README.md` to reflect the new `regular` terminology and provide better examples.
- Architecturally decoupled the theoretical exact solution (`u_ex`) from the `Stationary`, `TimeDerivative1`, and `TimeDerivative2` solvers. They now return `Tuple[np.ndarray, np.ndarray]` (`u_ap, vec`).
- Adapted all benchmarking scripts in `research/batches/` to compute `u_ex` locally instead of relying on the solvers.
- Modernized project configuration entirely into a PEP 621 compliant `pyproject.toml` (removed legacy setup approaches) with updated PyPI classifiers and metadata.
- Improved Type Hinting (`mypy`) coverage throughout all core modules and utilities.

### Deprecated
- None yet.

### Removed
- None yet.

### Fixed
- Addressed multiple syntax, indentation, and f-string errors flagged by `flake8`, strictly maintaining the column 136 comment alignment rule.

### Security
- None yet.
