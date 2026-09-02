# 📈 mGFD Library Changelog

All notable changes to the **mGFD** core Python library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [0.11.1] - 2026-09-02
### Added
- **Temporal Discretization Package (`mGFD.temporal`)**: Created dedicated `mGFD.temporal` package containing `estimate_cfl_dt` for calculating stable time steps $\Delta t$ and Courant CFL numbers automatically for 1st-order (parabolic/advective) and 2nd-order (hyperbolic) transient PDEs using characteristic average node spacing ($h_{avg}$) and propagation speed scales ($V_{adv} + \sqrt{\nu}$). `TimeDerivative1` and `TimeDerivative2` support automatic step count calculation (`t=None`) using a target `cfl` parameter (default 0.5) or custom `dt`.
- **Flexible Physical Time Domains (`t_span`)**: Added `t_span: Tuple[float, float] = (0.0, 1.0)` parameter to `TimeDerivative1` and `TimeDerivative2` allowing users to solve transient PDEs over any arbitrary physical time interval $[t_{start}, t_{end}]$ (e.g. `t_span=(0.0, 100.0)` or `t_span=(0.0, 86400.0)` seconds). Updated visualization functions `plot_transient` and `plot_transient_steps` in `mGFD.viz.graph` as well as `export_transient_vtk` in `mGFD.io.export_vtk` to automatically scale animation titles and ParaView PVD time steps to the physical domain $[t_{start}, t_{end}]$.
- **Independent Initial (`ic`) and Boundary (`bc`) Conditions**: Added optional `ic` (initial condition $u_0(x, y)$ at $t_{start}$) and `bc` (boundary condition $\phi(x, y, t)$ across boundary nodes) parameters to `TimeDerivative1` and `TimeDerivative2`. Falls back to `f` when omitted, guaranteeing 100% backward compatibility.
- **Independent Source / Forcing Term (`source`)**: Added support for non-homogeneous source terms $F_{source}(x, y, t)$ in `TimeDerivative1` ($\partial u/\partial t = L u + F_{source}$) and `TimeDerivative2` ($\partial^2 u/\partial t^2 = L u + F_{source}$). Supports scalars, 1D spatial arrays, 2D spatiotemporal matrices, and callables across both CPU and CUDA backends.
- **CFL Attributes in SolverResult**: `SolverResult` now encapsulates `dt`, `cfl`, and `t_steps` attributes while preserving 100% backward-compatible tuple unpacking (`u_ap, vec = result`).
- **High-Speed Rendering Dependency (`imageio-ffmpeg`)**: Integrated auto-detection of pre-compiled multi-threaded FFmpeg binaries in `mGFD.viz.graph` (`plot_transient`), reducing GIF compilation time per file from **4.19s down to 0.69s (~6x speedup)**.

### Changed
- **High-Performance CUDA Backend Optimization (~90x Speedup)**: Replaced CuPy's legacy `cp_factorized` wrapper with GPU-native `bicgstab` iterative sparse solves with warm-start initial state propagation (`x0 = u_ap_gpu[:, k-1]`) in `_backends/cuda/time_derivative1.py` and `time_derivative2.py`. Reduced CUDA solver execution time from **29.21s down to ~0.33s** for 500 time steps on NVIDIA GPUs.
- **High-Performance CPU Solver Optimization**: Accelerated callable evaluation (`bc`, `ic`, `g`, `source`) across CPU and CUDA solver backends (`time_derivative1.py`, `time_derivative2.py`) using vectorized spatial broadcasting and 2-arg/4-arg lambda fallback handlers, reducing boundary condition evaluation overhead by **~40x**. Re-used pre-allocated `RHS` memory buffers in implicit CPU time-stepping loops (`solve_cpu`) and eliminated zero-fill memory allocations in Numba stencil assembly (`gammas.py`).
- **Core Library Clean-Up**: Audited and purged all unused imports (`nnls`, `spilu`, `Point`, `STRtree`, `cm`, `npt`, `Callable`, `ExportVTK`, `Optional`, `iter_clouds`, `json`, `time`, `List`, `Any`) across 12 core `mGFD` modules (`stationary.py`, `gammas.py`, `neighbors.py`, `export_vtk.py`, `graph.py`, `core_utils.py`, `classification.py`, `poisson_generation.py`, `regular_generation.py`, `relaxation.py`, `visualization.py`, `cli.py`), keeping the core package 100% clean and zero-warning compliant.

### Fixed
- **CUDA Backend Missing `boun_idx` Definition**: Defined CPU-side integer index arrays `boun_idx` and `inne_idx` in `_backends/cuda/time_derivative1.py` and `_backends/cuda/time_derivative2.py` prior to boundary condition evaluation. Fixes `NameError: name 'boun_idx' is not defined` when evaluating 1D boundary functions on CUDA solvers.
- **Pyright / Mypy Type Annotations & Overload Resolutions**: Standardized type hints for `p`, `u`, `data`, and `triangles` across core visualization and utility modules. Resolves static type checker overload resolution errors.

---

## [0.11.0] - 2026-08-27
### Added
- **GPU Acceleration (CuPy)**: Fully integrated `device="cuda"` support for all PDE solvers (`Stationary`, `TimeDerivative1`, `TimeDerivative2`) via `cupyx.scipy.sparse`.
- **Krylov Preconditioning**: Added optional preconditioning parameter (`preconditioner="ilu"` or `"jacobi"`) to all solvers for improved convergence.
- **Matrix-Free Computing**: Implemented `matrix_free=True` mode utilizing Numba JIT-compiled on-the-fly matrix-vector multiplications via `scipy.sparse.linalg.LinearOperator`.
- **Adaptive Neighborhoods (Dynamic `nvec`)**: KDTree search algorithms in `neighbors.py` automatically build density-aware, condition-aware dynamic stencils using `np.linalg.cond`.

---

## [0.10.0] - 2026-08-26
### Added
- `mGFD.exceptions` custom exception hierarchy (`mGFDError`, `CloudShapeError`, `InputTypeError`, `DimensionMismatchError`, `OperatorFormatError`, `ParameterError`).
- `mGFD.solvers.results` module defining `SolverResult` dataclass to standardize PDE solver outputs.
- `mGFD.utils.adapters` for native Pandas DataFrames and Xarray DataArrays support.
