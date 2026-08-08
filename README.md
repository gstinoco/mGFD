# mGFD: Meshless Generalized Finite Differences :triangular_ruler:

<div align="center">

<img src="docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0;">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg)](https://github.com/gstinoco/mGFD) [![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/) [![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-013243.svg?logo=numpy)](https://numpy.org/) [![SciPy](https://img.shields.io/badge/SciPy-Numerical%20Algorithms-8CAAE6.svg?logo=scipy)](https://scipy.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A high-performance Python package for solving 2D PDEs on highly irregular domains using Generalized Finite Differences.**

</div>

## :star2: Overview

**mGFD** is a meshless computational engine designed to solve Partial Differential Equations (PDEs) directly on unstructured point clouds. Unlike traditional Finite Element or Finite Volume methods, mGFD requires no mesh generation, making it exceptionally powerful for complex, irregular geometries like geographic lakes, islands, or custom engineering domains.

### Features
- **:triangular_ruler: Pure Meshless Discretization:** Solves PDEs using only local neighbor stencils.
- **:zap: Stationary Solvers:** Out-of-the-box support for Poisson-type problems.
- **:fire: Transient Solvers:** First-order (Heat, Advection-Diffusion) and Second-order (Wave) time integrations.
- **:gear: Extensible API:** Easily define custom operators, forcing functions, and boundary conditions using standard Python callbacks.

---

## :package: Installation

mGFD requires Python 3.8+ and standard scientific libraries (`numpy`, `scipy`).

To install the package directly from GitHub:

```bash
pip install git+https://github.com/gstinoco/mGFD.git
```

For development or local testing, clone the repository and install it in editable mode:

```bash
git clone https://github.com/gstinoco/mGFD.git
cd mGFD
pip install -e .
```

---

## :rocket: Quick Start API

mGFD makes solving complex PDEs simple. Here is a minimal example of solving a **Stationary Equation** over a custom point cloud:

```python
import numpy as np
from mGFD import Stationary

# 1. Provide a point cloud `p` of shape (m, 3). 
# Columns are [x, y, flag], where flag=0 is interior, flag=0.5 is boundary.
p = np.array([
    [0.5, 0.5, 0.0],  # Interior point
    [0.0, 0.5, 0.5],  # Boundary point
    # ... your points here
])

# 2. Define the analytical boundary condition
phi = lambda x, y: 2 * np.exp(2*x + y)

# 3. Define the right-hand side forcing function
f = lambda x, y: 10 * np.exp(2*x + y)

# 4. Solve the equation!
u_ap, u_ex, vec = Stationary(p, phi, f)
```

The `mGFD` module exports three main solvers:
- `Stationary(p, phi, f)`: For time-independent equations.
- `TimeDerivative1(p, f, t, coef)`: For 1st-order transient PDEs (e.g. Heat).
- `TimeDerivative2(p, f, g, t, coef)`: For 2nd-order transient PDEs (e.g. Wave).

---

## :microscope: Research, Datasets & Benchmarks

Are you looking for the **mathematical formulation**, the **lake point-cloud datasets** (Caspio, Poopo, etc.), or the reproducible **benchmarking scripts**?

👉 **[Go to the Research Laboratory (`/research/README.md`)](./research/README.md)**

The research folder contains the complete academic suite, including historical experiment data, visualization scripts, VTK results, and theoretical explanations associated with our scientific publications.

---

## :memo: Citation & License

This project is open-sourced under the **MIT License**.

If you use this library or the core mathematical methodology in your research, please cite our reference paper:

```bibtex
@software{tinoco2024mGFD,
  title={mGFD: Meshless Generalized Finite Differences for Partial Differential Equations},
  author={Tinoco-Guerrero, Gerardo and Domínguez-Mota, Francisco Javier and Guzmán-Torres, José Alberto and Arias-Rojas, Heriberto},
  year={2024},
  institution={Universidad Michoacana de San Nicolás de Hidalgo},
  url={https://github.com/gstinoco/mGFD}
}
```

<div align="center">
<i>Built with ❤️ by the scientific computing community in Morelia, México.</i><br/>
<a href="https://github.com/gstinoco/mGFD/issues">Report a Bug</a> | <a href="mailto:gerardo.tinoco@umich.mx">Contact Author</a>
</div>
