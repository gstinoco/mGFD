# mGFD: Meshless Generalized Finite Differences 📐☁️

<div align="center">

<img src="docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0; border-radius: 8px;">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg?style=for-the-badge)](https://github.com/gstinoco/mGFD) 
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/downloads/) 
[![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-013243.svg?style=for-the-badge&logo=numpy)](https://numpy.org/) 
[![SciPy](https://img.shields.io/badge/SciPy-Numerical%20Algorithms-8CAAE6.svg?style=for-the-badge&logo=scipy&logoColor=white)](https://scipy.org/) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

**A high-performance Python ecosystem for Point Cloud Generation and PDE solving on highly irregular domains using Generalized Finite Differences.**

</div>

---

## 🌟 Overview

**mGFD** is a complete meshless computational suite. Traditional methods (like Finite Elements or Finite Volumes) require complex and restrictive mesh generation. **mGFD** completely bypasses this limitation by solving Partial Differential Equations (PDEs) directly on **unstructured point clouds**. 

This makes it exceptionally powerful for modeling physics in complex, real-world geometries such as natural lakes, islands, or custom engineering domains.

### 🚀 Key Features

* **☁️ Integrated Cloud Generator:** A powerful engine to automatically generate 2D point clouds from geographic contour boundaries (supports Poisson-Disk sampling, Lloyd relaxation, and grid-based methods).
* **📐 Pure Meshless Solvers:** Discretize and solve PDEs using only local neighbor stencils. No mesh required!
* **⚡ Stationary Solvers:** Out-of-the-box support for Poisson-type equations and generalized stationary problems.
* **🔥 Transient Solvers:** First-order (Heat, Advection-Diffusion) and Second-order (Wave) time integrations.
* **🧩 Modular Architecture:** Highly professional PEP-8 compliant sub-package architecture (`mGFD.solvers`, `mGFD.cloud_generator`, `mGFD.core`).

---

## 📦 Installation

**mGFD** relies on a robust scientific stack (`numpy`, `scipy`, `shapely`, `opencv-python-headless`).

To install the package from source (useful for development):

```bash
git clone https://github.com/gstinoco/mGFD.git
cd mGFD
pip install -e .
```

This will automatically install both the Python API and the `mgfd-cloud` command-line interface.

---

## 🛠️ Quick Start

The true power of `mGFD` lies in its ability to go from raw geographic contours to a solved PDE in just two steps.

### Step 1: Generate a Point Cloud

Use the built-in CLI tool to generate a point cloud from a CSV containing boundary coordinates:

```bash
# Generate a cloud using the Natural (Poisson-Disk) method
mgfd-cloud generate --input contours.csv --output my_cloud.csv --method natural --density 0.1
```

### Step 2: Solve a PDE

Use the Python API to load your cloud and run a meshless solver. Here is a minimal example for a **Stationary Equation**:

```python
import numpy as np
from mGFD.io.io import load_points
from mGFD import Stationary

# 1. Load the generated point cloud
p = load_points("my_cloud.csv")

# 2. Define the analytical boundary condition
phi = lambda x, y: np.exp(x + y)

# 3. Define the right-hand side forcing function (e.g. Laplacian of phi)
f_stat = lambda x, y: 2 * np.exp(x + y)

# 4. Define the differential operator [D, E, A, B, C, F]
# This corresponds to: 0*u_x + 0*u_y + 2*u_xx + 0*u_xy + 2*u_yy
L_stat = np.vstack([[0], [0], [2], [0], [2], [0]])

# 5. Solve the equation!
u_ap, u_ex, vec = Stationary(p, phi, f_stat, operator=L_stat)
```

The root `mGFD` module exposes the primary solvers for easy access:
* `Stationary(p, phi, f, ...)`: Time-independent equations.
* `TimeDerivative1(p, f, t, coef, ...)`: 1st-order transient PDEs (e.g., Heat).
* `TimeDerivative2(p, f, g, t, coef, ...)`: 2nd-order transient PDEs (e.g., Wave).

---

## 🔬 Research & Datasets

Looking for the **mathematical formulation**, **real-world geographic lake datasets** (Lake Patzcuaro, Caspian Sea, etc.), or reproducible **benchmarking scripts**?

👉 **[Explore the Research Laboratory (`/research/README.md`)](./research/README.md)**

The `research/` directory contains our complete academic suite, including experimental data, geographic boundary files, VTK results, and the exact theoretical explanations associated with our scientific publications.

---

## 📜 Citation & Credits

This project is open-sourced under the **MIT License**.

If you use this library or the core mathematical methodology in your research, please cite our reference paper:

```bibtex
@article{tinoco2025mgfd,
  title={mGFD: A meshless generalized finite difference method},
  author={Tinoco-Guerrero, Gerardo and Domínguez-Mota, Francisco Javier and Guzmán-Torres, José Alberto and Pedraza-Jiménez, Gabriela and Tinoco-Ruiz, José Gerardo},
  journal={Computers & Mathematics with Applications},
  volume={195},
  pages={396--418},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.camwa.2025.07.034}
}
```

<div align="center">
<br/>
<i>Built with ❤️ by the scientific computing community in Morelia, México.</i>
<br/><br/>
<b>Dr. Gerardo Tinoco-Guerrero</b><br/>
<b>Dr. Francisco Javier Domínguez-Mota</b><br/>
<b>Dr. José Alberto Guzmán-Torres</b><br/>
Universidad Michoacana de San Nicolás de Hidalgo<br/>
gerardo.tinoco@umich.mx
<br/><br/>
<a href="https://github.com/gstinoco/mGFD/issues">Report a Bug</a> | <a href="mailto:gerardo.tinoco@umich.mx">Contact Author</a>
</div>
