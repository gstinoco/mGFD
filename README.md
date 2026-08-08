# mGFD: Meshless Generalized Finite Differences 📐☁️

<div align="center">

<img src="docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0; border-radius: 8px;">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg?style=for-the-badge)](https://github.com/gstinoco/mGFD) 
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/downloads/) 
[![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-013243.svg?style=for-the-badge&logo=numpy)](https://numpy.org/) 
[![SciPy](https://img.shields.io/badge/SciPy-Numerical%20Algorithms-8CAAE6.svg?style=for-the-badge&logo=scipy&logoColor=white)](https://scipy.org/) 
[![Tests](https://img.shields.io/badge/Tests-Passing-brightgreen.svg?style=for-the-badge)](https://github.com/gstinoco/mGFD/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

**A high-performance Python ecosystem for Point Cloud Generation and PDE solving on highly irregular domains using Generalized Finite Differences.**

</div>

---

## 📑 Table of Contents
- [🌟 Overview](#-overview)
- [📦 Installation](#-installation)
- [🛠️ Quick Start](#️-quick-start)
- [🔬 Research & Datasets](#-research--datasets)
- [📜 Citation & Credits](#-citation--credits)

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

The `mgfd-cloud` CLI tool automates the process of converting geometric boundaries into valid computational domains.

<details open>
<summary><b>View available CLI commands</b></summary>

```bash
# General help
mgfd-cloud --help

# Generate a point cloud (Poisson-Disk)
mgfd-cloud generate --input contours.csv --output my_cloud.csv --method natural --density 0.1

# Reduce the density of an existing cloud silently
mgfd-cloud -q reduce --input high_density.csv --output low_density.csv --factor 0.5
```
</details>

### Step 2: Solve a PDE

Use the Python API to load your cloud and run a meshless solver.

<details open>
<summary><b>Stationary Equation (e.g. Poisson)</b></summary>

```python
import numpy as np
from mGFD.io.io import load_points
from mGFD import Stationary

# 1. Load the generated point cloud
p = load_points("my_cloud.csv")

# 2. Define the analytical boundary condition
phi = lambda x, y: np.exp(x + y)

# 3. Define the right-hand side forcing function
f_stat = lambda x, y: 2 * np.exp(x + y)

# 4. Define the differential operator [D, E, A, B, C, F]
L_stat = np.vstack([[0], [0], [2], [0], [2], [0]])

# 5. Solve the equation! (Verbose mode on by default in scripts, but off in the core library)
u_ap, u_ex, vec = Stationary(p, phi, f_stat, operator=L_stat, verbose=True)
```
</details>

<details>
<summary><b>First-Order Transient Equation (e.g. Heat)</b></summary>

```python
import numpy as np
from mGFD.io.io import load_points
from mGFD import TimeDerivative1

p = load_points("my_cloud.csv")
v, t1 = 0.01, 100

f_heat = lambda x, y, t, coef: np.exp(-2 * np.pi**2 * coef[0] * t) * np.cos(np.pi * x) * np.cos(np.pi * y)
L_heat = np.vstack([[0], [0], [2*v], [0], [2*v], [0]])

u_ap, u_ex, vec = TimeDerivative1(p, f_heat, t1, [v], operator=L_heat, implicit=True, lam=0.5, verbose=True)
```
</details>

<details>
<summary><b>Second-Order Transient Equation (e.g. Wave)</b></summary>

```python
import numpy as np
from mGFD.io.io import load_points
from mGFD import TimeDerivative2

p = load_points("my_cloud.csv")
c, t2 = 0.5, 50

f_wave = lambda x, y, t, coef: np.cos(np.sqrt(2) * np.pi * coef[0] * t) * np.sin(np.pi * x) * np.sin(np.pi * y)
g_wave = lambda x, y, t, coef: 0.0 * x  # Initial velocity
L_wave = np.vstack([[0], [0], [2*c**2], [0], [2*c**2], [0]])

u_ap, u_ex, vec = TimeDerivative2(p, f_wave, g_wave, t2, [c], operator=L_wave, implicit=True, lam=0.5, verbose=True)
```
</details>

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
<i>Developed for the advancement of meshless numerical methods and scientific computing.</i>
<br/><br/>
<b>Dr. Gerardo Tinoco-Guerrero</b><br/>
<b>Dr. Francisco Javier Domínguez-Mota</b><br/>
<b>Dr. José Alberto Guzmán-Torres</b><br/>
Universidad Michoacana de San Nicolás de Hidalgo<br/>
gerardo.tinoco@umich.mx
<br/><br/>
<a href="https://github.com/gstinoco/mGFD/issues">Report a Bug</a> | <a href="mailto:gerardo.tinoco@umich.mx">Contact Author</a>
</div>
