# mGFD: Meshless Generalized Finite Differences 📐☁️

<div align="center">

<img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0; border-radius: 8px;">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg?style=for-the-badge)](https://github.com/gstinoco/mGFD) 
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/downloads/) 
[![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-013243.svg?style=for-the-badge&logo=numpy)](https://numpy.org/) 
[![SciPy](https://img.shields.io/badge/SciPy-Numerical%20Algorithms-8CAAE6.svg?style=for-the-badge&logo=scipy&logoColor=white)](https://scipy.org/) 
[![CI Tests](https://img.shields.io/github/actions/workflow/status/gstinoco/mGFD/python-tests.yml?style=for-the-badge&label=CI%20Tests)](https://github.com/gstinoco/mGFD/actions/workflows/python-tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

**A high-performance Python ecosystem for Point Cloud Generation and PDE solving on highly irregular domains using Generalized Finite Differences.**

</div>

---

## 📑 Table of Contents
- [🌟 Overview & The Mathematical Core](#-overview--the-mathematical-core)
- [📦 Installation](#-installation)
- [🚀 End-to-End Workflow: Lake Pátzcuaro](#-end-to-end-workflow-lake-pátzcuaro)
- [📚 Python API Reference](#-python-api-reference)
- [🔬 Research & Datasets](#-research--datasets)
- [🏗️ Repository Structure](#️-repository-structure)
- [🤝 Contributing & Citation](#-contributing--citation)

---

## 🌟 Overview & The Mathematical Core

**mGFD** is a complete meshless computational suite. Traditional methods (like Finite Elements or Finite Volumes) require complex and restrictive mesh generation. **mGFD** completely bypasses this limitation by solving Partial Differential Equations (PDEs) directly on **unstructured point clouds**. 

This makes it exceptionally powerful for modeling physics in complex, real-world geometries such as natural lakes, islands, or custom engineering domains.

### 🧮 The Mathematical Core (El Operador Diferencial)

One of the most powerful features of **mGFD** is its ability to solve *almost any* second-order linear Partial Differential Equation by simply defining a 6-element mathematical vector. 

Based on the theoretical formulation (Tinoco-Guerrero et al., 2025), a general second-order linear PDE operator ($L u$) can be expressed as:

$$ L u = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u $$

Using this exact nomenclature, `mGFD` can solve three main families of problems out-of-the-box. You just construct a 6-element Python array `[D, E, A, B, C, F]` and pass it to the solver:

1. **Stationary Problems** ($L u = f(x, y)$)
   * **Poisson Equation**: $u_{xx} + u_{yy} = f(x, y)$
   * In this case, $A = 1$ and $C = 1$. The operator is:
     ```python
     # [D, E, 2*A, B, 2*C, F]
     L_poisson = np.vstack([[0], [0], [2], [0], [2], [0]])
     ```

2. **First-Order Transients** ($L u = u_t$)
   * **Advection-Diffusion Equation**: $\alpha_1 u_x + \alpha_2 u_y - \nu(u_{xx} + u_{yy}) = u_t$
   * Here $D = \alpha_1$, $E = \alpha_2$, $A = -\nu$, and $C = -\nu$. The operator is:
     ```python
     # [D, E, 2*A, B, 2*C, F]
     L_adv_diff = np.vstack([[a1], [a2], [-2*nu], [0], [-2*nu], [0]])
     ```

3. **Second-Order Transients** ($L u = u_{tt}$)
   * **Wave Equation**: $c^2 u_{xx} + c^2 u_{yy} = u_{tt}$
   * Here $A = c^2$ and $C = c^2$. The operator is:
     ```python
     # [D, E, 2*A, B, 2*C, F]
     L_wave = np.vstack([[0], [0], [2*c**2], [0], [2*c**2], [0]])
     ```

### 🚀 Key Features

* **☁️ Integrated Cloud Generator:** A powerful engine to automatically generate 2D point clouds from geographic contour boundaries.
* **📐 Pure Meshless Solvers:** Discretize and solve PDEs using only local neighbor stencils. No mesh required!
* **🏎️ High Performance:** Fully optimized with Numba JIT compilation and parallelized with LLVM multithreading and KDTree C++ backends.
* **🧩 Modular Architecture:** Highly professional PEP-8 compliant sub-package architecture.

---

## 📦 Installation

**mGFD** relies on a robust scientific stack (`numpy`, `scipy`, `shapely`, `opencv-python-headless`).

The easiest way to install the package is directly from PyPI:

```bash
pip install mGFD
```

To install the package from source (useful for development):

```bash
git clone https://github.com/gstinoco/mGFD.git
cd mGFD
pip install -e .
```

This will automatically install both the Python API and the `mgfd-cloud` command-line interface.

---

## 🚀 End-to-End Workflow: Lake Pátzcuaro

The true power of `mGFD` lies in its ability to go from a raw visual image to a solved PDE in a single, cohesive workflow. Here is a complete real-world scenario using **Lake Pátzcuaro**.

### Step 1: Visual Contour Extraction
Instead of relying on idealized shapes (like circles or squares), **mGFD** excels at real geographic domains. We start by identifying our domain from a satellite map.

Before generating the computational point cloud, we must extract the geometric contours of the lake into a CSV file (`Patzcuaro_contours.csv`). 

> [!TIP]
> **Visual Contour Extraction**
> The easiest way to generate this contour file is by uploading your image to our interactive [mGFD CloudGenerator Web Tool](https://malla.umich.mx/CloudGenerator/). You can manually trace the boundaries and export them directly to CSV!

<div align="center">
  <table>
    <tr>
      <td align="center"><b>Satellite Map</b></td>
      <td align="center"><b>Extracted Contours (CSV Plot)</b></td>
    </tr>
    <tr>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Patzcuaro.png" alt="Lake Patzcuaro Satellite" width="350" style="border-radius: 8px;"></td>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Patzcuaro_contours.png" alt="Lake Patzcuaro Contours" width="350" style="border-radius: 8px;"></td>
    </tr>
  </table>
</div>

### Step 2: Generate the Cloud (CLI)
Using the `mgfd-cloud` CLI, we fill the interior of these extracted contours with a computationally valid, Poisson-Disk distributed point cloud:

```bash
mgfd-cloud generate --input Patzcuaro_contours.csv --output Patzcuaro_cloud.csv --method natural --density 3
```

<div align="center">
  <img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Patzcuaro_cloud.png" alt="Lake Patzcuaro Point Cloud" width="350" style="border-radius: 8px;">
</div>

### Step 3: Python Solution (API)
Create a script named `solve_patzcuaro.py` to load your new cloud, apply the meshless solver for the Poisson equation, and render the results using the built-in visualization utilities:

```python
import numpy as np
from mGFD.io.io import load_points
from mGFD import Stationary
from mGFD.viz.graph import plot_stationary

# 1. Load the generated point cloud
p = load_points("Patzcuaro_cloud.csv")

# 2. Define the analytical boundary condition
phi = lambda x, y: np.exp(x + y)

# 3. Define the forcing function (RHS)
f_stat = lambda x, y: 2 * np.exp(x + y)

# 4. Define the differential operator for Poisson (Laplacian)
# Using [D, E, 2*A, B, 2*C, F]
L_stat = np.vstack([[0], [0], [2], [0], [2], [0]])

# 5. Execute the solver
print("Solving Poisson Equation for Lake Pátzcuaro...")
u_ap, vec = Stationary(p, phi, f_stat, operator=L_stat, verbose=True)

# 6. Visualize the numerical approximation
print("Rendering plot...")
plot_stationary(p, u_ap, save=False, title='Pátzcuaro Poisson Solution', verbose=True)
```

Run your script, and a Matplotlib 3D window will pop up showing the surface of your solution perfectly adapted to the irregular contours of the lake!

### Step 4: Visualize Results

When running batch solvers or keeping `save=True`, **mGFD** automatically exports high-quality PNGs (for stationary problems) and animated GIFs (for transient problems).

Here are the actual simulation outputs for Lake Pátzcuaro using the exact cloud generated above:

<div align="center">
  <table>
    <tr>
      <td align="center"><b>Poisson (Stationary)</b></td>
      <td align="center"><b>Heat (Transient)</b></td>
    </tr>
    <tr>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Poisson_Approximation.png" alt="Patzcuaro Poisson" width="400"></td>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Heat_Approximation.gif" alt="Patzcuaro Heat" width="400"></td>
    </tr>
    <tr>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Poisson_Approximation_top.png" alt="Patzcuaro Poisson" width="400"></td>
      <td><img src="https://raw.githubusercontent.com/gstinoco/mGFD/main/docs/Patzcuaro/Heat_Approximation_top.gif" alt="Patzcuaro Heat" width="400"></td>
    </tr>
  </table>
</div>

---

## 📚 Python API Reference

**mGFD** provides a clean and modular Python API located in the `mGFD` module. Below is an exhaustive manual of the core functions available for point cloud generation, I/O operations, PDE solving, and visualization.

### 1. Core Solvers (`mGFD.solvers`)

These functions are the mathematical heart of **mGFD**. They compute spatial derivatives using generalized finite differences over local neighborhoods.

#### `Stationary`
Solves stationary problems (e.g., Poisson equation) with Dirichlet boundary conditions.
```python
def Stationary(
    p: np.ndarray, 
    phi: Callable, 
    f: Callable, 
    operator: np.ndarray, 
    upwind: bool = False, 
    vec: Optional[np.ndarray] = None, 
    nvec: int = 12, 
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
```
*   `p` *(m x 3)*: Point cloud array with columns `[x, y, flag]`. Interior points have `flag=0`, boundary points have `flag=1` or `flag=2`.
*   `phi`: A Python lambda/function `phi(x, y)` representing the analytical Dirichlet boundary condition.
*   `f`: A Python lambda/function `f(x, y)` representing the right-hand side (forcing term).
*   `operator` *(6 x 1)*: Differential operator vector `[D, E, A, B, C, F]`.
*   `upwind`: Set to `True` for advection-dominated problems to automatically bias the neighbor search against the velocity field.
*   `vec`: Precomputed neighbor list. If `None`, it is computed automatically.
*   `nvec`: The target number of neighbors to find for each point (default: 12). If instability is detected, the solver may dynamically retry with expanded neighborhoods (up to 30).
*   **Returns**: `(u_ap, vec)`, where `u_ap` is the computed numerical solution and `vec` is the neighborhood list.

#### `TimeDerivative1`
Solves first-order-in-time problems (e.g., Heat and Advection-Diffusion equations).
```python
def TimeDerivative1(
    p: np.ndarray, 
    f: Callable, 
    t: int, 
    coef: List[float], 
    operator: np.ndarray, 
    implicit: bool = False, 
    lam: float = 0.5, 
    upwind: bool = False, 
    vec: Optional[np.ndarray] = None, 
    nvec: int = 12, 
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
```
*   `f`: A time-dependent function `f(x, y, T, coef)` representing the exact solution (used for boundary conditions during time steps).
*   `t`: The number of time steps. Transients use a normalized time grid `T = np.linspace(0, 1, t)`.
*   `coef`: A list of physical coefficients passed to `f` (e.g., `[viscosity]`).
*   `implicit`: If `True`, uses a fully implicit matrix solve (SciPy `bicgstab` is required). Highly recommended for stability.
*   `lam`: The interpolation factor $\lambda$ for time-stepping ($\lambda = 0.5$ implies a Crank-Nicolson scheme).

#### `TimeDerivative2`
Solves second-order-in-time problems (e.g., Wave equation).
```python
def TimeDerivative2(
    p: np.ndarray, 
    f: Callable, 
    g: Callable, 
    t: int, 
    coef: List[float], 
    operator: np.ndarray, 
    implicit: bool = False, 
    lam: float = 0.5, 
    upwind: bool = False, 
    vec: Optional[np.ndarray] = None, 
    nvec: int = 12, 
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
```
*   `g`: A function `g(x, y, T, coef)` representing the initial velocity of the system.

### 2. Cloud Generation (`mGFD.cloud_generator`)


#### Command-Line Interface (`mgfd-cloud`)
The `mGFD` installation includes a powerful CLI to rapidly generate and reduce computational domains.

```bash
mgfd-cloud generate -i INPUT -o OUTPUT -m {natural,regular} [-d DENSITY] [--inside-regions]

mgfd-cloud reduce -i INPUT -o OUTPUT -m MULTIPLIER
```
*   `generate`: Fills a boundary contour with points.
    *   `-i` / `--input`: Path to the boundary CSV file.
    *   `-o` / `--output`: Path for the generated cloud CSV.
    *   `-m` / `--method`: Use `natural` (Poisson-Disk) or `regular` (Grid).
    *   `-d` / `--density`: Density multiplier. Higher values create denser, tighter clouds.
    *   `--inside-regions`: If present, internal closed boundaries are treated as solid islands (holes in the domain).
*   `reduce`: Reduces the density of an existing cloud.
    *   `-m` / `--multiplier`: Scaling factor to space out points.

While the CLI (`mgfd-cloud`) is the easiest way to generate domains, you can invoke the generation directly via Python.

#### `generate_cloud_natural` / `generate_cloud_regular`
Generates an unstructured (Poisson-Disk) or regular (Grid) point cloud inside a contour.
```python
def generate_cloud_natural(
    csv_file: str, 
    output_file: str, 
    inside_regions: bool = False, 
    verbose: bool = False
)

def generate_cloud_regular(
    csv_file: str, 
    output_file: str, 
    inside_regions: bool = False, 
    verbose: bool = False
)
```
*   `csv_file`: Path to the input contour CSV file.
*   `output_file`: Path where the generated cloud CSV will be saved.
*   `inside_regions`: If `True`, interprets internal closed contours as solid islands (holes in the domain).

#### `reduce_points_by_region`
Reduces the point density of an existing cloud by a scaling factor.
```python
def reduce_points_by_region(
    input_csv: str, 
    output_csv: str, 
    multiplier: int = 2
) -> Optional[List[Dict[str, Any]]]
```
*   `multiplier`: A scaling factor for the distance between nodes. A multiplier of `2` significantly reduces the number of points.

### 3. I/O Operations (`mGFD.io.io`)

#### `load_points`
```python
def load_points(
    file_path: str, 
    verbose: bool = False
) -> np.ndarray
```
*   Reads a CSV file containing `x, y, flag` data and returns it as a standardized `(m, 3)` NumPy array ready for the solvers.

#### `export_stationary_vtk` / `export_transient_vtk` (`mGFD.io.export_vtk`)
Exports the numerical solutions directly to VTK format using PyVista, perfectly structured for advanced 3D analysis in **Paraview**.
```python
def export_stationary_vtk(
    p: np.ndarray, 
    u_ap: np.ndarray, 
    out_dir: str, 
    basename: str = "Stationary_Solution", 
    cloud_path: Optional[str] = None
) -> None

def export_transient_vtk(
    p: np.ndarray, 
    u_ap: np.ndarray, 
    t: int, 
    T: np.ndarray, 
    out_dir: str, 
    basename: str = "Transient_Solution", 
    cloud_path: Optional[str] = None
) -> None
```
*   `u_ap`: The computed numerical approximation array (or matrix for transients).
*   `out_dir`: Directory where the `.vtk` files will be saved.
*   `cloud_path`: (Optional) Passing the original cloud path allows `mGFD` to automatically locate and use the cached triangulation to build the PyVista mesh instantly.

### 4. Visualization Utilities (`mGFD.viz.graph`)

#### `plot_stationary`
```python
def plot_stationary(
    p: np.ndarray, 
    u: np.ndarray, 
    save: bool = False, 
    nom: str = '', 
    title: str = 'Solution', 
    verbose: bool = True
) -> None
```
*   Renders a 3D scatter surface of the solution over the point cloud. If `save=True`, it exports `.png` snapshots to the `nom` path.

#### `plot_transient`
```python
def plot_transient(
    p: np.ndarray, 
    u: np.ndarray, 
    save: bool = False, 
    nom: str = '', 
    title: str = 'Solution', 
    verbose: bool = True
) -> None
```
*   `u` *(m x t)*: A matrix containing the solution at each node for every time step.
*   Renders a 3D Matplotlib animation over time. If `save=True`, it exports an animated `.gif` to the `nom` path.

### 5. Core Utilities & Triangulations (`mGFD.core.utils`)

Because **mGFD** is a strictly *meshless* numerical method, the PDE solvers operate entirely on unstructured points and local distance-based neighborhoods (stencils). Meshes/triangles are **never** used during the numerical solution. 

However, rendering 3D surfaces in Matplotlib and exporting meshes to Paraview (VTK) mathematically requires a triangulation to color the faces between points. 

#### `get_valid_triangulation`
```python
def get_valid_triangulation(p: np.ndarray, nom: Optional[str] = None) -> Optional[np.ndarray]
```
*   **What it does:** Computes a constrained Delaunay triangulation of the point cloud, aggressively filtering out triangles that fall outside concave boundaries or inside holes (islands).
*   **Caching system:** Computing Delaunay simplices on large clouds is extremely expensive. When `mGFD` computes a valid triangulation, it automatically caches the resulting simplices into a `_triangulation.csv` file next to your original cloud. Subsequent calls to plotting or VTK exports will load this cache, saving massive amounts of computation time.

---

## 🔬 Research & Datasets

Looking for the **mathematical formulation**, **real-world geographic lake datasets** (Lake Patzcuaro, Caspian Sea, etc.), or reproducible **benchmarking scripts**?

👉 **[Explore the Research Laboratory (`/research/README.md`)](https://github.com/gstinoco/mGFD/tree/main/research)**

The `research/` directory contains our complete academic suite, including experimental data, geographic boundary files, VTK results, and the exact theoretical explanations associated with our scientific publications.

---

## 🏗️ Repository Structure

Understanding the layout of this project is helpful if you want to contribute or explore the underlying math:

```text
mGFD/
├── docs/                     # Media assets for the README (images, gifs)
├── src/mGFD/                 # The core Python library (uploaded to PyPI)
│   ├── cloud_generator/      # PDE domain creation (natural, regular, etc.)
│   ├── core/                 # Krylov solvers, KDTree neighbors, stencils
│   ├── io/                   # CSV and VTK import/export operations
│   ├── solvers/              # The high-level PDE solvers (Stationary, TimeDerivative)
│   └── viz/                  # Matplotlib 3D and 2D rendering
├── research/                 # The academic reproducible suite
│   ├── codes/                # Batch benchmarking scripts (e.g. main.py)
│   ├── data/                 # Geographic boundaries (CSV files)
│   ├── results/              # Output GIFs, VTKs, PNGs, and metrics
│   └── README.md             # Detailed mathematical foundations
└── tests/                    # CI/CD Unit and Integration tests
```

---

## 🤝 Contributing & Citation

### Contributing
We welcome contributions from the scientific and open-source community! 

1. **Fork** the repository.
2. **Clone** your fork locally and create a new feature branch (`git checkout -b feature-new-solver`).
3. **Develop** your feature, ensuring all code strictly adheres to our PEP-8 standard (especially the inline `#` alignment at column 136).
4. **Test** your code (we run a comprehensive CI suite via GitHub Actions).
5. **Open a Pull Request**!

For major changes, please open an issue first to discuss what you would like to change.

### Citation
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
