<div align="center">

# 🗂️ Benchmark Datasets

*Highly irregular, real-world geometric point clouds for numerical validation*

</div>

The `mGFD` research laboratory includes a comprehensive suite of geographic datasets. These domains represent highly irregular, real-world geometries that challenge traditional mesh-based methods. All point clouds were generated using the companion Poisson-Disk algorithm from the `mGFD CloudGenerator`.

---

## 1. 🌍 Geographic Regions

The laboratory provides point clouds for **20 distinct geographic regions** (primarily lakes from around the world). These regions test the method's ability to handle complex boundary curvatures, narrow bottlenecks, and multi-scale features.

| Region Name | Abbreviation | Location | Geometric Description |
|:---|:---:|:---|:---|
| **Alchichica** | `ALC` | 🇲🇽 Mexico | Maar crater lake testing pure radial symmetry and smooth curves. |
| **Baikal** | `BAI` | 🇷🇺 Russia | Deepest freshwater lake; tests elongated, high-aspect-ratio tectonic geometry. |
| **Balkhash** | `BAL` | 🇰🇿 Kazakhstan | Narrow, crescent-shaped lake testing flow through complex straits. |
| **Caspian** | `CAS` | 🌍 Eurasia | Largest enclosed inland body; vast multi-scale domain with smooth and jagged coasts. |
| **Catemaco** | `CAT` | 🇲🇽 Mexico | Volcanic lake featuring interior islands (challenging internal boundary holes). |
| **Chapala** | `CHA` | 🇲🇽 Mexico | Extensive lake with smooth, elongated boundaries. |
| **Cuitzeo** | `CUI` | 🇲🇽 Mexico | Highly irregular lake prone to severe fragmentation and narrow land bridges. |
| **Erie** | `ERI` | 🌎 North America | Great Lake with asymmetric borders and large-scale curvature. |
| **Huron** | `HUR` | 🌎 North America | Features massive bays and complex island clusters (e.g., Manitoulin). |
| **Ladoga** | `LAD` | 🇷🇺 Russia | Large lake with a highly fragmented, island-rich northern coast. |
| **Malawi** | `MAL` | 🌍 Africa | Long, narrow tectonic trench testing extreme length-to-width aspect ratios. |
| **Metztitlán** | `MET` | 🇲🇽 Mexico | Variable boundaries; the primary visual benchmark for unstructured clouds. |
| **Ness** | `NES` | 🏴󠁧󠁢󠁳󠁣󠁴󠁿 Scotland | Extremely narrow, straight glacial loch testing highly confined flow. |
| **Patzcuaro** | `PAT` | 🇲🇽 Mexico | C-shaped volcanic lake with prominent central islands. |
| **Poopo** | `POO` | 🇧🇴 Bolivia | Shallow saline lake with highly variable, vanishing borders. |
| **Santa Maria del Oro**| `SMO` | 🇲🇽 Mexico | Small, nearly perfectly circular crater lake. |
| **Titicaca** | `TIT` | 🌎 South America | Complex domain featuring deeply divided sub-basins connected by a strait. |
| **Yuriria** | `YUR` | 🇲🇽 Mexico | Artificial volcanic crater lake with a distinctive irregular perimeter. |
| **Zempoala** | `ZEM` | 🇲🇽 Mexico | Small, interconnected alpine lakes with intricate boundary merging. |
| **Zirahuen** | `ZIR` | 🇲🇽 Mexico | Deep mountain lake with a compact, teardrop-like shape. |

---

## 2. 📏 Density Scales

To evaluate numerical convergence and execution performance of the `mGFD` method, each geographic region is discretized across **4 calibrated node densities** (Scales 1 through 4) generated via independent Poisson-Disk Sampling (`generate_cloud_natural`) with calibrated density multipliers ($d = 0.10, 0.20, 0.30, 0.40$).

Inside `research/data/<Region>/`, you will find subdirectories numbered `1` through `4`:

> [!NOTE]
> **Scale Range:**
> *   `1/`: The **coarsest** point cloud representation (~500 to ~2,500 nodes; ideal for fast testing and parameter sweeps).
> *   `2/`: **Intermediate** resolution (~1,500 to ~9,500 nodes; default scale for graphical animation and MP4 rendering).
> *   `3/`: **Fine** resolution (~3,000 to ~21,000 nodes; high-fidelity physics).
> *   `4/`: **Ultra-dense** resolution (~5,000 to ~37,000 nodes; asymptotic convergence and GPU speedup benchmarking).
>
> *Note on Graphical Outputs:* To optimize repository storage and execution times during massive automated parameter sweeps, MP4 animations and PNG renderings are generated exclusively for **Scale 2** by default (`plot_scales: ["2"]` in `sweep_config.json`).

### 📊 Complete Discretization Census Across All 20 Geometries

The table below summarizes the exact node counts across all 80 benchmark point clouds (20 geographic regions $\times$ 4 scales):

| Geographic Region | Scale 1 (Coarse) | Scale 2 (Nominal) | Scale 3 (Fine) | Scale 4 (Dense) | Asymptotic Refinement Factor (4 / 1) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Alchichica** | 2,464 | 9,455 | 20,894 | 36,895 | **14.97x** |
| **Baikal** | 673 | 2,000 | 4,159 | 7,141 | **10.61x** |
| **Balkhash** | 525 | 1,487 | 2,797 | 4,680 | **8.91x** |
| **Caspio** | 806 | 2,775 | 5,889 | 10,277 | **12.75x** |
| **Catemaco** | 1,444 | 4,914 | 10,523 | 18,395 | **12.74x** |
| **Chapala** | 1,001 | 3,505 | 7,553 | 13,191 | **13.18x** |
| **Cuitzeo** | 1,151 | 4,013 | 8,675 | 15,135 | **13.15x** |
| **Erie** | 1,134 | 3,817 | 8,235 | 14,339 | **12.64x** |
| **Huron** | 905 | 2,816 | 5,747 | 9,819 | **10.85x** |
| **Ladoga** | 1,540 | 5,318 | 11,454 | 19,793 | **12.85x** |
| **Malawi** | 714 | 2,388 | 5,068 | 8,796 | **12.32x** |
| **Metztitlán** | 667 | 2,301 | 4,803 | 8,334 | **12.49x** |
| **Ness** | 492 | 1,512 | 3,046 | 5,252 | **10.67x** |
| **Patzcuaro** | 958 | 2,915 | 5,986 | 10,202 | **10.65x** |
| **Poopo** | 1,060 | 3,736 | 7,927 | 13,939 | **13.15x** |
| **Santa Maria del Oro** | 2,134 | 8,225 | 18,206 | 32,026 | **15.01x** |
| **Titicaca** | 733 | 2,181 | 4,339 | 7,326 | **9.99x** |
| **Yuriria** | 1,160 | 4,120 | 8,773 | 15,324 | **13.21x** |
| **Zempoala** | 909 | 3,200 | 6,701 | 11,695 | **12.87x** |
| **Zirahuen** | 1,469 | 5,577 | 12,118 | 21,260 | **14.47x** |

---

## 3. 🧩 Cloud Components

Within every scale directory (e.g., `research/data/Metztitlán/2/`), you will find the spatial domain files essential for the simulations.

### 📍 The Unstructured Cloud (`*_cloud.csv`)
This is the primary spatial domain. It contains the $m \times 3$ coordinate matrix where columns represent: `[x, y, flag]`. 
*   `flag = 0`: Interior computational nodes
*   `flag = 1`: Primary outer boundary nodes
*   `flag = 2`: Inner boundary nodes (e.g., islands or internal boundary holes)

### 🕸️ Dynamic In-Memory KD-Tree Stencils
In modern `mGFD` (v0.11+), neighbor stencils (`vec`) are computed dynamically in-memory via Numba JIT-accelerated spatial KDTrees (`Neighbors.compute_neighbors(p, nvec)` and `Neighbors.compute_upwind_neighbors(p, a, b, nvec)`). Dynamic stencil generation executes in **$< 0.05$ seconds** even for 37,000 nodes, completely eliminating the need for stale on-disk neighbor caches (`.csv`) and guaranteeing fresh, equation-specific stencils.

> [!WARNING]
> ### 🔺 Visualization Mesh (`*_triangulation.csv` / VTK)
> `mGFD` is strictly a **meshless** method. The mathematical solvers *never* use elements, cells, or connectivity matrices. However, plotting tools and external post-processors (ParaView) require surface elements to render continuous color maps. `SolverResult.export_vtk("solution.vtu")` automatically generates high-quality UnstructuredGrid VTK files.

---

## 4. 💻 Loading and Solving on the Datasets (Python)

Using the modern `mGFD.oop` interface, loading a geographic cloud, constructing the differential operator, solving the PDE, and exporting results is clean, expressive, and robust:

```python
import numpy as np
from mGFD.oop import Cloud, Dirichlet, Domain, PDE, Solver
from mGFD.io.io import load_points

# 1. Load the unstructured point cloud
cloud = Cloud.from_csv("research/data/Metztitlán/2/Metztitlán_cloud.csv")
print(f"Loaded {cloud.n_nodes} nodes (Interior: {cloud.n_interior}, Boundary: {cloud.n_boundary})")

# 2. Dynamic KD-Tree stencil generation in-memory (<0.05s)
cloud.compute_neighbors(nvec=16)

# 3. Define Dirichlet boundary condition matching exact solution u = 2 * exp(2x + y)
bc = Dirichlet(lambda x, y: 2.0 * np.exp(2.0 * x + y))
domain = cloud.set_boundary(bc)

# 4. Define the Poisson PDE: u_xx + u_yy = 10 * exp(2x + y)
# Operator format: [D, E, 2*A, B, 2*C, F] -> [0, 0, 2, 0, 2, 0]
poisson_pde = PDE(
    operator=[0, 0, 2, 0, 2, 0],
    f=lambda x, y: 10.0 * np.exp(2.0 * x + y)
)

# 5. Solve on CPU (or pass device="cuda" for GPU acceleration)
solver = Solver(domain=domain, pde=poisson_pde, device="cpu")
result = solver.solve()

# 6. Evaluate accuracy and export to ParaView
u_exact = 2.0 * np.exp(2.0 * cloud.nodes[:, 0] + cloud.nodes[:, 1])
rmse = np.sqrt(np.mean((result.solution - u_exact) ** 2))
print(f"Solve completed in {result.execution_time:.4f}s | RMSE: {rmse:.3e}")

# Export for high-end rendering in ParaView
result.export_vtk("research/results/Metztitlan_Poisson_Scale2.vtu")
```

