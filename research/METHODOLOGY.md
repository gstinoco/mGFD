<div align="center">

# 🧮 Mathematical Methodology: mGFD

*Rigorous mathematical derivation of the meshless Generalized Finite Difference method*

</div>

The meshless Generalized Finite Difference (mGFD) method is a flexible, highly stable numerical scheme for solving Partial Differential Equations (PDEs) directly on unstructured point clouds, bypassing the complex mesh-generation phase required by traditional FEM or FVM.

<br/>

> [!NOTE]
> **Core Novelty (2025):** Unlike other Generalized Finite Difference Methods that require empirical spatial weighting functions (e.g., $w = 1/r^3$ or Gaussian functions) to ensure solvability, the mGFD method enforces consistency through a direct optimization formulation using the Moore-Penrose pseudoinverse. This elegantly minimizes the Euclidean norm of the weights, eliminating the need for arbitrary weight functions and user-defined tuning parameters.

---

## 1. 📐 The General Partial Differential Operator

A general second-order linear partial differential operator $L$ over a scalar field $u(x, y)$ can be expressed as:

<div align="center">

$$ L u = A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u $$

</div>

### 🧩 Operator Components

| Coefficient | Corresponding Term | Physical Interpretation (Example) |
| :---: | :--- | :--- |
| **$A, B, C$** | $u_{xx}, u_{xy}, u_{yy}$ | Diffusion / Dispersive components |
| **$D, E$** | $u_x, u_y$ | Advection / Convective velocity field |
| **$F$** | $u$ | Reaction / Source / Sink term |

In the mGFD scheme, we consider a non-uniform distribution of nodes $P = (x, y)$ as shown below.

<div align="center">
  <img src="docs/images/methodology_nodes.png" alt="Arbitrary distribution of nodes" width="400" style="background-color: #ffffff; padding: 10px; border-radius: 8px; margin: 15px 0; box-shadow: 0 4px 6px rgba(0,0,0,0.08);">
  <br/>
  <sub><em>Fig 1. Arbitrary distribution of nodes.</em></sub>
</div>
<br/>

For an arbitrary central node $p_i = (x_i, y_i)$, the continuous operator $L u$ is approximated discretely using a linear combination of function values at its local neighborhood:

> $$ L_0 u|_{p_i} = \Gamma_{i,0} u(p_i) + \sum_{j=1}^{n_i} \Gamma_{i,j} u(q_{i,j}) $$

**Where:**
*   $q_{i,j} = (\hat{x}_j, \hat{y}_j)$ represents the $j$-th neighbor of $p_i$.
*   $n_i$ is the number of neighbors around $p_i$.
*   $\Gamma_{i,j}$ are the **stencil weights** (or gammas) that must be computed to satisfy consistency.

---

## 2. ⚖️ The Consistency Condition

A finite difference scheme $L_0$ is consistent with the continuous operator $L$ if the local truncation error $\tau|_{p_i}$ tends to zero as the distance between the neighboring nodes and the central node $p_i$ goes to zero.

### 🔍 Taylor Series Expansion

To enforce this, we first expand the function value at each neighboring node $u(q_{i,j})$ using a second-order **Taylor series** around the central node $p_i$:

$$ u(q_{i,j}) = u(p_i) + \Delta x_{i,j} u_x|_{p_i} + \Delta y_{i,j} u_y|_{p_i} + \frac{(\Delta x_{i,j})^2}{2} u_{xx}|_{p_i} + \Delta x_{i,j} \Delta y_{i,j} u_{xy}|_{p_i} + \frac{(\Delta y_{i,j})^2}{2} u_{yy}|_{p_i} + \mathcal{O}(h^3) $$

*(where $\Delta x_{i,j} = \hat{x}_j - x_i$ and $\Delta y_{i,j} = \hat{y}_j - y_i$)*

Substituting this expansion back into the discrete operator $L_0 u|_{p_i}$ and grouping the terms yields a set of algebraic constraints. For the truncation error to vanish perfectly, we generate the following system:

```math
    \begin{aligned}
        A(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} \frac{(\Delta x_{i,j})^2}{2} &= 0 \\
        B(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} \Delta x_{i,j} \Delta y_{i,j} &= 0 \\
        C(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} \frac{(\Delta y_{i,j})^2}{2} &= 0 \\
        D(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} \Delta x_{i,j} &= 0 \\
        E(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} \Delta y_{i,j} &= 0 \\
        F(p_i) - \sum_{j=0}^{n_i} \Gamma_{i,j} &= 0
    \end{aligned}
```

---

## 3. 🎯 Linear System & Optimization

The consistency conditions defined above can be assembled into a purely geometric linear system $M \Gamma = \beta$ that dictates the values of the neighbor weights $\Gamma_{i,j}$ (for $j \geq 1$):

<div align="center">

```math
    \begin{pmatrix}
        \Delta x_{i,1} & \dots & \Delta x_{i,n_i} \\
        \Delta y_{i,1} & \dots & \Delta y_{i,n_i} \\
        (\Delta x_{i,1})^2 & \dots & (\Delta x_{i,n_i})^2 \\
        \Delta x_{i,1} \Delta y_{i,1} & \dots & \Delta x_{i,n_i} \Delta y_{i,n_i} \\
        (\Delta y_{i,1})^2 & \dots & (\Delta y_{i,n_i})^2
    \end{pmatrix}
    \begin{pmatrix}
        \Gamma_{i,1} \\
        \Gamma_{i,2} \\
        \vdots \\
        \Gamma_{i,n_i}
    \end{pmatrix}
        =
    \begin{pmatrix}
        D(p_i) \\
        E(p_i) \\
        2A(p_i) \\
        B(p_i) \\
        2C(p_i)
    \end{pmatrix}
```

</div>

<div align="center">
  <img src="docs/images/methodology_neighbors.png" alt="Selection of the neighbor nodes" width="400" style="background-color: #ffffff; padding: 10px; border-radius: 8px; margin: 15px 0; box-shadow: 0 4px 6px rgba(0,0,0,0.08);">
  <br/>
  <sub><em>Fig 2. Selection of the neighbor nodes within an influence radius.</em></sub>
</div>
<br/>

Because the number of neighbors $n_i$ is typically chosen to be strictly greater than the number of derivative constraints (5 equations for a 2D second-order problem), the system is underdetermined.

> [!IMPORTANT]
> **The Moore-Penrose Pseudoinverse:** The mGFD method elegantly resolves this system by computing the pseudoinverse $M^+$. Mathematically, this automatically finds the unique solution that **minimizes the Euclidean norm** of the weights ($||\Gamma||_2$), guaranteeing optimal numerical stability.

<div align="center">

$$ \Gamma = M^+ \beta $$

</div>

Finally, the central weight $\Gamma_{i,0}$ is recovered explicitly from the isolated reaction term constraint:

$$ \Gamma_{i,0} = F(p_i) - \sum_{j=1}^{n_i} \Gamma_{i,j} $$

### 🎛️ Stencil Conditioning & Distance-Weighted Pseudoinverse

On dense or irregular point clouds where nodal spacing satisfies $h \ll 1$, powers of $\Delta x$ and $\Delta y$ can lead to severe numerical ill-conditioning in the Vandermonde-like matrix $M$. 

To guarantee unconditional numerical stability:
1. **Local Coordinate Scaling ($h_i$):** Local offsets are normalized by the maximum neighborhood radius:
   $$ h_i = \max_{j=1,\dots,n_i} \sqrt{(\Delta x_{i,j})^2 + (\Delta y_{i,j})^2}, \quad \Delta \bar{x}_{i,j} = \frac{\Delta x_{i,j}}{h_i}, \quad \Delta \bar{y}_{i,j} = \frac{\Delta y_{i,j}}{h_i} $$
   This bounds the condition number of the scaled matrix $\kappa(\bar{M}) \le 30$, completely eliminating floating-point precision loss. The continuous weights $\Gamma$ are recovered by inverse powers of $h_i$.

2. **Distance-Weighted Least Squares:** To prioritize immediate neighbors and enforce smooth truncation decay, quadratic distance weighting $W = \text{diag}\left( \frac{1}{\bar{r}_j^2 + \epsilon} \right)$ is incorporated, solving the weighted least-squares problem:
   $$ \min_{\Gamma} \| W^{1/2} (M \Gamma - \beta) \|_2 $$

### 🚪 Enforcement of Boundary Conditions in the Global System

For a point cloud with $m$ total nodes ($m_{\text{int}}$ interior nodes and $m_{\text{bnd}}$ boundary nodes), the discrete operator is assembled into a global sparse matrix $K \in \mathbb{R}^{m \times m}$ (Compressed Sparse Row format, `scipy.sparse.csr_matrix` on CPU and `cupyx.scipy.sparse.csr_matrix` on CUDA GPU).

Boundary conditions are enforced row-by-row on the global system:

* **Dirichlet Boundary Conditions ($u|_{\partial \Omega} = g(x, y)$):**
  For each boundary node $p_i \in \partial \Omega$:
  $$ K_{i, i} = 1, \quad K_{i, j} = 0 \quad (\forall j \ne i), \quad b_i = g(p_i) $$
  This uncouples the boundary degrees of freedom, preserving the exact prescribed values without perturbing the interior stencil equations.

* **Neumann Boundary Conditions ($\frac{\partial u}{\partial \hat{n}}|_{\partial \Omega} = h(x, y)$):**
  Given the unit outward normal vector $\hat{n} = (n_x, n_y)$ at boundary node $p_i$, the directional derivative is approximated via a first-order mGFD stencil with continuous operator $L_{\hat{n}} = n_x \frac{\partial}{\partial x} + n_y \frac{\partial}{\partial y}$:
  $$ \sum_{j=0}^{n_i} \Gamma_{i,j}^{(\hat{n})} u(q_{i,j}) = h(p_i) $$
  Row $i$ of the global matrix $K$ is populated directly with the directional weights $\Gamma_{i,j}^{(\hat{n})}$.

---

## 4. 🌊 The Upwind mGFD Scheme

When solving **advection-dominated problems** (such as the Advection-Diffusion equation with high Péclet numbers), standard central stencils can induce severe numerical instabilities and spurious oscillations.

The mGFD method mitigates this through a rigorously defined **Upwind Stencil** geometric strategy:

| Step | Action | Description |
| :---: | :--- | :--- |
| **1** | **Flow Evaluation** | The velocity vector $\vec{v} = (v_x, v_y)$ is evaluated directly at the central node $p_i$. |
| **2** | **Upstream Selection** | An imaginary line perpendicular to $\vec{v}$ is drawn through $p_i$. Only the nodes $q_{i,j}$ located *upstream* (in the opposite direction of the flow) are selected. |

<div align="center">
  <img src="docs/images/methodology_upwind.png" alt="Upwind neighbor selection" width="420" style="background-color: #ffffff; padding: 10px; border-radius: 8px; margin: 15px 0; box-shadow: 0 4px 6px rgba(0,0,0,0.08);">
  <br/>
  <sub><em>Fig 3. Selection of the neighbor nodes for an upwind scheme.</em></sub>
</div>
<br/>
| **3** | **Information Propagation** | This guarantees that relevant flow information propagates physically and stably in the advective direction. |

> [!TIP]
> By enforcing the consistency conditions exclusively on the upwind nodes, mGFD preserves rigorous mathematical stability without requiring artificial viscosity, stabilizing parameters, or ad-hoc upwind blending factors.

---

## 5. ⏳ Time Discretization & Stability (Transient Schemes)

For time-dependent PDEs, the mGFD method couples the spatial discretization with unconditionally stable, high-order finite difference time-stepping schemes. 

### A. Parabolic Equations (First-Order in Time: Heat & Advection-Diffusion)
For general first-order transient PDEs:

<div align="center">

$$ \frac{\partial u}{\partial t} = L u + F(x, y, t) $$

</div>

The time derivative is discretized using an implicit $\theta$-weighted scheme ($\theta \in [0, 1]$):

<div align="center">

$$ \frac{u^{k+1} - u^k}{\Delta t} = \theta K u^{k+1} + (1 - \theta) K u^k + F^{k+\theta} $$

</div>

Grouped into a linear system for the future state $u^{k+1}$:

<div align="center">

$$ (I - \theta \Delta t K) u^{k+1} = (I + (1 - \theta) \Delta t K) u^k + \Delta t F^{k+\theta} $$

</div>

> [!TIP]
> **Crank-Nicolson ($\theta = 0.5$):** Setting $\theta = 0.5$ yields the second-order Crank-Nicolson scheme ($\mathcal{O}(\Delta t^2)$). Because the scheme is A-stable, it bypasses the prohibitive explicit parabolic Courant stability restriction ($\Delta t \propto h^2$), allowing time steps to scale linearly with characteristic spacing ($\Delta t \propto h_{\text{char}}$) and dramatically accelerating transient solves.

---

### B. Hyperbolic Equations & Conservative Laplacian Symmetrization (Second-Order in Time: Wave)
For second-order hyperbolic systems:

<div align="center">

$$ \frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) + F(x, y, t) $$

</div>

#### The Asymmetric Energy-Pumping Problem
On unstructured point clouds, the $k$-nearest-neighbor graph contains **18% to 24% unidirectional edges** ($j$ is a neighbor of $i$, but $i$ is not a neighbor of $j$). Consequently, the Taylor spatial matrix $K$ is non-symmetric:

$$ K = K_{\text{sym}} + K_{\text{skew}}, \quad K_{\text{skew}} = \frac{1}{2}(K - K^T) \ne 0 $$

In an oscillatory, energy-conserving hyperbolic PDE, the skew-symmetric component acts as an unphysical energy pump:

$$ \frac{dE}{dt} = u_t^T K_{\text{skew}} u \ne 0 $$

This produces eigenvalues with non-zero imaginary parts ($\text{Im}(\lambda) \ne 0$). Under centered time integration, the amplification factor roots exit the unit circle ($|r| > 1.0$), resulting in catastrophic exponential divergence on dense meshes.

#### The Conservative Symmetrization Solution
To enforce strict conservation of energy and unconditional stability, `mGFD` implements off-diagonal edge symmetrization with exact row-sum diagonal reconstruction:

1. **Edge Symmetrization**:
   $$ W = \max(0, K_{\text{off}}), \quad W_{\text{sym}} = \frac{1}{2}(W + W^T) $$
2. **Conservative Diagonal Preservation**:
   $$ D_{ii} = -\sum_{j \ne i} W_{\text{sym}, ij} \quad (+ c_{\text{react}} \text{ if reaction present}) $$
   $$ K_{\text{cons}} = W_{\text{sym}} + \text{diag}(D) $$

**Guaranteed Mathematical Properties:**
- $K_{\text{cons}} = K_{\text{cons}}^T$ strictly symmetric.
- $x^T K_{\text{cons}} x \le 0$ unconditionally negative semi-definite.
- All eigenvalues are purely real and non-positive: $\text{Im}(\lambda) \equiv 0, \; \lambda \le 0$.
- The Newmark-$\beta$ amplification roots satisfy $|r| \equiv 1$ unconditionally, guaranteeing 100% stability across all irregular geometries.

#### Temporal Schemes: Newmark-$\beta$ and Hilber-Hughes-Taylor ($\alpha$)
The second-order system is integrated in time using the generalized Newmark-$\beta$ family:

<div align="center">

$$ u^{k+1} = u^k + \Delta t v^k + \Delta t^2 \left[ \left(\frac{1}{2} - \beta\right) a^k + \beta a^{k+1} \right] $$
$$ v^{k+1} = v^k + \Delta t \left[ (1 - \gamma) a^k + \gamma a^{k+1} \right] $$

</div>

For high-frequency noise dissipation on highly irregular domains, `mGFD` supports **Hilber-Hughes-Taylor ($\alpha$)** numerical dissipation ($\alpha \in [-0.333, 0.0]$, $\gamma = 0.5 - \alpha$, $\beta = 0.25(1 - \alpha)^2$) alongside physical velocity damping ($\eta u_t$), damping unresolvable grid-scale oscillations while preserving physical low-frequency wave modes.

---

### C. Adaptive Characteristic CFL Time-Stepping
To guarantee Courant stability across diverse point cloud scales without manual tuning, `mGFD` incorporates an automatic CFL estimator (`estimate_cfl_dt`):

<div align="center">

$$ \Delta t_{\text{max}} = \text{CFL} \cdot \frac{h_{\text{char}}}{V_{\text{adv}} + \sqrt{\nu}} \quad (\text{Order 1 Parabolic/Advective}) $$
$$ \Delta t_{\text{max}} = 0.60 \cdot \text{CFL} \cdot \frac{h_{\text{char}}}{c\sqrt{2}} \quad (\text{Order 2 Hyperbolic Wave}) $$

</div>

where the characteristic node spacing $h_{\text{char}} = \max(h_{\text{min}}, 0.25 h_{\text{avg}})$ accounts for local density variations near irregular lake boundaries.

---

## 6. 🧱 Discontinuous Coefficients (Multilayer Interfaces)

The mGFD method natively handles multilayer domains where physical coefficients (e.g., diffusivity $\nu$) present sharp discontinuities across internal interfaces.

To enforce the continuity of the normal flux across an interface boundary, the method equates the directional derivatives between the two adjoining layers:

<div align="center">

$$ \nu_1 \frac{\partial u}{\partial \hat{n}}(p_i) = \nu_2 \frac{\partial u}{\partial \hat{n}}(p_i) $$

</div>

<div align="center">
  <img src="docs/images/methodology_ghost.png" alt="Ghost node selection" width="380" style="background-color: #ffffff; padding: 10px; border-radius: 8px; margin: 15px 0; box-shadow: 0 4px 6px rgba(0,0,0,0.08);">
  <br/>
  <sub><em>Fig 4. Ghost node projection across an internal interface.</em></sub>
</div>
<br/>

To resolve this numerically, **ghost nodes** are dynamically projected across the interface into the adjacent domain along the normal vector $\hat{n}$. The standard mGFD spatial discretization is then applied to the interface balance equation, establishing a coupled relationship between the true boundary nodes and the ghost nodes, ensuring continuous flux without requiring adaptive meshing.
