"""
Example 02: Solving the Transient Heat Equation using mGFD (OOP Interface)

Overview:
    This tutorial demonstrates how to solve the Heat Equation (First-order in time):
        u_t - k * (u_xx + u_yy) = F(x, y, t)
    On an irregular 2D domain (Star), using the modern OOP Architecture of mGFD.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

Date:
    September, 2026.
Last Modification:
    September, 2026.
"""

## Library importation.
import csv, os                                                                                                                          # Standard OS and CSV interfaces.
import numpy as np                                                                                                                      # Core numerical array operations.
import mGFD as mgfd                                                                                                                     # Import mGFD library.

def main() -> None:                                                                                                                     # Main execution routine.
    print("=================================================================================")                                          # Separator log.
    print("                 Example 02: Transient Heat PDE (OOP API)                       ")                                           # Title log.
    print("=================================================================================\n")                                        # Separator log.

    # 1. Define geometry contour & generate point cloud via mGFD.Cloud
    print("Step 1: Building irregular point cloud geometry...")                                                                         # Log step 1.
    star_contour = [(0.5, 1.0), (0.65, 0.65), (1.0, 0.5), (0.65, 0.35), (0.5, 0.0), (0.35, 0.35), (0.0, 0.5), (0.35, 0.65)]             # Define star contour vertices.
    contour_file = 'star_contour.csv'                                                                                                   # Contour filename.
    cloud_file   = 'star_cloud.csv'                                                                                                     # Output cloud filename.

    with open(contour_file, 'w', newline='') as f:                                                                                      # Write CSV contour.
        writer = csv.writer(f)                                                                                                          # Initialize CSV writer.
        writer.writerow(['x', 'y'])                                                                                                     # Write header.
        for pt in star_contour: writer.writerow(pt)                                                                                     # Write vertices.

    cloud = mgfd.Cloud.generate_natural(contour_file, cloud_file, cloud_size=0.03, save=False, show=False)                              # High-level point cloud generation.
    print(f"Generated {cloud}: {cloud.num_nodes} nodes.")                                                                               # Log cloud summary.

    # 2. Set Dirichlet boundary conditions and construct Domain
    print("\nStep 2: Binding Dirichlet boundary condition and domain...")                                                               # Log step 2.
    k_diff  = 0.08                                                                                                                      # Thermal diffusivity constant.
    bc_func = lambda x, y, t=0, coef=None: np.exp(-2 * np.pi**2 * k_diff * t) * np.sin(np.pi * x) * np.sin(np.pi * y)                   # Spatiotemporal boundary condition.
    domain  = cloud.set_boundary(mgfd.Dirichlet(bc_func))                                                                               # Bind Dirichlet boundary condition.

    # 3. Define Physics & Instantiate Solver
    print("\nStep 3: Formulating Heat PDE physics and solver...")                                                                       # Log step 3.
    ic_func     = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)                                                                    # Initial condition u(x, y, 0).
    heat_source = lambda x, y, t=0, coef=None: 10.0 * np.sin(4 * np.pi * t) * np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))               # Pulsed heat laser source F(x, y, t).
    
    pde         = mgfd.HeatEquation(k=k_diff, ic=ic_func, source=heat_source)                                                           # Formulate Heat PDE physics.
    solver      = mgfd.Solver(domain, pde, nvec=15, cfl=0.5, verbose=True)                                                              # Instantiate high-level Solver.

    # 4. Solve over t_span=(0, 2.0) and Plot
    print("\nStep 4: Solving Heat PDE over t_span=(0.0, 2.0)...")                                                                       # Log step 4.
    t_span = (0.0, 2.0)                                                                                                                 # Physical time span domain.
    result = solver.solve(t_span=t_span)                                                                                                # Execute transient solver.
    print(f"Transient solution computed in {result.compute_time:.4f} seconds (dt = {result.dt:.6f}, t_steps = {result.t_steps})!")      # Log metrics.

    result.plot(save=False, show=True, filename='02_heat_result', title="mGFD Solution: Heat Equation", t_span=t_span)                  # Render solution plot.

    # Clean up temporary CSV files
    if os.path.exists(contour_file): os.remove(contour_file)                                                                            # Clean contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                # Clean cloud file.

if __name__ == '__main__':                                                                                                              # Entry point guard.
    main()                                                                                                                              # Run main.
