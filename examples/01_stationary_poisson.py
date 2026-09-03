"""
Example 01: Solving the Stationary Poisson Equation using mGFD (OOP Interface)

Overview:
    This tutorial demonstrates how to solve the classic Poisson Equation:
        u_xx + u_yy = f(x, y)
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
    print("                Example 01: Stationary Poisson PDE (OOP API)                     ")                                          # Title log.
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
    bc_func = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)                                                                        # Analytical boundary condition.
    domain  = cloud.set_boundary(mgfd.Dirichlet(bc_func))                                                                               # Bind Dirichlet boundary condition.

    # 3. Define Physics & Instantiate Solver
    print("\nStep 3: Formulating Poisson PDE physics and solver...")                                                                    # Log step 3.
    force_func = lambda x, y: -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)                                                   # Forcing term function.
    pde        = mgfd.PoissonEquation(source=force_func)                                                                                # Formulate Poisson PDE.
    solver     = mgfd.Solver(domain, pde, nvec=15, verbose=True)                                                                        # Instantiate high-level Solver.

    # 4. Solve and Plot
    print("\nStep 4: Solving PDE and rendering results...")                                                                             # Log step 4.
    result = solver.solve()                                                                                                             # Execute solver.
    print(f"Solution completed in {result.compute_time:.4f} seconds!")                                                                  # Log compute time.

    result.plot(save=False, show=True, filename='01_poisson_result', title="mGFD Solution: Poisson on Star Domain")                     # Render solution plot.

    # Clean up temporary CSV files
    if os.path.exists(contour_file): os.remove(contour_file)                                                                            # Clean contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                # Clean cloud file.

if __name__ == '__main__':                                                                                                              # Entry point guard.
    main()                                                                                                                              # Run main.
