r"""
Example 03: Solving the Wave Equation using mGFD (OOP Interface)

Overview:
    This tutorial demonstrates how to solve a hyperbolic PDE (the Wave Equation):
        u_tt + \eta u_t = c^2 (u_xx + u_yy) + F(x, y, t)
    On a star domain, using the modern OOP Architecture of mGFD with HHT-\alpha and damping.

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
    print("                 Example 03: Wave Equation PDE (OOP API)                        ")                                           # Title log.
    print("=================================================================================\n")                                        # Separator log.

    # 1. Define geometry contour & generate point cloud via mGFD.Cloud
    print("Step 1: Building irregular point cloud geometry...")                                                                         # Log step 1.
    theta  = np.linspace(0, 2 * np.pi, 150)                                                                                             # Parametric angular coordinates.
    r_star = 0.5 + 0.15 * np.sin(5 * theta)                                                                                             # Polar radius of star domain.
    star_contour = [(0.5 + r_star[i] * np.cos(theta[i]), 0.5 + r_star[i] * np.sin(theta[i])) for i in range(len(theta))]                # Compute (x, y) contour boundary.

    contour_file = 'star_contour_wave.csv'                                                                                              # Contour filename.
    cloud_file   = 'star_cloud_wave.csv'                                                                                                # Output cloud filename.

    with open(contour_file, 'w', newline='') as f:                                                                                      # Write CSV contour.
        writer = csv.writer(f)                                                                                                          # Initialize CSV writer.
        writer.writerow(['x', 'y'])                                                                                                     # Write header.
        for pt in star_contour: writer.writerow(pt)                                                                                     # Write vertices.

    cloud = mgfd.Cloud.generate_natural(contour_file, cloud_file, cloud_size=0.015, save=False, show=False)                             # High-level point cloud generation.
    print(f"Generated {cloud}: {cloud.num_nodes} nodes.")                                                                               # Log cloud summary.

    # 2. Set Dirichlet boundary conditions and construct Domain
    print("\nStep 2: Binding Dirichlet boundary condition and domain...")                                                               # Log step 2.
    domain = cloud.set_boundary(mgfd.Dirichlet(0.0))                                                                                    # Bind zero Dirichlet boundary condition.

    # 3. Define Physics & Instantiate Solver
    print("\nStep 3: Formulating Wave PDE physics with HHT-alpha and velocity damping...")                                              # Log step 3.
    ic_func = lambda x, y: np.exp(-100 * ((x - 0.5)**2 + (y - 0.5)**2))                                                                 # Initial position Gaussian bump.
    c_wave  = 0.8                                                                                                                       # Wave speed propagation.
    
    pde     = mgfd.WaveEquation(c=c_wave, damping=0.05, alpha=-0.15, ic=ic_func, g=0.0)                                                 # Formulate Wave PDE physics.
    solver  = mgfd.Solver(domain, pde, nvec=16, cfl=0.25, verbose=True)                                                                 # Instantiate high-level Solver.

    # 4. Solve over t_span=(0, 3.0) and Plot
    print("\nStep 4: Solving Wave PDE over t_span=(0.0, 3.0)...")                                                                       # Log step 4.
    t_span = (0.0, 3.0)                                                                                                                 # Physical time span domain.
    result = solver.solve(t_span=t_span)                                                                                                # Execute transient solver.
    print(f"Wave solution computed in {result.compute_time:.4f} seconds (dt = {result.dt:.6f}, t_steps = {result.t_steps})!")           # Log metrics.

    result.export_vtk("wave_solution.vtu")                                                                                              # Export to ParaView VTK file.
    print("Exported VTK file: wave_solution.vtu")                                                                                       # Log VTK export confirmation.

    result.plot(save=False, show=True, filename='03_wave_result', title="mGFD Solution: Wave Equation", t_span=t_span)                  # Render solution plot.

    # Clean up temporary CSV files
    if os.path.exists(contour_file): os.remove(contour_file)                                                                            # Clean contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                # Clean cloud file.
    if os.path.exists("wave_solution.vtu"): os.remove("wave_solution.vtu")                                                              # Clean VTK file.

if __name__ == '__main__':                                                                                                              # Entry point guard.
    main()                                                                                                                              # Run main.
