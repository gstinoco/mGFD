"""
Example 04: High Performance & Optimizers — GPU CUDA and CPU Solving using mGFD (OOP Interface)

Overview:
    This tutorial demonstrates how to use mGFD's advanced optimizers to solve highly complex,
    large-scale PDEs efficiently. We simulate an Advection-Diffusion Equation:
        u_t = v(u_xx + u_yy) - (v_x u_x + v_y u_y) + F(x, y, t)
    on a Kite-shaped domain using the modern OOP Architecture of mGFD, comparing CPU vs CUDA backends.

Public API:
    main                                    Main execution routine for the High Performance optimizers tutorial.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx
    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México.
        Aula CIMNE-Morelia. México.
        SIIIA-MATH: Soluciones de Ingeniería. México.

    Based on the theoretical concepts presented in:
        "mGFD: A meshless generalized finite difference method",
        Gerardo Tinoco-Guerrero, Francisco Javier Domínguez-Mota, José Alberto Guzmán-Torres, 
        Gabriela Pedraza-Jiménez, José Gerardo Tinoco-Ruiz,
        Computers & Mathematics with Applications, Volume 195 (2025) 396-418.
        https://doi.org/10.1016/j.camwa.2025.07.034

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
    """
    main
    Executes Example 04 tutorial comparing CPU and GPU CUDA solvers for advection-diffusion.

    Input:
        None

    Output:
        None
    """
    print("=================================================================================")                                          # Separator log.
    print("            Example 04: Advection-Diffusion & High Performance (OOP API)         ")                                          # Title log.
    print("=================================================================================\n")                                        # Separator log.

    # 1. Define geometry contour & generate point cloud via mGFD.Cloud
    print("Step 1: Building irregular point cloud geometry...")                                                                         # Log step 1.
    kite_contour = [
        (1.0, 0.5), (0.85, 0.7), (0.5, 1.0), (0.15, 0.7),
        (0.0, 0.5), (0.15, 0.3), (0.5, 0.0), (0.85, 0.3)
    ]                                                                                                                                   # Define Kite boundary vertices.

    contour_file = 'kite_contour.csv'                                                                                                   # Contour filename.
    cloud_file   = 'kite_cloud.csv'                                                                                                     # Output cloud filename.

    with open(contour_file, 'w', newline='') as f:                                                                                      # Write CSV contour.
        writer = csv.writer(f)                                                                                                          # Initialize CSV writer.
        writer.writerow(['x', 'y'])                                                                                                     # Write header.
        for pt in kite_contour: writer.writerow(pt)                                                                                     # Write vertices.

    cloud = mgfd.Cloud.generate_natural(contour_file, cloud_file, cloud_size=0.015, save=False, show=False)                             # High-level point cloud generation.
    print(f"Generated {cloud}: {cloud.num_nodes} nodes.")                                                                               # Log cloud summary.

    # 2. Set Dirichlet boundary conditions and construct Domain
    print("\nStep 2: Binding Dirichlet boundary condition and domain...")                                                               # Log step 2.
    domain = cloud.set_boundary(mgfd.Dirichlet(0.0))                                                                                    # Bind zero Dirichlet boundary condition.

    # 3. Define Physics & Instantiate Solvers
    print("\nStep 3: Formulating generalized Advection-Diffusion PDE physics...")                                                       # Log step 3.
    v_diff, v_x, v_y = 0.01, 0.15, 0.15                                                                                                 # Physical diffusion and velocity components.
    ic_func          = lambda x, y: np.exp(-30 * ((x - 0.25)**2 + (y - 0.25)**2))                                                       # Initial injection Gaussian blob.
    emitter_source   = lambda x, y, t=0, coef=None: 3.0 * np.exp(-30 * ((x - 0.25)**2 + (y - 0.25)**2))                                 # Continuous chimney emitter source F(x, y, t).
    
    # For Adv-Diff: u_t = -vx*u_x - vy*u_y + k*(u_xx + u_yy).
    # The operator vector [D, E, A, B, C, F] is [-vx, -vy, 2k, 0, 2k, 0] (since A=2k*u_xx, C=2k*u_yy).
    custom_operator = [-v_x, -v_y, 2.0 * v_diff, 0.0, 2.0 * v_diff, 0.0]                                                                # Custom Advection-Diffusion Operator.
    pde             = mgfd.PDE(operator=custom_operator, source=emitter_source, ic=ic_func, order=1)                                    # Instantiate fully generalized 1st order PDE.

    t_span          = (0.0, 1.5)                                                                                                        # Physical time span domain.

    # ------------------------------------------------------------------------------------------------------------------------------ #
    # Strategy A: CPU Solver
    # ------------------------------------------------------------------------------------------------------------------------------ #
    print("\n--- Strategy A: CPU Execution (Pre-factorized Direct Solver) ---")                                                         # Log Strategy A header.
    solver_cpu = mgfd.Solver(domain, pde, device="cpu", cfl=0.5, implicit=True, upwind=False, verbose=True)                             # Instantiate CPU Solver.
    result_cpu = solver_cpu.solve(t_span=t_span)                                                                                        # Execute CPU solver.
    print(f"Strategy A completed in {result_cpu.compute_time:.4f} seconds!")                                                            # Log CPU duration.

    # ------------------------------------------------------------------------------------------------------------------------------ #
    # Strategy B: GPU CUDA Solver
    # ------------------------------------------------------------------------------------------------------------------------------ #
    print("\n--- Strategy B: GPU Execution (CUDA Accelerated) ---")                                                                     # Log Strategy B header.
    try:                                                                                                                                # Attempt CUDA solver instantiation.
        solver_cuda = mgfd.Solver(domain, pde, device="cuda", cfl=0.5, implicit=True, upwind=False, verbose=True)                       # Instantiate CUDA Solver.
        result_cuda = solver_cuda.solve(t_span=t_span)                                                                                  # Execute CUDA solver.
        print(f"Strategy B completed in {result_cuda.compute_time:.4f} seconds!")                                                       # Log CUDA duration.
    except ImportError:                                                                                                                 # Catch missing CuPy or CUDA hardware.
        print("CuPy is not installed or CUDA GPU is not available. Skipping Strategy B.")                                               # Soft fail notice.

    # 4. Render Solution Plot
    print("\nStep 4: Rendering Advection-Diffusion animation...")                                                                       # Log step 4.
    result_cpu.plot(save=False, show=True, filename='04_optimizers_result', title="mGFD Solution: Advection-Diffusion", t_span=t_span)  # Render solution plot.

    # Clean up temporary CSV files
    if os.path.exists(contour_file): os.remove(contour_file)                                                                            # Clean contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                # Clean cloud file.

if __name__ == '__main__':                                                                                                              # Entry point guard.
    main()                                                                                                                              # Run main.
