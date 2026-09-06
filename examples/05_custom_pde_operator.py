"""
Example 05: Arbitrary Custom PDE — Solving Custom Physics using mGFD (OOP API)

Overview:
    This tutorial demonstrates the ultimate generalized power of mGFD.
    Instead of using standard helper classes (Poisson, Heat, Wave), we will formulate
    a completely arbitrary Anisotropic Advection-Diffusion-Reaction Equation:
        u_t = -v_x u_x - v_y u_y + k_x u_xx + k_y u_yy + k_xy u_xy - R u
    
    Using the base `mgfd.PDE` class, we can directly inject any mathematical operator
    vector to solve novel physics without altering the core library.

Public API:
    main                                    Main execution routine for the custom PDE tutorial.

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
    Executes Example 05 tutorial solving an arbitrary anisotropic advection-diffusion-reaction PDE.

    Input:
        None

    Output:
        None
    """
    print("=================================================================================")                                          # Separator log.
    print("                Example 05: Arbitrary Custom PDE (OOP API)                       ")                                          # Title log.
    print("=================================================================================\n")                                        # Separator log.

    # 1. Define geometry contour & generate point cloud via mGFD.Cloud
    print("Step 1: Building irregular point cloud geometry...")                                                                         # Log step 1.
    hexagon_contour = [
        (0.866, 0.5), (0.433, 1.0), (-0.433, 1.0), 
        (-0.866, 0.5), (-0.433, 0.0), (0.433, 0.0)
    ]                                                                                                                                   # Define Hexagon boundary vertices.

    contour_file = 'hex_contour.csv'                                                                                                    # Contour filename.
    cloud_file   = 'hex_cloud.csv'                                                                                                      # Output cloud filename.

    with open(contour_file, 'w', newline='') as f:                                                                                      # Write CSV contour.
        writer = csv.writer(f)                                                                                                          # Initialize CSV writer.
        writer.writerow(['x', 'y'])                                                                                                     # Write header.
        for pt in hexagon_contour: writer.writerow(pt)                                                                                  # Write vertices.

    cloud = mgfd.Cloud.generate_natural(contour_file, cloud_file, cloud_size=0.03, save=False, show=False)                              # High-level point cloud generation.
    print(f"Generated {cloud}: {cloud.num_nodes} nodes.")                                                                               # Log cloud summary.

    # 2. Set Dirichlet boundary conditions and construct Domain
    print("\nStep 2: Binding Dirichlet boundary condition and domain...")                                                               # Log step 2.
    domain = cloud.set_boundary(mgfd.Dirichlet(0.0))                                                                                    # Bind zero Dirichlet boundary condition.

    # 3. Define Custom Physics & Instantiate Solver
    print("\nStep 3: Formulating completely arbitrary anisotropic PDE physics...")                                                      # Log step 3.
    
    # Let's invent a weird physics equation (Anisotropic Diffusion + Advection + Reaction):
    # u_t = -0.5*u_x - 0.2*u_y + 0.1*u_xx + 0.01*u_yy + 0.05*u_xy - 2.0*u + F(x,y,t)
    # The operator vector [D, E, A, B, C, F] maps to [u_x, u_y, (1/2)u_xx, u_xy, (1/2)u_yy, 1]!
    # So to get 0.1*u_xx, we set A = 0.2. To get 0.01*u_yy, we set C = 0.02. To get 0.05*u_xy, B = 0.05.
    
    custom_operator = [
        -0.5,                                                                                                                           # D (Advection X)
        -0.2,                                                                                                                           # E (Advection Y)
         0.2,                                                                                                                           # A (Diffusion XX)
         0.05,                                                                                                                          # B (Mixed XY)
         0.02,                                                                                                                          # C (Diffusion YY)
        -2.0                                                                                                                            # F (Reaction Term)
    ]                                                                                                                                   # Custom operator array.
    
    ic_func     = lambda x, y: np.exp(-50 * ((x)**2 + (y - 0.5)**2))                                                                    # Initial condition blob.
    source_func = lambda x, y, t=0, coef=None: 10.0 * np.sin(10 * np.pi * t) * np.exp(-100 * (x**2 + (y-0.5)**2))                       # Pulsating source term.
    
    # We use the base `PDE` class and explicitly tell it this is a First Order (transient) PDE.
    pde    = mgfd.PDE(operator=custom_operator, ic=ic_func, source=source_func, order=1)                                                # Instantiate arbitrary PDE.
    
    solver = mgfd.Solver(domain, pde, nvec=15, cfl=0.3, verbose=True)                                                                   # Instantiate high-level Solver.

    # 4. Solve and Plot
    print("\nStep 4: Solving PDE and rendering results...")                                                                             # Log step 4.
    result = solver.solve(t_span=(0.0, 1.0))                                                                                            # Execute solver.
    
    print("\nStep 5: Visualizing custom PDE dynamics...")                                                                               # Log step 5.
    result.plot(save=False, show=True, title="Arbitrary Custom PDE (mGFD Generalized)")                                                 # Plot 3D interactive surface.

    # Clean up temporary CSV files
    if os.path.exists(contour_file): os.remove(contour_file)                                                                            # Clean contour file.
    if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                # Clean cloud file.

if __name__ == '__main__':                                                                                                              # Main execution block.
    main()                                                                                                                              # Call main routine.
