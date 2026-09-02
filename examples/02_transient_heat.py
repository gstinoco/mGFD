"""
Example 02: Solving the Transient Heat Equation using mGFD

Overview:
    This tutorial demonstrates how to solve the Heat Equation (First-order in time):
        u_t - v(u_xx + u_yy) = F(x, y, t)
    On an irregular 2D domain (Star), using the Meshless Generalized Finite Difference (mGFD) method.
    Highlights independent spatiotemporal source terms F(x, y, t), custom physical time spans t_span=(0, 2.0),
    and automatic CFL time step estimation.

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
    May, 2024.
Last Modification:
    September, 2026.
"""

import os                                                                                                                               # Standard OS imports.
import csv                                                                                                                              # Standard CSV imports.
import numpy as np                                                                                                                      # Numpy for arrays.
import pandas as pd                                                                                                                     # Pandas for dataframes.

from mGFD import TimeDerivative1                                                                                                        # Import the first-order time derivative solver.
from mGFD.viz.graph import plot_transient                                                                                               # Import the transient plotting tool.
from mGFD.cloud_generator.core.generator import generate_cloud_natural                                                                  # Import the point cloud generator.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Problem Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining the Problem Physics...")                                                                                        # Log step 1.

# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# For the Heat Equation (u_t = v(u_xx + u_yy) + F), A = v, C = v.
def create_heat_operator(v):                                                                                                            # Function to generate the correct operator.
    return np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                            # Returns exactly [D, E, 2A, B, 2C, F].

v = 0.08                                                                                                                                 # Define thermal diffusivity constant.
L = create_heat_operator(v)                                                                                                             # Compute operator matrix.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 2. Build the Point Cloud using mGFD's Cloud Generator
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 2: Building the irregular point cloud...")                                                                                  # Log step 2.

star_contour = [
    (0.5, 1.0), (0.65, 0.65), (1.0, 0.5), (0.65, 0.35), 
    (0.5, 0.0), (0.35, 0.35), (0.0, 0.5), (0.35, 0.65)
]                                                                                                                                       # Define coordinates of the irregular 4-pointed star.

contour_file = 'star_contour.csv'                                                                                                       # Define contour filename.
with open(contour_file, 'w', newline='') as f:                                                                                          # Open file for writing.
    writer = csv.writer(f)                                                                                                              # Initialize CSV writer.
    writer.writerow(['x', 'y'])                                                                                                         # Write CSV header.
    for pt in star_contour:                                                                                                             # Iterate over the star coordinates.
        writer.writerow(pt)                                                                                                             # Write each point to the file.

cloud_file = 'star_cloud.csv'                                                                                                           # Define output cloud filename.
print(f"Generating Natural Cloud on Star domain (spacing = 0.03)...")                                                                   # Output progress to user.

generate_cloud_natural(contour_file, cloud_file, cloud_size=0.03, save=False)                                                           # Call the mGFD Poisson Disk generator.

df         = pd.read_csv(cloud_file)                                                                                                    # Read the generated cloud into a DataFrame.
df['flag'] = df['classification'].map({'boundary': 1, 'interior': 0})                                                                   # Map the generator's string labels to mGFD numeric flags.
p          = df[['x', 'y', 'flag']].to_numpy(dtype=np.float64)                                                                          # Convert the DataFrame directly to a Numpy array.

print(f"Successfully generated point cloud with {len(p)} nodes.")                                                                       # Print the total number of nodes generated.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 3. Create Initial, Boundary, and Independent Pulsing Source Term Functions
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Creating initial, boundary, and independent heat source term functions...")                                              # Log step 3.

def initial_condition(x, y):                                                                                                            # Initial state function at t=0.
    return np.sin(np.pi * x) * np.sin(np.pi * y)                                                                                        # Initial heat profile.

def boundary_condition(x, y, t, coef):                                                                                                  # Spatiotemporal boundary function.
    diff_v = coef[0]                                                                                                                    # Diffusivity coefficient.
    return np.exp(-2 * np.pi**2 * diff_v * t) * np.sin(np.pi * x) * np.sin(np.pi * y)                                                   # Zero boundary cooling profile.

def heat_source(x, y, t, coef):                                                                                                         # Independent pulsed heat laser source term F(x, y, t).
    r2 = (x - 0.5)**2 + (y - 0.5)**2                                                                                                    # Distance squared from laser center.
    return 10.0 * np.sin(4 * np.pi * t) * np.exp(-50 * r2)                                                                              # Pulsing laser source at center.

print("Prepared initial, boundary, and independent heat source F(x, y, t) functions.")                                                  # Report successful data formatting.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Solve the PDE over a Custom Time Span t_span=(0, 2.0) with Automatic CFL dt
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Solving the Heat equation with mGFD (t_span=(0, 2.0), CFL Auto-Dt)...")                                                 # Log step 4.

# We specify custom physical time span t_span=(0, 2.0) and pass independent source term F(x, y, t).
result       = TimeDerivative1(                                                                                                         # Call first-order transient solver.
    p, boundary_condition, ic=initial_condition, source=heat_source,                                                                    # Pass physical callables & independent source.
    t_span=(0.0, 2.0), coef=[v], operator=L, nvec=15, implicit=True, lam=0.5, cfl=0.5, verbose=True                                     # Custom t_span & automatic CFL step count.
)                                                                                                                                       # Execute solver call.

u_approx     = result.solution                                                                                                          # Extract the computed spatiotemporal solution array.
compute_time = result.compute_time                                                                                                      # Extract the computation elapsed time.

print(f"Transient solution computed in {compute_time:.4f} seconds (dt = {result.dt:.6f}, t_steps = {result.t_steps})!")                 # Report total solving time and dt.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 5. Visualize the Result using mGFD's native renderer
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 5: Rendering the animation...")                                                                                             # Log step 5.
plot_transient(p, u_approx, save=False, nom='02_heat_result', title=f"mGFD Solution: Heat Equation with Source (t_span=[0, 2])")        # Render and save animation.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.
