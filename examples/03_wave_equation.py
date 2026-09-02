"""
Example 03: Solving the Wave Equation using mGFD

Overview:
    This tutorial demonstrates how to solve a hyperbolic PDE (the Wave Equation):
        u_tt = c^2 (u_xx + u_yy) + F(x, y, t)
    On a circular/star domain, using the Meshless Generalized Finite Difference (mGFD) method.
    Highlights independent acoustic driving source terms F(x, y, t), custom physical time spans t_span=(0, 3.0),
    and energy-conserving implicit Newmark-beta.

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

from mGFD import TimeDerivative2                                                                                                        # Import the second-order time derivative solver.
from mGFD.viz.graph import plot_transient                                                                                               # Import the transient plotting tool.
from mGFD.cloud_generator.core.generator import generate_cloud_regular                                                                 # Import the regular point cloud generator.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Problem Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining the Problem Physics...")                                                                                        # Log step 1.

# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# For the Wave Equation (u_tt = c^2(u_xx + u_yy) + F), A = c^2, C = c^2.
def create_wave_operator(c):                                                                                                            # Function to generate the correct operator.
    return np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                                                      # Returns exactly [D, E, 2A, B, 2C, F].

c = 0.8                                                                                                                                 # Define wave speed propagation constant.
L = create_wave_operator(c)                                                                                                             # Compute operator matrix.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 2. Build the Point Cloud
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 2: Building the irregular point cloud...")                                                                                  # Log step 2.

theta = np.linspace(0, 2 * np.pi, 150)                                                                                                  # Parametric angular coordinates.
r_star = 0.5 + 0.15 * np.sin(5 * theta)                                                                                                 # Polar radius of star domain.
star_contour = [(0.5 + r_star[i] * np.cos(theta[i]), 0.5 + r_star[i] * np.sin(theta[i])) for i in range(len(theta))]                    # Compute (x, y) contour boundary.

contour_file = 'star_contour_wave.csv'                                                                                                  # Define contour filename.
with open(contour_file, 'w', newline='') as f:                                                                                          # Open file for writing.
    writer = csv.writer(f)                                                                                                              # Initialize CSV writer.
    writer.writerow(['x', 'y'])                                                                                                         # Write CSV header.
    for pt in star_contour:                                                                                                             # Iterate over the star coordinates.
        writer.writerow(pt)                                                                                                             # Write each point to the file.

cloud_file = 'star_cloud_wave.csv'                                                                                                      # Define output cloud filename.
print(f"Generating Regular Cloud on Star domain (spacing = 0.02)...")                                                                   # Output progress to user.
generate_cloud_regular(contour_file, cloud_file, cloud_size=0.02, save=False)                                                           # Call the mGFD point cloud generator.

df = pd.read_csv(cloud_file)                                                                                                            # Read the generated cloud into a DataFrame.
df['flag'] = df['classification'].map({'boundary': 1, 'interior': 0})                                                                   # Map the generator's string labels to mGFD numeric flags.
p = df[['x', 'y', 'flag']].to_numpy(dtype=np.float64)                                                                                   # Convert the DataFrame directly to a Numpy array.

print(f"Successfully generated point cloud with {len(p)} nodes.")                                                                       # Print the total number of nodes generated.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 3. Create Initial Position, Initial Velocity, Boundary, and Acoustic Source Functions
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Creating initial position, velocity, boundary, and acoustic driving source functions...")                                # Log step 3.

def initial_position(x, y):                                                                                                             # Initial Gaussian pulse centered at (0.5, 0.5).
    r2 = (x - 0.5)**2 + (y - 0.5)**2                                                                                                    # Distance squared from center.
    return np.exp(-100 * r2)                                                                                                            # Gaussian bump.

def initial_velocity(x, y):                                                                                                             # Initial velocity u_t(x, y, 0).
    return 0.0                                                                                                                          # Starts from rest.

def boundary_condition(x, y, t, coef):                                                                                                  # Fixed Dirichlet boundary condition.
    return 0.0                                                                                                                          # Rigid boundary walls.

def acoustic_source(x, y, t, coef):                                                                                                     # Independent oscillating acoustic source term F(x, y, t).
    r2 = (x - 0.5)**2 + (y - 0.5)**2                                                                                                    # Distance squared from center speaker.
    return 5.0 * np.cos(6 * np.pi * t) * np.exp(-80 * r2)                                                                               # Oscillating acoustic wave driver.

print("Prepared initial pulse, zero velocity, boundary, and acoustic driver source F(x, y, t) functions.")                              # Log success.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Solve the Wave Equation over Custom Physical Time Span t_span=(0, 3.0)
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Solving the Wave equation with mGFD (t_span=(0, 3.0), CFL Auto-Dt)...")                                                 # Log step 4.

# We pass physical time domain t_span=(0, 3.0), independent source term acoustic_source, and implicit Newmark-beta.
result = TimeDerivative2(                                                                                                               # Call second-order transient solver.
    p, boundary_condition, g=initial_velocity, ic=initial_position, source=acoustic_source,                                             # Pass physical callables & independent source.
    t_span=(0.0, 3.0), coef=[c], operator=L, nvec=12, implicit=True, lam=0.25, cfl=0.5, verbose=True                                    # Physical t_span & automatic CFL step count.
)                                                                                                                                       # Execute solver call.

u_approx     = result.solution                                                                                                          # Extract the computed spatiotemporal solution array.
u_approx_vis = np.clip(u_approx, -2.0, 2.0)                                                                                             # Clip focal singularity for visualization.
compute_time = result.compute_time                                                                                                      # Extract the computation elapsed time.

print(f"Transient solution computed in {compute_time:.4f} seconds (dt = {result.dt:.6f}, steps = {result.t_steps})!")                    # Report total solving time.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 5. Visualize the Result using mGFD's native renderer
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 5: Rendering the animation...")                                                                                             # Log step 5.
plot_transient(p, u_approx_vis, save=False, nom='03_wave_result', title=f"mGFD Solution: Wave Equation with Source (t_span=[0, 3])")     # Render and save animation.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.