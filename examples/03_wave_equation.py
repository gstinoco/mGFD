"""
Example 03: Solving the Wave Equation using mGFD

Overview:
    This tutorial demonstrates how to solve a hyperbolic PDE (the Wave Equation):\n        u_tt = c^2 (u_xx + u_yy)\n    On a circular domain, using the Meshless Generalized Finite Difference (mGFD) method.

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
    August, 2026.
Last Modification:
    August, 2026.
"""

import os                                                                                                                               # Standard OS imports.
import csv                                                                                                                              # Standard CSV imports.
import numpy as np                                                                                                                      # Numpy for arrays.
import pandas as pd                                                                                                                     # Pandas for dataframes.

from PIL import ExifTags

from mGFD import TimeDerivative2                                                                                                        # Import the second-order time derivative solver.
from mGFD.viz.graph import plot_transient                                                                                               # Import the transient plotting tool.
from mGFD.cloud_generator.core.generator import generate_cloud_natural                                                                  # Import the point cloud generator.
from mGFD.cloud_generator.core.generator import generate_cloud_regular

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Problem Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining the Problem Physics...")                                                                                        # Log step 1.
# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# For the Wave Equation (u_tt = c^2(u_xx + u_yy)), A = c^2, C = c^2.
def create_wave_operator(c):                                                                                                            # Function to generate the correct operator.
    return np.vstack([[0], [0], [2 * c**2], [0], [2 * c**2], [0]])                                                                      # Returns exactly [D, E, 2A, B, 2C, F].

c = 0.8                                                                                                                                 # Define wave speed propagation constant.
L = create_wave_operator(c)                                                                                                             # Compute operator matrix.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 2. Build the Point Cloud
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 2: Building the irregular point cloud...")                                                                                  # Log step 2.
# We will solve the Wave Equation on a beautiful STAR-shaped irregular domain to show off mGFD.
# We use time_steps = 400 to show the reflections before long-term numerical dispersion degrades the unstructured simulation.
theta = np.linspace(0, 2 * np.pi, 150)
r_star = 0.5 + 0.15 * np.sin(5 * theta)
star_contour = [(0.5 + r_star[i] * np.cos(theta[i]), 0.5 + r_star[i] * np.sin(theta[i])) for i in range(len(theta))]

contour_file = 'star_contour_wave.csv'                                                                                                  # Define contour filename.
with open(contour_file, 'w', newline='') as f:                                                                                          # Open file for writing.
    writer = csv.writer(f)                                                                                                              # Initialize CSV writer.
    writer.writerow(['x', 'y'])                                                                                                         # Write CSV header.
    for pt in star_contour:                                                                                                             # Iterate over the star coordinates.
        writer.writerow(pt)                                                                                                             # Write each point to the file.

cloud_file = 'star_cloud_wave.csv'                                                                                                      # Define output cloud filename.
print(f"Generating Regular Cloud on Star domain (spacing = 0.02)...")                                                                   # Output progress to user.
generate_cloud_regular(contour_file, cloud_file, cloud_size=0.02, save=True)                                                            # Call the mGFD Poisson Disk generator.

df = pd.read_csv(cloud_file)                                                                                                            # Read the generated cloud into a DataFrame.
df['flag'] = df['classification'].map({'boundary': 1, 'interior': 0})                                                                   # Map the generator's string labels to mGFD numeric flags.
p = df[['x', 'y', 'flag']].to_numpy(dtype=np.float64)                                                                                   # Convert the DataFrame directly to a Numpy array.

X = p[:, 0]                                                                                                                             # Extract the X coordinates vector.
Y = p[:, 1]                                                                                                                             # Extract the Y coordinates vector.
print(f"Successfully generated point cloud with {len(p)} nodes.")                                                                       # Print the total number of nodes generated.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 3. Create Boundary & Initial Conditions Data (Using Pandas DataFrame)
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Creating boundary and initial conditions using Pandas...")                                                               # Log step 3.
# For the Wave equation (TimeDerivative2), we need to specify:
# 1. The primary values matrix 'f' (initial position + boundary conditions)
# 2. The secondary values matrix 'g' (initial velocity + boundary conditions velocity)

time_steps = 4000                                                                                                                       # Total number of time steps.
f_arr = np.zeros((len(p), time_steps))                                                                                                  # Initialize the primary forcing matrix array.
for k in range(time_steps):                                                                                                             # Iterate over each time step.
    if k == 0:                                                                                                                          # Check if evaluating the initial condition step.
        r2 = (X - 0.5)**2 + (Y - 0.5)**2                                                                                                # Calculate squared distance from center.
        f_arr[:, 0] = np.exp(-100 * r2)                                                                                                 # Apply Gaussian pulse (ripple) at the center of the domain.
    else:                                                                                                                               # Evaluating boundary conditions for subsequent steps.
        f_arr[:, k] = 0.0                                                                                                               # Fix boundaries to 0 to simulate rigid walls.

g_arr = np.zeros(len(p))                                                                                                                # Initial velocity is perfectly zero everywhere.

force_f = pd.DataFrame(f_arr)                                                                                                           # Convert 'f' to a Pandas DataFrame.
force_g = pd.Series(g_arr.tolist())                                                                                                     # Convert 'g' to a Pandas Series.
print("Prepared initial conditions and velocities.")                                                                                    # Log success.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Solve the Equation
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Solving the Wave equation with mGFD...")                                                                                 # Log step 4.

# The TimeDerivative2 solver calculates the second-order temporal evolution.
# Notice we are passing force_f and force_g concurrently to the solver.
result = TimeDerivative2(p, force_f, force_g, time_steps, [c], operator=L, nvec=12, verbose=True, implicit=True, lam=0.25)              # Execute the solver with unconditionally stable energy-conserving Newmark-beta.

u_approx = result.solution                                                                                                              # Extract the computed spatiotemporal solution array.
u_approx_vis = np.clip(u_approx, -2.0, 2.0)                                                                                             # Clip the focal singularity for visualization to prevent Z-axis squashing.
compute_time = result.compute_time                                                                                                      # Extract the computation elapsed time.
print(f"Transient solution computed in {compute_time:.4f} seconds!")                                                                    # Report total solving time.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 5. Visualize the Result using mGFD's native renderer
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 5: Rendering the animation...")                                                                                             # Log step 5.
plot_transient(p, u_approx_vis, save=False, nom='03_wave_result', title=f"mGFD Solution: Wave Equation (t={time_steps})")               # Render and save animation.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.
if os.path.exists(cloud_file.replace('.csv', '.png')): os.remove(cloud_file.replace('.csv', '.png'))                                    # Clean up generated cloud PNG.
if os.path.exists(cloud_file.replace('.csv', '.svg')): os.remove(cloud_file.replace('.csv', '.svg'))                                    # Clean up generated cloud EPS.