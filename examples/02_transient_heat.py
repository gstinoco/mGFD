"""
Example 02: Solving the Transient Heat Equation using mGFD

Overview:
    This tutorial demonstrates how to solve the classic Heat Equation (First-order in time):\n        u_t - v(u_xx + u_yy) = f(x, y, t)\n    On an irregular 2D domain (Star), using the Meshless Generalized Finite Difference (mGFD) method.

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

from mGFD import TimeDerivative1                                                                                                        # Import the first-order time derivative solver.
from mGFD.viz.graph import plot_transient                                                                                               # Import the transient plotting tool.
from mGFD.cloud_generator.core.generator import generate_cloud_natural                                                                  # Import the point cloud generator.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Problem Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining the Problem Physics...")                                                                                        # Log step 1.

# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# For the Heat Equation (u_t = v(u_xx + u_yy) + f), A = v, C = v.
def create_heat_operator(v):                                                                                                            # Function to generate the correct operator.
    return np.vstack([[0], [0], [2 * v], [0], [2 * v], [0]])                                                                            # Returns exactly [D, E, 2A, B, 2C, F].

v = 0.08                                                                                                                                 # Define thermal diffusivity constant.
L = create_heat_operator(v)                                                                                                             # Compute operator matrix.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 2. Build the Point Cloud using mGFD's Cloud Generator
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 2: Building the irregular point cloud...")                                                                                  # Log step 2.

# To showcase the power of Meshless methods, we will solve the PDE on an irregular star-shaped domain.
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

# In mGFD, flag = 0 indicates an interior node, flag = 1 indicates a boundary node.
df         = pd.read_csv(cloud_file)                                                                                                    # Read the generated cloud into a DataFrame.
df['flag'] = df['classification'].map({'boundary': 1, 'interior': 0})                                                                   # Map the generator's string labels to mGFD numeric flags.
p          = df[['x', 'y', 'flag']].to_numpy(dtype=np.float64)                                                                          # Convert the DataFrame directly to a Numpy array.

print(f"Successfully generated point cloud with {len(p)} nodes.")                                                                       # Print the total number of nodes generated.


# ------------------------------------------------------------------------------------------------------------------------------ #
# 3. Create Boundary & Initial Conditions Data (Using Pandas DataFrame)
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Creating boundary and initial conditions using Pandas...")                                                               # Log step 3.

# The TimeDerivative1 solver expects a spatiotemporal matrix 'f' for its right-hand side.
# Rows = spatial nodes, Columns = time steps.
time_steps = 200                                                                                                                        # Define number of time steps.
T_arr      = np.linspace(0, 1, time_steps)                                                                                              # Create an array representing time from t=0 to t=1.
f_arr      = np.zeros((len(p), time_steps))                                                                                             # Initialize the forcing matrix array.

for k in range(time_steps):                                                                                                             # Iterate over each time step.
    t_k         = T_arr[k]                                                                                                              # Extract the current time step value.
    f_arr[:, k] = np.exp(-2 * np.pi**2 * v * t_k) * np.sin(np.pi * p[:, 0]) * np.sin(np.pi * p[:, 1])                                   # Evaluate analytical solution.

force_df   = pd.DataFrame(f_arr)                                                                                                        # Convert the numpy matrix into a Pandas DataFrame for mGFD natively.

print("Prepared initial/boundary conditions as a Pandas DataFrame.")                                                                    # Report successful data formatting.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Solve the PDE
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Solving the Heat equation with mGFD...")                                                                                 # Log step 4.

# We pass the cloud (p), the spatiotemporal DataFrame (force_df), and the number of steps.
result       = TimeDerivative1(p, force_df, time_steps, [v], operator=L, nvec=15, implicit=True, lam=0.5, verbose=True)                 # Execute the transient solver.

u_approx     = result.solution                                                                                                          # Extract the computed spatiotemporal solution array.
compute_time = result.compute_time                                                                                                      # Extract the computation elapsed time.

print(f"Transient solution computed in {compute_time:.4f} seconds!")                                                                    # Report total solving time.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 5. Visualize the Result using mGFD's native renderer
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 5: Rendering the animation...")                                                                                             # Log step 5.

# The plot_transient module computes the Delaunay triangulation and generates an animated visualization 
plot_transient(p, u_approx, save=False, nom='02_heat_result', title=f"mGFD Solution: Heat Equation (t={time_steps})")                   # Render and save animation.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.
