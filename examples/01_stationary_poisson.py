"""
Example 01: Solving the Stationary Poisson Equation using mGFD

Overview:
    This tutorial demonstrates how to solve the classic Poisson Equation:\n        u_xx + u_yy = f(x, y)\n    On an irregular 2D domain (Star), using the Meshless Generalized Finite Difference (mGFD) method.

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

from mGFD import Stationary                                                                                                             # Import the stationary solver.
from mGFD.viz.graph import plot_stationary                                                                                              # Import the stationary plotting tool.
from mGFD.cloud_generator.core.generator import generate_cloud_natural                                                                  # Import the point cloud generator.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Problem Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining the Problem Physics...")                                                                                        # Log step 1.

# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# For Poisson (u_xx + u_yy = f), A=1, C=1.
# The mGFD library requires a 6-element operator vector defined exactly as: [D, E, 2A, B, 2C, F]
L = np.vstack([[0], [0], [2], [0], [2], [0]])                                                                                           # Construct the operator vector.

def force(x, y):                                                                                                                        # Define the forcing term function.
    return -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)                                                                      # Evaluate the analytical force.

def boundary(x, y):                                                                                                                     # Define the boundary condition function.
    return np.sin(np.pi * x) * np.sin(np.pi * y)                                                                                        # Analytical solution used for boundaries.

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
# 3. Solve the PDE
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Solving the Poisson equation with mGFD...")                                                                              # Log step 3.

# The Stationary solver handles everything: KDTree neighbor search, sparse matrix assembly, and the linear solve.
# We pass the cloud (p), the boundary function (boundary), the forcing function (force), and the operator (L).
result       = Stationary(p, boundary, force, operator=L, nvec=15, verbose=True)                                                        # Execute the stationary solver.

# The solver returns a SolverResult object containing the approximate solution and the neighbor map.
u_approx     = result.solution                                                                                                          # Extract the computed solution values.
compute_time = result.compute_time                                                                                                      # Extract the computation elapsed time.
print(f"Solution found in {compute_time:.4f} seconds!")                                                                                 # Report total solving time.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Visualize the Result using mGFD's native renderer
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Rendering the output...")                                                                                                # Log step 4.
# mGFD includes an advanced 3D visualization module that automatically computes Delaunay triangulations!
plot_stationary(p, u_approx, save=False, nom='01_poisson_result', title="mGFD Solution: Poisson on Star Domain")                         # Render and save plot.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.
