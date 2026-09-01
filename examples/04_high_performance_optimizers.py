"""
Example 04: High Performance & Optimizers using mGFD

Overview:
    This tutorial demonstrates how to use mGFD's advanced optimizers to solve highly complex,\n    large-scale PDEs efficiently. We simulate an Advection-Diffusion Equation:\n        u_t = v(u_xx + u_yy) - (v_x u_x + v_y u_y)\n    We highlight Upwind Schemes, Memory-efficient Solvers, and Iterative Preconditioners on an irregular Kite domain.

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
from mGFD.exceptions import ParameterError                                                                                              # Import exception handling.
from mGFD.cloud_generator.core.generator import generate_cloud_natural                                                                  # Import the point cloud generator.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 1. Define the Advection-Diffusion Physics
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 1: Defining Advection-Diffusion Physics...")                                                                                # Log step 1.
# The PDE Operator: L u = A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u
# Advection-Diffusion: u_t = v(u_xx + u_yy) - (v_x u_x + v_y u_y) + f
# A = v, C = v, D = -v_x, E = -v_y.
def create_adv_diff_operator(v, v_x, v_y):                                                                                              # Function to generate the correct operator.
    return np.vstack([[-v_x], [-v_y], [2 * v], [0], [2 * v], [0]])                                                                      # Returns exactly [D, E, 2A, B, 2C, F].

v = 0.05                                                                                                                                # Diffusion coefficient.
v_x, v_y = 0.2, 0.2                                                                                                                     # Advection velocity components (strong wind blowing diagonally).
L = create_adv_diff_operator(v, v_x, v_y)                                                                                               # Compute operator matrix.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 2. Build a Dense Point Cloud (Irregular Kite Shape)
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 2: Generating the irregular point cloud...")                                                                                # Log step 2.

# We will use an irregular Kite geometry.
kite_contour = [
    (1.0, 0.5), (0.85, 0.7), (0.5, 1.0), (0.15, 0.7),
    (0.0, 0.5), (0.15, 0.3), (0.5, 0.0), (0.85, 0.3)
]                                                                                                                                       # Define Kite boundary.

contour_file = 'kite_contour.csv'                                                                                                       # Define contour filename.
with open(contour_file, 'w', newline='') as f:                                                                                          # Open file for writing.
    writer = csv.writer(f)                                                                                                              # Initialize CSV writer.
    writer.writerow(['x', 'y'])                                                                                                         # Write CSV header.
    for pt in kite_contour:                                                                                                             # Iterate over the Kite coordinates.
        writer.writerow(pt)                                                                                                             # Write each point to the file.

cloud_file = 'kite_cloud.csv'                                                                                                           # Define output cloud filename.
# We use a very dense spacing to demonstrate large-scale PDE solving.
generate_cloud_natural(contour_file, cloud_file, cloud_size=0.015, save=False)                                                          # Call generator.

df = pd.read_csv(cloud_file)                                                                                                            # Read the generated cloud into a DataFrame.
df['flag'] = df['classification'].map({'boundary': 1, 'interior': 0})                                                                   # Map the generator's string labels to mGFD numeric flags.
p = df[['x', 'y', 'flag']].to_numpy(dtype=np.float64)                                                                                   # Convert the DataFrame directly to a Numpy array.

X = p[:, 0]                                                                                                                             # Extract the X coordinates vector.
Y = p[:, 1]                                                                                                                             # Extract the Y coordinates vector.
print(f"Generated dense point cloud with {len(p)} nodes.")                                                                              # Print the total number of nodes generated.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 3. Initial and Boundary Conditions
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 3: Creating boundary and initial conditions...")                                                                            # Log step 3.
time_steps = 500                                                                                                                        # Define number of time steps.
T_arr = np.linspace(0, 1, time_steps)                                                                                                   # Create time array.

f_arr = np.zeros((len(p), time_steps))                                                                                                  # Initialize forcing array.
r2 = (X - 0.25)**2 + (Y - 0.25)**2                                                                                                      # Calculate distance from origin injection point.
f_arr[:, 0] = np.exp(-100 * r2)                                                                                                         # Inject a high-concentration pollutant blob as the initial condition.

for k in range(1, time_steps):                                                                                                          # Iterate over each time step.
    f_arr[:, k] = 0.0                                                                                                                   # Apply zero Dirichlet boundary conditions (free outflow).

force_df = pd.DataFrame(f_arr)                                                                                                          # Convert the numpy matrix into a Pandas DataFrame.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 4. Solvers & Optimizers Comparison
# ------------------------------------------------------------------------------------------------------------------------------ #
print("Step 4: Solvers and Optimizers Comparison")                                                                                      # Log step 4.
print("\n--- Strategy A: CPU Execution (Pre-factorized Direct Solver) ---")                                                             # Separator log.

# Strategy A uses CPU direct sparse pre-factorized LU decomposition which provides maximum precision and speed.
# We use `upwind=True` which biases the neighbor search in the opposite direction of the velocity vector (v_x, v_y),
# heavily stabilizing the convective terms.
try:                                                                                                                                    # Attempt CPU execution.
    result_A = TimeDerivative1(                                                                                                         # Call the time derivative solver.
        p, force_df, time_steps, coef=[v, v_x, v_y], operator=L,                                                                        # Pass fundamental math params.
        nvec=15,                                                                                                                        # Specify neighbor count.
        implicit=True,                                                                                                                  # Guarantee numerical stability with fully implicit time stepping.
        upwind=True,                                                                                                                    # Apply upwind stabilization for convective fluxes.
        device="cpu",                                                                                                                   # High-performance CPU direct solver.
        verbose=True                                                                                                                    # Output solver logs.
    )                                                                                                                                   # End of solver call.
    print(f"Strategy A completed in {result_A.compute_time:.4f} seconds.")                                                              # Log Strategy A completion time.
except Exception as e:                                                                                                                  # Catch any numerical divergence issues.
    print(f"Strategy A failed: {e}")                                                                                                    # Log Strategy A failure reason.

print("\n--- Strategy B: GPU Accelerated (CUDA) ---")                                                                                   # Separator log.
# Strategy B uses GPU acceleration via CuPy pre-factorized direct sparse solver.
try:                                                                                                                                    # Attempt GPU execution via CuPy.
    result_B = TimeDerivative1(                                                                                                         # Call the time derivative solver.
        p, force_df, time_steps, coef=[v, v_x, v_y], operator=L,                                                                        # Pass fundamental math params.
        nvec=15,                                                                                                                        # Specify neighbor count.
        implicit=True,                                                                                                                  # Guarantee numerical stability.
        upwind=True,                                                                                                                    # Apply upwind stabilization.
        device='cuda',                                                                                                                  # Exploit GPU parallelism via CuPy!
        verbose=True                                                                                                                    # Output solver logs.
    )                                                                                                                                   # End of solver call.
    print(f"Strategy B completed in {result_B.compute_time:.4f} seconds.")                                                              # Log Strategy B completion time.
except ImportError:                                                                                                                     # Catch missing CuPy library on non-CUDA systems.
    print("CuPy is not installed or no NVIDIA GPU detected. Skipping Strategy B.")                                                      # Soft fail and notify user.
except ParameterError as e:                                                                                                             # Catch mGFD strict parameter validation errors.
    print(f"mGFD Configuration Error: {e}")                                                                                             # Log configuration error.

# ------------------------------------------------------------------------------------------------------------------------------ #
# 5. Visualize Strategy A
# ------------------------------------------------------------------------------------------------------------------------------ #
print("\nStep 5: Rendering Advection-Diffusion animation...")                                                                           # Log step 5.
plot_transient(p, result_A.solution, save=True, nom='04_optimizers_result', title="mGFD Solution: Advection-Diffusion")                 # Render and save plot.

if os.path.exists(contour_file): os.remove(contour_file)                                                                                # Clean up temporary contour CSV.
if os.path.exists(cloud_file): os.remove(cloud_file)                                                                                    # Clean up temporary cloud CSV.
