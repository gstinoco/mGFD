"""
Meshless Generalized Finite Differences to solve Advection-Diffusion Equation on irregular regions.

This script contains all the requirements to:
    - Read the files with the data of all the regions in the Data folder.
    - State the conditions for the problem.
    - Solve the problem using a meshless Generalized Finite Difference approach.
    - Save the results.
    - Plot the results.

All the codes presented below were developed by:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

With the funding of:
    National Council of Science and Technology, CONACyT (Consejo Nacional de Ciencia y Tecnología, CONACyT). México.
    Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
    Aula CIMNE-Morelia. México

Date:
    May, 2024.

Last Modification:
    May 2024
"""

# Library importation
import os
import re
import numpy as np
import Scripts.Graph as Graph
import Scripts.Errors as Errors
from mGFD import TimeDerivative1

## Create a dictionary to get all the regions in da Data folder.
def group_files_by_region(files):
    pattern = re.compile(r'^(.*?)(_p\.csv|_tt\.csv)$')                                      # Look for the files ending in "_p.csv" and "_tt.csv"
    regions = {}                                                                            # Dictionary for the regions.
    for file in files:                                                                      # For each of the files in clouds.
        match = pattern.match(file)                                                         # Check for the match of the pattern with the file.
        if match:                                                                           # If is a match.
            region, suffix = match.groups()                                                 # Get the name of the region.
            if region not in regions:                                                       # If the file haven't been added.
                regions[region] = {}                                                        # Create an empty entry.
            regions[region][suffix] = file                                                  # Add the file to the regions.
    return regions                                                                          # Return the regions dictionary.

## Process the regions and compute the solutions.
def process_region(region, files, data_path, results_path, save):
    print(f'Working on region: {region}')
    if '_p.csv' in files and '_tt.csv' in files:                                            # Check the existence of points and triangles.
        p_file_path  = os.path.join(data_path, files['_p.csv'])                             # Get the file path for the points.
        tt_file_path = os.path.join(data_path, files['_tt.csv'])                            # Get the file path fot the triangles.

        p  = np.genfromtxt(p_file_path,  delimiter = ',', skip_header = 0)                  # Load the coordinates of the points.
        tt = np.genfromtxt(tt_file_path, delimiter = ',', skip_header = 0)                  # Load the triangles correspondence.

        u_ap, u_ex, vec = TimeDerivative1(p, f, t, [v, a, b], operator = L, triangulation = False, tt = [], implicit = False, lam = 0.5)
                                                                                            # Compute the numerical solution.

        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)                                     # Compute the error.
        print(f'\tError: {np.mean(er)}')                                                    # Print the mean of the error.

        if save:                                                                            # If we are going to save.
            os.makedirs(os.path.join(results_path, 'Advection-Diffusion', region), exist_ok = True)
                                                                                            # Ensure the directory exists.
            error_path = os.path.join(results_path, 'Advection-Diffusion', region, 'Error.txt')
                                                                                            # Set the name of the file for the error.
            with open(error_path, 'w') as file:                                             # Create the file.
                file.write(str(np.mean(er)))                                                # Save the error.

            computed_solution_path = os.path.join(results_path, 'Advection-Diffusion', region, 'Computed Solution.csv')
                                                                                            # Set the name of the file for the computed solution.
            np.savetxt(computed_solution_path, u_ap, delimiter = ',', fmt = '%d')           # Save the computed solution.

            theoretical_solution_path = os.path.join(results_path, 'Advection-Diffusion', region, 'Theoretical Solution.csv')
                                                                                            # Set the name of the file for the theoretical solution.
            np.savetxt(theoretical_solution_path, u_ex, delimiter = ',', fmt = '%d')        # Save the theoretical solution.

            plot_path = os.path.join(results_path, 'Advection-Diffusion', region, 'Solution')
                                                                                            # Set the name for the resulting graphs.
            Graph.Cloud_Transient_Steps(p, tt, u_ap, u_ex, nom = plot_path)                 # Save the resulting graphs.

        plot_path = os.path.join(results_path, 'Advection-Diffusion', region, 'Solution.mp4')
                                                                                            # Set the name for the resulting video.
        Graph.Cloud_Transient(p, tt, u_ap, u_ex, save = Save, nom = plot_path)              # Save the resulting video.

# Read the files with the data of all the regions in the Data folder.
## Define the paths for the data and the results for unstructured clouds.
data_clouds    = 'Data/Clouds/'                                                             # Folder with the data of the regions.
results_clouds = 'Results/Clouds/'                                                          # Folder to save the results.

## Define the paths for the data and the results for unstructured clouds with holes.
data_holes     = 'Data/Holes/'                                                              # Folder with the data of the regions.
results_holes  = 'Results/Holes/'                                                           # Folder to save the results.

## Create lists of the data files.
clouds = os.listdir(data_clouds)                                                            # List for the clouds.
holes  = os.listdir(data_holes)                                                             # List for the clouds with holes.

# Group the files by regions.
regions_c = group_files_by_region(clouds)                                                   # Create a dictionary for all the regions in Clouds.
regions_h = group_files_by_region(holes)                                                    # Create a dictionary for all the regions in Holes.

# State the conditions for the problem.
## Problem Parameters
v = 0.1                                                                                     # Diffusion coefficient.
a = 0.3                                                                                     # Transport velocity on the x direction.
b = 0.2                                                                                     # Transport velocity on the y direction.
t = 2000                                                                                    # Number of time-steps.

## Functions for the problem.
f = lambda x, y, t, coef: (1/(4*t+1))*np.exp(-(x-coef[1]*t-0.5)**2/(coef[0]*(4*t+1)) - (y-coef[2]*t-0.5)**2/(coef[0]*(4*t+1)))

## Operator L = [C, D, A, B, C, F]
L = np.vstack([[-a], [-b], [2*v], [0], [2*v], [0]])                                         # Operator coefficients for Au_{xx} + Bu_{xy} + Cu_{yy} + Du_{x} + Eu_{y} + Fu

# Should I save the results?
Save = True                                                                                 # Choose wether the results must be saved.

# Solve the problem using a meshless Generalized Finite Difference approach.
# Solve in clouds.
print('Processing Clouds of points.')
for region, files in regions_c.items():                                                     # For each of the regions.
    process_region(region, files, data_clouds, results_clouds, Save)                        # Process the region.

# Solve in clouds with holes.
print('Processing Clouds of points with Holes.')
for region, files in regions_h.items():                                                     # For each of the regions.
    process_region(region, files, data_holes, results_holes, Save)                          # Process the region.