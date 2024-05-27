import os
import re
import numpy as np
import Scripts.Graph as Graph
import Scripts.Errors as Errors
from mGFD import Stationary

# Define the paths for the data and the results.
data_path    = 'Data/Clouds/'
results_path = 'Results/Clouds/'

# Create a list of the data files.
files = os.listdir(data_path)

# A dictionary to group the files by region.
regions = {}

# Find the files of the regions.
pattern = re.compile(r'^(.*?)(_p\.csv|_tt\.csv)$')

# Group the files by regions.
for file in files:
    match = pattern.match(file)
    if match:
        region, suffix = match.groups()
        if region not in regions:
            regions[region] = {}
        regions[region][suffix] = file

# Should I save the results?
Save = True

# Load data files for each region.
counter = 1
for region, files in regions.items():
    print('Working on region: ' + region + '. (' + str(counter) + '/' + str(len(regions)) + ').')
    counter +=1
    if '_p.csv' in files and '_tt.csv' in files:
        p_file_path = os.path.join(data_path, files['_p.csv'])
        tt_file_path = os.path.join(data_path, files['_tt.csv'])

        ## The region is loaded into p and tt.
        p  = np.genfromtxt(p_file_path, delimiter=',', skip_header = 0)
        tt = np.genfromtxt(tt_file_path, delimiter=',', skip_header = 0)

        print('\tPoisson Problem')
        ## Problem parameters
        phi = lambda x, y: 2*np.exp(2*x + y)
        f   = lambda x, y: 10*np.exp(2*x + y)

        ## Operator
        L = np.vstack([[0], [0], [2], [0], [2], [0]])   

        ## Problem Solving
        u_ap, u_ex, vec = Stationary(p, phi, f, operator = L, triangulation = False, tt = None)

        ##Compute the Error
        er = Errors.Cloud_Stationary(p, vec, u_ap, u_ex)
        print('\t\tError: ', np.mean(er))

        if Save:
            ## Save the Error and the solution
            nom = os.path.join(results_path, 'Poisson', region, 'Error.txt')
            with open(nom, 'w') as file:
                file.write(str(np.mean(er)))

            nom = os.path.join(results_path, 'Poisson', region, 'Computed Solution.csv')
            np.savetxt(nom, u_ap, delimiter=',', fmt='%d')
            
            nom = os.path.join(results_path, 'Poisson', region, 'Theoretical Solution.csv')
            np.savetxt(nom, u_ex, delimiter=',', fmt='%d')

        ## Plot the solution
        nom = os.path.join(results_path, 'Poisson', region, 'Solution.png')
        Graph.Cloud_Stationary(p, tt, u_ap, u_ex, save = Save, nom = nom)