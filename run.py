import os
import re
import numpy as np
import Scripts.Graph as Graph
import Scripts.Errors as Errors
from mGFD import Stationary
from mGFD import TimeDerivative1
from mGFD import TimeDerivative2

# Define the paths for the data and the results.
data_path    = 'Data/Holes/'
results_path = 'Results/Holes/'

# Create a list of the data files.
files = os.listdir(data_path)

# A dictionary to group the files by region.
regions = {}

# Find the files of the regions.
pattern = re.compile(r'^(.*?)(_p\.csv|_tt\.csv)$')

# Group the file by regions.
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

        #########################  Poisson Equation #########################
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
        #################################################################
         
        ######################### Heat Equation #########################
        print('\tHeat Problem')
        ## Problem Parameters
        v       = 0.2
        t       = 2000                                                                      # Number of time-steps.
        f = lambda x, y, t, coef: np.exp(-2*np.pi**2*coef[0]*t)*np.cos(np.pi*x)*np.cos(np.pi*y)

        ## Operator
        L = np.vstack([[0], [0], [2*v], [0], [2*v], [0]])                                   # L = [D, E, A, B, C, F] from operator Lu = Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + F_u

        ## Problem Solving
        u_ap, u_ex, vec = TimeDerivative1(p, f, t, [v], operator = L, triangulation = False, tt = [], implicit = False, lam = 0.5)

        ##Compute the Error
        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
        print('\t\tError: ', np.mean(er))

        if Save:
            ## Save the Error and the solution
            nom = os.path.join(results_path, 'Heat', region, 'Error.txt')
            with open(nom, 'w') as file:
                file.write(str(np.mean(er)))

            nom = os.path.join(results_path, 'Heat', region, 'Computed Solution.csv')
            np.savetxt(nom, u_ap, delimiter=',', fmt='%d')
            
            nom = os.path.join(results_path, 'Heat', region, 'Theoretical Solution.csv')
            np.savetxt(nom, u_ex, delimiter=',', fmt='%d')

            nom = os.path.join(results_path, 'Heat', region, 'Solution')
            Graph.Cloud_Transient_Steps(p, tt, u_ap, u_ex, nom = nom)

        ## Plot the solution
        nom = os.path.join(results_path, 'Heat', region, 'Solution.mp4')
        Graph.Cloud_Transient(p, tt, u_ap, u_ex, save = Save, nom = nom)
        #################################################################

        #########################  Advection-Diffusion Equation #########################
        print('\tAdvection-Diffusion Problem')
        ## Problem Parameters
        v    = 0.1
        a    = 0.3
        b    = 0.2
        t    = 2000                                                                      # Number of time-steps.
        f = lambda x, y, t, coef: (1/(4*t+1))*np.exp(-(x-coef[1]*t-0.5)**2/(coef[0]*(4*t+1)) - (y-coef[2]*t-0.5)**2/(coef[0]*(4*t+1)))

        ## Operator
        L = np.vstack([[-a], [-b], [2*v], [0], [2*v], [0]])                              # L = [D, E, A, B, C, F] from operator Lu = Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + F_u

        ## Problem Solving
        u_ap, u_ex, vec = TimeDerivative1(p, f, t, [v, a, b], operator = L, triangulation = False, tt = [], implicit = False, lam = 0.5)

        ##Compute the Error
        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
        print('\t\tError: ', np.mean(er))

        if Save:
            ## Save the Error and the solution
            nom = os.path.join(results_path, 'Advection-Diffusion', region, 'Error.txt')
            with open(nom, 'w') as file:
                file.write(str(np.mean(er)))

            nom = os.path.join(results_path, 'Advection-Diffusion', region, 'Computed Solution.csv')
            np.savetxt(nom, u_ap, delimiter=',', fmt='%d')
            
            nom = os.path.join(results_path, 'Advection-Diffusion', region, 'Theoretical Solution.csv')
            np.savetxt(nom, u_ex, delimiter=',', fmt='%d')

            nom = os.path.join(results_path, 'Advection-Diffusion', region, 'Solution')
            Graph.Cloud_Transient_Steps(p, tt, u_ap, u_ex, nom = nom)

        ## Plot the solution
        nom = os.path.join(results_path, 'Advection-Diffusion', region, 'Solution.mp4')
        Graph.Cloud_Transient(p, tt, u_ap, u_ex, save = Save, nom = nom)
        #################################################################

        #########################  Wave Equation #########################
        print('\tWave Problem')
        ## Problem Parameters
        c       = np.sqrt(1/2)                                                              # Wave coefficient.
        t       = 2000                                                                      # Number of time-steps.
        f = lambda x, y, t, coef: np.cos(np.pi*t)*np.sin(np.pi*(x+y))                          # f = \cos{\pi t}\sin{\pi(x + y)}
        g = lambda x, y, t, coef: -np.sin(np.pi*t)*np.sin(np.pi*(x+y))                         # g = -\pi\sin{\pi t}\sin{\pi(x + y)}

        ## Operator
        L = np.vstack([[0], [0], [2*c**2], [0], [2*c**2], [0]])                             # L = [D, E, A, B, C, F] from operator Lu = Au_xx + Bu_xy + Cu_yy + Du_x + Eu_y + F_u

        ## Problem Solving
        u_ap, u_ex, vec = TimeDerivative2(p, f, g, t, [c], operator = L, triangulation = False, tt = None, implicit = False, lam = 1)

        ##Compute the Error
        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
        print('\t\tError: ', np.mean(er))

        if Save:
            ## Save the Error and the solution
            nom = os.path.join(results_path, 'Wave', region, 'Error.txt')
            with open(nom, 'w') as file:
                file.write(str(np.mean(er)))

            nom = os.path.join(results_path, 'Wave', region, 'Computed Solution.csv')
            np.savetxt(nom, u_ap, delimiter=',', fmt='%d')
            
            nom = os.path.join(results_path, 'Wave', region, 'Theoretical Solution.csv')
            np.savetxt(nom, u_ex, delimiter=',', fmt='%d')

            nom = os.path.join(results_path, 'Wave', region, 'Solution')
            Graph.Cloud_Transient_Steps(p, tt, u_ap, u_ex, nom = nom)

        ## Plot the solution
        nom = os.path.join(results_path, 'Wave', region, 'Solution.mp4')
        Graph.Cloud_Transient(p, tt, u_ap, u_ex, save = Save, nom = nom)
        #################################################################