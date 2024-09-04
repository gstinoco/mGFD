"""
All the codes presented below were developed by:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

With the funding of:
    National Council of Humanities, Sciences and Technologies, CONAHCyT (Consejo Nacional de Humanidades, Ciencias y Tecnologías, CONAHCyT). México.
    Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
    Aula CIMNE-Morelia. México

Date:
    May, 2024.

Last Modification:
    May, 2024.
"""

## Library importation.
import numpy as np
from scipy.spatial import KDTree
from joblib import Parallel, delayed

def Triangulation(p, tt, nvec):
    """
    Triangulation
    Function to find the neighbor nodes in a triangulation.
    
    Input:
        p           m x 2           double          Array with the coordinates of the nodes.
        tt          n x 3           double          Array with the correspondence of the n triangles.
        nvec                        integer         Maximum number of neighbors.
    
    Output:
        vec         m x nvec        double          Array with matching neighbors of each node.
    """

    ## Variable initialization.
    m   = len(p[:, 0])                                                              # The size if the triangulation is obtained.
    vec = np.zeros([m, nvec], dtype=int)-1                                          # The array for the neighbors is initialized.

    ## Neighbor search.
    for i in np.arange(m):                                                          # For each of the nodes.
        kn    = np.argwhere(tt == i)                                                # Search in which triangles the node appears.
        vec2  = np.setdiff1d(tt[kn[:, 0]], i)                                       # Neighbors are stored inside vec2.
        vec2  = np.vstack([vec2])                                                   # Convert vec2 to a column.
        nvec2 = sum(vec2[0, :] != -1)                                               # The number of neighbors of the node is calculated.
        nnvec = np.minimum(nvec, nvec2)                                             # The real number of neighbors.
        for j in np.arange(nnvec):                                                  # For each of the nodes.
            vec[i, j] = vec2[0, j]                                                  # Neighbors are saved.
    return vec

def Cloud(p, nvec):
    """
    Cloud
    Function to find the neighbor nodes in a cloud of points.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and a flag for the boundary.
        nvec                        integer         Maximum number of neighbors.
    
    Output:
        vec         m x nvec        ndarray         Array with matching neighbors of each node.
    """

    ## Delta computation.
    dist = find_distances(p, mode = 2)

    ## Neighbor search.
    vec = find_neighbors(p, dist, nvec, mode = 3)

    return vec

def CloudAdv(p, a, b, nvec):
    """
    CloudAdv
    Function to find the neighbor nodes in a cloud of points; it considers an Upwind scheme
    
    Input:
        p               ndarray         Array with the coordinates of the nodes and a flag for the boundary.
        nvec            int             Maximum number of neighbors.
    
    Output:
        vec             ndarray         Array with matching neighbors of each node.
    """

    ## Delta computation.
    dist = find_distances(p, mode = 2)

    ## Neighbor search.
    vec = find_neighbors_adv(p, dist, a, b, nvec)

    return vec

def find_distances(p, mode = 2):
    """
    find_distances
    Function to find the distances between all the give nodes.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and a flag for the boundary.
        mode                        integer         Choose the way to compute the distances:
                                                    1: brute force
                                                    2: optimized (default)
    
    Output:
        dist                        float           The maximum distance between two consecutive nodes.
    """

    ## Variable initialization.
    m    = len(p[:, 0])                                                             # The size if the triangulation is obtained.

    if mode == 1:
        ## Brute Force
        dmin = np.zeros([m, 1]) + 10                                                # dmin initialization with a "big" value.
        for i in np.arange(m):                                                      # For each of the nodes.
            x    = p[i, 0]                                                          # x coordinate of the central node.
            y    = p[i, 1]                                                          # y coordinate of the central node.
            for j in np.arange(m):                                                  # For all the nodes.
                if i != j:                                                          # If the the node is different to the central one.
                    x1 = p[j, 0]                                                    # x coordinate of the possible neighbor.
                    y1 = p[j, 1]                                                    # y coordinate of the possible neighbor.
                    d  = np.sqrt((x - x1)**2 + (y - y1)**2)                         # Distance from the possible neighbor to the central node.
                    dmin[i] = min(dmin[i], d)                                       # Look for the distance to the closest node.
        dist = (3/2)*np.max(dmin)                                                   # The distance is the maximum distance between two consecutive nodes.

    if mode == 2:
        ## Optimized.
        p_expanded    = np.expand_dims(p, axis=1)                                   # Expand p to use vectorized operations.
        differences   = p_expanded - p                                              # The distance between each node and all the others is computed.
        distances     = np.sum(differences**2, axis=2)                              # The sum x^2 + y^2 is perform to compute the Euclidean norm.
        np.fill_diagonal(distances, np.inf)                                         # Distances to the self node are state as infinity and not zero.
        min_distances = np.sqrt(np.min(distances, axis=1))                          # Look for the distance to the closest node.
        dist          = (3/2)*np.max(min_distances)                                 # The distance is the maximum distance between two consecutive nodes.
    
    return dist

def find_neighbors(p, dist, nvec, mode = 2):
    """
    find_neighbors
    Function to find all the neighbors of a node withing a given distance.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and a flag for the boundary.
        dist                        float           Radius distance to look for neighbors.
        nvec                        integer         Maximum number of neighbors.
        mode                        integer         Choose the way to compute the distances:
                                                    1: brute force
                                                    2: optimized (default)
    
    Output:
        vec         m x nvec        ndarray         Array with matching neighbors of each node.
    """

    ## Variable initialization.
    m    = len(p[:, 0])                                                             # The size if the triangulation is obtained.
    vec  = np.zeros([m, nvec], dtype = int) - 1                                     # The array for the neighbors is initialized.

    if mode == 1:
        ## Brute Force
        for i in np.arange(m):                                                      # For each of the nodes.
            x, y  = p[i, 0], p[i, 1]                                                # x, y coordinates of the central node.
            temp_neighbors = []                                                     # Create an empty array for neighbors.
            for j in np.arange(m):                                                  # For all the interior nodes.
                if i != j:                                                          # Check that we are not working with the central node.
                    x1, y1 = p[j, 0], p[j, 1]                                       # x, y coordinates of the possible neighbor.
                    d = np.sqrt((x - x1)**2 + (y - y1)**2)                          # Distance from the possible neighbor to the central node.
                    if d < dist:                                                    # If the distance is smaller or equal to the tolerance distance.
                        temp_neighbors.append((d, j))                               # Store the neighbor node.
            temp_neighbors.sort()                                                   # Sort the neighbors by distance.
            for idx, (d, j) in enumerate(temp_neighbors[:nvec]):                    # For each stored neighbor.
                vec[i, idx] = j                                                     # Store the neighbor.
    
    elif mode == 2:
        ## Optimized
        dx = np.expand_dims(p[:, 0 ], 1) - np.expand_dims(p[:, 0], 0)                # Compute dx between all the nodes.
        dy = np.expand_dims(p[:, 1 ], 1) - np.expand_dims(p[:, 1], 0)                # Compute dy between all the nodes.
        radius = np.sqrt(dx**2 + dy**2)                                              # Get the distances from each node to all the others.
        
        for i in range(m):                                                           # For each node.
            neighbors = np.where((radius[i, :] < dist) & (np.arange(m) != i))[0]     # The neighbors are all the nodes within the radius.

            if len(neighbors) > 0:                                                   # If there are more neighbors than the requested.
                neighbors = neighbors[np.argsort(radius[i, neighbors])][:nvec]       # The neighbors are sorted by distance and only the closest nvec neighbors remains.
                vec[i, :len(neighbors)] = neighbors                                  # The nvec neighbors are stored.
            else:
                neighbors = neighbors[np.argsort(radius[i, neighbors])]              # The neighbors are sorted by distance and only the closest nvec neighbors remains.
                vec[i, :len(neighbors)] = neighbors                                  # The neighbors are stored.
    
    elif mode == 3:
        # KDTree
        tree = KDTree(p[:, :2])                                                     # Create a KDTree using the first two columns of p (x and y coordinates).
        for i in range(m):                                                          # For each node.
            distances, indices = tree.query(p[i, :2], k = nvec + 1, distance_upper_bound = dist)
            valid_indices = indices[distances < dist]                               # Filter out invalid distances.
            valid_indices = valid_indices[valid_indices != i]                       # Remove self from neighbors.
            vec[i, :min(len(valid_indices), nvec)] = valid_indices[:nvec]           # Store the neighbors.

    return vec

def find_neighbors_adv(p, dist, a, b, nvec, n_jobs = -1):
    """
    find_neighbors
    Function to find all the neighbors of a node within a given distance, with parallelization.
    
    Input:
        p                   ndarray         Array with the coordinates of the nodes and a flag for the boundary.
        dist                float           Radius distance to look for neighbors.
        nvec                int             Maximum number of neighbors.
        mode                int             Choose the way to compute the distances:
                                                    1: brute force
                                                    2: optimized (default)
                                                    3: KDTree
    
    Output:
        vec                 ndarray         Array with matching neighbors of each node.
    """
    m   = len(p[:, 0])                                                              # The size if the triangulation is obtained.
    vec = np.zeros([m, nvec], dtype=int) - 1                                        # The array for the neighbors is initialized.

    # Brute Force.
    vec_rows = Parallel(n_jobs=n_jobs)(delayed(find_neighbors_brute_force)(i, p, a, b, dist, nvec) for i in range(m))
    vec      = np.array(vec_rows)

    return vec

def find_neighbors_brute_force(i, p, a, b, dist, nvec):
    """
    Helper function for finding neighbors using brute force method.
    """
    x, y = p[i, 0], p[i, 1]                                                         # x, y coordinates of the central node.
    temp_neighbors = []                                                             # Create an empty array for neighbors.
    for j in np.arange(len(p)):                                                     # For all the nodes.
        if i != j:                                                                  # Check that we are not working with the central node.
            x1, y1 = p[j, 0], p[j, 1]                                               # x, y coordinates of the possible neighbor.
            xt, yt = x1 - x, y1 - y                                                 # Change the "origin" of x and y.
            di = np.sign(np.dot([xt, yt],[a, b]))                                   # Check the direction between (xt, yt) and (a,b).
            d  = np.sqrt((x - x1)**2 + (y - y1)**2)                                 # Distance from the possible neighbor to the central node.
            if di == -1 and d < dist:                                               # If the neighbors are in the correct position and the distance is smaller or equal to the tolerance distance.
                temp_neighbors.append((d, j))                                       # Store the neighbor node.
    temp_neighbors.sort()                                                           # Sort the neighbors by distance.
    neighbors = np.array([j for _, j in temp_neighbors[:nvec]])                     # Extract neighbor indices.
    vec_row   = np.zeros(nvec, dtype=int) - 1                                       # Initialize row with -1.
    vec_row[:len(neighbors)] = neighbors                                            # Store the neighbors.
    return vec_row