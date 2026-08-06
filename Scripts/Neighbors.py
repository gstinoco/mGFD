"""
Neighbors — Neighbor search utilities for 2D point clouds

Overview:
    Helper routines to build neighbor lists (vec) for a 2D point cloud used by the mGFD solvers.
    Two main workflows are provided:
    - Isotropic neighbor search within a radius distance (Cloud / find_neighbors)
    - Upwind-biased neighbor selection for advection-dominated problems (CloudAdv / find_neighbors_adv)

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag].
    vec     (m, nvec) ndarray[int]
            Neighbor list. Each row contains neighbor indices; unused slots are padded with -1.
    nvec    int
            Maximum number of neighbors stored per node.
    dist    float
            Search radius used to collect candidate neighbors.
    a, b    float
            Advection direction components for upwind selection (vector (a, b)).

Public API:
    Cloud                   Convenience wrapper: compute dist and then build vec (KDTree mode).
    CloudAdv                Convenience wrapper: compute dist and then build upwind vec (KDTree-based).
    find_distances          Compute a characteristic spacing and convert it to a search radius.
    find_neighbors          Build isotropic neighbor list (brute force / vectorized / KDTree).
    find_neighbors_adv      Build upwind-biased neighbor list using KDTree candidates.
    find_neighbors_brute_force  Helper used by older parallel approaches (kept for compatibility).

Dependencies:
    SciPy is required for KDTree-based modes. joblib is imported for compatibility with older
    parallel implementations, even if the current implementation runs serially.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI (Secretaria de Ciencia, Humanidades, Tecnología e Innovación). México.
        Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
        Aula CIMNE-Morelia. México
        SIIIA-MATH: Soluciones de Ingeniería. México

    Date:
        May, 2024.
    Last Modification:
        August, 2026.
"""
## Library importation.
import numpy as np                                                                                      # Core numerical operations.
from scipy.spatial import KDTree                                                                        # Spatial indexing for fast neighbor queries.




def Cloud(p, nvec):
    """
    Cloud
    Convenience function to build the neighbor list vec for a point cloud.
    
    This wrapper computes a characteristic search radius dist from the point distribution (KDTree-based)
    and then calls find_neighbors() in KDTree mode.
    
    Input:
        p           m x 3           ndarray         Point cloud [x, y, flag].
        nvec                        int             Maximum number of neighbors per node.
    
    Output:
        vec         m x nvec        ndarray[int]    Neighbor indices (global) for each node (padded with -1).
    """

    m    = int(p.shape[0])                                                                              # Total number of nodes.
    dist = find_distances(p, mode = 3)                                                                  # Search radius from nearest-neighbor spacing.
    vec  = find_neighbors(p, dist, nvec, mode = 3)                                                      # Build neighbor list using KDTree query.

    return vec                                                                                          # Return global neighbor list.

def CloudAdv(p, a, b, nvec):                                                                            # Compute an upwind-biased neighbor list for advection problems.
    """
    CloudAdv
    Convenience function to build an upwind-biased neighbor list for a point cloud.
    
    This wrapper computes a characteristic search radius dist and then calls find_neighbors_adv(), which
    selects neighbors preferentially in the upwind direction defined by (a, b).
    
    Input:
        p               ndarray         Point cloud [x, y, flag].
        a, b            float           Advection direction components.
        nvec            int             Maximum number of neighbors per node.
    
    Output:
        vec             ndarray[int]    Neighbor indices (global) for each node (padded with -1).
    """

    m    = int(p.shape[0])                                                                              # Total number of nodes.
    dist = find_distances(p, mode = 3)                                                                  # Search radius from nearest-neighbor spacing.
    vec  = find_neighbors_adv(p, dist, a, b, nvec)                                                      # Build neighbor list using upwind KDTree query.

    return vec                                                                                          # Return global neighbor list.

def find_distances(p, mode = 3):                                                                        # Estimate a search radius dist from point spacing.
    """
    find_distances
    Estimate a characteristic spacing from the point cloud and convert it to a search radius dist.
    
    The output dist is defined as (3/2) times the maximum nearest-neighbor distance, computed according
    to the selected mode. This provides a radius that tends to include a small neighborhood around each
    point, while remaining robust to moderate non-uniformity.
    
    Input:
        p           m x 3           ndarray         Point cloud [x, y, flag] (only x,y are used).
        mode                        int             Distance estimation mode:
                                                    1: brute force O(m^2)
                                                    2: vectorized O(m^2) (memory heavy)
                                                    3: KDTree-based (memory efficient)
    
    Output:
        dist                        float           Search radius used to collect candidate neighbors.
    """

    ## Variable initialization.
    m = int(len(p[:, 0]))                                                                               # Total number of nodes.

    if mode == 1:                                                                                       # Brute-force nearest-neighbor scan.
        ## Brute Force
        dmin = np.zeros([m, 1]) + 10                                                                    # Initialize with a large distance.
        for i in np.arange(m):                                                                          # For each node i.
            x = p[i, 0]                                                                                 # x coordinate of node i.
            y = p[i, 1]                                                                                 # y coordinate of node i.
            for j in np.arange(m):                                                                      # For each candidate node j.
                if i != j:                                                                              # Skip self.
                    x1 = p[j, 0]                                                                        # x coordinate of node j.
                    y1 = p[j, 1]                                                                        # y coordinate of node j.
                    d  = np.sqrt((x - x1)**2 + (y - y1)**2)                                             # Euclidean distance ||p_i - p_j||.
                    dmin[i] = min(dmin[i], d)                                                           # Track the nearest neighbor distance for i.
        dist = (3 / 2) * np.max(dmin)                                                                   # Convert spacing to a search radius.

    if mode == 2:                                                                                       # Vectorized O(m^2) method (memory heavy).
        ## Optimized.
        xy          = p[:, :2]                                                                          # Use only (x, y) coordinates for distance computation.
        p_expanded  = np.expand_dims(xy, axis = 1)                                                      # Expand xy for broadcasting.
        differences = p_expanded - xy                                                                   # Pairwise differences p_i - p_j.
        distances   = np.sum(differences**2, axis = 2)                                                  # Squared distances.
        np.fill_diagonal(distances, np.inf)                                                             # Ignore self distances.
        min_distances = np.sqrt(np.min(distances, axis = 1))                                            # Nearest-neighbor distance for each node.
        dist = (3 / 2) * np.max(min_distances)                                                          # Convert spacing to a search radius.

    if mode == 3:                                                                                       # KDTree-based nearest-neighbor query.
        ## KDTree-based (memory efficient).
        tree = KDTree(p[:, :2])                                                                         # Build KDTree on (x, y) coordinates.
        distances, _ = tree.query(p[:, :2], k = 2)                                                      # Query self + nearest neighbor.
        min_distances = distances[:, 1]                                                                 # Keep nearest-neighbor distances (skip self at k=0).
        dist = (3 / 2) * float(np.max(min_distances))                                                   # Convert spacing to a search radius.
    
    return dist                                                                                         # Return radius distance.

def find_neighbors(p, dist, nvec, mode = 2):                                                            # Build an isotropic neighbor list within a radius.
    """
    find_neighbors
    Find up to nvec neighbors for each point within a radius distance dist.
    
    Each row vec[i, :] contains indices of nearby points (excluding i). Unused slots are padded with -1.
    
    Input:
        p           m x 3           ndarray         Point cloud [x, y, flag] (only x,y are used).
        dist                        float           Search radius.
        nvec                        int             Maximum number of neighbors.
        mode                        int             Search mode:
                                                    1: brute force O(m^2)
                                                    2: vectorized O(m^2) (memory heavy)
                                                    3: KDTree-based (memory efficient)
    
    Output:
        vec         m x nvec        ndarray[int]    Neighbor indices per node (padded with -1).
    """

    ## Variable initialization.
    m   = int(len(p[:, 0]))                                                                             # Total number of nodes.
    vec = np.zeros([m, nvec], dtype = int) - 1                                                          # Initialize neighbor matrix with padding -1.

    if mode == 1:                                                                                       # Brute force radius scan.
        ## Brute Force
        for i in np.arange(m):                                                                          # For each node i.
            x, y = p[i, 0], p[i, 1]                                                                     # Coordinates of node i.
            temp_neighbors = []                                                                         # Accumulator for (distance, index) pairs.
            for j in np.arange(m):                                                                      # For each candidate node j.
                if i != j:                                                                              # Skip self.
                    x1, y1 = p[j, 0], p[j, 1]                                                           # Coordinates of node j.
                    d = np.sqrt((x - x1)**2 + (y - y1)**2)                                              # Distance from i to j.
                    if d < dist:                                                                        # Keep only candidates inside the radius.
                        temp_neighbors.append((d, j))                                                   # Store candidate neighbor and its distance.
            temp_neighbors.sort()                                                                       # Sort by increasing distance.
            for idx, (d, j) in enumerate(temp_neighbors[:nvec]):                                        # Keep at most the closest nvec neighbors.
                vec[i, idx] = j                                                                         # Store neighbor index.
    
    elif mode == 2:                                                                                     # Vectorized all-pairs distances (memory heavy).
        ## Optimized
        dx     = np.expand_dims(p[:, 0], 1) - np.expand_dims(p[:, 0], 0)                                # Pairwise dx.
        dy     = np.expand_dims(p[:, 1], 1) - np.expand_dims(p[:, 1], 0)                                # Pairwise dy.
        radius = np.sqrt(dx**2 + dy**2)                                                                 # Pairwise Euclidean distances.
        
        for i in range(m):                                                                              # For each node i.
            neighbors = np.where((radius[i, :] < dist) & (np.arange(m) != i))[0]                        # Candidate neighbors inside radius, excluding self.

            if len(neighbors) > 0:                                                                      # If there is at least one neighbor candidate.
                neighbors = neighbors[np.argsort(radius[i, neighbors])][:nvec]                          # Sort by distance and keep closest nvec.
                vec[i, :len(neighbors)] = neighbors                                                     # Store found neighbors.
            else:                                                                                       # No candidates inside radius.
                neighbors = neighbors[np.argsort(radius[i, neighbors])]                                 # Keep as empty (no-op, preserved for legacy structure).
                vec[i, :len(neighbors)] = neighbors                                                     # Store nothing.
    
    elif mode == 3:                                                                                     # KDTree radius query per node.
        # KDTree
        tree = KDTree(p[:, :2])                                                                         # Build KDTree on (x, y) coordinates.
        for i in range(m):                                                                              # For each node i.
            distances, indices = tree.query(p[i, :2], k = nvec + 1, distance_upper_bound = dist)        # Query up to nvec+1 neighbors inside radius.
            valid_indices = indices[distances < dist]                                                   # Keep only neighbors with valid distances.
            valid_indices = valid_indices[valid_indices != i]                                           # Remove self from the neighbor list.
            vec[i, :min(len(valid_indices), nvec)] = valid_indices[:nvec]                               # Store up to nvec neighbors.

    return vec                                                                                          # Return neighbor list.

def find_neighbors_adv(p, dist, a, b, nvec, n_jobs = -1):                                               # Build an upwind-biased neighbor list.
    """
    find_neighbors_adv
    Find up to nvec neighbors per node with an upwind preference for direction (a, b).
    
    The procedure collects a set of KDTree candidate neighbors for each node, sorts them by distance,
    then splits them into upwind (dx*a + dy*b <= tol) and downwind groups. Upwind candidates are taken
    first, and remaining slots are filled with downwind candidates.
    
    Input:
        p                   ndarray         Point cloud [x, y, flag] (only x,y are used).
        dist                float           Search radius.
        a, b                float           Advection direction components.
        nvec                int             Maximum number of neighbors.
        n_jobs              int             Reserved for legacy parallel implementations (unused here).
    
    Output:
        vec                 ndarray[int]    Neighbor indices per node (padded with -1).
    """
    m   = int(len(p[:, 0]))                                                                             # Total number of nodes.
    vec = np.zeros([m, nvec], dtype = int) - 1                                                          # Initialize neighbor matrix with padding -1.

    a    = float(a)                                                                                     # Normalize a to float.
    b    = float(b)                                                                                     # Normalize b to float.
    nvec = int(nvec)                                                                                    # Normalize neighbor count to int.
    if m == 0 or nvec <= 0:                                                                             # Handle empty input or invalid neighbor count.
        return vec                                                                                      # Return empty neighbor matrix.

    speed = float(np.hypot(a, b))                                                                       # Magnitude of advection direction.
    if (not np.isfinite(speed)) or speed == 0.0:                                                        # Fallback when direction is invalid.
        dist0 = find_distances(p, mode = 3)                                                             # Use isotropic distance estimate.
        return find_neighbors(p, dist0, nvec, mode = 3)                                                 # Return isotropic KDTree neighbors.

    k = max(nvec + 1, 128, 40 * nvec)                                                                   # Candidate pool size for KDTree query.
    k = int(min(m, k))                                                                                  # Clamp to m to avoid invalid k.

    tree = KDTree(p[:, :2])                                                                             # Build KDTree on (x, y) coordinates.
    distances, indices = tree.query(p[:, :2], k = k)                                                    # Query k nearest candidates for all points.
    if k == 1:                                                                                          # Ensure 2D shapes for consistent slicing.
        distances = distances.reshape(-1, 1)                                                            # Reshape distances to (m, 1).
        indices   = indices.reshape(-1, 1)                                                              # Reshape indices to (m, 1).

    tol = 1e-14 * speed                                                                                 # Tolerance to decide upwind vs downwind.
    for i in range(m):                                                                                  # Loop over nodes.
        cand   = indices[i, 1:]                                                                         # Candidate neighbor indices (skip self at position 0).
        dist_c = distances[i, 1:]                                                                       # Candidate neighbor distances.
        good   = (cand >= 0) & (cand < m) & np.isfinite(dist_c)                                         # Filter invalid indices and infinite distances.
        cand   = cand[good].astype(np.int32, copy = False)                                              # Keep only valid candidates.
        dist_c = dist_c[good]                                                                           # Keep aligned distances for sorting.
        if cand.size == 0:                                                                              # No valid candidates.
            continue                                                                                    # Leave row as -1 padding.

        order  = np.lexsort((cand, dist_c))                                                             # Sort primarily by distance, then by index.
        cand   = cand[order]                                                                            # Apply sorting to candidate indices.
        dist_c = dist_c[order]                                                                          # Apply sorting to candidate distances.

        dx  = p[cand, 0] - p[i, 0]                                                                      # Candidate dx offsets.
        dy  = p[cand, 1] - p[i, 1]                                                                      # Candidate dy offsets.
        dot = dx * a + dy * b                                                                           # Directional dot product to classify upwind/downwind.
        up  = cand[dot <= tol]                                                                          # Upwind candidates (including near-orthogonal).
        dn  = cand[dot > tol]                                                                           # Downwind candidates.

        if up.size >= nvec:                                                                             # If upwind set is large enough.
            vec[i, :] = up[:nvec]                                                                       # Take the closest nvec upwind neighbors.
        else:                                                                                           # Otherwise, fill remaining slots with downwind.
            fill = np.concatenate([up, dn])                                                             # Concatenate upwind-first ordering.
            take = min(nvec, int(fill.size))                                                            # Limit to available candidates.
            vec[i, :take] = fill[:take]                                                                 # Store selected neighbors.

    return vec                                                                                          # Return upwind neighbor list.

def find_neighbors_brute_force(i, p, a, b, dist, nvec):                                                 # Brute-force upwind neighbor selection for a single node.
    """
    find_neighbors_brute_force
    Helper function for finding upwind neighbors for a single node using brute force.
    
    This function scans all nodes j != i, keeps those within dist and with upwind direction
    (sign(dot([xt, yt], [a, b])) == -1), sorts by distance, and returns a padded vec row.
    
    Input:
        i                           int             Central node index.
        p                           ndarray         Point cloud [x, y, flag].
        a, b                        float           Advection direction components.
        dist                        float           Search radius.
        nvec                        int             Maximum number of neighbors.
    
    Output:
        vec_row                      ndarray[int]   Length-nvec array of neighbor indices (padded with -1).
    """
    x, y = p[i, 0], p[i, 1]                                                                             # Coordinates of the central node i.
    temp_neighbors = []                                                                                 # Accumulator for (distance, index) pairs.
    for j in np.arange(len(p)):                                                                         # Loop over all nodes j.
        if i != j:                                                                                      # Skip self.
            x1, y1 = p[j, 0], p[j, 1]                                                                   # Coordinates of candidate node j.
            xt, yt = x1 - x, y1 - y                                                                     # Offset vector from i to j.
            di = np.sign(np.dot([xt, yt], [a, b]))                                                      # Determine if candidate is upwind/downwind.
            d  = np.sqrt((x - x1)**2 + (y - y1)**2)                                                     # Distance from i to j.
            if di == -1 and d < dist:                                                                   # Keep only upwind neighbors inside radius.
                temp_neighbors.append((d, j))                                                           # Store candidate neighbor.
    temp_neighbors.sort()                                                                               # Sort by increasing distance.
    neighbors = np.array([j for _, j in temp_neighbors[:nvec]])                                         # Extract up to nvec neighbor indices.
    vec_row   = np.zeros(nvec, dtype = int) - 1                                                         # Initialize padded row with -1.
    vec_row[:len(neighbors)] = neighbors                                                                # Store selected neighbors.
    return vec_row                                                                                      # Return one vec row.
