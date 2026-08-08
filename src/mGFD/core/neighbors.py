"""
Neighbors — Neighbor search utilities for 2D point clouds

Overview:
    Helper routines to build neighbor lists (vec) for a 2D point cloud used by the mGFD solvers.
    Three main workflows are provided:
    - Isotropic neighbor search within a radius distance (compute_neighbors)
    - Quadrant-balanced neighbor selection to stabilize stencils (compute_balanced_neighbors)
    - Upwind-biased neighbor selection for advection-dominated problems (compute_upwind_neighbors)

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
    compute_neighbors                   Convenience wrapper: computes radius and builds quadrant-balanced vec.
    compute_upwind_neighbors            Convenience wrapper: computes radius and builds upwind vec.
    find_distances                      Compute a characteristic spacing and convert it to a search radius.
    find_neighbors                      Build isotropic neighbor list using KDTree.
    find_neighbors_balanced             Build quadrant-balanced neighbor list for stencil stability.
    find_neighbors_adv                  Build upwind-biased neighbor list using KDTree candidates.

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
    August, 2026.
"""
## Library importation.
import numpy as np                                                                                                  # Core numerical operations.

from scipy.spatial import KDTree                                                                                    # Spatial indexing for fast neighbor queries.

def compute_neighbors(p, nvec):
    """
    compute_neighbors
    Convenience function to build the neighbor list vec for a point cloud.

    This wrapper computes a characteristic search radius dist from the point distribution (KDTree-based)
    and then calls find_neighbors().

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        nvec                        int             Maximum number of neighbors for each node.

    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    m    = int(p.shape[0])                                                                                          # Total number of nodes.
    dist = find_distances(p)                                                                                        # Search radius from nearest-neighbor spacing.

    if nvec > 12:                                                                                                   # Scale dist if asking for expanded stencil.
        dist = dist * np.sqrt(nvec / 12.0)                                                                          # Area scales with r^2, so r scales with sqrt(n).

    vec  = find_neighbors_balanced(p, dist, nvec)                                                                   # Build balanced neighbor list using KDTree candidates.

    return vec                                                                                                      # Return global neighbor list.

def compute_upwind_neighbors(p, a, b, nvec):                                                                        # Compute an upwind-biased neighbor list for advection problems.
    """
    compute_upwind_neighbors
    Convenience function to build an upwind-biased neighbor list for a point cloud.

    This wrapper computes a characteristic search radius dist and then calls find_neighbors_adv(), which
    selects neighbors preferentially in the upwind direction defined by (a, b).

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        a, b                        float           Advection direction components.
        nvec                        int             Maximum number of neighbors for each node.

    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    m    = int(p.shape[0])                                                                                          # Total number of nodes.
    dist = find_distances(p)                                                                                        # Search radius from nearest-neighbor spacing.

    if nvec > 12:                                                                                                   # Scale dist if asking for expanded stencil.
        dist = dist * np.sqrt(nvec / 12.0)                                                                          # Area scales with r^2, so r scales with sqrt(n).

    vec  = find_neighbors_adv(p, dist, a, b, nvec)                                                                  # Build neighbor list using upwind KDTree query.

    return vec                                                                                                      # Return global neighbor list.

def find_distances(p):                                                                                              # Estimate a search radius dist from point spacing.
    """
    find_distances
    Estimate a characteristic spacing from the point cloud and convert it to a search radius dist.

    The output dist is defined as (3/2) times the maximum nearest-neighbor distance.
    This provides a radius that tends to include a small neighborhood around each
    point, while remaining robust to moderate non-uniformity.

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.

    Output:
        dist                        float           Search radius used to collect candidate neighbors.
    """
    ## KDTree-based (memory efficient).
    tree          = KDTree(p[:, :2])                                                                                # Build KDTree on (x, y) coordinates.
    distances, _  = tree.query(p[:, :2], k = 2)                                                                     # Query self + nearest neighbor.
    min_distances = distances[:, 1]                                                                                 # Keep nearest-neighbor distances (skip self at k=0).
    dist          = (3 / 2) * float(np.max(min_distances))                                                          # Convert spacing to a search radius.

    return dist                                                                                                     # Return radius distance.

def find_neighbors(p, dist, nvec):                                                                                  # Build an isotropic neighbor list within a radius.
    """
    find_neighbors
    Find up to nvec neighbors for each point within a radius distance dist.

    Each row vec[i, :] contains indices of nearby points (excluding i). Unused slots are padded with -1.

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        dist                        float           Search radius.
        nvec                        int             Maximum number of neighbors.

    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    ## Variable initialization.
    m   = int(len(p[:, 0]))                                                                                         # Total number of nodes.
    vec = np.zeros([m, nvec], dtype = int) - 1                                                                      # Initialize neighbor matrix with padding -1.

    # KDTree radius query per node.
    tree = KDTree(p[:, :2])                                                                                         # Build KDTree on (x, y) coordinates.
    for i in range(m):                                                                                              # For each node i.
        distances, indices = tree.query(p[i, :2], k = nvec + 1, distance_upper_bound = dist)                        # Query up to nvec+1 neighbors inside radius.
        valid_indices      = indices[distances < dist]                                                              # Keep only neighbors with valid distances.
        valid_indices      = valid_indices[valid_indices != i]                                                      # Remove self from the neighbor list.
        vec[i, :min(len(valid_indices), nvec)] = valid_indices[:nvec]                                               # Store up to nvec neighbors.

    return vec                                                                                                      # Return neighbor list.

def find_neighbors_balanced(p, dist, nvec):                                                                         # Build a quadrant-balanced neighbor list.
    """
    find_neighbors_balanced
    Find up to nvec neighbors for each point, ensuring a balanced spatial distribution by
    attempting to select an equal number of neighbors from each of the 4 geometric quadrants.

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        dist                        float           Search radius.
        nvec                        int             Maximum number of neighbors.

    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    m   = int(len(p[:, 0]))                                                                                         # Total number of nodes.
    vec = np.zeros([m, nvec], dtype = int) - 1                                                                      # Initialize neighbor matrix with padding -1.

    nvec = int(nvec)
    if m == 0 or nvec <= 0:
        return vec

    # Query a large candidate pool to ensure we find points in all quadrants if possible.
    k_cand = min(m, 10 * nvec)                                                                                      # Expand candidate search pool.
    tree = KDTree(p[:, :2])                                                                                         # Build KDTree on (x, y).
    distances, indices = tree.query(p[:, :2], k = k_cand)                                                           # Query candidates.

    if k_cand == 1:
        distances = distances.reshape(-1, 1)
        indices   = indices.reshape(-1, 1)

    target_per_quad = int(np.ceil(nvec / 4.0))                                                                      # Target number of points per quadrant.

    for i in range(m):
        cand   = indices[i, 1:]                                                                                     # Candidate indices (skip self).
        dist_c = distances[i, 1:]                                                                                   # Candidate distances.

        # Filter valid candidates within a slightly expanded radius to allow filling empty quadrants.
        good   = (cand >= 0) & (cand < m) & np.isfinite(dist_c) & (dist_c <= 1.5 * dist)
        cand   = cand[good].astype(np.int32, copy = False)
        dist_c = dist_c[good]

        if cand.size == 0:
            continue

        dx = p[cand, 0] - p[i, 0]                                                                                   # DX relative to central node.
        dy = p[cand, 1] - p[i, 1]                                                                                   # DY relative to central node.

        # Define quadrants
        q1 = (dx >= 0) & (dy >= 0)                                                                                  # NE quadrant.
        q2 = (dx < 0)  & (dy >= 0)                                                                                  # NW quadrant.
        q3 = (dx < 0)  & (dy < 0)                                                                                   # SW quadrant.
        q4 = (dx >= 0) & (dy < 0)                                                                                   # SE quadrant.

        selected = []
        for q in [q1, q2, q3, q4]:                                                                                  # Extract balanced points.
            q_cand = cand[q]
            selected.extend(q_cand[:target_per_quad])

        # If we couldn't fill the stencil purely with balanced points, fill with closest remaining.
        if len(selected) < nvec:
            remaining = [c for c in cand if c not in selected]
            needed = nvec - len(selected)
            selected.extend(remaining[:needed])

        selected = selected[:nvec]                                                                                  # Truncate if we overshot.
        take = len(selected)
        vec[i, :take] = selected                                                                                    # Save selected neighbors.

    return vec                                                                                                      # Return balanced neighbor list.

def find_neighbors_adv(p, dist, a, b, nvec):                                                                        # Build an upwind-biased neighbor list.
    """
    find_neighbors_adv
    Find up to nvec neighbors per node with an upwind preference for direction (a, b).

    The procedure collects a set of KDTree candidate neighbors for each node, sorts them by distance,
    then splits them into upwind (dx*a + dy*b <= tol) and downwind groups. Upwind candidates are taken
    first, and remaining slots are filled with downwind candidates.

    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        dist                        float           Search radius.
        a, b                        float           Advection direction components.
        nvec                        int             Maximum number of neighbors.

    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.  
    """
    m   = int(len(p[:, 0]))                                                                                         # Total number of nodes.
    vec = np.zeros([m, nvec], dtype = int) - 1                                                                      # Initialize neighbor matrix with padding -1.

    a    = float(a)                                                                                                 # Normalize a to float.
    b    = float(b)                                                                                                 # Normalize b to float.
    nvec = int(nvec)                                                                                                # Normalize neighbor count to int.

    if m == 0 or nvec <= 0:                                                                                         # Handle empty input or invalid neighbor count.
        return vec                                                                                                  # Return empty neighbor matrix.

    speed = float(np.hypot(a, b))                                                                                   # Magnitude of advection direction.

    if (not np.isfinite(speed)) or speed == 0.0:                                                                    # Fallback when direction is invalid.
        dist0 = find_distances(p)                                                                                   # Use isotropic distance estimate.
        return find_neighbors(p, dist0, nvec)                                                                       # Return isotropic KDTree neighbors.

    k = max(nvec + 1, 128, 40 * nvec)                                                                               # Candidate pool size for KDTree query.
    k = int(min(m, k))                                                                                              # Clamp to m to avoid invalid k.

    tree               = KDTree(p[:, :2])                                                                           # Build KDTree on (x, y) coordinates.
    distances, indices = tree.query(p[:, :2], k = k)                                                                # Query k nearest candidates for all points.

    if k == 1:                                                                                                      # Ensure 2D shapes for consistent slicing.
        distances = distances.reshape(-1, 1)                                                                        # Reshape distances to (m, 1).
        indices   = indices.reshape(-1, 1)                                                                          # Reshape indices to (m, 1).

    tol = 1e-14 * speed                                                                                             # Tolerance to decide upwind vs downwind.

    for i in range(m):                                                                                              # Loop over nodes.
        cand   = indices[i, 1:]                                                                                     # Candidate neighbor indices (skip self at position 0).
        dist_c = distances[i, 1:]                                                                                   # Candidate neighbor distances.
        good   = (cand >= 0) & (cand < m) & np.isfinite(dist_c)                                                     # Filter invalid indices and infinite distances.
        cand   = cand[good].astype(np.int32, copy = False)                                                          # Keep only valid candidates.
        dist_c = dist_c[good]                                                                                       # Keep aligned distances for sorting.

        if cand.size == 0:                                                                                          # No valid candidates.
            continue                                                                                                # Leave row as -1 padding.

        order  = np.lexsort((cand, dist_c))                                                                         # Sort primarily by distance, then by index.
        cand   = cand[order]                                                                                        # Apply sorting to candidate indices.
        dist_c = dist_c[order]                                                                                      # Apply sorting to candidate distances.

        dx     = p[cand, 0] - p[i, 0]                                                                               # Candidate dx offsets.
        dy     = p[cand, 1] - p[i, 1]                                                                               # Candidate dy offsets.
        dot    = dx * a + dy * b                                                                                    # Directional dot product to classify upwind/downwind.
        up     = cand[dot <= tol]                                                                                   # Upwind candidates (including near-orthogonal).
        dn     = cand[dot > tol]                                                                                    # Downwind candidates.

        if up.size >= nvec:                                                                                         # If upwind set is large enough.
            vec[i, :] = up[:nvec]                                                                                   # Take the closest nvec upwind neighbors.
        else:                                                                                                       # Otherwise, fill remaining slots with downwind.
            fill          = np.concatenate([up, dn])                                                                # Concatenate upwind-first ordering.
            take          = min(nvec, int(fill.size))                                                               # Limit to available candidates.
            vec[i, :take] = fill[:take]                                                                             # Store selected neighbors.

    return vec                                                                                                      # Return upwind neighbor list.

