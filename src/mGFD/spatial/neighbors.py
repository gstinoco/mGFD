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
    compute_mesh_spacing                Compute characteristic spatial node spacing (h_min, h_avg) for a 2D point cloud.
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
import numpy as np                                                                                                                      # Core numerical operations.
import numba as nb                                                                                                                      # JIT compilation.

from scipy.spatial import KDTree                                                                                                        # Spatial indexing for fast neighbor queries.
from typing import Optional, Tuple                                                                                                      # Type hinting.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Assign @nb.njit(cache.
def _find_neighbors_jit(indices: np.ndarray, distances: np.ndarray, dist: float, m: int, nvec: int) -> np.ndarray:
    """
    _find_neighbors_jit
    Numba JIT-compiled helper to extract valid neighbors from KDTree candidate pool.
    
    Input:
        indices     m x k           ndarray         Array with neighbor indices from KDTree query.
        distances   m x k           ndarray         Array with neighbor distances from KDTree query.
        dist                        float           Search radius.
        m                           int             Total number of nodes.
        nvec                        int             Maximum number of neighbors per node.
        
    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.
    """
    vec = np.zeros((m, nvec), dtype=np.int64) - 1                                                                                       # Initialize neighbor matrix with padding -1.
    
    for i in nb.prange(m):                                                                                                              # type: ignore
        idx_count = 0                                                                                                                   # Counter for valid neighbors.
        
        for j in range(indices.shape[1]):                                                                                               # Iterate over queried neighbors.
            idx = indices[i, j]                                                                                                         # Fetch neighbor index.
            d   = distances[i, j]                                                                                                       # Fetch neighbor distance.
            
            if idx >= 0 and idx < m and idx != i and d < dist:                                                                          # Valid neighbor condition.
                vec[i, idx_count] = idx                                                                                                 # Store neighbor.
                idx_count        += 1                                                                                                   # Increment counter.
                
                if idx_count >= nvec:                                                                                                   # Stop if we reached nvec.
                    break                                                                                                               # Exit neighbor loop.
    return vec                                                                                                                          # Return populated neighbor list.

@nb.njit(cache=True, fastmath=True)                                                                                                     # Assign @nb.njit(cache.
def _check_stencil_condition_jit(dx: np.ndarray, dy: np.ndarray) -> float:
    """
    _check_stencil_condition_jit
    Numba JIT-compiled helper to compute the condition number of the local geometric matrix.
    
    Input:
        dx          n x 1           ndarray         x-offsets of the selected neighbors.
        dy          n x 1           ndarray         y-offsets of the selected neighbors.
        
    Output:
        cond                        float           Condition number of the M^T M geometric matrix.
    """
    # 1. Variable initialization
    n = len(dx)                                                                                                                         # Get the number of provided neighbors.
    if n < 5:                                                                                                                           # Minimum neighbors for 2D Taylor polynomial.
        return np.inf                                                                                                                   # Return infinite condition number.
        
    a = 0.0                                                                                                                             # Initialization of the (0,0) covariance element.
    b = 0.0                                                                                                                             # Initialization of the (0,1) covariance element.
    c = 0.0                                                                                                                             # Initialization of the (1,1) covariance element.
    
    # 2. Covariance matrix computation
    for j in range(n):                                                                                                                  # Iterate over all local neighbors.
        a += dx[j] * dx[j]                                                                                                              # Accumulate xx cross term.
        b += dx[j] * dy[j]                                                                                                              # Accumulate xy cross term.
        c += dy[j] * dy[j]                                                                                                              # Accumulate yy cross term.
        
    T = a + c                                                                                                                           # Trace of the 2x2 covariance matrix.
    D = a * c - b * b                                                                                                                   # Determinant of the 2x2 covariance matrix.
    
    if D <= 0.0:                                                                                                                        # If determinant is negative or zero.
        return np.inf                                                                                                                   # Return infinite condition number.
        
    # 3. Eigenvalues and condition number
    sqrt_disc = np.sqrt(max(0.0, T * T - 4 * D))                                                                                        # Square root of the discriminant.
    L1        = (T + sqrt_disc) / 2.0                                                                                                   # Largest eigenvalue computation.
    L2        = (T - sqrt_disc) / 2.0                                                                                                   # Smallest eigenvalue computation.
    
    if L2 <= 1e-14 * L1:                                                                                                                # Check for extremely ill-conditioned matrix.
        return np.inf                                                                                                                   # Return infinite condition number.
        
    return L1 / L2                                                                                                                      # Return geometric condition number.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Assign @nb.njit(cache.
def _find_neighbors_balanced_jit(p: np.ndarray, indices: np.ndarray, distances: np.ndarray, dist: float, m: int, nvec: int, target_per_quad: int) -> np.ndarray:
    """
    _find_neighbors_balanced_jit
    Numba JIT-compiled helper to extract quadrant-balanced and condition-aware neighbors.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        indices     m x k           ndarray         Array with neighbor indices from KDTree query.
        distances   m x k           ndarray         Array with neighbor distances from KDTree query.
        dist                        float           Search radius.
        m                           int             Total number of nodes.
        nvec                        int             Maximum number of neighbors per node.
        target_per_quad             int             Number of neighbors to attempt extracting per quadrant.
        
    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.
    """
    vec      = np.zeros((m, nvec), dtype=np.int64) - 1                                                                                  # Initialize neighbor matrix with padding -1.
    max_cand = indices.shape[1]                                                                                                         # Maximum number of candidates per node.
    
    for i in nb.prange(m):                                                                                                              # type: ignore
        cand       = np.zeros(max_cand, dtype=np.int64)                                                                                 # Allocate array for valid candidates.
        dist_c     = np.zeros(max_cand, dtype=np.float64)                                                                               # Allocate array for candidate distances.
        cand_count = 0                                                                                                                  # Counter for valid candidates.
        
        for j in range(1, max_cand):                                                                                                    # Skip self at j=0.
            idx = indices[i, j]                                                                                                         # Fetch candidate index.
            d   = distances[i, j]                                                                                                       # Fetch candidate distance.
            
            if idx >= 0 and idx < m and d <= 5.0 * dist:                                                                                # Generous distance constraint for stabilization.
                cand[cand_count]   = idx                                                                                                # Store candidate index.
                dist_c[cand_count] = d                                                                                                  # Store candidate distance.
                cand_count        += 1                                                                                                  # Increment candidate count.
        
        if cand_count > 0:                                                                                                              # Process if candidates exist.
            q1                     = np.zeros(cand_count, dtype=np.int64)                                                               # Array for NE quadrant.
            q2                     = np.zeros(cand_count, dtype=np.int64)                                                               # Array for NW quadrant.
            q3                     = np.zeros(cand_count, dtype=np.int64)                                                               # Array for SW quadrant.
            q4                     = np.zeros(cand_count, dtype=np.int64)                                                               # Array for SE quadrant.
            q1_c, q2_c, q3_c, q4_c = 0, 0, 0, 0                                                                                         # Counters for quadrants.
            
            for j in range(cand_count):                                                                                                 # Sort candidates into quadrants.
                idx = cand[j]                                                                                                           # Candidate global index.
                dx  = p[idx, 0] - p[i, 0]                                                                                               # DX relative to node i.
                dy  = p[idx, 1] - p[i, 1]                                                                                               # DY relative to node i.
                
                if dx >= 0 and dy >= 0:                                                                                                 # NE quadrant.
                    q1[q1_c] = idx                                                                                                      # Assign to NE.
                    q1_c    += 1                                                                                                        # Increment NE counter.
                elif dx < 0 and dy >= 0:                                                                                                # NW quadrant.
                    q2[q2_c] = idx                                                                                                      # Assign to NW.
                    q2_c    += 1                                                                                                        # Increment NW counter.
                elif dx < 0 and dy < 0:                                                                                                 # SW quadrant.
                    q3[q3_c] = idx                                                                                                      # Assign to SW.
                    q3_c    += 1                                                                                                        # Increment SW counter.
                else:                                                                                                                   # SE quadrant.
                    q4[q4_c] = idx                                                                                                      # Assign to SE.
                    q4_c += 1                                                                                                           # Increment SE counter.
            
            selected = np.zeros(nvec, dtype=np.int64) - 1                                                                               # Array to store selected neighbors.
            sel_c    = 0                                                                                                                # Counter for selected neighbors.
            q_idx    = np.zeros(4, dtype=np.int64)                                                                                      # Index trackers for the 4 quadrants.
            q_lens   = np.array([q1_c, q2_c, q3_c, q4_c], dtype=np.int64)                                                               # Number of points in each quadrant.
            
            condition = 1e6                                                                                                             # Initialize condition to prevent Numba undefined behavior.
            dx_arr    = np.zeros(nvec, dtype=np.float64)                                                                                # Preallocate array for DX.
            dy_arr    = np.zeros(nvec, dtype=np.float64)                                                                                # Preallocate array for DY.
            
            while sel_c < nvec:                                                                                                         # Until we hit maximum allowed neighbors.
                added_any = False                                                                                                       # Flag to check if we added a new point in this round.
                
                for quad in range(4):                                                                                                   # Round-robin over quadrants.
                    if q_idx[quad] < q_lens[quad]:                                                                                      # If quadrant has available points.
                        if quad == 0:                                                                                                   # Quadrant NE.
                            idx = q1[q_idx[quad]]                                                                                       # Fetch candidate index from NE.
                        elif quad == 1:                                                                                                 # Quadrant NW.
                            idx = q2[q_idx[quad]]                                                                                       # Fetch candidate index from NW.
                        elif quad == 2:                                                                                                 # Quadrant SW.
                            idx = q3[q_idx[quad]]                                                                                       # Fetch candidate index from SW.
                        else:                                                                                                           # Quadrant SE.
                            idx = q4[q_idx[quad]]                                                                                       # Fetch candidate index from SE.
                            
                        selected[sel_c] = idx                                                                                           # Add point to selected list.
                        sel_c          += 1                                                                                             # Increment selected counter.
                        q_idx[quad]    += 1                                                                                             # Advance quadrant index.
                        added_any       = True                                                                                          # Mark as added.
                        
                        if sel_c >= nvec:                                                                                               # If we hit the absolute maximum.
                            break                                                                                                       # Exit quadrant loop.
                            
                        if sel_c >= 9:                                                                                                  # Minimum neighbors required for stability.
                            for k in range(sel_c):                                                                                      # Collect relative coordinates.
                                dx_arr[k] = p[selected[k], 0] - p[i, 0]                                                                 # DX computation.
                                dy_arr[k] = p[selected[k], 1] - p[i, 1]                                                                 # DY computation.
                            condition = _check_stencil_condition_jit(dx_arr[:sel_c], dy_arr[:sel_c])                                    # Check 2x2 spatial condition.
                            if condition < 100.0:                                                                                       # If spatial condition is well distributed.
                                break                                                                                                   # Stop adding neighbors dynamically!
                                
                if sel_c >= 9 and condition < 100.0:                                                                                    # Check again after breaking the inner loop.
                    break                                                                                                               # Stop adding neighbors dynamically!
                    
                if not added_any:                                                                                                       # If no more candidates can be added from ANY quadrant.
                    break                                                                                                               # Break condition loop.
            
            # 3. Transfer to final vec array
            for j in range(sel_c):                                                                                                      # Transfer to vec.
                vec[i, j] = selected[j]                                                                                                 # Save selected neighbors.
    
    return vec                                                                                                                          # Return dynamically balanced neighbor list.

@nb.njit(cache=True, fastmath=True, parallel=True)                                                                                      # Assign @nb.njit(cache.
def _find_neighbors_adv_jit(p: np.ndarray, indices: np.ndarray, distances: np.ndarray, a: float, b: float, tol: float, m: int, nvec: int) -> np.ndarray:
    """
    _find_neighbors_adv_jit
    Numba JIT-compiled helper to extract upwind-biased neighbors from KDTree candidate pool.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        indices     m x k           ndarray         Array with neighbor indices from KDTree query.
        distances   m x k           ndarray         Array with neighbor distances from KDTree query.
        a, b                        float           Advection direction components.
        tol                         float           Tolerance for upwind dot product.
        m                           int             Total number of nodes.
        nvec                        int             Maximum number of neighbors per node.
        
    Output:
        vec         m x nvec        ndarray         Array with the correspondence of the neighbors of each node.
    """
    vec      = np.zeros((m, nvec), dtype=np.int64) - 1                                                                                  # Initialize neighbor matrix with padding -1.
    max_cand = indices.shape[1]                                                                                                         # Maximum number of candidates per node.
    
    for i in nb.prange(m):                                                                                                              # type: ignore
        cand       = np.zeros(max_cand, dtype=np.int64)                                                                                 # Allocate array for valid candidates.
        dist_c     = np.zeros(max_cand, dtype=np.float64)                                                                               # Allocate array for candidate distances.
        cand_count = 0                                                                                                                  # Counter for valid candidates.
        
        for j in range(1, max_cand):                                                                                                    # Skip self at j=0.
            idx = indices[i, j]                                                                                                         # Fetch candidate index.
            d   = distances[i, j]                                                                                                       # Fetch candidate distance.
            
            if idx >= 0 and idx < m and d < 1e10:                                                                                       # Check valid index constraint.
                cand[cand_count]   = idx                                                                                                # Store candidate index.
                dist_c[cand_count] = d                                                                                                  # Store candidate distance.
                cand_count        += 1                                                                                                  # Increment candidate count.
        
        if cand_count > 0:                                                                                                              # Process if candidates exist.
            for j in range(cand_count):                                                                                                 # Simple bubble sort by distance.
                for k in range(0, cand_count - j - 1):                                                                                  # Bubble pass.
                    if dist_c[k] > dist_c[k+1]:                                                                                         # If out of order.
                        tmp_d       = dist_c[k]                                                                                         # Swap distance.
                        dist_c[k]   = dist_c[k+1]                                                                                       # Swap distance.
                        dist_c[k+1] = tmp_d                                                                                             # Swap distance.
                        tmp_idx     = cand[k]                                                                                           # Swap index.
                        cand[k]     = cand[k+1]                                                                                         # Swap index.
                        cand[k+1]   = tmp_idx                                                                                           # Swap index.
            
            up         = np.zeros(cand_count, dtype=np.int64)                                                                           # Upwind candidates array.
            dn         = np.zeros(cand_count, dtype=np.int64)                                                                           # Downwind candidates array.
            up_c, dn_c = 0, 0                                                                                                           # Upwind and downwind counters.
            
            for j in range(cand_count):                                                                                                 # Classify by upwind direction.
                idx = cand[j]                                                                                                           # Candidate index.
                dx  = p[idx, 0] - p[i, 0]                                                                                               # DX offset.
                dy  = p[idx, 1] - p[i, 1]                                                                                               # DY offset.
                
                if dx * a + dy * b <= tol:                                                                                              # Upwind condition.
                    up[up_c]  = idx                                                                                                     # Add to upwind list.
                    up_c     += 1                                                                                                       # Increment upwind count.
                else:                                                                                                                   # Downwind condition.
                    dn[dn_c]  = idx                                                                                                     # Add to downwind list.
                    dn_c     += 1                                                                                                       # Increment downwind count.
            
            if up_c >= nvec:                                                                                                            # Enough upwind candidates.
                for j in range(nvec):                                                                                                   # Take top nvec upwind.
                    vec[i, j] = up[j]                                                                                                   # Store in vec.
            else:                                                                                                                       # Not enough upwind candidates.
                for j in range(up_c):                                                                                                   # Take all upwind.
                    vec[i, j] = up[j]                                                                                                   # Store in vec.
                
                rem     = nvec - up_c                                                                                                   # Remaining slots.
                take_dn = min(rem, dn_c)                                                                                                # How many downwind to take.
                
                for j in range(take_dn):                                                                                                # Take downwind.
                    vec[i, up_c + j] = dn[j]                                                                                            # Store in vec.
    
    return vec                                                                                                                          # Return upwind neighbor list.

def compute_neighbors(p: np.ndarray, nvec: int) -> np.ndarray:
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
    # 1. Variable initialization
    m    = int(p.shape[0])                                                                                                              # Total number of nodes.
    dist = find_distances(p)                                                                                                            # Search radius from nearest-neighbor spacing.

    # 2. Search radius scaling
    if nvec > 12:                                                                                                                       # Scale dist if asking for expanded stencil.
        dist = dist * np.sqrt(nvec / 12.0)                                                                                              # Area scales with r^2, so r scales with sqrt(n).

    # 3. Neighbor search
    vec  = find_neighbors_balanced(p, dist, nvec)                                                                                       # Build balanced neighbor list using KDTree candidates.

    return vec                                                                                                                          # Return global neighbor list.

def compute_upwind_neighbors(p: np.ndarray, a: float, b: float, nvec: int) -> np.ndarray:
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
    m    = int(p.shape[0])                                                                                                              # Total number of nodes.
    dist = find_distances(p)                                                                                                            # Search radius from nearest-neighbor spacing.

    # 2. Search radius scaling
    if nvec > 12:                                                                                                                       # Scale dist if asking for expanded stencil.
        dist = dist * np.sqrt(nvec / 12.0)                                                                                              # Area scales with r^2, so r scales with sqrt(n).

    # 3. Upwind neighbor search
    vec  = find_neighbors_adv(p, dist, a, b, nvec)                                                                                      # Build neighbor list using upwind KDTree query.

    return vec                                                                                                                          # Return global neighbor list.

def find_distances(p: np.ndarray) -> float:
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
    # 1. KDTree-based search
    tree          = KDTree(p[:, :2])                                                                                                    # Build KDTree on (x, y) coordinates.
    distances, _  = tree.query(p[:, :2], k = 2, workers=-1)                                                                             # Query self + nearest neighbor.
    min_distances = distances[:, 1]                                                                                                     # Keep nearest-neighbor distances (skip self at k=0).
    dist          = (3 / 2) * float(np.max(min_distances))                                                                              # Convert spacing to a search radius.

    return dist                                                                                                                         # Return radius distance.

def find_neighbors(p: np.ndarray, dist: float, nvec: int) -> np.ndarray:
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
    # 1. Variable initialization
    m    = int(p.shape[0])                                                                                                              # Total number of nodes.
    
    if m == 0 or nvec <= 0:                                                                                                             # Edge case check.
        return np.zeros((m, nvec), dtype=np.int64) - 1                                                                                  # Return empty matrix.
        
    # 2. KDTree vectorized query
    tree               = KDTree(p[:, :2])                                                                                               # Build KDTree on (x, y) coordinates.
    distances, indices = tree.query(p[:, :2], k=nvec+1, distance_upper_bound=dist, workers=-1)                                          # Bulk query all nodes.
    
    # 3. JIT-accelerated filtering and assignment
    vec = _find_neighbors_jit(indices, distances, dist, m, nvec)                                                                        # Delegate array manipulation to Numba.
    
    return vec                                                                                                                          # Return neighbor list.

def find_neighbors_balanced(p: np.ndarray, dist: float, nvec: int) -> np.ndarray:
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
    # 1. Variable initialization
    m    = int(p.shape[0])                                                                                                              # Total number of nodes.
    
    if m == 0 or nvec <= 0:                                                                                                             # Edge case check.
        return np.zeros((m, nvec), dtype=np.int64) - 1                                                                                  # Return empty matrix.
        
    # 2. KDTree vectorized candidate pool query
    k_cand             = min(m, 10 * nvec)                                                                                              # Expand candidate search pool.
    tree               = KDTree(p[:, :2])                                                                                               # Build KDTree on (x, y).
    distances, indices = tree.query(p[:, :2], k=k_cand, workers=-1)                                                                     # Bulk query all candidates.
    
    if k_cand == 1:                                                                                                                     # Handle edge case.
        distances = distances.reshape(-1, 1)                                                                                            # Reshape distances.
        indices   = indices.reshape(-1, 1)                                                                                              # Reshape indices.
        
    target_per_quad = min(3, int(np.ceil(nvec / 4.0)))                                                                                  # Target number of points per quadrant (max 3 for base seed).
    
    # 3. JIT-accelerated quadrant-balanced selection
    vec = _find_neighbors_balanced_jit(np.asarray(p, dtype=np.float64), indices, distances, dist, m, nvec, target_per_quad)             # Delegate balancing to Numba.
    
    return vec                                                                                                                          # Return balanced neighbor list.

def find_neighbors_adv(p: np.ndarray, dist: float, a: float, b: float, nvec: int) -> np.ndarray:
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
    # 1. Variable initialization
    m    = int(p.shape[0])                                                                                                              # Total number of nodes.
    
    if m == 0 or nvec <= 0:                                                                                                             # Handle empty input or invalid neighbor count.
        return np.zeros((m, nvec), dtype=np.int64) - 1                                                                                  # Return empty neighbor matrix.

    # 2. Advection direction magnitude
    speed = float(np.hypot(a, b))                                                                                                       # Magnitude of advection direction.
    
    if not np.isfinite(speed) or speed == 0.0:                                                                                          # Fallback when direction is invalid.
        dist0 = find_distances(p)                                                                                                       # Use isotropic distance estimate.
        return find_neighbors(p, dist0, nvec)                                                                                           # Return isotropic KDTree neighbors.

    # 3. KDTree vectorized candidate query
    k                  = min(m, max(nvec + 1, 128, 40 * nvec))                                                                          # Candidate pool size for KDTree query.
    tree               = KDTree(p[:, :2])                                                                                               # Build KDTree on (x, y) coordinates.
    distances, indices = tree.query(p[:, :2], k=k, workers=-1)                                                                          # Query candidates.

    if k == 1:                                                                                                                          # Ensure 2D shapes for consistent slicing.
        distances = distances.reshape(-1, 1)                                                                                            # Reshape distances.
        indices   = indices.reshape(-1, 1)                                                                                              # Reshape indices.

    tol = 1e-14 * speed                                                                                                                 # Tolerance to decide upwind vs downwind.

    # 4. JIT-accelerated upwind filtering and selection
    vec = _find_neighbors_adv_jit(np.asarray(p, dtype=np.float64), indices, distances, a, b, tol, m, nvec)                              # Delegate filtering and assignment to Numba.

    return vec                                                                                                                          # Return upwind neighbor list.

def compute_mesh_spacing(p: np.ndarray, vec: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """
    compute_mesh_spacing
    Compute characteristic spatial node spacing (h_min, h_avg) for a 2D point cloud.

    Input:
        p           m x 3           ndarray         Point cloud array [x, y, flag].
        vec         m x nvec        ndarray         Cached neighbor matrix (optional).

    Output:
        h_min                       float           Minimum spatial step between neighbor nodes.
        h_avg                       float           Average spatial step between neighbor nodes.
    """
    tree         = KDTree(p[:, :2])                                                                                                     # Build KDTree on (x, y).
    distances, _ = tree.query(p[:, :2], k=2, workers=-1)                                                                                # Query self + nearest neighbor.
    min_dists    = distances[:, 1]                                                                                                      # Distance to nearest neighbor for each node.
    h_min        = float(np.min(min_dists))                                                                                             # Minimum node spacing.
    h_avg        = float(np.mean(min_dists))                                                                                            # Average node spacing.
    return h_min, h_avg                                                                                                                 # Return spatial metrics.

