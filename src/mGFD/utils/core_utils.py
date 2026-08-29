"""
Utils — General utility functions for mGFD
    
Overview:
    This module provides generic helper functions for point cloud processing,
    geometry analysis, and general utilities used across the mGFD codebase.
    
Public API:
    poly_area                   Compute polygon area by the shoelace formula.
    get_valid_triangulation     Compute a Delaunay triangulation respecting boundaries.
    get_aspect_and_bounds       Calculate physical bounds and aspect ratio for a point cloud.
    
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

## Library importation.
import os                                                                                                                               # OS interfaces for file/directory paths.
import numpy as np                                                                                                                      # Core numerical operations.

from scipy.spatial import Delaunay                                                                                                      # Delaunay triangulation algorithm.
from typing import Callable, Optional, Tuple, List                                                                                      # Type hinting.

from mGFD.spatial.neighbors import find_distances                                                                                       # Distance estimator.

def poly_area(x: np.ndarray, y: np.ndarray) -> float:
    """
    poly_area
    Compute the area of a polygon defined by its ordered vertices (x, y).
    
    This function uses the shoelace formula and returns the absolute area.
    
    Input:
        x           (n,)            ndarray         x-coordinates of polygon vertices.
        y           (n,)            ndarray         y-coordinates of polygon vertices.
    
    Output:
        area                        float           Area of the polygon.
    
    Notes:
        The vertices must be ordered along the polygon boundary (clockwise or counterclockwise).
        If the polygon is self-intersecting or the vertex order is arbitrary, the result may be incorrect.
    """
    # 1. Area computation
    x    = np.asarray(x, dtype = float)                                                                                                 # Normalize x to ndarray.
    y    = np.asarray(y, dtype = float)                                                                                                 # Normalize y to ndarray.
    area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))                                                            # Shoelace formula (absolute area).
    
    return float(area)                                                                                                                  # Return scalar area.


def get_valid_triangulation(p: np.ndarray, nom: Optional[str] = None) -> Optional[np.ndarray]:
    """
    get_valid_triangulation
    Compute a valid triangulation for the point cloud p by filtering out Delaunay
    triangles whose baricenters lie outside the water domain (defined in *_contours.csv).
    
    It processes regions (islands/disconnected bodies) independently to prevent triangles
    from crossing between them, and uses a local cache file to optimize performance.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
        nom                         str             Base filename or path to identify the dataset and scale.
    
    Output:
        global_triangles            ndarray         Array of shape (n_triangles, 3) containing valid Delaunay simplices, or None if failed.
    """
    # 1. Variable initialization
    cache_path = None                                                                                                                   # Initialize cache_path variable.
    
    # 2. Cache parsing
    if nom:                                                                                                                             # Check if a valid filename was provided to attempt caching.
        parts          = nom.split(os.sep)                                                                                              # Split path to parse dataset details.
        dataset        = None                                                                                                           # Initialize dataset variable.
        scale          = None                                                                                                           # Initialize scale variable.
        workspace_root = None                                                                                                           # Initialize workspace root variable.
    
        if len(parts) >= 4:                                                                                                             # Check if path has enough components.
            scale   = parts[-2]                                                                                                         # Extract scale name.
            dataset = parts[-3]                                                                                                         # Extract dataset name.
    
            temp    = nom                                                                                                               # Initialize temp path for workspace search.
            
            for _ in range(10):                                                                                                         # Climb up directory tree.
                temp = os.path.dirname(temp)                                                                                            # Go up one directory level.
                
                if os.path.exists(os.path.join(temp, 'Data')):                                                                          # Check if Data folder exists.
                    workspace_root = temp                                                                                               # Set workspace root.
                    break                                                                                                               # Stop climbing directory tree.
    
        if workspace_root and dataset and scale:                                                                                        # Check if path info was successfully parsed.
            cloud_file = f"{dataset}_cloud.csv"                                                                                         # Format cloud filename.
            cloud_path = os.path.join(workspace_root, 'Data', dataset, scale, cloud_file)                                               # Construct cloud path.
            
            if not os.path.exists(cloud_path):                                                                                          # Check if cloud file exists.
                lake_scale_dir = os.path.join(workspace_root, 'Data', dataset, scale)                                                   # Construct fallback directory.
                
                if os.path.exists(lake_scale_dir):                                                                                      # Check if fallback directory exists.
                    
                    for f in os.listdir(lake_scale_dir):                                                                                # Iterate over fallback directory files.
                        
                        if f.endswith('.csv') and '_neighbors' not in f and '_triangulation' not in f:                                  # Filter valid cloud files.
                            cloud_path = os.path.join(lake_scale_dir, f)                                                                # Construct cloud path.
                            break                                                                                                       # Stop climbing directory tree.
            
            if cloud_path and os.path.exists(cloud_path):                                                                               # Check if valid cloud path was found.
                cache_path = cloud_path.replace('.csv', '_triangulation.csv')                                                           # Set cache path based on cloud path.
    
        # 3. Cache loading
        if cache_path and os.path.exists(cache_path):                                                                                   # Check if cache file exists.
            try:                                                                                                                        # Attempt operation.
                triangles = np.loadtxt(cache_path, dtype = np.int32, delimiter = ',')                                                   # Load triangulation from cache.
                
                if triangles.ndim == 1:                                                                                                 # Check if triangles is 1D.
                    triangles = triangles.reshape(1, -1)                                                                                # Reshape triangles to 2D.
                
                return triangles                                                                                                        # Return cached triangulation.
            
            except Exception:                                                                                                           # Ignore exceptions.
                pass                                                                                                                    # Do nothing.

    # 4. Coordinates extraction
    x                = p[:, 0]                                                                                                          # Extract X coordinates directly from p.
    y                = p[:, 1]                                                                                                          # Extract Y coordinates directly from p.

    dist             = find_distances(np.column_stack([x, y, np.zeros_like(x)]))                                                        # Estimate typical point spacing.
    threshold        = 1.3 * dist                                                                                                       # Max allowed edge length to preserve concavities perfectly.

    global_triangles = []                                                                                                               # Initialize empty list for valid triangles.
    xy               = np.column_stack([x, y])                                                                                          # Create 2D coordinate array.
    
    # 5. Triangulation
    try:                                                                                                                                # Attempt operation.
        tri       = Delaunay(xy)                                                                                                        # Compute Delaunay triangulation.
        simplices = tri.simplices                                                                                                       # Extract Delaunay simplices.
        
        for t in simplices:                                                                                                             # Iterate over all triangles.
            pts = xy[t]                                                                                                                 # Get vertices of current triangle.
            L1  = np.linalg.norm(pts[0] - pts[1])                                                                                       # Calculate edge length.
            L2  = np.linalg.norm(pts[1] - pts[2])                                                                                       # Calculate edge length.
            L3  = np.linalg.norm(pts[2] - pts[0])                                                                                       # Calculate edge length.
            if max(L1, L2, L3) < threshold:                                                                                             # Filter out large boundary triangles.
                global_triangles.append(t)                                                                                              # Store valid triangle.
    except Exception:                                                                                                                   # Ignore exceptions.
        pass                                                                                                                            # Do nothing.

    # 6. Return and cache saving
    if len(global_triangles) > 0:                                                                                                       # Check if any valid triangles were found.
        global_triangles_arr = np.array(global_triangles, dtype = np.int32)                                                             # Convert list to NumPy array.
        
        if cache_path:                                                                                                                  # Check if cache path is defined.
            try:                                                                                                                        # Attempt operation.
                np.savetxt(cache_path, global_triangles_arr, delimiter = ',', fmt = '%d')                                               # Save triangulation to cache.
            except Exception:                                                                                                           # Ignore exceptions.
                pass                                                                                                                    # Do nothing.
        
        return global_triangles_arr                                                                                                     # Return computed triangulation.

    return None                                                                                                                         # Return None.


def get_aspect_and_bounds(p: np.ndarray) -> Tuple[Tuple[float, float, float], List[float], List[float]]:
    """
    get_aspect_and_bounds
    Extract physical aspect ratio and bounding box from a point cloud.
    
    This helper is used to determine isometric scaling and layout parameters
    for 3D plotting and geometry processing.
    
    Input:
        p           m x 3           ndarray         Array with the coordinates of the nodes and the boundary flag.
    
    Output:
        box_aspect                  Tuple           Tuple (x_range, y_range, z_scale) representing physical proportions.
        x_bounds                    List            List [x_min, x_max] of global X bounds.
        y_bounds                    List            List [y_min, y_max] of global Y bounds.
    """
    # 1. Variable initialization
    x_min, x_max     = p[:, 0].min(), p[:, 0].max()                                                                                     # Extract global X bounds.
    y_min, y_max     = p[:, 1].min(), p[:, 1].max()                                                                                     # Extract global Y bounds.
    x_range, y_range = x_max - x_min, y_max - y_min                                                                                     # Calculate physical X and Y ranges.
    max_range        = max(x_range, y_range) if max(x_range, y_range) > 0 else 1.0                                                      # Determine maximum range for isometric scaling.
    box_aspect       = (x_range, y_range, 0.4 * max_range)                                                                              # Set 3D box aspect ratio to preserve geometry.
    
    return box_aspect, [x_min, x_max], [y_min, y_max]                                                                                   # Return layout parameters.
