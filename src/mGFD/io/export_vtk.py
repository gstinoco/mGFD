"""
ExportVTK — VTK export utilities for mGFD

Overview:
    This module provides functions to export point cloud simulation results into VTK formats (.vtp, .pvd)
    for high-performance, interactive 3D visualization in ParaView.
    It replaces the slow and disk-heavy CSV saving method.

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.
    u_ap    (m,) or (m, t) ndarray
            Approximate numerical solution.
    u_ap    (m,) or (m, t) ndarray
            Approximate numerical solution.

Public API:
    export_stationary_vtk     Export a single-step (stationary) result to a .vtp file.
    export_transient_vtk      Export a multi-step (transient) result to a .pvd series of .vtp files.

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
import logging                                                                                                                          # Standard logging module.
import numpy as np                                                                                                                      # Core numerical operations.
import pyvista as pv                                                                                                                    # 3D visualization utilities.

from typing import Callable, Optional, Tuple, List                                                                                      # Type hinting.

from mGFD.utils.core_utils import get_valid_triangulation                                                                               # Geometry utilities.

logger = logging.getLogger(__name__)                                                                                                    # Module level logger.

def _create_mesh(p: np.ndarray, cloud_path: Optional[str] = None) -> pv.PolyData:
    """
    _create_mesh
    Helper function to construct a PyVista PolyData mesh from point coordinates and triangles.
    
    Input:
        p           m x 3           ndarray         Point cloud coordinates and boundary flags.
        cloud_path                  str             (Optional) Path to external triangulation source.
        
    Output:
        mesh                        pv.PolyData     The constructed PyVista mesh with bound point data.
    """
    # 1. 3D point creation
    points    = np.column_stack([p[:, 0], p[:, 1], np.zeros(p.shape[0])])                                                               # Create 3D points with Z=0.
    
    # 2. Triangulation generation
    triangles = get_valid_triangulation(p, cloud_path)                                                                                  # Retrieve optimal triangulation.
        
    # 3. Mesh construction
    if triangles is not None:                                                                                                           # If triangulation is available.
        faces        = np.empty((triangles.shape[0], 4), dtype = np.int32)                                                              # Preallocate array for faces.
        faces[:, 0]  = 3                                                                                                                # Set padding indicating 3 points per face.
        faces[:, 1:] = triangles                                                                                                        # Set vertices indices.
        mesh         = pv.PolyData(points, faces)                                                                                       # type: ignore
    else:                                                                                                                               # Fallback to point cloud if triangulation fails.
        mesh         = pv.PolyData(points)                                                                                              # Create PyVista mesh as point cloud.
        
    mesh.point_data['Flag'] = p[:, 2]                                                                                                   # Store boundary flag in the mesh.
    
    return mesh                                                                                                                         # Return created mesh.

def export_stationary_vtk(p: np.ndarray, u_ap: np.ndarray, out_dir: str, basename: str = "Stationary_Solution", cloud_path: Optional[str] = None, verbose: bool = True) -> None:
    """
    export_stationary_vtk
    Export a stationary PDE solution to a VTK (.vtp) file.
    
    This function computes the absolute error, constructs a 3D PolyData mesh (with Z=0),
    and attaches the approximate solution, exact solution, and error as point data.
    
    Input:
        p           m x 3           ndarray         Point cloud coordinates and boundary flags.
        u_ap        m               ndarray         Approximate solution values.
        out_dir                     str             Directory where the VTK file will be saved.
        basename                    str             Base filename for the output file.
        cloud_path                  str             (Optional) Path to the point cloud CSV, used to load cached boundary-aware triangulations.
        verbose                     bool            (Optional) If True, prints status messages to standard output.
    """
    # 1. Initialization and mesh creation
    os.makedirs(out_dir, exist_ok = True)                                                                                               # Ensure output directory exists.
    mesh = _create_mesh(p, cloud_path = cloud_path)                                                                                     # Create base mesh for the geometry.
    
    # 2. Data attachment
    mesh.point_data['U_ap']           = u_ap                                                                                            # Store approximate solution.
    
    # 3. VTK saving
    filepath = os.path.join(out_dir, f"{basename}.vtp")                                                                                 # Assemble full output path.
    mesh.save(filepath)                                                                                                                 # Save VTP file to disk.
    
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f"\tSaved VTK to {filepath}")                                                                                       # Report successful save.

def export_transient_vtk(p: np.ndarray, u_ap: np.ndarray, t: int, T: np.ndarray, out_dir: str, basename: str = "Transient_Solution", cloud_path: Optional[str] = None, verbose: bool = True) -> None:
    """
    export_transient_vtk
    Export a transient PDE solution to a time-series VTK (.pvd + .vtp) format.
    
    This function constructs a PyVista mesh for each saved time step and links them together
    using a PVD (ParaView Data) file, which allows ParaView to play the sequence as an animation.
    
    Input:
        p           m x 3           ndarray         Point cloud coordinates and boundary flags.
        u_ap        m x t           ndarray         Approximate solution values over time.
        t                           int             Total number of time steps.
        T           t               ndarray         Array of physical time values.
        out_dir                     str             Directory where the VTK files will be saved.
        basename                    str             Base filename for the output PVD and VTP files.
        cloud_path                  str             (Optional) Path to the point cloud CSV, used to load cached boundary-aware triangulations.
        verbose                     bool            (Optional) If True, prints status messages to standard output.
    """
    # 1. Initialization and mesh creation
    out_dir = os.path.join(out_dir, 'VTK')                                                                                              # Append VTK subdirectory to avoid cluttering.
    os.makedirs(out_dir, exist_ok = True)                                                                                               # Ensure output directory exists.
    mesh = _create_mesh(p, cloud_path = cloud_path)                                                                                     # Create base mesh for the geometry.
    
    # 2. PVD Initialization
    pvd_content = [                                                                                                                     # Initialize PVD XML structure.
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '  <Collection>'
    ]
    
    step         = max(1, t // 500)                                                                                                     # Calculate step size to prevent massive storage.
    frames_saved = 0                                                                                                                    # Initialize frame counter.
    
    # 3. Frame iteration
    for k in range(0, t, step):                                                                                                         # Iterate over chosen time steps.
        time_val = float(T[k])                                                                                                          # Physical time value.
        
        mesh.point_data['U_ap']           = u_ap[:, k]                                                                                  # Store approximate solution at time k.
        
        filename = f"{basename}_{k:05d}.vtp"                                                                                            # Frame filename.
        filepath = os.path.join(out_dir, filename)                                                                                      # Assemble full path for the frame.
        mesh.save(filepath)                                                                                                             # Save individual VTP frame.
        
        pvd_content.append(f'    <DataSet timestep="{time_val}" group="" part="0" file="{filename}"/>')                                 # Register frame in PVD file.
        frames_saved += 1                                                                                                               # Increment frame counter.
        
    # 4. PVD finalization
    pvd_content.append('  </Collection>')                                                                                               # Close Collection XML tag.
    pvd_content.append('</VTKFile>')                                                                                                    # Close VTKFile XML tag.
    
    # 5. PVD saving
    pvd_filepath = os.path.join(out_dir, f"{basename}.pvd")                                                                             # Assemble PVD file path.
    
    with open(pvd_filepath, "w") as f:                                                                                                  # Open PVD file for writing.
        f.write("\n".join(pvd_content))                                                                                                 # Write all lines to PVD.
        
    if verbose:                                                                                                                         # Check if verbosity is enabled.
        logger.info(f"\tSaved VTK series ({frames_saved} frames) to {pvd_filepath}")                                                    # Report successful save.
