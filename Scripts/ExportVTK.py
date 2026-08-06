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
    u_ex    (m,) or (m, t) ndarray
            Exact reference solution.

Public API:
    export_stationary_vtk     Export a single-step (stationary) result to a .vtp file.
    export_transient_vtk      Export a multi-step (transient) result to a .pvd series of .vtp files.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco Guerrero
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
import os
import numpy as np
import pyvista as pv
from scipy.spatial import Delaunay
from Scripts.Neighbors import find_distances

def _get_triangulation(p):
    """
    Compute a valid 2D Delaunay triangulation for the point cloud, filtering out large boundary triangles.
    """
    x = p[:, 0]
    y = p[:, 1]
    
    # Estimate point spacing to set a threshold for edge lengths
    pts_3d = np.column_stack([x, y, np.zeros_like(x)])
    dist = find_distances(pts_3d, mode=3)
    mean_spacing = np.mean(dist)
    threshold = 2.5 * mean_spacing

    xy = np.column_stack([x, y])
    try:
        tri = Delaunay(xy)
        simplices = tri.simplices
        
        valid_triangles = []
        for t in simplices:
            pts = xy[t]
            L1 = np.linalg.norm(pts[0] - pts[1])
            L2 = np.linalg.norm(pts[1] - pts[2])
            L3 = np.linalg.norm(pts[2] - pts[0])
            if max(L1, L2, L3) < threshold:
                valid_triangles.append(t)
                
        if valid_triangles:
            return np.array(valid_triangles, dtype=np.int32)
    except Exception:
        pass
    
    return None

def _create_mesh(p):
    """
    Create a PyVista PolyData object representing the domain mesh/cloud.
    """
    points = np.column_stack([p[:, 0], p[:, 1], np.zeros(p.shape[0])])
    triangles = _get_triangulation(p)
    
    if triangles is not None:
        # PyVista PolyData faces require padding with the number of points per face (3 for triangles)
        faces = np.empty((triangles.shape[0], 4), dtype=np.int32)
        faces[:, 0] = 3
        faces[:, 1:] = triangles
        mesh = pv.PolyData(points, faces)
    else:
        # Fallback to point cloud if triangulation fails
        mesh = pv.PolyData(points)
        
    mesh.point_data['Flag'] = p[:, 2]
    return mesh

def export_stationary_vtk(p, u_ap, u_ex, out_dir, basename="Stationary_Solution"):
    """
    Export stationary problem results to a VTK (.vtp) file.
    """
    os.makedirs(out_dir, exist_ok=True)
    mesh = _create_mesh(p)
    
    mesh.point_data['U_ap'] = u_ap
    mesh.point_data['U_ex'] = u_ex
    mesh.point_data['Absolute_Error'] = np.abs(u_ap - u_ex)
    
    filepath = os.path.join(out_dir, f"{basename}.vtp")
    mesh.save(filepath)
    print(f"\tSaved VTK to {filepath}")

def export_transient_vtk(p, u_ap, u_ex, t, T, out_dir, basename="Transient_Solution"):
    """
    Export transient problem results to a VTK time series (.pvd + .vtp files).
    """
    os.makedirs(out_dir, exist_ok=True)
    mesh = _create_mesh(p)
    
    pvd_content = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '  <Collection>'
    ]
    
    # Optional: Only save every N steps if there are too many (e.g. max 500 frames)
    step = max(1, t // 500)
    
    frames_saved = 0
    for k in range(0, t, step):
        time_val = float(T[k])
        
        # Update point data for this time step
        mesh.point_data['U_ap'] = u_ap[:, k]
        mesh.point_data['U_ex'] = u_ex[:, k]
        mesh.point_data['Absolute_Error'] = np.abs(u_ap[:, k] - u_ex[:, k])
        
        # Save VTP frame
        filename = f"{basename}_{k:05d}.vtp"
        filepath = os.path.join(out_dir, filename)
        mesh.save(filepath)
        
        pvd_content.append(f'    <DataSet timestep="{time_val}" group="" part="0" file="{filename}"/>')
        frames_saved += 1
        
    pvd_content.append('  </Collection>')
    pvd_content.append('</VTKFile>')
    
    pvd_filepath = os.path.join(out_dir, f"{basename}.pvd")
    with open(pvd_filepath, "w") as f:
        f.write("\n".join(pvd_content))
        
    print(f"\tSaved VTK series ({frames_saved} frames) to {pvd_filepath}")
