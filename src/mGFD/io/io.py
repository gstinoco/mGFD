"""
IO — Data and neighbor-list I/O helpers

Overview:
    Small utilities to:
    - Load a point cloud from CSV into the (m, 3) convention used by mGFD (x, y, flag)

Data conventions:
    p       (m, 3) ndarray
            Point cloud with columns [x, y, flag]. flag = 0 for interior; flag = 1/2 for boundary.

Public API:
    load_points         Load a point cloud from CSV into the (m, 3) format used by the solvers.

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
"""
import numpy as np

def load_points(file_path, verbose=False):
    """
    load_points
    Load a point cloud from a CSV file into the (m, 3) format used by mGFD.
    
    Supported formats:
        - CSV with header columns including x, y, and optionally flag
        - CSV with header columns x, y, and optionally classification (mapped to flag)
        - Headerless CSV with 2 columns (x, y) or 3 columns (x, y, flag)
    
    Input:
        file_path                    str            Path to a CSV file.
        verbose                      bool           If True, prints confirmation of loaded points.
    
    Output:
        p                            ndarray         Array with columns [x, y, flag] (float, float, int).
    """
    try:
        data = np.genfromtxt(file_path, delimiter=',', names=True, dtype=None, encoding='utf-8')
        if data.dtype.names and 'x' in data.dtype.names and 'y' in data.dtype.names:
            x = np.asarray(data['x'], dtype=float)
            y = np.asarray(data['y'], dtype=float)
            if 'flag' in data.dtype.names:
                flag = np.asarray(data['flag'], dtype=int)
            elif 'classification' in data.dtype.names:
                cls  = np.asarray(data['classification']).astype(str)
                cls  = np.char.lower(cls)
                flag = np.zeros(cls.shape, dtype=int)
                flag[np.isin(cls, ['boundary'])] = 1
                flag[np.isin(cls, ['hole', 'interior_boundary', 'inner_boundary'])] = 2
            else:
                flag = np.zeros(x.shape, dtype=int)
            if verbose:
                print(f"Loaded {x.shape[0]} points from {file_path}")
            return np.column_stack([x, y, flag])
    except Exception:
        pass

    p = np.genfromtxt(file_path, delimiter=',', skip_header=0)
    if p.ndim == 1:
        p = np.vstack([p])
    if p.shape[1] == 2:
        p = np.column_stack([p, np.zeros(p.shape[0])])
    
    if verbose:
        print(f"Loaded {p.shape[0]} points from {file_path}")
    
    return p
