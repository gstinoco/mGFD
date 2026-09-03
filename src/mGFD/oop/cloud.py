"""
mGFD.oop.cloud — Point Cloud Abstraction for mGFD OOP Interface

Overview:
    This module provides the Cloud class, an object-oriented representation of 2D meshless point clouds.
    Supports loading from CSV, creating from NumPy arrays, generating via Poisson/Regular algorithms,
    computing neighbor stencils, and binding boundary conditions.

Public API:
    Cloud       OOP Point Cloud class encapsulating coordinates, flags, and neighbors.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

    With the funding of:
        Secretary of Science, Humanities, Technology and Innovation, SECIHTI. México.
        Coordination of Scientific Research, CIC-UMSNH. México.
        Aula CIMNE-Morelia. México.
        SIIIA-MATH: Soluciones de Ingeniería. México.

Date:
    September, 2026.
"""

## Library importation.
import numpy as np                                                                                                      # Core numerical operations.
from typing import Union, Optional, Tuple, Dict, Any, Callable, TYPE_CHECKING                           # Type hinting interfaces.

if TYPE_CHECKING:                                                                                       # Static type checker block.
    from mGFD.oop.domain import Domain                                                                  # Forward declaration for Domain.

from mGFD.io.io import load_points                                                                                      # Point cloud loading utility.
from mGFD.spatial.neighbors import compute_neighbors, compute_upwind_neighbors                                          # Neighbor computation routines.
from mGFD.cloud_generator.core.generator import generate_cloud_natural, generate_cloud_regular                          # Point cloud generators.
from mGFD.cloud_generator.viz.visualization import create_visualization                                                 # Cloud visualization tool.
from mGFD.oop.boundary import BoundaryCondition, Dirichlet                                                              # Boundary condition abstractions.

class Cloud:                                                                                                            # OOP Point Cloud class.
    """
    mGFD Point Cloud abstraction.
    
    Encapsulates node coordinates (x, y), node flags (0=interior, 1/2=boundary), and spatial neighbor list (vec).
    Provides factory constructors (.from_csv, .from_array, .generate_natural, .generate_regular), neighbor calculation,
    and boundary binding (.set_boundary).
    """
    def __init__(self, points: np.ndarray, neighbors: Optional[np.ndarray] = None) -> None:                             # Constructor.
        """
        __init__
        Initialize Cloud with an (m, 3) NumPy array [x, y, flag] and an optional neighbor matrix.
        
        Input:
            points      m x 3 ndarray           Point cloud [x, y, flag].
            neighbors   m x nvec ndarray[int]   Optional precomputed neighbor list.
        """
        from mGFD.exceptions import CloudShapeError                                                                     # Import CloudShapeError.
        if not isinstance(points, np.ndarray) or points.ndim != 2:                                                      # Validate 2D numpy array.
            raise CloudShapeError("Point cloud 'p' must be a 2D numpy array with shape (m, 3)")                         # Raise CloudShapeError.
        self.p = np.asarray(points, dtype=np.float64)                                                                   # Store point coordinates array.
        self.neighbors = neighbors                                                                                      # Store neighbor matrix if provided.
        self._bc: Optional[BoundaryCondition] = None                                                                    # Boundary condition holder.

    @classmethod                                                                                                        # Execute statement.
    def from_csv(cls, filepath: str) -> "Cloud":                                                                        # Factory constructor from CSV.
        """
        from_csv
        Load point cloud from a CSV file.
        
        Input:
            filepath    str         Path to CSV file containing [x, y, flag] or [x, y, classification].
            
        Output:
            cloud       Cloud       Instantiated Cloud object.
        """
        p = load_points(filepath)                                                                                       # Load points using mGFD.io.
        return cls(p)                                                                                                   # Return instantiated Cloud object.

    @classmethod                                                                                                        # Execute statement.
    def from_array(cls, array: np.ndarray) -> "Cloud":                                                                  # Factory constructor from array.
        """
        from_array
        Instantiate Cloud from a NumPy array [x, y, flag].
        
        Input:
            array       ndarray     m x 3 array [x, y, flag].
            
        Output:
            cloud       Cloud       Instantiated Cloud object.
        """
        return cls(array)                                                                                               # Return instantiated Cloud object.

    @classmethod                                                                                                        # Execute statement.
    def generate_natural(cls, csv_contour: str, output_csv: str = "generated_cloud.csv",
                         cloud_size: Optional[float] = None, density_multiplier: float = 1.0,                           # Assign cloud_size: Optional[float].
                         inside_regions: bool = False, save: bool = False, show: bool = False) -> "Cloud":              # Factory constructor via Natural Poisson Disk.
        """
        generate_natural
        Generate an irregular point cloud using Poisson Disk Sampling from a boundary contour CSV.
        
        Input:
            csv_contour         str         Path to boundary contour CSV file.
            output_csv          str         Target CSV path to write generated cloud.
            cloud_size          float       Target point spacing.
            density_multiplier  float       Density scaling multiplier.
            inside_regions      bool        Multi-region hole processing flag.
            save                bool        Save PNG/SVG visualization to disk.
            show                bool        Show interactive Matplotlib plot window.
            
        Output:
            cloud               Cloud       Instantiated Cloud object.
        """
        generate_cloud_natural(csv_contour, output_csv, inside_regions=inside_regions,                                  # Assign generate_cloud_natural(csv_contour, output_csv, inside_regions.
                               cloud_size=cloud_size, density_multiplier=density_multiplier, save=save, show=show)      # Execute Poisson Disk generator.
        p = load_points(output_csv)                                                                                     # Load generated points array.
        return cls(p)                                                                                                   # Return instantiated Cloud object.

    @classmethod                                                                                                        # Execute statement.
    def generate_regular(cls, csv_contour: str, output_csv: str = "generated_cloud.csv",
                        cloud_size: Optional[float] = None, density_multiplier: float = 1.0,                            # Assign cloud_size: Optional[float].
                        inside_regions: bool = False, save: bool = False, show: bool = False) -> "Cloud":               # Factory constructor via Regular Grid.
        """
        generate_regular
        Generate a regular grid point cloud from a boundary contour CSV.
        
        Input:
            csv_contour         str         Path to boundary contour CSV file.
            output_csv          str         Target CSV path to write generated cloud.
            cloud_size          float       Target point spacing.
            density_multiplier  float       Density scaling multiplier.
            inside_regions      bool        Multi-region hole processing flag.
            save                bool        Save PNG/SVG visualization to disk.
            show                bool        Show interactive Matplotlib plot window.
            
        Output:
            cloud               Cloud       Instantiated Cloud object.
        """
        generate_cloud_regular(csv_contour, output_csv, inside_regions=inside_regions,                                  # Assign generate_cloud_regular(csv_contour, output_csv, inside_regions.
                               cloud_size=cloud_size, density_multiplier=density_multiplier, save=save, show=show)      # Execute Regular Grid generator.
        p = load_points(output_csv)                                                                                     # Load generated points array.
        return cls(p)                                                                                                   # Return instantiated Cloud object.

    @property                                                                                                           # Execute statement.
    def num_nodes(self) -> int:                                                                                         # Number of nodes property.
        """Returns total number of nodes in point cloud."""
        return int(self.p.shape[0])                                                                                     # Return node count.

    @property                                                                                                           # Execute statement.
    def interior_mask(self) -> np.ndarray:                                                                              # Interior mask property.
        """Returns boolean mask array for interior nodes (flag == 0)."""
        return self.p[:, 2] == 0                                                                                        # Return boolean interior array.

    @property                                                                                                           # Execute statement.
    def boundary_mask(self) -> np.ndarray:                                                                              # Boundary mask property.
        """Returns boolean mask array for boundary nodes (flag != 0)."""
        return self.p[:, 2] != 0                                                                                        # Return boolean boundary array.

    @property                                                                                                           # Execute statement.
    def points(self) -> np.ndarray:                                                                                     # Points property.
        """Returns raw (m, 3) ndarray [x, y, flag]."""
        return self.p                                                                                                   # Return points array.

    def compute_neighbors(self, nvec: int = 15, upwind: bool = False, a: float = 0.0, b: float = 0.0) -> np.ndarray:    # Compute neighbor list.
        """
        compute_neighbors
        Computes spatial neighbor list (vec) and caches it on the instance.
        
        Input:
            nvec        int     Number of neighbors per node.
            upwind      bool    Use upwind-biased neighbor selection.
            a, b        float   Advection velocity components (if upwind=True).
            
        Output:
            vec         ndarray Neighbor list array.
        """
        if upwind and (a != 0.0 or b != 0.0):                                                                           # Check if upwind search requested.
            self.neighbors = compute_upwind_neighbors(self.p, a, b, nvec)                                               # Compute upwind neighbors.
        else:                                                                                                           # Isotropic search branch.
            self.neighbors = compute_neighbors(self.p, nvec=nvec)                                                       # Compute KDTree neighbors.
        return self.neighbors                                                                                           # Return cached neighbor matrix.

    def set_boundary(self, bc_spec: Union[BoundaryCondition, float, Callable, np.ndarray]) -> "Domain":                 # Bind boundary condition.
        """
        set_boundary
        Binds a boundary condition specification to the cloud and returns a Domain object.
        
        Input:
            bc_spec     Union[BoundaryCondition, float, Callable, ndarray]  Boundary condition value/class.
            
        Output:
            domain      Domain                                              Constructed Domain object.
        """
        from mGFD.oop.domain import Domain                                                                              # Lazy import to avoid circular import.
        if isinstance(bc_spec, BoundaryCondition):                                                                      # If already BoundaryCondition object.
            self._bc = bc_spec                                                                                          # Assign directly.
        else:                                                                                                           # Otherwise wrap in Dirichlet.
            self._bc = Dirichlet(bc_spec)                                                                               # Wrap in Dirichlet boundary condition.
        return Domain(self, self._bc)                                                                                   # Return Domain object.

    def plot(self, save: bool = False, show: bool = True, filename: str = "cloud_plot") -> bool:                        # Point cloud plot method.
        """
        plot
        Render a 2D scatter plot of the point cloud.
        
        Input:
            save        bool    Save PNG/SVG file to disk.
            show        bool    Display interactive plot window.
            filename    str     Base filename for saving.
            
        Output:
            success     bool    True if plot succeeded.
        """
        regions = [1] * len(self.p)                                                                                     # Single region mapping list.
        classifications = ['boundary' if f != 0 else 'interior' for f in self.p[:, 2]]                                  # Map flags to classification strings.
        return create_visualization(self.p[:, :2], regions, filename, classifications, save=save, show=show)            # Call visualization generator.

    def __repr__(self) -> str:                                                                                          # String representation.
        return f"Cloud(nodes={self.num_nodes}, interior={np.sum(self.interior_mask)}, boundary={np.sum(self.boundary_mask)})"# Return formatted summary string.
