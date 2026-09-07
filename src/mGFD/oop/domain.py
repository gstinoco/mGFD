"""
mGFD.oop.domain — Domain Abstraction for mGFD OOP Interface

Overview:
    This module provides the Domain class, combining a Cloud instance with its boundary conditions.

Public API:
    Domain      Combines point cloud topology and boundary condition specifications.

Credits:
    All the codes presented below were developed by:
        Dr. Gerardo Tinoco-Guerrero
        Dr. Francisco Javier Domínguez-Mota
        Dr. José Alberto Guzmán-Torres
        Universidad Michoacana de San Nicolás de Hidalgo
        gerardo.tinoco@umich.mx

Date:
    September, 2026.
"""

## Library importation.
from typing import Any                                                                                                                  # Typing support.
from mGFD.oop.boundary import BoundaryCondition, Dirichlet                                                                              # Boundary condition abstractions.

class Domain:                                                                                                                           # Domain class definition.
    """
    mGFD Physical Domain representation.
    
    Acts as the primary wrapper tying together the spatial discretization (Cloud) and the mathematical 
    boundary constraints (BoundaryCondition). The solver consumes a Domain object to understand 
    where points are located spatially and how the borders of that space are mathematically bounded.
    """
    def __init__(self, cloud: Any, boundary: Any = None) -> None:                                                                       # Constructor.
        """
        __init__
        Initialize Domain with a Cloud and a BoundaryCondition.
        
        Input:
            cloud       Cloud               Point cloud object.
            boundary    BoundaryCondition   Boundary condition object (defaults to Dirichlet(0.0)).
        """
        self.cloud = cloud                                                                                                              # Store point cloud instance.
        self.boundary: BoundaryCondition                                                                                                # Declare boundary condition type.
        if boundary is None:                                                                                                            # If boundary condition not provided.
            self.boundary = Dirichlet(0.0)                                                                                              # Default to Dirichlet(0.0).
        elif isinstance(boundary, BoundaryCondition):                                                                                   # If BoundaryCondition instance.
            self.boundary = boundary                                                                                                    # Store directly.
        else:                                                                                                                           # Otherwise wrap in Dirichlet.
            self.boundary = Dirichlet(boundary)                                                                                         # Wrap in Dirichlet.

        if hasattr(self.cloud, 'p') and self.cloud.p is not None:                                                                       # Attach cloud points to boundary wrapper for indexing.
            self.boundary.points = self.cloud.p                                                                                         # Store cloud points reference.

    @property                                                                                                                           # Execute statement.
    def points(self) -> Any:                                                                                                            # Points property.
        """Returns the underlying (m, 3) points array [x, y, flag]."""
        return self.cloud.p                                                                                                             # Return points array.

    def __repr__(self) -> str:                                                                                                          # String representation.
        return f"Domain(cloud={self.cloud}, boundary={self.boundary})"                                                                  # Return summary string.
