import numpy as np
from mGFD.cloud_generator.core.point_generation.boundary import generate_boundary_points

def test_boundary_generation():
    contour = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    cloud_size = 0.5
    
    p = generate_boundary_points(contour, cloud_size)
    
    assert isinstance(p, np.ndarray)
    assert p.shape[1] == 2
    # A square of 1x1 with cloud_size 0.5 has points at 0, 0.5, 1 on each side
    # Corners: 4, midpoints: 4, total = 8 boundary points.
    assert p.shape[0] >= 4
