"""
CloudGenerator.utils — Utilities

Overview:
    This module provides essential geometric utilities and helper functions used across
    the Cloud Generation system. It encapsulates core mathematical logic to ensure
    consistency and reusability.

Data conventions:
    None.

Public API:
    calculate_cloud_size                    Calculate adaptive cloud size based on region geometry characteristics.
    create_closed_contour                   Ensure a contour is properly closed for geometric operations.
    calculate_dynamic_boundary_refinement   Calculate dynamic boundary refinement based on cloud density.

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
    March, 2026.
Last Modification:
    August, 2026.
"""

## Library importation.
import numpy as np                                                                                                                      # Numerical arrays and mathematical operations.

from typing import List, Tuple, Optional                                                                                                # Import type hints.

CLOUD_FACTORS = {                                                                                                                       # Global configuration parameters.
    "adaptive_factor": 0.15,                                                                                                            # Factor to scale the cloud size.
    "default_cloud_size": 0.05,                                                                                                         # Default cloud size if calculation fails.
    "boundary_refinement": 0.02                                                                                                         # Default boundary refinement distance.
}

def calculate_cloud_size(region_points: List[Tuple[float, float]]) -> float:
    """
    Calculate adaptive cloud size based on region geometry characteristics.
    
    This function analyzes the spacing between consecutive boundary points to determine
    an optimal cloud size for point generation. The cloud size is calculated as a
    fraction of the average distance between boundary points, ensuring appropriate
    point density for the given geometry.
    
    Input:
        region_points       List[Tuple[float, float]]       List of (x, y) coordinate tuples defining the region boundary.
    
    Output:
        cloud_size          float                           Calculated cloud size (spacing between points) based on geometry analysis.
                                                            Returns default cloud size if calculation fails or insufficient points.
    """
    try:                                                                                                                                # Start execution block.
        if len(region_points) < 2:                                                                                                      # Check if enough points are provided.
            return CLOUD_FACTORS["default_cloud_size"]                                                                                  # Return default if not.
        
        # 1. Calculate distances between boundary points
        distances = []                                                                                                                  # List to hold distances.
        for i in range(len(region_points) - 1):                                                                                         # Iterate through adjacent points.
            x1, y1   = region_points[i]                                                                                                 # Extract first point coordinates.
            x2, y2   = region_points[i + 1]                                                                                             # Extract second point coordinates.
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                                                                             # Calculate Euclidean distance.
            distances.append(distance)                                                                                                  # Store the calculated distance.
        
        # 2. Compute cloud size
        if distances:                                                                                                                   # Ensure distances list is not empty.
            avg_distance = np.mean(distances)                                                                                           # Compute average distance.
            cloud_size   = avg_distance * CLOUD_FACTORS["adaptive_factor"]                                                              # Apply adaptive factor.
            return float(cloud_size)                                                                                                    # Return computed cloud size.
        else:                                                                                                                           # If no distances were calculated.
            return CLOUD_FACTORS["default_cloud_size"]                                                                                  # Return default cloud size.
            
    except Exception as e:                                                                                                              # If execution fails.
        return CLOUD_FACTORS["default_cloud_size"]                                                                                      # Return default cloud size.

def create_closed_contour(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Ensure a contour is properly closed for geometric operations.
    
    This function validates that a contour forms a closed polygon by checking if the
    first and last points are identical. If not, it adds the first point at the end
    to create a closed contour, which is required for proper polygon operations.
    
    Input:
        points              List[Tuple[float, float]]       List of (x, y) coordinate tuples defining the contour.
    
    Output:
        points              List[Tuple[float, float]]       Closed contour with first point repeated at the end if necessary.
    """
    if len(points) < 3:                                                                                                                 # A contour must have at least 3 points.
        raise ValueError("At least 3 points are needed to create a contour")                                                            # Raise exception otherwise.
    
    if points[0] != points[-1]:                                                                                                         # Check if first and last points do not match.
        return points + [points[0]]                                                                                                     # Append the first point to the end.
    else:                                                                                                                               # If already closed.
        return points                                                                                                                   # Return as is.

def calculate_dynamic_boundary_refinement(points: np.ndarray, cloud_size: Optional[float] = None) -> float:
    """
    Calculate dynamic boundary refinement based on cloud density.
    
    This function calculates an appropriate boundary tolerance based on:
    1. Average distance between nearest neighbors (point density)
    2. Cloud size parameter if provided
    3. Domain size for scaling
    
    Input:
        points              m x 2 ndarray   Array of point coordinates.
        cloud_size          float           Cloud size parameter used in generation.
    
    Output:
        dynamic_refinement  float           Dynamic boundary refinement value.
    """
    if len(points) < 2:                                                                                                                 # Ensure enough points are present.
        return CLOUD_FACTORS["boundary_refinement"]                                                                                     # Fallback to default.
    
    try:                                                                                                                                # Start execution block.
        
        # 1. Convert to numpy array for calculations
        points_array = np.array(points, dtype=float)                                                                                    # Ensure points are in a NumPy array.
        
        # 2. Calculate domain dimensions
        domain_width  = np.max(points_array[:, 0]) - np.min(points_array[:, 0])                                                         # Calculate maximum width in x.
        domain_height = np.max(points_array[:, 1]) - np.min(points_array[:, 1])                                                         # Calculate maximum height in y.
        domain_area   = domain_width * domain_height                                                                                    # Compute bounding box area.
        
        # 3. Method 1: Based on average nearest neighbor distance
        sample_size    = min(500, len(points_array))                                                                                    # Sample a subset of points for efficiency.
        sample_indices = np.random.choice(len(points_array), sample_size, replace=False)                                                # Randomly select indices.
        sample_points  = points_array[sample_indices]                                                                                   # Extract the sampled points.
        
        diff = sample_points[:, np.newaxis, :] - sample_points[np.newaxis, :, :]                                                        # Calculate pairwise coordinate differences.
        distances = np.sqrt(np.sum(diff**2, axis=-1))                                                                                   # Compute Euclidean distance matrix.
        
        np.fill_diagonal(distances, np.inf)                                                                                             # Exclude self-distances from minimum calculation.
        nearest_distances = np.min(distances, axis=1)                                                                                   # Find distance to closest neighbor for each point.
        
        avg_nearest_distance = np.mean(nearest_distances)                                                                               # Calculate average nearest neighbor distance.
        
        # 4. Method 2: Based on point density
        point_density            = len(points_array) / domain_area if domain_area > 0 else 1                                            # Compute points per unit area.
        density_based_refinement = 1.0 / np.sqrt(point_density) if point_density > 0 else 0.01                                          # Estimate inter-point distance based on density.
        
        # 5. Method 3: Based on cloud_size if provided
        cloud_size_based_refinement = cloud_size * 0.5 if cloud_size else None                                                          # Estimate boundary threshold based on cloud size.
        
        # 6. Combine methods with weights
        refinement_candidates = []                                                                                                      # List to hold the refinement estimates.
        
        refinement_candidates.append(avg_nearest_distance * 0.3)                                                                        # Primary: Average nearest neighbor distance.
        refinement_candidates.append(density_based_refinement * 0.1)                                                                    # Secondary: Density-based calculation.
        
        if cloud_size_based_refinement:                                                                                                 # If cloud size estimate exists.
            refinement_candidates.append(float(cloud_size) * 0.05 if cloud_size else 0.0)                                               # Tertiary: Cloud size based (if available).
        
        # 7. Use a weighted approach or median
        if cloud_size:                                                                                                                  # Favor cloud size if provided.
            dynamic_refinement = cloud_size * 0.05                                                                                      # Set refinement relative to cloud size.
        else:                                                                                                                           # Otherwise use statistical approaches.
            dynamic_refinement = np.median(refinement_candidates)                                                                       # Fallback to median of candidates.
        
        # 8. Apply reasonable bounds
        min_refinement = 0.0001                                                                                                         # Minimum acceptable threshold.
        if cloud_size:                                                                                                                  # Adjust minimum if cloud size is given.
            min_refinement = max(min_refinement, cloud_size * 0.005)                                                                    # Do not drop below a fraction of cloud size.
            
        max_refinement = min(domain_width, domain_height) * 0.02                                                                        # Maximum acceptable threshold relative to domain.
        
        if cloud_size:                                                                                                                  # Adjust maximum if cloud size is given.
            max_refinement = max(max_refinement, cloud_size * 1.0)                                                                      # Cap max relative to cloud size.
        
        dynamic_refinement = np.clip(dynamic_refinement, min_refinement, max_refinement)                                                # Clamp the value within acceptable bounds.

        return dynamic_refinement                                                                                                       # Return the computed refinement.
        
    except Exception as e:                                                                                                              # Fallback for errors during calculation.
        return CLOUD_FACTORS["boundary_refinement"]                                                                                     # Return default tolerance.
