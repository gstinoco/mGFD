"""
Exceptions — Custom exception hierarchy for the mGFD package

Overview:
    This module defines the custom exception hierarchy for the mGFD package.
    Using specific exceptions allows developers to catch and handle errors
    more granularly (e.g. malformed point clouds vs mismatched array sizes)
    rather than relying on generic ValueErrors or TypeErrors.

Public API:
    mGFDError                   Base exception class.
    CloudShapeError             Raised for invalid point cloud matrices.
    InputTypeError              Raised for unsupported input data types.
    DimensionMismatchError      Raised when matrix dimensions do not align.
    OperatorFormatError         Raised for invalid differential operator weights.
    ParameterError              Raised for strictly invalid numerical parameters.

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

class mGFDError(Exception):
    """
    Base exception class for all custom errors in the mGFD package.
    """
    pass                                                                                                                                # Base marker, no specific logic inside.

class CloudShapeError(mGFDError):
    """
    Raised when the point cloud 'p' does not conform to the expected format.
    Specifically, it must be a 2D numpy array with 3 columns: [x, y, flag].
    """
    pass                                                                                                                                # Marker for geometry-related errors.

class InputTypeError(mGFDError):
    """
    Raised when an input parameter (like a forcing term or boundary condition)
    is of an unsupported type. Supported types generally include Callables,
    NumPy ndarrays, floats, or ints.
    """
    pass                                                                                                                                # Marker for unsupported parameter type errors.

class DimensionMismatchError(mGFDError):
    """
    Raised when the dimensions of provided empirical arrays (e.g. data vectors)
    do not match the shape of the point cloud or the time discretization grid.
    """
    pass                                                                                                                                # Marker for algebraic dimension mismatches.

class OperatorFormatError(mGFDError):
    """
    Raised when the provided differential operator vector 'L' does not have
    the required number of coefficients (typically at least 5 for the 2D PDE).
    """
    pass                                                                                                                                # Marker for malformed PDE operator inputs.

class ParameterError(mGFDError):
    """
    Raised when a numerical parameter is strictly invalid 
    (e.g., negative number of time steps, invalid iteration count).
    """
    pass                                                                                                                                # Marker for invalid numerical bounds in arguments.
