"""
CloudGenerator — Point Cloud Generation System for mGFD

Overview:
    This package contains the core modular components for the Cloud Generation system.
    It is organized into specialized modules to separate concerns and improve maintainability.

Data conventions:
    None.

Public API:
    utils               Geometric utilities and helper functions.
    data_processing     Input data handling and validation.
    point_generation    Core point generation algorithms (Regular, Natural).
    classification      Node classification (boundary vs. interior).
    visualization       Point cloud rendering and plotting.
    export              Data output and verification.
    reduction           Post-generation point cloud decimation.

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
    May, 2024.
Last Modification:
    August, 2026.
"""
