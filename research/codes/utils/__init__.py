"""
Utils — Research Laboratory Utilities

Overview:
    This package aggregates and exposes common utility functions used across the mGFD
    research benchmarking laboratory, including neighbor caching, batch runner orchestration,
    area-weighted error metrics, and sweep visualization.

Public API:
    project_root                Resolve absolute root directory of the repository.
    neighbors_path              Generate standard filesystem paths for neighbor cache files.
    load_neighbors              Load cached neighbor lists from CSV files.
    save_neighbors              Save computed neighbor indices to CSV cache.
    import_module_from_file     Dynamically load and execute a Python runner module.
    run_job                     Worker execution wrapper for process pool parallelization.
    iter_clouds                 Generator yielding point cloud paths across datasets and scales.
    run_batch_suite             Execute a batch processing routine across cloud datasets.
    save_metrics                Record and serialize execution time and numerical error metrics.
    generate_sweep_summary      Aggregate individual batch JSON metric files into a master CSV.
    compute_rmse_stationary     Compute area-weighted RMSE for a stationary PDE solution snapshot.
    compute_rmse_transient      Compute area-weighted RMSE across time steps for transient solutions.
    Compute_Metrics_Stationary  Calculate comprehensive stationary numerical error metrics.
    Compute_Metrics_Transient   Calculate comprehensive transient numerical error metrics.
    find_latest_summary         Locate the most recently generated parameter sweep summary CSV.
    plot_sweep_results          Generate scaling and device comparison plots from sweep summary.
    extract_device              Extract target hardware backend device from configuration string.
    extract_config_label        Construct human-readable labels from configuration identifiers.

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
    May, 2024.
Last Modification:
    September, 2026.
"""

## Library importation.
from .batch_utils import (                                                                                                              # Import batch execution utilities.
    project_root,                                                                                                                       # Root locator.
    neighbors_path,                                                                                                                     # Neighbor path resolver.
    load_neighbors,                                                                                                                     # Neighbor cache loader.
    save_neighbors,                                                                                                                     # Neighbor cache saver.
    import_module_from_file,                                                                                                            # Dynamic module importer.
    run_job,                                                                                                                            # Multiprocessing job wrapper.
    iter_clouds,                                                                                                                        # Point cloud generator iterator.
    run_batch_suite,                                                                                                                    # Batch runner executor.
    save_metrics,                                                                                                                       # Metric writer.
    generate_sweep_summary,                                                                                                             # Sweep aggregator.
)                                                                                                                                       # End batch_utils imports.

from .metrics import (                                                                                                                  # Import numerical metric utilities.
    compute_rmse_stationary,                                                                                                            # Stationary RMSE calculator.
    compute_rmse_transient,                                                                                                             # Transient RMSE calculator.
    Compute_Metrics_Stationary,                                                                                                         # Stationary metrics dictionary.
    Compute_Metrics_Transient,                                                                                                          # Transient metrics dictionary.
)                                                                                                                                       # End metrics imports.

from .plot_sweep import (                                                                                                               # Import visualization utilities.
    find_latest_summary,                                                                                                                # Summary file finder.
    plot_sweep_results,                                                                                                                 # Sweep result visualizer.
    extract_device,                                                                                                                     # Device label extractor.
    extract_config_label,                                                                                                               # Configuration label builder.
)                                                                                                                                       # End plot_sweep imports.

__all__ = [                                                                                                                             # Public symbol table.
    "project_root",                                                                                                                     # Project root locator.
    "neighbors_path",                                                                                                                   # Neighbor path resolver.
    "load_neighbors",                                                                                                                   # Neighbor cache loader.
    "save_neighbors",                                                                                                                   # Neighbor cache saver.
    "import_module_from_file",                                                                                                          # Dynamic module importer.
    "run_job",                                                                                                                          # Multiprocessing job wrapper.
    "iter_clouds",                                                                                                                      # Point cloud generator iterator.
    "run_batch_suite",                                                                                                                  # Batch runner executor.
    "save_metrics",                                                                                                                     # Metric writer.
    "generate_sweep_summary",                                                                                                           # Sweep aggregator.
    "compute_rmse_stationary",                                                                                                          # Stationary RMSE calculator.
    "compute_rmse_transient",                                                                                                           # Transient RMSE calculator.
    "Compute_Metrics_Stationary",                                                                                                       # Stationary metrics dictionary.
    "Compute_Metrics_Transient",                                                                                                        # Transient metrics dictionary.
    "find_latest_summary",                                                                                                              # Summary file finder.
    "plot_sweep_results",                                                                                                               # Sweep result visualizer.
    "extract_device",                                                                                                                   # Device label extractor.
    "extract_config_label",                                                                                                             # Configuration label builder.
]                                                                                                                                       # End symbol table.
