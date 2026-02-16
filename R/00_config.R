# =============================================================================
# 00_config.R — Configurable parameters for FTIR–Raman particle matching
# =============================================================================

#' Create default configuration
#'
#' Returns a list of all pipeline parameters. Modify this function or override
#' individual values in main.R to adjust the pipeline behavior.
#'
#' @param ftir_path Path to FTIR Excel file
#' @param raman_path Path to Raman Excel file
#' @param output_dir Directory for output files
#' @return Named list of configuration parameters
make_config <- function(ftir_path  = NULL,
                        raman_path = NULL,
                        output_dir = "output") {

  list(
    # --- File paths ---
    ftir_path   = ftir_path,
    raman_path  = raman_path,
    output_dir  = output_dir,

    # --- Sheet names ---
    ftir_sheet  = "Long_Table",
    raman_sheet = "Long_Table",

    # --- Pre-filtering ---
    raman_hqi_threshold   = 70,    # minimum HQI score for Raman particles
    ftir_quality_threshold = 0,    # minimum AAU score for FTIR (0 = keep all)
    min_particle_size_um  = 0,     # minimum Feret size in µm (0 = keep all)

    # --- Coordinate normalization ---
    normalize_scale = FALSE,  # scale both clouds to unit variance (usually not needed)

    # --- RANSAC alignment ---
    ransac_coarse_step_deg = 2,     # rotation grid step for coarse search (degrees)
    ransac_n_iterations    = 1000,  # number of RANSAC iterations for refinement
    ransac_min_samples     = 3,     # minimum point pairs per RANSAC sample
    ransac_inlier_dist_um  = 100,   # distance threshold to count as inlier (µm)
    ransac_allow_mirror    = TRUE,  # also search over reflections

    # --- ICP refinement ---
    icp_max_iterations      = 50,   # maximum ICP iterations
    icp_convergence_thresh  = 0.01, # stop when RMS improvement < this (µm)
    icp_max_pair_dist_um    = 200,  # max distance for ICP correspondences (µm)

    # --- Particle matching ---
    match_dist_threshold_um = 50,   # max distance between matched particles (µm)
    match_size_weight       = 0,    # weight for size similarity (0 = spatial only)
    match_size_metric       = "feret_max_um",  # which size metric to compare

    # --- Diagnostics ---
    plot_width  = 10,   # plot width in inches
    plot_height = 8     # plot height in inches
  )
}
