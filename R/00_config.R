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

    # --- Material-based alignment anchors ---
    # Only particles matching these patterns (case-insensitive regex) are used
    # for spatial alignment. This ensures the transform is driven by particles
    # that genuinely overlap between instruments.
    align_ftir_materials = c("PET", "Polypro"),
    align_raman_materials = c("Polyethylene terephtalate", "Polypropylene"),
    align_raman_min_size_um = 20,  # exclude Raman particles below FTIR detection limit

    # --- Landmark-first alignment (Tier 1) ---
    # Large particles and fibers are used as high-confidence anchors.
    # If enough landmarks match with tight residuals, the expensive full
    # RANSAC grid search (Tier 2) is skipped entirely.
    landmark_min_size_um       = 100,  # particles >= this size are landmarks
    landmark_fiber_aspect_ratio = 3.0, # major/minor >= this → fiber (landmark even if smaller)
    landmark_fiber_min_size_um  = 100, # minimum size for fiber landmarks
    landmark_min_count          = 4,   # need at least this many landmarks per dataset
    landmark_confidence_min_inlier_ratio = 0.5,  # >= 50% landmarks must be inliers
    landmark_confidence_max_residual_um  = 50,   # mean residual must be below this (µm)
    landmark_skip_full_ransac  = TRUE, # if TRUE, skip Tier 2 when landmarks are confident

    # --- RANSAC alignment ---
    ransac_coarse_step_deg = 1,     # rotation grid step for coarse search (degrees)
    ransac_n_iterations    = 2000,  # number of RANSAC iterations for refinement
    ransac_min_samples     = 3,     # minimum point pairs per RANSAC sample
    ransac_inlier_dist_um  = 200,   # distance threshold to count as inlier (µm)
    ransac_allow_mirror    = TRUE,  # also search over reflections

    # --- ICP refinement ---
    icp_max_iterations      = 100,  # maximum ICP iterations
    icp_convergence_thresh  = 0.01, # stop when RMS improvement < this (µm)
    icp_max_pair_dist_um    = 500,  # max distance for ICP correspondences (µm)
    icp_reciprocal          = TRUE, # only keep mutual nearest-neighbor pairs
    icp_trim_pct            = 0.10, # discard worst 10% of pairs each iteration

    # --- Particle matching ---
    match_dist_threshold_um = 100,  # max distance between matched particles (µm)
    match_size_weight       = 0,    # weight for size similarity (0 = spatial only)
    match_size_metric       = "feret_max_um",  # which size metric to compare

    # --- Material equivalence mapping (for agreement scoring) ---
    # User-defined mapping of FTIR material names → canonical types.
    # Will be populated with actual mapping next session.
    # Format: named list where names are canonical types, values are regex
    # patterns matching instrument-specific material names.
    material_map_ftir  = NULL,  # e.g., list(PET = "PET", PP = "Polypro|PP", ...)
    material_map_raman = NULL,  # e.g., list(PET = "Polyethylene terephtalate|PET", ...)

    # --- Diagnostics ---
    plot_width  = 10,   # plot width in inches
    plot_height = 8     # plot height in inches
  )
}
