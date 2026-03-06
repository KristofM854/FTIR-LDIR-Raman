# =============================================================================
# 00_config.R — Configurable parameters for multi-instrument particle matching
# =============================================================================

#' Create default configuration
#' 
#' # ---- Python configuration (must run before any reticulate use) ----

Sys.setenv(RETICULATE_USE_UV = "0")  # prevents uv auto-install attempts

if (requireNamespace("reticulate", quietly = TRUE)) {
  reticulate::use_python(
    "C:/Program Files/Python314/python.exe",
    required = TRUE
  )
}
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
    min_particle_size_um  = 0,     # minimum Feret size in um (0 = keep all)

    # --- Coordinate normalization ---
    normalize_scale = FALSE,  # scale both clouds to unit variance (usually not needed)

    # --- Material-based alignment anchors ---
    align_ftir_materials = c("PET", "Polypro"),
    align_raman_materials = c("Polyethylene terephtalate", "Polypropylene"),
    align_raman_min_size_um = 20,  # exclude Raman particles below FTIR detection limit

    # --- LDIR-specific settings ---
    align_ldir_materials = c("Polyethylene terephthalate", "Polypropylene",
                             "Polycarbonate"),
    ldir_quality_threshold = 0.6,  # Agilent quality score (0-1)
    ldir_scan_diameter_um  = 13000, # 13mm filter diameter
    ldir_flip_y_for_alignment = FALSE, # map_pixels_to_um_circle already inverts y;
                                       # Raman y is also upward — no second flip needed
    ldir_rotate_deg_for_alignment = -90, # clockwise rotation applied before alignment;
                                          # corrects instrument export convention vs Raman.
                                          # must be one of: 0, 90, -90, 180

    # Fixed µm-per-pixel scale for the Raman microscope image.
    # NULL = fall back to particle-extent method (may cause systematic viewer drift).
    # Set to the instrument-specific value, e.g. raman_um_per_px = 2.5
    raman_um_per_px = NULL,

    # --- Descriptor RANSAC (optional Tier 2 replacement) ---
    # Set TRUE to use descriptor-based RANSAC instead of the coarse-grid material
    # RANSAC for LDIR→Raman alignment.  Backward-compatible default: FALSE.
    ldir_use_descriptor_ransac = FALSE,

    # Transform guardrail thresholds (applied to every alignment path):
    icp_min_scale        = 0.5,   # scale below this → WARN (likely degenerate)
    icp_max_scale        = 2.0,   # scale above this → WARN (likely degenerate)
    icp_max_rotation_deg = 90,    # |rotation| above this → WARN (likely spurious)

    # Enable/disable scan-circle detection for LDIR image coordinate mapping.
    # TRUE  (default): detect the circular scan area and use it for µm calibration.
    # FALSE: skip circle detection entirely and map the full image to scan bounds.
    #        Use this when the LDIR image already covers the full scan field without
    #        a visible circular crop (e.g. tiled mosaic exports).
    ldir_use_circle_detection = TRUE,

    # Manual scan-circle override (pixels). NULL = auto-detect (recommended).
    # Set when detect_ldir_scan_circle() fails and prints the hard-stop message:
    #   config$ldir_circle_manual = list(cx = 1000, cy = 1000, r = 950)
    ldir_circle_manual = NULL,

    # LDIR export format: "auto" (detect circle vs mosaic), "circular" (require circle),
    # "mosaic" (force full-image mapping, skip circle detection).
    ldir_export_format = "auto",

    # Named explicit landmark correspondences: LDIR particle_id → Raman particle_id
    # Example: c("A3" = "A3", "MP_11" = "Raman_190")
    # NULL = no explicit map, fall back to size-based landmark RANSAC
    ldir_landmark_map = NULL,

    # When TRUE and ldir_landmark_map is set: use Procrustes as the FINAL
    # transform (do not let ICP override it). ICP still runs for diagnostics.
    ldir_procrustes_lock = TRUE,

    # Particle IDs to trace stage-by-stage in debug mode
    debug_trace_ids = c("A3", "MP_11"),

    # --- Debug mode ---
    debug = FALSE,  # set TRUE for debug artifacts

    # --- Landmark-first alignment (Tier 1) ---
    landmark_min_size_um       = 100,
    landmark_fiber_aspect_ratio = 3.0,
    landmark_fiber_min_size_um  = 100,
    landmark_min_count          = 4,
    landmark_confidence_min_inlier_ratio = 0.5,
    landmark_confidence_max_residual_um  = 50,
    landmark_skip_full_ransac  = TRUE,

    # --- RANSAC alignment ---
    ransac_coarse_step_deg = 1,
    ransac_n_iterations    = 2000,
    ransac_min_samples     = 3,
    ransac_inlier_dist_um  = 200,
    ransac_allow_mirror    = TRUE,

    # --- ICP refinement ---
    icp_max_iterations      = 100,
    icp_convergence_thresh  = 0.01,
    icp_max_pair_dist_um    = 500,
    icp_reciprocal          = TRUE,
    icp_trim_pct            = 0.10,
    icp_elongation_downweight = TRUE,
    icp_elongation_alpha    = 0.5,

    # --- Particle matching ---
    match_method               = "hungarian",  # "hungarian" or "greedy"
    match_dist_threshold_um    = 100,
    match_adaptive_dist_factor = 0.15,
    match_size_weight          = 0.2,    # (used by greedy only)
    match_size_metric          = "feret_max_um",

    # Hungarian cost function weights
    match_lambda_area          = 0.3,   # weight for |log(area ratio)| penalty
    match_lambda_feret         = 0.3,   # weight for |log(feret ratio)| penalty
    match_lambda_aspect        = 0.1,   # weight for aspect ratio difference

    # Ambiguity detection
    ambiguity_radius_um        = 50,    # radius for counting nearby candidates

    # --- Material equivalence mapping ---
    material_map_ftir  = NULL,
    material_map_raman = NULL,
    material_map_ldir  = NULL,

    # --- Diagnostics ---
    plot_width  = 10,
    plot_height = 8
  )
}
