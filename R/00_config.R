# =============================================================================
# 00_config.R — Configurable parameters for multi-instrument particle matching
# =============================================================================

#' Create default configuration
#' 
#' # ---- Python configuration (must run before any reticulate use) ----
#
# The LDIR image backend is optional: if Python or its packages are missing,
# the pipeline falls back to an R-based segmentation (see README). Python is
# therefore never hard-required here — a hardcoded, mandatory interpreter path
# would break the pipeline on any machine that does not have that exact path.
#
# BUT we must still PIN a real interpreter when one exists. If we leave the
# choice to reticulate, it auto-provisions an ephemeral 'r-reticulate'
# virtualenv — which on locked-down/offline machines tries to download `uv`
# (SSL error) and points at a python.exe that does not exist, so a perfectly
# good Python install is wrongly reported as "not installed" and the LDIR
# detector silently falls back to (weaker) R extraction.
#
# Interpreter priority:
#   1. RETICULATE_PYTHON environment variable — explicit override. Set this
#      (e.g. in ~/.Renviron) if your interpreter is not auto-detected:
#        RETICULATE_PYTHON=/usr/bin/python3          (Linux/macOS)
#        RETICULATE_PYTHON=C:/Program Files/Python314/python.exe   (Windows)
#   2. A known interpreter that actually exists on THIS machine (list below).
#   3. python3 / python on PATH.
# Each candidate is used only when the file exists, so listing a Windows path
# is harmless on Linux/macOS. Nothing is forced with required = TRUE.

Sys.setenv(RETICULATE_USE_UV = "0")  # prevents uv auto-install attempts

if (requireNamespace("reticulate", quietly = TRUE)) {
  .py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (!nzchar(.py)) {
    .cands <- c(
      "C:/Program Files/Python314/python.exe",
      "C:/Program Files/Python313/python.exe",
      "C:/Program Files/Python312/python.exe",
      "C:/Program Files/Python311/python.exe",
      unname(Sys.which("python3")),
      unname(Sys.which("python"))
    )
    .cands <- .cands[nzchar(.cands) & file.exists(.cands)]
    if (length(.cands) > 0) .py <- .cands[[1]]
    rm(.cands)
  }
  if (nzchar(.py)) {
    # Pin it so reticulate uses THIS interpreter and does not build a venv /
    # download uv. required = FALSE: a bad value degrades to the R fallback
    # instead of aborting the pipeline at source() time.
    Sys.setenv(RETICULATE_PYTHON = .py)
    try(reticulate::use_python(.py, required = FALSE), silent = TRUE)
  }
  rm(.py)
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
    ldir_quality_threshold = 0,    # Agilent quality score (0-1). 0 = keep ALL
                                   # ingested particles — no quality pre-filter.
                                   # Filter interactively in the viewer instead;
                                   # alignment/agreement apply their own
                                   # anchor-quality criteria independently.

    # ldir_hqi_unknown_threshold: particles whose Agilent quality (HQI) is below
    #   this value are RELABELLED as "unknown" (never dropped), mirroring the
    #   LDIR software's own display rule.  The raw identification is preserved in
    #   the identification_raw column.  The relabelled material flows into the
    #   viewer AND the cross-instrument agreement analysis.  Set to NULL/0 to
    #   disable relabelling and keep every raw identification.
    ldir_hqi_unknown_threshold = 0.85,
    ldir_scan_diameter_um  = 13000, # 13mm filter diameter
    ldir_flip_y_for_alignment = FALSE, # map_pixels_to_um_circle already inverts y;
                                       # Raman y is also upward — no second flip needed
    ldir_rotate_deg_for_alignment = -90, # clockwise rotation applied before alignment;
                                          # corrects instrument export convention vs Raman.
                                          # must be one of: 0, 90, -90, 180
    ldir_image_width_um = NULL,    # physical width (µm) of the LDIR image export.
                                   # NULL = assume the image spans the full scan
                                   # circle (ldir_scan_diameter_um). SET THIS when
                                   # the export covers only the deposit region —
                                   # otherwise every LDIR coordinate is inflated.
                                   # tools/diagnose_ldir_rotation.R measures it:
                                   # width = best_scale x ldir_scan_diameter_um.
                                   # PER-DATASET value, like raman_image_*.

    # --- Raman microscope image placement (WITec metadata) ---
    # Enter the values EXACTLY as WITec's Particle Scout shows them —
    # Width, Height, Center X, Center Y (all µm), including a negative
    # Center Y if that is what the panel reports.  WITec gives the center in
    # its video/image frame (Y pointing DOWN) while the particle export
    # ("Visual Center Point Y") is stage-frame (Y UP); the viewer resolves
    # the Y-axis convention automatically by checking which interpretation
    # contains the run's particles (raman_image_extent_from_config()).
    # When all four are set, the image is placed at exact physical bounds
    # regardless of the uploaded image's pixel resolution (resize-invariant).
    # NULL = fall back to particle-extent method.

    # NOTE: these values are PER-DATASET — read them from WITec's Particle
    # Scout for each new Raman scan and update them here, or the viewer
    # falls back to heuristic (particle-bbox) placement for that run.
    raman_image_width_um    = 12471.2,
    raman_image_height_um   = 12313.0,
    raman_image_center_x_um = 1487.6,
    raman_image_center_y_um = -4394.0,
    # raman_image_width_um    = NULL,
    # raman_image_height_um   = NULL,
    # raman_image_center_x_um = NULL,
    # raman_image_center_y_um = NULL,
    # Fixed µm-per-pixel scale for the Raman microscope image.
    # Used only as Priority 2 fallback when the four fields above are NULL.
    # NULL = fall back to TIFF DPI auto-detection or particle-extent method.
    # raman_um_per_px = 0.38,
    raman_um_per_px = NULL,
    # --- Descriptor RANSAC (optional Tier 2 replacement) ---
    # Set TRUE to use descriptor-based RANSAC instead of the coarse-grid material
    # RANSAC for LDIR→Raman alignment.  Backward-compatible default: FALSE.
    ldir_use_descriptor_ransac = FALSE,

    # Robust global registration for LDIR->Raman alignment (rotation x scale
    # sweep with translation voting + one-to-one inliers). Runs alongside the
    # coarse RANSAC and the transform with more inliers wins — fixes the
    # sparse-anchor case where RANSAC locks onto a poor local optimum
    # (observed: 5 inliers where 23 are achievable). TRUE = enabled.
    ldir_use_global_register = TRUE,

    # Allow a reflection (mirror) in the LDIR->Raman transform. FALSE by
    # default: both instruments image the same filter from the same side, so
    # there is no physical mirror between them. Allowing one lets the aligner
    # pick a spurious reflected optimum on sparse data (matching particles to
    # the wrong neighbours). Only set TRUE if a dataset genuinely needs it.
    ldir_allow_reflection = FALSE,

    # TPS local refinement of LDIR->Raman alignment: after the global
    # similarity + match, warp LDIR coordinates with a regularized thin-plate
    # spline fit to the residual displacement at the confident matches, then
    # re-match. Recovers peripheral particles left unmatched by non-rigid
    # distortion the global transform can't capture; kept only if it increases
    # matches. No-op under ldir_force_complete_match.
    ldir_tps_refine       = TRUE,
    ldir_tps_lambda       = 0.5,   # spline smoothing (normalized units);
                                   # larger = smoother/more conservative
    ldir_tps_min_controls = 6,     # min confident matches to attempt TPS

    # Transform guardrail thresholds (applied to every alignment path):
    icp_min_scale        = 0.5,   # scale below this → WARN (likely degenerate)
    icp_max_scale        = 2.0,   # scale above this → WARN (likely degenerate)
    icp_max_rotation_deg = 90,    # |rotation| above this → WARN (likely spurious)

    # Enable/disable scan-circle detection for LDIR image coordinate mapping.
    # TRUE  (default): detect the circular scan area and use it for µm calibration.
    # FALSE: skip circle detection entirely and map the full image to scan bounds.
    #        Use this when the LDIR image already covers the full scan field without
    #        a visible circular crop (e.g. tiled mosaic exports).
    ldir_use_circle_detection = FALSE,

    # Manual scan-circle override (pixels). NULL = auto-detect (recommended).
    # Set when detect_ldir_scan_circle() fails and prints the hard-stop message:
    #   config$ldir_circle_manual = list(cx = 1000, cy = 1000, r = 950)
    ldir_circle_manual = NULL,

    # LDIR export format: "auto" (detect circle vs mosaic), "circular" (require circle),
    # "mosaic" (force full-image mapping, skip circle detection).
    ldir_export_format = "auto",

    # LDIR particle detection tuning:
    # overshoot_factor: detect this many times more particles than expected from Excel.
    #   Values > 1 (e.g. 1.3) deliberately over-extract so that the Hungarian matcher
    #   has more candidates to choose from, improving match rate.
    ldir_overshoot_factor = 1.3,

    # merge_fibers: when TRUE, merge nearby elongated (high aspect-ratio) particles
    #   that are likely fragmented fiber detections back into single particles.
    ldir_merge_fibers = TRUE,

    # closing_radius: morphological closing radius (pixels) applied to the binary
    #   mask before connected-component labeling. Bridges small gaps in fiber
    #   detections. Set to 0 to disable.
    ldir_closing_radius = 2L,

    # match_threshold: maximum Hungarian assignment cost to accept a coordinate join.
    #   Higher values allow more lenient size-mismatch tolerance.
    #   Only used when ldir_force_coord_match = FALSE (the default is TRUE).
    ldir_match_threshold = 2.0,

    # --- Coordinate-join cost weights (join_ldir_coords) ---
    # The base cost is: log_area + 0.5*log_feret + w_ar*ar_term + w_rank*rank_term
    #
    # ldir_join_weight_ar: weight for the normalized aspect-ratio difference term.
    #   0 = disable. Higher values penalise shape-order swaps more strongly.
    ldir_join_weight_ar = 0.3,

    # ldir_join_weight_rank: weight for the normalized rank-consistency penalty.
    #   LDIR numbers particles largest-first; a large rank discrepancy between
    #   an Excel row and an image blob indicates a probable size-order swap.
    #   Set high (2.0+) so rank dominates over size when the two disagree —
    #   size terms are tanh-compressed to [0,1] so rank can always compete.
    #   0 = disable.
    ldir_join_weight_rank = 2.0,

    # ldir_join_weight_shape: weight for the shape-fingerprint term. When the
    #   image blobs carry rotation/scale-invariant shape descriptors
    #   (eccentricity, circularity, solidity) — as the processed-image extractor
    #   produces — each is rank-matched against the same Excel column and the
    #   normalized rank differences are averaged.  Because the processed overlay
    #   is the machine's own segmentation, this fingerprint disambiguates
    #   particles of near-identical size that the size/rank terms alone swap.
    #   The term is inert when either side lacks the descriptors (e.g. the
    #   optical-image path), so it only affects processed-image runs. 0 = disable.
    ldir_join_weight_shape = 1.0,

    # ldir_join_blob_keep_factor: before matching, sort image blobs by area
    #   (descending) and keep only the top ceiling(n_excel * factor) blobs.
    #   Removes spurious small fragments that shift the size rank ordering and
    #   cause large particles to be mismatched. Inf = keep all blobs.
    ldir_join_blob_keep_factor = 1.1,

    # ldir_join_confidence_threshold: base cost below which a pair is eligible
    #   to be locked in the confidence-first pass (Pass 1). Set 0 to skip Pass 1
    #   and fall back to a single global Hungarian solve.
    ldir_join_confidence_threshold = 0.5,

    # ldir_join_confidence_margin: uniqueness ratio for Pass 1 locking.
    #   The second-best competing Excel row for a candidate blob must cost at
    #   least this many times the best cost before the match is locked.
    #   Higher values = stricter uniqueness requirement.
    ldir_join_confidence_margin = 2.0,

    # ldir_join_rank_first_frac: fraction of Excel particles (from the top, i.e.
    #   the largest) that are matched in a rank-only mini-Hungarian (Pass 0)
    #   before the full cost is applied.  For the largest particles the
    #   LDIR-guaranteed descending ordering is more reliable than the image-
    #   derived sizes, so anchoring these pairs first prevents cascade failures.
    #   Set 0 to disable Pass 0.
    ldir_join_rank_first_frac = 0.25,

    # ldir_join_merge_dist_um: if non-NULL and > 0, image blobs whose centroids
    #   are within this distance (in µm) are fused into a single pseudo-particle
    #   before matching (Union-Find clustering).  Useful when the image
    #   segmenter splits a large particle into several fragments.
    #   NULL = disabled.
    ldir_join_merge_dist_um = NULL,

    # --- LDIR processed (analyzed) particle-overlay image ---
    # The optical (white-light) LDIR image over-sizes large particle blobs by
    # ~3x in linear dimension (optical halo/scattering), which breaks size-based
    # matching for large particles.  The LDIR software can additionally export an
    # "analyzed" overlay image where the SAME particles are drawn as solid
    # coloured blobs (green, blue, …) on a pure-black background; those blobs are
    # the machine's own segmentation and their sizes match the Excel data.  When
    # such an image is present alongside the optical image, the pipeline uses it
    # (via extract_ldir_processed_image_coords) as the image_df fed to
    # join_ldir_coords, so the matcher receives correctly-scaled blobs.

    # ldir_processed_image_suffix: filename suffix to search for when
    #   looking for the LDIR software's analyzed particle overlay image.
    #   Set to NULL to disable processed-image extraction entirely.
    ldir_processed_image_suffix = "_analyzed",

    # ldir_processed_image_min_brightness: per-pixel RGB sum threshold
    #   below which a pixel is treated as background (pure black = 0).
    ldir_processed_image_min_brightness = 30L,

    # ldir_min_blob_area_px: minimum connected-component area (pixels) for a
    #   processed-image blob to be kept.  Smaller components are discarded as
    #   noise / rendering artifacts.
    ldir_min_blob_area_px = 5L,

    # ldir_image_scale_um_per_px: µm-per-pixel scale applied to processed-image
    #   pixel measurements.  The processed overlay shares the optical image's
    #   pixel dimensions and scale, so this is the same µm/px as the optical
    #   calibration.  NULL = 1.0 (measurements stay in pixel units; size-based
    #   matching still works because only relative sizes/ranks drive the join).
    ldir_image_scale_um_per_px = NULL,

    # When TRUE, every Excel particle is assigned an image coordinate regardless
    # of size-match cost — no joins are rejected. Use the coord_match_cost slider
    # in the Shiny viewer to post-hoc filter poor-quality coordinate assignments.
    ldir_force_coord_match = TRUE,

    # When TRUE, Hungarian matching forces a 1-to-1 assignment for every LDIR
    # particle regardless of spatial distance — no pairs are rejected.
    #
    # DO NOT ENABLE. Forcing every particle removes the spatial gate, so the
    # Hungarian minimises the GLOBAL SUM of pairing costs across all forced
    # pairs. That lets a partner-less particle "rob" a good particle's correct
    # match through a multi-particle cascade: e.g. a particle sitting 22 µm from
    # its true partner was reassigned to one 2049 µm away so a partner-less
    # neighbour could take the 22 µm one (its own best option was ~1092 µm).
    # Post-hoc filtering (the viewer match-gate slider) CANNOT repair this — the
    # corrupted assignment has already given the partner away.
    #
    # FALSE (correct): the LDIR distance gate (match_dist_threshold_ldir_um,
    # with adaptive expansion for large particles) is a hard match-time
    # constraint. Only within-gate pairs are formed; particles with no plausible
    # partner are reported unmatched instead of stealing someone else's.
    ldir_force_complete_match = FALSE,

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

    # --- Determinism ---
    # Seed for the stochastic samplers in RANSAC / global registration. Fixing
    # it makes alignment (and therefore match counts) reproducible run-to-run.
    # The seeding is RNG-safe: each aligner saves and restores the global RNG
    # state, so this does not perturb randomness elsewhere in the pipeline.
    align_seed = 1L,

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
    match_dist_threshold_um    = 100,   # FTIR<->Raman gate (fine, dense)
    # Dedicated, looser gate for any pairing that involves LDIR. LDIR
    # centroids are coarser (circle calibration, large particles) so genuine
    # matches sit at larger residuals — on real runs the 90th percentile of
    # true LDIR<->Raman match distance is ~220 um, which the 100 um gate cuts.
    # Measure headroom for a given run with tools/diagnose_matching.R.
    # NULL = fall back to match_dist_threshold_um.
    match_dist_threshold_ldir_um = 250,
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
