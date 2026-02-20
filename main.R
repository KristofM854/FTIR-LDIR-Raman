# =============================================================================
# main.R — Multi-instrument particle matching pipeline
# =============================================================================
#
# Orchestrates the full pipeline for aligning and matching particles
# detected by FTIR, Raman, and LDIR microspectroscopy on the same filter.
#
# Input modes:
#   Mode 1 (Hardcoded): Set ftir_file / raman_file / ldir_file before sourcing
#   Mode 2 (Interactive): Set input_mode <- "interactive" to get file pickers
#
# Usage:
#   source("main.R")
#
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Setup: load packages and source modules
# ---------------------------------------------------------------------------

required_packages <- c("readxl", "ggplot2", "RANN")

missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace,
                                          logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

library(ggplot2)

# Source all modules (relative to project root)
source("R/utils.R")
source("R/00_config.R")
source("R/00b_file_input.R")
source("R/01_ingest.R")
source("R/01b_ingest_image.R")
source("R/utils_python.R")
source("R/01c_ingest_ldir.R")
source("R/02_prefilter.R")
source("R/03_normalize.R")
source("R/03b_landmark_align.R")
source("R/04_ransac.R")
source("R/05_transform.R")
source("R/06_icp_refine.R")
source("R/07_match.R")
source("R/08_agreement.R")
source("R/08b_material_map.R")
source("R/09_diagnostics.R")
source("R/10_export.R")

# ---------------------------------------------------------------------------
# 1. File input and configuration
# ---------------------------------------------------------------------------

# --- Input mode selection ---
# Set input_mode <- "interactive" before sourcing for file picker dialogs.
# Default: "hardcoded" — uses ftir_file / raman_file / ldir_file variables.
if (!exists("input_mode")) input_mode <- "hardcoded"

if (input_mode == "interactive") {
  # ------ Mode 2: Interactive file picker ------
  file_manifest <- collect_files_interactive()
  grouped       <- group_files_by_instrument(file_manifest)

  ftir_file      <- grouped$FTIR$tabular
  ftir_image     <- grouped$FTIR$image
  raman_file     <- grouped$Raman$tabular
  ldir_file      <- grouped$LDIR$tabular
  ldir_image     <- grouped$LDIR$image

} else {
  # ------ Mode 1: Hardcoded paths ------
  # Defaults for test data (set these before sourcing, or leave for defaults)
  if (!exists("ftir_file"))  ftir_file  <- NULL
  if (!exists("raman_file")) raman_file <- NULL
  if (!exists("ldir_file"))  ldir_file  <- NULL
  if (!exists("ftir_image")) ftir_image <- NULL
  if (!exists("ldir_image")) ldir_image <- NULL

  # Prompt for mandatory files if not set
  if (is.null(ftir_file)) {
    if (interactive()) {
      message("Please select the FTIR data file (.csv or .xlsx) ...")
      ftir_file <- file.choose()
    } else {
      stop("ftir_file must be set before running in non-interactive mode")
    }
  }
  if (is.null(raman_file)) {
    if (interactive()) {
      message("Please select the Raman data file (.csv or .xlsx) ...")
      raman_file <- file.choose()
    } else {
      stop("raman_file must be set before running in non-interactive mode")
    }
  }
  # ldir_file is optional — NULL means skip LDIR
}

config <- make_config(
  ftir_path  = ftir_file,
  raman_path = raman_file,
  output_dir = "output"
)

# Store LDIR path in config
config$ldir_path  <- if (exists("ldir_file"))  ldir_file  else NULL
config$ftir_image <- if (exists("ftir_image")) ftir_image else NULL
config$ldir_image <- if (exists("ldir_image")) ldir_image else NULL

# Create a timestamped run subfolder (output/YYYY-MM-DD_1, _2, ...)
config$output_dir <- make_run_dir(config$output_dir)

# Override any defaults as needed:
# config$raman_hqi_threshold     <- 75
# config$match_dist_threshold_um <- 25
# config$ransac_allow_mirror     <- TRUE

# ---------------------------------------------------------------------------
# 2. Data ingestion
# ---------------------------------------------------------------------------

log_message(strrep("=", 60))
log_message("Multi-Instrument Particle Matching Pipeline")
log_message(strrep("=", 60))

# --- FTIR ---
ftir_raw  <- ingest_ftir(config$ftir_path, sheet = config$ftir_sheet)

# --- Raman ---
raman_raw <- ingest_raman(config$raman_path, sheet = config$raman_sheet)

# --- LDIR (optional) ---
ldir_raw <- NULL
has_ldir <- !is.null(config$ldir_path) && nzchar(config$ldir_path)
if (has_ldir) {
  ldir_raw <- ingest_ldir(config$ldir_path)
  log_message("LDIR data loaded: ", nrow(ldir_raw), " particles (no coordinates)")
} else {
  log_message("LDIR: not provided — skipping LDIR analysis")
}

# --- FTIR image-based coordinate extraction (if image provided) ---
ftir_scan_bounds <- NULL
if (!is.null(config$ftir_image) && nzchar(config$ftir_image)) {
  log_message(strrep("-", 50))
  log_message("FTIR image-based particle extraction")

  # Compute scan bounds from image dimensions and 25 µm grid step.
  # The PerkinElmer Spotlight renders ~6 image pixels per grid cell.
  ftir_img_raw <- png::readPNG(config$ftir_image)
  ftir_grid_nx <- round((ncol(ftir_img_raw) + 1) / 6)
  ftir_grid_ny <- round((nrow(ftir_img_raw) + 1) / 6)
  ftir_scan_bounds <- list(
    x_min = 0,
    x_max = ftir_grid_nx * 25,
    y_min = 0,
    y_max = ftir_grid_ny * 25
  )
  rm(ftir_img_raw)
  log_message("  Scan bounds: [0, ", ftir_scan_bounds$x_max, "] x [0, ",
              ftir_scan_bounds$y_max, "] µm")

  ftir_from_image <- extract_particles_from_image(
    config$ftir_image,
    scan_bounds     = ftir_scan_bounds,
    adaptive_radius = 75L,
    adaptive_offset = 0.02,
    min_pixels      = 4,
    min_size_um     = 20,
    expected_count  = nrow(ftir_raw),
    instrument      = "FTIR"
  )

  log_message("Image extraction: ", nrow(ftir_from_image),
              " particles from FTIR image")
}

# ---------------------------------------------------------------------------
# 3. Pre-filtering
#
# Strategy: use only PLASTIC particles for spatial alignment (they are the
# ones that genuinely overlap between instruments). Then apply quality
# filters for the final matching & agreement steps.
# ---------------------------------------------------------------------------

# Minimal filtering (just remove invalid coords)
ftir_clean  <- prefilter_ftir(ftir_raw, min_quality = 0, min_size_um = 0)
raman_clean <- prefilter_raman(raman_raw, min_hqi = 0, min_size_um = 0)

# Extract PLASTIC particles from FTIR for alignment anchoring
ftir_plastic_mask <- grepl(
  paste(config$align_ftir_materials, collapse = "|"),
  ftir_clean$material, ignore.case = TRUE
)
ftir_for_align <- ftir_clean[ftir_plastic_mask, ]
log_message("FTIR plastic anchor particles: ", nrow(ftir_for_align),
            " (materials: ", paste(unique(ftir_for_align$material), collapse = ", "), ")")

# Filter Raman alignment targets by material (if configured)
if (!is.null(config$align_raman_materials)) {
  raman_plastic_mask <- grepl(
    paste(config$align_raman_materials, collapse = "|"),
    raman_clean$material, ignore.case = TRUE
  )
  raman_for_align <- raman_clean[raman_plastic_mask, ]
} else {
  raman_for_align <- raman_clean
}

# Remove Raman particles below FTIR detection limit
min_size <- config$align_raman_min_size_um
if (!is.null(min_size) && min_size > 0 && any(!is.na(raman_for_align$feret_max_um))) {
  size_mask <- is.na(raman_for_align$feret_max_um) | raman_for_align$feret_max_um >= min_size
  n_before <- nrow(raman_for_align)
  raman_for_align <- raman_for_align[size_mask, ]
  log_message("Raman alignment: removed ", n_before - nrow(raman_for_align),
              " particles < ", min_size, " um")
}
# Apply HQI threshold to material alignment anchors (only reliably
# identified particles should drive material-based alignment)
if (!is.null(config$raman_hqi_threshold) && config$raman_hqi_threshold > 0 &&
    any(!is.na(raman_for_align$quality))) {
  hqi_mask <- !is.na(raman_for_align$quality) &
              raman_for_align$quality >= config$raman_hqi_threshold
  n_before_hqi <- nrow(raman_for_align)
  raman_for_align <- raman_for_align[hqi_mask, ]
  log_message("Raman alignment: HQI filter (>= ", config$raman_hqi_threshold,
              "): kept ", nrow(raman_for_align), " of ", n_before_hqi)
}
log_message("Raman alignment target particles: ", nrow(raman_for_align),
            " (material + >= ", min_size, " um + HQI >= ",
            config$raman_hqi_threshold, ")")

# Quality-filtered sets for matching (applied after alignment)
# FTIR: apply quality filter as usual
ftir_for_match <- prefilter_ftir(
  ftir_raw,
  min_quality = config$ftir_quality_threshold,
  min_size_um = config$min_particle_size_um
)
# Raman: use ALL particles for spatial matching (no HQI filter here).
# HQI filtering is applied only inside analyze_agreement() so that
# spatial matching has maximum coverage while agreement scoring
# only considers particles with reliable Raman identifications.
raman_for_match <- prefilter_raman(
  raman_raw,
  min_hqi     = 0,
  min_size_um = config$min_particle_size_um
)
log_message("Raman for matching: ", nrow(raman_for_match),
            " particles (all HQI, spatial only; HQI filter applied in agreement)")

# ---------------------------------------------------------------------------
# 4. Coordinate normalization
#
# Compute centroids from the PLASTIC alignment subsets so that centering
# reflects the particle population that actually overlaps.
# ---------------------------------------------------------------------------

norm_result <- normalize_coordinates(
  ftir_for_align, raman_for_align,
  normalize_scale = config$normalize_scale
)

ftir_norm_align  <- norm_result$ftir
raman_norm_align <- norm_result$raman

# Also normalize ALL clean particles using the SAME centroids (for ICP & diagnostics)
ftir_clean$x_norm  <- ftir_clean$x_um - norm_result$ftir_centroid[1]
ftir_clean$y_norm  <- ftir_clean$y_um - norm_result$ftir_centroid[2]
raman_clean$x_norm <- raman_clean$x_um - norm_result$raman_centroid[1]
raman_clean$y_norm <- raman_clean$y_um - norm_result$raman_centroid[2]

# Build spatial transform set: all Raman >= 20 um (visible to FTIR).
# This set is used for ALL spatial transform steps (landmarks, ICP)
# regardless of material or HQI — only size matters for geometry.
min_size_spatial <- config$align_raman_min_size_um
raman_for_transform <- raman_clean
if (!is.null(min_size_spatial) && min_size_spatial > 0 &&
    any(!is.na(raman_clean$feret_max_um))) {
  raman_for_transform <- raman_clean[
    is.na(raman_clean$feret_max_um) | raman_clean$feret_max_um >= min_size_spatial, ]
}
log_message("Raman for spatial transform (>= ", min_size_spatial, " um): ",
            nrow(raman_for_transform), " particles")
log_message("ICP refinement sets: FTIR ", nrow(ftir_clean),
            ", Raman ", nrow(raman_for_transform))

# ---------------------------------------------------------------------------
# 5. Tiered alignment (FTIR ↔ Raman)
#
# Tier 1 — Landmark alignment: use large particles & fibers to quickly
#   determine the spatial transform. If confident, skip Tier 2.
# Tier 2 — Full RANSAC: exhaustive grid search on material-filtered anchors.
#   Only runs if Tier 1 was not confident enough.
# ICP refinement always runs to polish the transform.
# ---------------------------------------------------------------------------

# --- Tier 1: Landmark alignment ---
# Use size-filtered Raman (>= 20 um) — landmarks are selected by size inside
landmark_result <- landmark_align(ftir_clean, raman_for_transform, config)

use_landmark_transform <- landmark_result$confident && config$landmark_skip_full_ransac

if (use_landmark_transform) {
  log_message("Using landmark transform (Tier 1) — skipping full RANSAC")
  alignment_transform <- landmark_result$transform
  alignment_method    <- "landmark"
} else {
  # --- Tier 2: Full material-based RANSAC ---
  log_message(strrep("-", 50))
  log_message("Tier 2: Full RANSAC alignment (material-based anchors)")

  ransac_result <- ransac_align(ftir_norm_align, raman_norm_align, config)

  log_message("RANSAC transform: scale = ", round(ransac_result$params$scale, 4),
              ", rotation = ", round(ransac_result$params$rotation_deg, 2), " deg",
              ", reflected = ", ransac_result$params$reflected,
              ", inliers = ", ransac_result$n_inliers)

  alignment_transform <- ransac_result$transform
  alignment_method    <- "ransac"
}

# --- ICP refinement (always runs to polish the transform) ---
icp_result <- icp_refine(ftir_clean, raman_for_transform, alignment_transform, config)

log_message("ICP refined transform: scale = ", round(icp_result$params$scale, 4),
            ", rotation = ", round(icp_result$params$rotation_deg, 2), " deg",
            ", converged = ", icp_result$converged)

# ---------------------------------------------------------------------------
# 7. Apply transform to particles for matching
#
# Normalize using the SAME centroids from step 4 (plastic-based),
# then apply the ICP-refined transform.
# ---------------------------------------------------------------------------

# Normalize filtered FTIR using alignment centroids
ftir_for_match$x_norm <- ftir_for_match$x_um - norm_result$ftir_centroid[1]
ftir_for_match$y_norm <- ftir_for_match$y_um - norm_result$ftir_centroid[2]

# Normalize filtered Raman using alignment centroids
raman_for_match$x_norm <- raman_for_match$x_um - norm_result$raman_centroid[1]
raman_for_match$y_norm <- raman_for_match$y_um - norm_result$raman_centroid[2]

# Apply ICP-refined transform to filtered FTIR
ftir_aligned <- apply_ftir_transform(ftir_for_match, icp_result$transform)

# Also align the full FTIR set for diagnostics
ftir_aligned_all <- apply_ftir_transform(ftir_clean, icp_result$transform)

# ---------------------------------------------------------------------------
# 8. Particle matching (FTIR ↔ Raman, spatial)
# ---------------------------------------------------------------------------

match_result <- match_particles(ftir_aligned, raman_for_match, config)

# ---------------------------------------------------------------------------
# 9. Agreement analysis (FTIR ↔ Raman)
# ---------------------------------------------------------------------------

agreement <- analyze_agreement(match_result, config)

# ---------------------------------------------------------------------------
# 10. TPS assessment (systematic distortion check)
# ---------------------------------------------------------------------------

tps_assessment <- assess_tps_need(match_result)
log_message("TPS assessment: ", tps_assessment$message)
if (isTRUE(tps_assessment$recommend_tps)) {
  log_message("  RECOMMENDATION: Consider thin-plate spline warp for improved alignment",
              level = "WARN")
}

# ---------------------------------------------------------------------------
# 11. Composite matching (1:many) for unmatched FTIR
# ---------------------------------------------------------------------------

composites <- data.frame()
if (nrow(match_result$unmatched_ftir) > 0 && nrow(match_result$unmatched_raman) > 0) {
  # Ensure unmatched FTIR particles have aligned coordinates
  unmatched_ftir_src <- match_result$unmatched_ftir
  if (!"x_aligned" %in% names(unmatched_ftir_src)) {
    unmatched_ftir_src$x_norm <- unmatched_ftir_src$x_um - norm_result$ftir_centroid[1]
    unmatched_ftir_src$y_norm <- unmatched_ftir_src$y_um - norm_result$ftir_centroid[2]
    tf <- apply_transform_points(unmatched_ftir_src$x_norm,
                                 unmatched_ftir_src$y_norm,
                                 icp_result$transform)
    unmatched_ftir_src$x_aligned <- tf$x_transformed
    unmatched_ftir_src$y_aligned <- tf$y_transformed
  }

  composites <- find_composite_matches(
    unmatched_ftir_src, match_result$unmatched_raman, config
  )

  if (nrow(composites) > 0) {
    log_message("Composite matches found: ", nrow(composites),
                " FTIR particles matched to multiple Raman fragments")
  }
}

# ---------------------------------------------------------------------------
# 12. LDIR spatial pipeline (image → coordinates → alignment → matching)
# ---------------------------------------------------------------------------

ldir_results <- NULL
if (has_ldir && !is.null(ldir_raw)) {
  log_message(strrep("-", 50))
  log_message("LDIR Spatial Pipeline")

  # 12a. Pre-filter LDIR particles
  ldir_clean <- prefilter_ldir(
    ldir_raw,
    min_quality = config$ldir_quality_threshold,
    min_size_um = config$min_particle_size_um
  )

  # 12b. Extract coordinates from LDIR image (if available)
  ldir_with_coords <- ldir_clean
  has_ldir_coords  <- FALSE

  if (!is.null(config$ldir_image) && nzchar(config$ldir_image)) {
    log_message("Extracting LDIR coordinates from companion image")

    # Compute scan bounds from configured scan diameter (circular filter)
    ldir_scan_diam <- config$ldir_scan_diameter_um
    if (is.null(ldir_scan_diam)) ldir_scan_diam <- 13000
    ldir_scan_bounds <- list(
      x_min = 0, x_max = ldir_scan_diam,
      y_min = 0, y_max = ldir_scan_diam
    )

    ldir_image_particles <- extract_ldir_image_coords(
      config$ldir_image,
      scan_bounds    = ldir_scan_bounds,
      expected_count = nrow(ldir_clean)
    )

    # Save raw image-extracted coordinates (before join) for Shiny viewer
    ldir_image_extracted <- ldir_image_particles

    # Join image coordinates with Excel data via size-based Hungarian matching
    ldir_with_coords <- join_ldir_coords(ldir_clean, ldir_image_particles)

    # Validate join quality via scan-order correlation
    scan_order <- validate_ldir_scan_order(ldir_with_coords)
    log_message("  Scan order validation: ", scan_order$message)

    n_with_coords <- sum(!is.na(ldir_with_coords$x_um))
    has_ldir_coords <- n_with_coords >= 10
    log_message("  LDIR particles with coordinates: ", n_with_coords,
                " of ", nrow(ldir_with_coords))
  }

  # 12c–j. Spatial alignment & matching (only if coordinates available)
  ldir_raman_match     <- NULL
  ldir_ftir_match      <- NULL
  ldir_raman_agreement <- NULL
  ldir_icp             <- NULL
  ldir_aligned         <- NULL
  triplets             <- data.frame()

  if (has_ldir_coords) {
    log_message("LDIR spatial matching enabled")

    # 12c. Normalize LDIR coordinates (center using LDIR's own centroid)
    ldir_valid    <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
    ldir_centroid <- c(mean(ldir_valid$x_um), mean(ldir_valid$y_um))
    ldir_with_coords$x_norm <- ldir_with_coords$x_um - ldir_centroid[1]
    ldir_with_coords$y_norm <- ldir_with_coords$y_um - ldir_centroid[2]

    # 12d. Tiered LDIR→Raman alignment
    #
    # Tier 1 — Landmark alignment: use large particles & fibers
    #   (same strategy as FTIR→Raman). Landmarks are instrument-agnostic
    #   geometric features that are visible across all instruments.
    # Tier 2 — Material-based RANSAC: PET/PP/PC anchors (fallback).
    ldir_landmark_result <- tryCatch({
      landmark_align(ldir_valid, raman_for_transform, config, src_label = "LDIR")
    }, error = function(e) {
      log_message("  LDIR landmark alignment failed: ", e$message, level = "WARN")
      list(success = FALSE, confident = FALSE,
           n_ftir_landmarks = 0, n_raman_landmarks = 0)
    })

    use_ldir_landmark <- ldir_landmark_result$confident &&
                         config$landmark_skip_full_ransac

    if (use_ldir_landmark) {
      log_message("  Using LDIR landmark transform (Tier 1) — skipping material RANSAC")
      ldir_ransac <- list(
        transform = ldir_landmark_result$transform,
        params    = ldir_landmark_result$params,
        n_inliers = ldir_landmark_result$n_inliers
      )
    } else {
      # Tier 2: Material-based RANSAC (fallback)
      log_message(strrep("-", 50))
      log_message("  LDIR Tier 2: Material-based RANSAC alignment")
      ldir_for_align <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
      if (!is.null(config$align_ldir_materials)) {
        ldir_mat_mask <- grepl(
          paste(config$align_ldir_materials, collapse = "|"),
          ldir_for_align$material, ignore.case = TRUE
        )
        if (sum(ldir_mat_mask) >= 4) {
          ldir_for_align <- ldir_for_align[ldir_mat_mask, ]
        } else {
          log_message("  Insufficient LDIR anchor materials (",
                      sum(ldir_mat_mask),
                      ") — using all LDIR particles for alignment")
        }
      }
      log_message("  LDIR alignment anchors: ", nrow(ldir_for_align), " particles")

      # 12e. LDIR→Raman alignment (RANSAC + ICP)
      # Raman reference: use the same normalized set from step 4
      ldir_ransac <- tryCatch({
        ransac_align(ldir_for_align, raman_norm_align, config)
      }, error = function(e) {
        log_message("  LDIR RANSAC failed: ", e$message, level = "WARN")
        NULL
      })
    }

    if (!is.null(ldir_ransac)) {
      log_message("  LDIR-Raman RANSAC: scale=",
                  round(ldir_ransac$params$scale, 4),
                  ", rot=", round(ldir_ransac$params$rotation_deg, 2), " deg",
                  ", reflected=", ldir_ransac$params$reflected,
                  ", inliers=", ldir_ransac$n_inliers)

      # ICP refinement on all LDIR particles with coordinates
      ldir_for_icp <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
      ldir_icp <- tryCatch({
        icp_refine(ldir_for_icp, raman_for_transform,
                   ldir_ransac$transform, config)
      }, error = function(e) {
        log_message("  LDIR ICP failed: ", e$message,
                    " — using RANSAC transform", level = "WARN")
        list(
          transform    = ldir_ransac$transform,
          params       = ldir_ransac$params,
          converged    = FALSE,
          n_iterations = 0,
          rms_history  = numeric(0)
        )
      })

      log_message("  LDIR-Raman ICP: scale=",
                  round(ldir_icp$params$scale, 4),
                  ", rot=", round(ldir_icp$params$rotation_deg, 2), " deg",
                  ", converged=", ldir_icp$converged)

      # Quality check: warn if LDIR ICP alignment is substantially worse than FTIR
      ldir_final_rms <- if (length(ldir_icp$rms_history) > 0)
        tail(ldir_icp$rms_history, 1) else NA_real_
      ftir_final_rms  <- if (length(icp_result$rms_history) > 0)
        tail(icp_result$rms_history, 1) else NA_real_

      if (!is.na(ldir_final_rms) && ldir_final_rms > 100) {
        log_message("  LDIR alignment quality: POOR (ICP RMS = ",
                    round(ldir_final_rms, 1), " µm). Possible causes:",
                    level = "WARN")
        log_message("    1. LDIR scan area (ldir_scan_diameter_um=",
                    config$ldir_scan_diameter_um, " µm) may not match actual scan",
                    level = "WARN")
        log_message("    2. Too few anchor material particles for robust RANSAC",
                    level = "WARN")
        log_message("    3. Systematic coordinate flip or offset in image extraction",
                    level = "WARN")
        log_message("    Recommendation: verify ldir_scan_diameter_um in config ",
                    "and check LDIR overlay plot (plots/ldir_overlay.png)",
                    level = "WARN")
      } else if (!is.na(ldir_final_rms) && !is.na(ftir_final_rms) &&
                 ldir_final_rms > 3 * ftir_final_rms) {
        log_message("  LDIR alignment quality: MARGINAL (RMS ",
                    round(ldir_final_rms, 1), " µm vs FTIR ",
                    round(ftir_final_rms, 1), " µm)", level = "WARN")
      }

      # 12f. Apply transform to all LDIR particles with coordinates
      ldir_aligned <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
      ldir_tf <- apply_transform_points(
        ldir_aligned$x_norm, ldir_aligned$y_norm, ldir_icp$transform
      )
      ldir_aligned$x_aligned <- ldir_tf$x_transformed
      ldir_aligned$y_aligned <- ldir_tf$y_transformed

      # 12g. LDIR↔Raman matching
      ldir_raman_match <- match_particles(
        ldir_aligned, raman_for_match, config,
        src_label = "ldir", ref_label = "raman"
      )

      # 12h. LDIR↔Raman agreement analysis
      ldir_raman_agreement <- analyze_agreement(
        ldir_raman_match, config,
        instrument_a = "LDIR", instrument_b = "Raman"
      )

      # 12i. LDIR↔FTIR matching (both already in Raman coordinate frame)
      # Set up aligned FTIR as reference: matcher expects ref with x_norm/y_norm
      ftir_as_ref <- ftir_aligned
      ftir_as_ref$x_norm <- ftir_as_ref$x_aligned
      ftir_as_ref$y_norm <- ftir_as_ref$y_aligned
      ldir_ftir_match <- match_particles(
        ldir_aligned, ftir_as_ref, config,
        src_label = "ldir", ref_label = "ftir"
      )
      log_message("  LDIR-FTIR direct matching: ",
                  ldir_ftir_match$match_stats$n_matched, " pairs")

      # 12j. Three-way triplets (FTIR↔Raman ∩ LDIR↔Raman via Raman ID)
      if (nrow(match_result$matched) > 0 &&
          nrow(ldir_raman_match$matched) > 0) {
        ftir_raman_pairs <- data.frame(
          raman_particle_id   = match_result$matched$raman_particle_id,
          ftir_particle_id    = match_result$matched$ftir_particle_id,
          ftir_material       = match_result$matched$ftir_material,
          ftir_raman_distance = match_result$matched$match_distance,
          stringsAsFactors = FALSE
        )
        ldir_raman_pairs <- data.frame(
          raman_particle_id   = ldir_raman_match$matched$raman_particle_id,
          ldir_particle_id    = ldir_raman_match$matched$ldir_particle_id,
          ldir_material       = ldir_raman_match$matched$ldir_material,
          ldir_raman_distance = ldir_raman_match$matched$match_distance,
          stringsAsFactors = FALSE
        )

        triplets <- merge(ftir_raman_pairs, ldir_raman_pairs,
                          by = "raman_particle_id")

        # Add Raman material
        raman_mat_idx <- match(triplets$raman_particle_id,
                               raman_for_match$particle_id)
        triplets$raman_material <- raman_for_match$material[raman_mat_idx]

        # Compute per-triplet material agreement quality (0–3 instruments agreeing)
        # Uses polymer family classification so minor name differences still score
        if (nrow(triplets) > 0) {
          ftir_fam  <- classify_family_vec(triplets$ftir_material)
          raman_fam <- classify_family_vec(triplets$raman_material)
          ldir_fam  <- classify_family_vec(triplets$ldir_material)

          # Count how many instrument families agree with Raman
          fr_agree   <- ftir_fam == raman_fam & ftir_fam != "Unknown"
          lr_agree   <- ldir_fam == raman_fam & ldir_fam != "Unknown"
          fl_agree   <- ftir_fam == ldir_fam  & ftir_fam != "Unknown"
          n_agree    <- as.integer(fr_agree) + as.integer(lr_agree) +
                        as.integer(fl_agree)
          # n_agree ranges 0–3 (3 = all three pairwise family comparisons agree)
          triplets$n_instrument_agreement <- n_agree
          triplets$ftir_family  <- ftir_fam
          triplets$raman_family <- raman_fam
          triplets$ldir_family  <- ldir_fam
          triplets$material_consensus <- ifelse(
            n_agree == 3, "Full agreement",
            ifelse(fr_agree & lr_agree, "All three agree",
            ifelse(fr_agree, "FTIR+Raman agree",
            ifelse(lr_agree, "LDIR+Raman agree",
            ifelse(fl_agree, "FTIR+LDIR agree",
            "No agreement"))))
          )

          n_full <- sum(n_agree == 3, na.rm = TRUE)
          n_partial <- sum(n_agree >= 2 & n_agree < 3, na.rm = TRUE)
          log_message("  Three-way triplets: ", nrow(triplets),
                      " total (", n_full, " full agreement, ",
                      n_partial, " partial, ",
                      nrow(triplets) - n_full - n_partial, " no agreement)")
        } else {
          log_message("  Three-way triplets: 0 particles matched")
        }
      }

    } else {
      log_message("  LDIR spatial alignment failed — ",
                  "material comparison only", level = "WARN")
    }

  } else {
    log_message("  LDIR: no spatial coordinates — material comparison only")
  }

  # Non-spatial material distribution (always available)
  ldir_mats <- table(ldir_clean$material)
  log_message("  LDIR material distribution (top 5):")
  for (m in names(head(sort(ldir_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", ldir_mats[m])
  }

  ldir_results <- list(
    ldir_clean           = ldir_clean,
    ldir_with_coords     = ldir_with_coords,
    ldir_aligned         = ldir_aligned,
    ldir_raman_match     = ldir_raman_match,
    ldir_ftir_match      = ldir_ftir_match,
    ldir_raman_agreement = ldir_raman_agreement,
    ldir_icp             = ldir_icp,
    ldir_landmark        = if (exists("ldir_landmark_result")) ldir_landmark_result else NULL,
    ldir_alignment_method = if (exists("use_ldir_landmark") && use_ldir_landmark) "landmark" else "ransac",
    triplets             = triplets,
    has_coords           = has_ldir_coords,
    ldir_material_dist   = ldir_mats,
    ldir_image_extracted = if (exists("ldir_image_extracted")) ldir_image_extracted else NULL
  )
}

# ---------------------------------------------------------------------------
# 13. Diagnostics (use full datasets for overlay, matched pairs for detail)
# ---------------------------------------------------------------------------

diagnostics <- generate_diagnostics(
  ftir_aligned_all, raman_clean,
  match_result, icp_result, agreement,
  config
)

# Add LDIR diagnostics if spatial matching was performed
if (!is.null(ldir_results) && !is.null(ldir_results$ldir_aligned) &&
    !is.null(ldir_results$ldir_raman_match)) {
  ldir_diag <- generate_ldir_diagnostics(
    ldir_results$ldir_aligned, raman_clean,
    ldir_results$ldir_raman_match, ldir_results$ldir_raman_agreement
  )
  diagnostics <- c(diagnostics, ldir_diag)
  log_message("  Added ", length(ldir_diag), " LDIR diagnostic plots")
}

# ---------------------------------------------------------------------------
# 14. Export
# ---------------------------------------------------------------------------

export_results(
  match_result, agreement, diagnostics,
  icp_result, norm_result, config,
  ftir_scan_bounds = ftir_scan_bounds,
  ldir_results     = ldir_results
)

# Export composite matches if found
if (nrow(composites) > 0) {
  write.csv(composites,
            file.path(config$output_dir, "composite_matches.csv"),
            row.names = FALSE)
  log_message("  Wrote composite_matches.csv (", nrow(composites), " composites)")
}

# Export TPS assessment
if (!is.null(tps_assessment$quadrant_residuals)) {
  write.csv(tps_assessment$quadrant_residuals,
            file.path(config$output_dir, "tps_quadrant_residuals.csv"),
            row.names = FALSE)
  tps_lines <- c(
    "# TPS (Thin-Plate Spline) Assessment",
    paste0("# Generated: ", Sys.time()),
    "",
    paste0("systematic_distortion: ", tps_assessment$systematic_distortion),
    paste0("recommend_tps:         ", tps_assessment$recommend_tps),
    paste0("dx_range_um:           ", tps_assessment$dx_range_um),
    paste0("dy_range_um:           ", tps_assessment$dy_range_um),
    paste0("assessment:            ", tps_assessment$message)
  )
  writeLines(tps_lines, file.path(config$output_dir, "tps_assessment.txt"))
}

# ---------------------------------------------------------------------------
# 15. Summary
# ---------------------------------------------------------------------------

log_message(strrep("=", 60))
log_message("Pipeline complete!")
log_message(strrep("=", 60))
log_message("  Alignment method:  ", alignment_method,
            if (alignment_method == "landmark")
              paste0(" (", landmark_result$n_inliers, " landmark inliers, ",
                     round(landmark_result$mean_residual, 1), " µm mean residual)")
            else "")
log_message("  FTIR landmarks:    ", landmark_result$n_ftir_landmarks,
            " (confident: ", landmark_result$confident, ")")
log_message("  FTIR plastic anchors (alignment): ", nrow(ftir_for_align))
log_message("  Raman targets (alignment):       ", nrow(raman_for_align))
log_message("  FTIR particles (matching):       ", nrow(ftir_for_match))
log_message("  Raman particles (matching):      ", nrow(raman_for_match))
log_message("  Matched pairs:      ", match_result$match_stats$n_matched)
log_message("  Unmatched FTIR:     ", match_result$match_stats$n_unmatched_ftir)
log_message("  Unmatched Raman:    ", match_result$match_stats$n_unmatched_raman)
log_message("  FTIR match rate:    ",
            round(match_result$match_stats$match_rate_ftir * 100, 1), "%")
if (!is.na(agreement$agreement_rate)) {
  log_message("  Material agreement: ", round(agreement$agreement_rate * 100, 1), "%")
}
if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
  tr <- agreement$tiered_rates
  log_message("  Tiered agreement:  Exact ", tr$exact_pct, "%, Family+ ",
              tr$family_or_better_pct, "%")
}
if (nrow(composites) > 0) {
  log_message("  Composite matches: ", nrow(composites))
}
log_message("  TPS assessment:    ", tps_assessment$message)
if (has_ldir && !is.null(ldir_results)) {
  log_message("  --- LDIR ---")
  log_message("  LDIR particles:    ", nrow(ldir_results$ldir_clean))
  log_message("  LDIR coordinates:  ",
              if (ldir_results$has_coords) "spatial (from image)" else "none")
  if (!is.null(ldir_results$ldir_landmark)) {
    log_message("  LDIR alignment:    ", ldir_results$ldir_alignment_method,
                " (landmarks: ", ldir_results$ldir_landmark$n_ftir_landmarks,
                ", confident: ", ldir_results$ldir_landmark$confident, ")")
  }
  if (!is.null(ldir_results$ldir_raman_match)) {
    log_message("  LDIR-Raman matched: ",
                ldir_results$ldir_raman_match$match_stats$n_matched)
  }
  if (!is.null(ldir_results$ldir_ftir_match)) {
    log_message("  LDIR-FTIR matched:  ",
                ldir_results$ldir_ftir_match$match_stats$n_matched)
  }
  if (nrow(ldir_results$triplets) > 0) {
    log_message("  Three-way triplets: ", nrow(ldir_results$triplets))
  }
  if (!is.null(ldir_results$ldir_raman_agreement) &&
      !is.null(ldir_results$ldir_raman_agreement$tiered_rates) &&
      ldir_results$ldir_raman_agreement$tiered_rates$n_total > 0) {
    ltr <- ldir_results$ldir_raman_agreement$tiered_rates
    log_message("  LDIR-Raman agreement: Exact ", ltr$exact_pct,
                "%, Family+ ", ltr$family_or_better_pct, "%")
  }
}
log_message("  Results in: ", config$output_dir)
