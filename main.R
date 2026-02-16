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
if (!is.null(config$ftir_image) && nzchar(config$ftir_image)) {
  log_message(strrep("-", 50))
  log_message("FTIR image-based particle extraction")

  # Estimate scan bounds from tabular FTIR data
  ftir_scan_bounds <- list(
    x_min = 0,
    x_max = max(ftir_raw$x_um, na.rm = TRUE) * 1.05,
    y_min = 0,
    y_max = max(ftir_raw$y_um, na.rm = TRUE) * 1.05
  )

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
# 10. LDIR comparison (non-spatial, material composition)
# ---------------------------------------------------------------------------

ldir_comparison <- NULL
if (has_ldir && !is.null(ldir_raw)) {
  log_message(strrep("-", 50))
  log_message("LDIR comparison (non-spatial: material composition)")

  ldir_clean <- prefilter_ldir(ldir_raw, min_quality = 0.6, min_size_um = 0)

  # Compare material distributions across all three instruments
  ftir_mats  <- table(ftir_clean$material)
  raman_mats <- table(raman_for_match$material)
  ldir_mats  <- table(ldir_clean$material)

  log_message("  FTIR material distribution (top 5):")
  for (m in names(head(sort(ftir_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", ftir_mats[m])
  }
  log_message("  Raman material distribution (top 5):")
  for (m in names(head(sort(raman_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", raman_mats[m])
  }
  log_message("  LDIR material distribution (top 5):")
  for (m in names(head(sort(ldir_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", ldir_mats[m])
  }

  ldir_comparison <- list(
    ldir_particles = ldir_clean,
    ftir_material_dist  = ftir_mats,
    raman_material_dist = raman_mats,
    ldir_material_dist  = ldir_mats,
    n_ldir = nrow(ldir_clean)
  )
}

# ---------------------------------------------------------------------------
# 11. Diagnostics (use full datasets for overlay, matched pairs for detail)
# ---------------------------------------------------------------------------

diagnostics <- generate_diagnostics(
  ftir_aligned_all, raman_clean,
  match_result, icp_result, agreement,
  config
)

# ---------------------------------------------------------------------------
# 12. Export
# ---------------------------------------------------------------------------

export_results(
  match_result, agreement, diagnostics,
  icp_result, norm_result, config
)

# Export LDIR comparison if available
if (!is.null(ldir_comparison)) {
  ldir_out <- file.path(config$output_dir, "ldir_material_distribution.csv")
  ldir_df <- data.frame(
    material = names(ldir_comparison$ldir_material_dist),
    count    = as.integer(ldir_comparison$ldir_material_dist),
    stringsAsFactors = FALSE
  )
  write.csv(ldir_df, ldir_out, row.names = FALSE)
  log_message("  Wrote LDIR material distribution to: ", basename(ldir_out))
}

# ---------------------------------------------------------------------------
# 13. Summary
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
if (has_ldir) {
  log_message("  LDIR particles:     ", ldir_comparison$n_ldir)
}
log_message("  Results in: ", config$output_dir)
