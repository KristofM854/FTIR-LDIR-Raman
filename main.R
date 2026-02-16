# =============================================================================
# main.R — FTIR–Raman particle matching pipeline
# =============================================================================
#
# Orchestrates the full pipeline for aligning and matching particles
# detected by FTIR and Raman microspectroscopy on the same physical filter.
#
# Usage:
#   1. Set ftir_path and raman_path below to your paired data files
#   2. Adjust config parameters as needed
#   3. Source this script:  source("main.R")
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
source("R/01_ingest.R")
source("R/02_prefilter.R")
source("R/03_normalize.R")
source("R/03b_landmark_align.R")
source("R/04_ransac.R")
source("R/05_transform.R")
source("R/06_icp_refine.R")
source("R/07_match.R")
source("R/08_agreement.R")
source("R/09_diagnostics.R")
source("R/10_export.R")

# ---------------------------------------------------------------------------
# 1. Configuration
# ---------------------------------------------------------------------------

# File paths — leave as NULL to get a file-picker dialog at runtime
ftir_file  <- NULL   # e.g. "Comparstic Spotlight F2Ba Au 240926.csv"
raman_file <- NULL   # e.g. "Comparstic Raman F2Ba Au IAEA 240930.csv"

if (is.null(ftir_file)) {
  message("Please select the FTIR data file (.csv or .xlsx) ...")
  ftir_file <- file.choose()
}
if (is.null(raman_file)) {
  message("Please select the Raman data file (.csv or .xlsx) ...")
  raman_file <- file.choose()
}

config <- make_config(
  ftir_path  = ftir_file,
  raman_path = raman_file,
  output_dir = "output"
)

# Override any defaults as needed:
# config$raman_hqi_threshold     <- 75
# config$match_dist_threshold_um <- 25
# config$ransac_allow_mirror     <- TRUE

# ---------------------------------------------------------------------------
# 2. Data ingestion
# ---------------------------------------------------------------------------

log_message(strrep("=", 60))
log_message("FTIR-Raman Particle Matching Pipeline")
log_message(strrep("=", 60))

ftir_raw  <- ingest_ftir(config$ftir_path, sheet = config$ftir_sheet)
raman_raw <- ingest_raman(config$raman_path, sheet = config$raman_sheet)

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
log_message("Raman alignment target particles: ", nrow(raman_for_align))

# Quality-filtered sets for matching (applied after alignment)
ftir_for_match <- prefilter_ftir(
  ftir_raw,
  min_quality = config$ftir_quality_threshold,
  min_size_um = config$min_particle_size_um
)
raman_for_match <- prefilter_raman(
  raman_raw,
  min_hqi     = config$raman_hqi_threshold,
  min_size_um = config$min_particle_size_um
)

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

# Build ICP sets: all Raman particles >= size threshold (visible to FTIR)
min_size_icp <- config$align_raman_min_size_um
raman_for_icp <- raman_clean
if (!is.null(min_size_icp) && min_size_icp > 0 && any(!is.na(raman_clean$feret_max_um))) {
  raman_for_icp <- raman_clean[
    is.na(raman_clean$feret_max_um) | raman_clean$feret_max_um >= min_size_icp, ]
}
log_message("ICP refinement sets: FTIR ", nrow(ftir_clean),
            ", Raman ", nrow(raman_for_icp), " (>= ", min_size_icp, " um)")

# ---------------------------------------------------------------------------
# 5. Tiered alignment
#
# Tier 1 — Landmark alignment: use large particles & fibers to quickly
#   determine the spatial transform. If confident, skip Tier 2.
# Tier 2 — Full RANSAC: exhaustive grid search on material-filtered anchors.
#   Only runs if Tier 1 was not confident enough.
# ICP refinement always runs to polish the transform.
# ---------------------------------------------------------------------------

# --- Tier 1: Landmark alignment ---
landmark_result <- landmark_align(ftir_clean, raman_clean, config)

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
icp_result <- icp_refine(ftir_clean, raman_for_icp, alignment_transform, config)

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
# 8. Particle matching (quality-filtered datasets only)
# ---------------------------------------------------------------------------

match_result <- match_particles(ftir_aligned, raman_for_match, config)

# ---------------------------------------------------------------------------
# 9. Agreement analysis
# ---------------------------------------------------------------------------

agreement <- analyze_agreement(match_result)

# ---------------------------------------------------------------------------
# 10. Diagnostics (use full datasets for overlay, matched pairs for detail)
# ---------------------------------------------------------------------------

diagnostics <- generate_diagnostics(
  ftir_aligned_all, raman_clean,
  match_result, icp_result, agreement,
  config
)

# ---------------------------------------------------------------------------
# 11. Export
# ---------------------------------------------------------------------------

export_results(
  match_result, agreement, diagnostics,
  icp_result, norm_result, config
)

# ---------------------------------------------------------------------------
# 12. Summary
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
log_message("  Results in: ", config$output_dir)
