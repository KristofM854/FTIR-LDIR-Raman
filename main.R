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

# Set your file paths here
config <- make_config(
  ftir_path  = "Comparstic Spotlight F2Ba Au 240926.csv",
  raman_path = "Comparstic Raman F2Ba Au IAEA 240930.csv",
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
# Strategy: use ALL particles (with valid coords) for spatial alignment,
# then apply quality filters only for the final matching & agreement steps.
# This maximizes the number of spatial reference points for alignment.
# ---------------------------------------------------------------------------

# Minimal filtering for alignment (just remove invalid coords)
ftir_for_align  <- prefilter_ftir(ftir_raw, min_quality = 0, min_size_um = 0)
raman_for_align <- prefilter_raman(raman_raw, min_hqi = 0, min_size_um = 0)

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
# 4. Coordinate normalization (using ALL particles for robust centering)
# ---------------------------------------------------------------------------

norm_result <- normalize_coordinates(
  ftir_for_align, raman_for_align,
  normalize_scale = config$normalize_scale
)

ftir_norm_all  <- norm_result$ftir
raman_norm_all <- norm_result$raman

# ---------------------------------------------------------------------------
# 5. RANSAC alignment estimation (using ALL particles)
# ---------------------------------------------------------------------------

ransac_result <- ransac_align(ftir_norm_all, raman_norm_all, config)

log_message("RANSAC transform: scale = ", round(ransac_result$params$scale, 4),
            ", rotation = ", round(ransac_result$params$rotation_deg, 2), " deg",
            ", reflected = ", ransac_result$params$reflected,
            ", inliers = ", ransac_result$n_inliers)

# ---------------------------------------------------------------------------
# 6. ICP refinement (using ALL particles)
# ---------------------------------------------------------------------------

icp_result <- icp_refine(ftir_norm_all, raman_norm_all, ransac_result$transform, config)

log_message("ICP refined transform: scale = ", round(icp_result$params$scale, 4),
            ", rotation = ", round(icp_result$params$rotation_deg, 2), " deg",
            ", converged = ", icp_result$converged)

# ---------------------------------------------------------------------------
# 7. Apply transform to QUALITY-FILTERED particles for matching
#
# Normalize the filtered datasets using the SAME centroids from step 4,
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
ftir_aligned_all <- apply_ftir_transform(ftir_norm_all, icp_result$transform)

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
  ftir_aligned_all, raman_norm_all,
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
log_message("  FTIR particles (alignment): ", nrow(ftir_for_align))
log_message("  Raman particles (alignment):", nrow(raman_for_align))
log_message("  FTIR particles (matching):  ", nrow(ftir_for_match))
log_message("  Raman particles (matching): ", nrow(raman_for_match))
log_message("  Matched pairs:      ", match_result$match_stats$n_matched)
log_message("  Unmatched FTIR:     ", match_result$match_stats$n_unmatched_ftir)
log_message("  Unmatched Raman:    ", match_result$match_stats$n_unmatched_raman)
log_message("  FTIR match rate:    ",
            round(match_result$match_stats$match_rate_ftir * 100, 1), "%")
if (!is.na(agreement$agreement_rate)) {
  log_message("  Material agreement: ", round(agreement$agreement_rate * 100, 1), "%")
}
log_message("  Results in: ", config$output_dir)
