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

required_packages <- c("readxl", "dplyr", "ggplot2", "RANN")

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
  ftir_path  = "processed_data_ftir_20260211_102404.xlsx",
  raman_path = "processed_data_raman_20260211_120938.xlsx",
  output_dir = "output"
)

# Override any defaults as needed:
# config$raman_hqi_threshold    <- 75
# config$match_dist_threshold_um <- 25
# config$ransac_allow_mirror     <- TRUE

# ---------------------------------------------------------------------------
# 2. Data ingestion
# ---------------------------------------------------------------------------

log_message("=" |> strrep(60))
log_message("FTIR-Raman Particle Matching Pipeline")
log_message("=" |> strrep(60))

ftir_raw  <- ingest_ftir(config$ftir_path, sheet = config$ftir_sheet)
raman_raw <- ingest_raman(config$raman_path, sheet = config$raman_sheet)

# ---------------------------------------------------------------------------
# 3. Pre-filtering
# ---------------------------------------------------------------------------

ftir_clean <- prefilter_ftir(
  ftir_raw,
  min_quality = config$ftir_quality_threshold,
  min_size_um = config$min_particle_size_um
)

raman_clean <- prefilter_raman(
  raman_raw,
  min_hqi     = config$raman_hqi_threshold,
  min_size_um = config$min_particle_size_um
)

# ---------------------------------------------------------------------------
# 4. Coordinate normalization
# ---------------------------------------------------------------------------

norm_result <- normalize_coordinates(
  ftir_clean, raman_clean,
  normalize_scale = config$normalize_scale
)

ftir_norm  <- norm_result$ftir
raman_norm <- norm_result$raman

# ---------------------------------------------------------------------------
# 5. RANSAC alignment estimation
# ---------------------------------------------------------------------------

ransac_result <- ransac_align(ftir_norm, raman_norm, config)

log_message("RANSAC transform: scale = ", round(ransac_result$params$scale, 4),
            ", rotation = ", round(ransac_result$params$rotation_deg, 2), "°",
            ", reflected = ", ransac_result$params$reflected,
            ", inliers = ", ransac_result$n_inliers)

# ---------------------------------------------------------------------------
# 6. Apply RANSAC transform to FTIR
# ---------------------------------------------------------------------------

ftir_aligned <- apply_ftir_transform(ftir_norm, ransac_result$transform)

# ---------------------------------------------------------------------------
# 7. ICP refinement
# ---------------------------------------------------------------------------

icp_result <- icp_refine(ftir_norm, raman_norm, ransac_result$transform, config)

log_message("ICP refined transform: scale = ", round(icp_result$params$scale, 4),
            ", rotation = ", round(icp_result$params$rotation_deg, 2), "°",
            ", converged = ", icp_result$converged)

# Update FTIR alignment with ICP-refined transform
ftir_aligned <- apply_ftir_transform(ftir_norm, icp_result$transform)

# ---------------------------------------------------------------------------
# 8. Particle matching
# ---------------------------------------------------------------------------

match_result <- match_particles(ftir_aligned, raman_norm, config)

# ---------------------------------------------------------------------------
# 9. Agreement analysis
# ---------------------------------------------------------------------------

agreement <- analyze_agreement(match_result)

# ---------------------------------------------------------------------------
# 10. Diagnostics
# ---------------------------------------------------------------------------

diagnostics <- generate_diagnostics(
  ftir_aligned, raman_norm,
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

log_message("=" |> strrep(60))
log_message("Pipeline complete!")
log_message("=" |> strrep(60))
log_message("  Matched pairs:      ", match_result$match_stats$n_matched)
log_message("  Unmatched FTIR:     ", match_result$match_stats$n_unmatched_ftir)
log_message("  Unmatched Raman:    ", match_result$match_stats$n_unmatched_raman)
log_message("  FTIR match rate:    ",
            round(match_result$match_stats$match_rate_ftir * 100, 1), "%")
if (!is.na(agreement$agreement_rate)) {
  log_message("  Material agreement: ", round(agreement$agreement_rate * 100, 1), "%")
}
log_message("  Results in: ", config$output_dir)
