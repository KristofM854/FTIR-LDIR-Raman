# =============================================================================
# load.R — load the pipeline's packages and modules
# =============================================================================
#
#   source("R/load.R"); load_pipeline()      # from the repository root
#
# then build a config and call run_pipeline(config) (R/run_pipeline.R).
# Sourcing this file only defines the two objects below.
# =============================================================================

# clue is required, not optional: without it match_particles() would switch
# from optimal (Hungarian) to greedy matching and change the results.
pipeline_required_packages <- c("readxl", "ggplot2", "RANN", "ggrepel", "clue",
                                "jsonlite", "png")

pipeline_modules <- c(
  "read_bmp.R",
  "measure_raman_placement.R",
  "tps_refine.R",
  "utils.R",
  "particle_contract.R",
  "00_config.R",
  "00b_file_input.R",
  "01_ingest.R",
  "01_readers.R",
  "01b_ingest_image.R",
  "ftir_image_placement.R",
  "utils_python.R",
  "01c_ingest_ldir.R",
  "02_prefilter.R",
  "03_normalize.R",
  "03b_landmark_align.R",
  "03c_procrustes_align.R",
  "align_helpers.R",
  "04_ransac.R",
  "05_transform.R",
  "06_icp_refine.R",
  "07_match.R",
  "rate_intervals.R",
  "08_agreement.R",
  "08b_material_map.R",
  "09_diagnostics.R",
  "10_alignment_diagnostics.R",
  "10_export.R",
  "align_and_match.R",
  "detection_rates.R",
  "ldir_join_check.R",
  "ldir_stage.R",
  "run_pipeline.R"
)

#' Check the required packages, attach ggplot2 and source every module.
#'
#' @param root Repository root
#' @param envir Environment to define the functions in
load_pipeline <- function(root = ".", envir = globalenv()) {
  missing <- pipeline_required_packages[!vapply(pipeline_required_packages,
                                                requireNamespace, logical(1),
                                                quietly = TRUE)]
  if (length(missing) > 0)
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c(",
         paste0('"', missing, '"', collapse = ", "), "))", call. = FALSE)
  library(ggplot2)
  for (m in pipeline_modules) sys.source(file.path(root, "R", m), envir = envir)
  invisible(pipeline_modules)
}
