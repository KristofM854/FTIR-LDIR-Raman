# =============================================================================
# global.R — Load pipeline output and prepare data for Shiny
# =============================================================================

library(shiny)
library(ggplot2)
library(png)

# ---------------------------------------------------------------------------
# Locate the most recent pipeline run
# ---------------------------------------------------------------------------
find_latest_run <- function(output_dir = file.path("..", "output")) {
  runs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  if (length(runs) == 0) return(NULL)
  runs <- runs[order(file.mtime(runs), decreasing = TRUE)]
  runs[1]
}

# ---------------------------------------------------------------------------
# Load pipeline CSVs
# ---------------------------------------------------------------------------
load_run_data <- function(run_dir) {
  data <- list(run_dir = run_dir)

  paths <- list(
    matched   = "matched_particles.csv",
    unmatched_ftir  = "unmatched_ftir.csv",
    unmatched_raman = "unmatched_raman.csv",
    agreement = "agreement_summary.csv"
  )

  for (nm in names(paths)) {
    fp <- file.path(run_dir, paths[[nm]])
    if (file.exists(fp)) data[[nm]] <- read.csv(fp, stringsAsFactors = FALSE)
  }

  for (txt in c("match_statistics.txt", "transform_params.txt")) {
    fp <- file.path(run_dir, txt)
    if (file.exists(fp)) data[[gsub("\\.txt$", "", gsub("_", "_", txt))]] <- readLines(fp, warn = FALSE)
  }

  data
}

# ---------------------------------------------------------------------------
# Build per-instrument data frames from pipeline output
# ---------------------------------------------------------------------------
build_instrument_dfs <- function(data) {
  result <- list()

  # --- FTIR ---
  ftir_parts <- list()
  if (!is.null(data$matched) && nrow(data$matched) > 0) {
    m <- data$matched
    ftir_parts[[1]] <- data.frame(
      particle_id  = m$ftir_particle_id,
      x = m$ftir_x_aligned, y = m$ftir_y_aligned,
      x_orig = m$ftir_x_um, y_orig = m$ftir_y_um,
      area_um2 = m$ftir_area_um2,
      major_um = m$ftir_major_um, minor_um = m$ftir_minor_um,
      feret_max = m$ftir_feret_max_um,
      material = m$ftir_material, quality = m$ftir_quality,
      match_status = "matched", match_id = m$match_id,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_ftir) && nrow(data$unmatched_ftir) > 0) {
    u <- data$unmatched_ftir
    ftir_parts[[length(ftir_parts) + 1]] <- data.frame(
      particle_id = u$particle_id,
      x = u$x_aligned, y = u$y_aligned,
      x_orig = u$x_um, y_orig = u$y_um,
      area_um2 = u$area_um2,
      major_um = u$major_um, minor_um = u$minor_um,
      feret_max = u$feret_max_um,
      material = u$material, quality = u$quality,
      match_status = "unmatched", match_id = NA_integer_,
      stringsAsFactors = FALSE)
  }
  result$ftir <- do.call(rbind, ftir_parts)

  # --- Raman ---
  raman_parts <- list()
  if (!is.null(data$matched) && nrow(data$matched) > 0) {
    m <- data$matched
    raman_parts[[1]] <- data.frame(
      particle_id = m$raman_particle_id,
      x = m$raman_x_um, y = m$raman_y_um,
      x_orig = m$raman_x_um, y_orig = m$raman_y_um,
      area_um2 = m$raman_area_um2,
      major_um = m$raman_major_um, minor_um = m$raman_minor_um,
      feret_max = m$raman_feret_max_um,
      material = m$raman_material, quality = m$raman_quality,
      match_status = "matched", match_id = m$match_id,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_raman) && nrow(data$unmatched_raman) > 0) {
    u <- data$unmatched_raman
    raman_parts[[length(raman_parts) + 1]] <- data.frame(
      particle_id = u$particle_id,
      x = u$x_um, y = u$y_um,
      x_orig = u$x_um, y_orig = u$y_um,
      area_um2 = u$area_um2,
      major_um = u$major_um, minor_um = u$minor_um,
      feret_max = u$feret_max_um,
      material = u$material, quality = u$quality,
      match_status = "unmatched", match_id = NA_integer_,
      stringsAsFactors = FALSE)
  }
  result$raman <- do.call(rbind, raman_parts)

  result
}

# ---------------------------------------------------------------------------
# Coordinate bounds
# ---------------------------------------------------------------------------
compute_bounds <- function(...) {
  dfs <- list(...)
  all_x <- unlist(lapply(dfs, function(d) d$x))
  all_y <- unlist(lapply(dfs, function(d) d$y))
  pad <- 200
  list(
    x = c(min(all_x, na.rm = TRUE) - pad, max(all_x, na.rm = TRUE) + pad),
    y = c(min(all_y, na.rm = TRUE) - pad, max(all_y, na.rm = TRUE) + pad)
  )
}

# ---------------------------------------------------------------------------
# Load a PNG image as a raster for ggplot annotation
# ---------------------------------------------------------------------------
load_image_raster <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  png::readPNG(path)
}

# ---------------------------------------------------------------------------
# Find nearest particle to a hover coordinate
# ---------------------------------------------------------------------------
nearest_particle <- function(df, hover_x, hover_y, max_dist_frac = 0.02) {
  if (is.null(df) || nrow(df) == 0 || is.null(hover_x) || is.null(hover_y)) return(NULL)
  dx <- df$x - hover_x
  dy <- df$y - hover_y
  dists <- sqrt(dx^2 + dy^2)
  idx <- which.min(dists)
  # Threshold: max_dist_frac of the plot range
  x_range <- diff(range(df$x, na.rm = TRUE))
  y_range <- diff(range(df$y, na.rm = TRUE))
  threshold <- max(x_range, y_range) * max_dist_frac
  if (dists[idx] > threshold) return(NULL)
  df[idx, , drop = FALSE]
}
