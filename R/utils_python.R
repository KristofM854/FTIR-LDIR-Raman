# =============================================================================
# utils_python.R — Reticulate bridge for Python particle detector
# =============================================================================

.python_detector_loaded <- FALSE

#' Initialize the Python environment for particle detection
#'
#' Sources the particle_detector.py module via reticulate.
#' Safe to call multiple times (idempotent).
#'
#' @return Invisible TRUE if successful, FALSE if not.
setup_python_detector <- function() {
  if (.python_detector_loaded) return(invisible(TRUE))

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    log_message("  [WARN] reticulate package not available — Python detector disabled")
    return(invisible(FALSE))
  }

  # Verify required Python packages
  required <- c("numpy", "scipy", "PIL")
  for (pkg in required) {
    if (!reticulate::py_module_available(pkg)) {
      log_message("  [WARN] Python package '", pkg, "' not found — ",
                  "install with: pip install ",
                  ifelse(pkg == "PIL", "Pillow", pkg))
      return(invisible(FALSE))
    }
  }

  # Find detector script relative to project root
  detector_path <- file.path("inst", "python", "particle_detector.py")
  if (!file.exists(detector_path)) {
    # Try from shiny_app working directory
    detector_path <- file.path("..", "inst", "python", "particle_detector.py")
  }

  if (!file.exists(detector_path)) {
    log_message("  [WARN] particle_detector.py not found")
    return(invisible(FALSE))
  }

  tryCatch({
    reticulate::source_python(detector_path, envir = globalenv())
    .python_detector_loaded <<- TRUE
    log_message("  Python particle detector loaded successfully")
    invisible(TRUE)
  }, error = function(e) {
    log_message("  [WARN] Failed to load Python detector: ", conditionMessage(e))
    invisible(FALSE)
  })
}


#' Run particle detection on an LDIR mosaic image using Python backend
#'
#' Wrapper around the Python run_full_pipeline() function.
#' Converts Python dict-of-lists to a proper R data.frame with
#' physical (µm) coordinates.
#'
#' @param image_path Character. Path to the LDIR mosaic PNG file.
#' @param scan_bounds List with xmin, xmax, ymin, ymax in µm (for coordinate conversion).
#' @param expected_count Integer. Target particle count for auto-tuning threshold.
#'   If NULL or 0, uses default threshold.
#' @param grid_rows,grid_cols Integer. Tile grid dimensions (default 4x4).
#' @param bg_sigma Numeric. Gaussian sigma for background estimation.
#' @param threshold Numeric. Detection threshold (used only if expected_count is NULL).
#' @param min_area Integer. Minimum particle area in pixels.
#'
#' @return A data.frame with columns: particle_id, x_um, y_um, area_px,
#'   feret_max_um, equivalent_diameter_px, aspect_ratio, tile_row, tile_col, etc.
#'   Returns NULL on failure.
detect_particles_python <- function(image_path, scan_bounds = NULL,
                                     expected_count = NULL,
                                     grid_rows = 4L, grid_cols = 4L,
                                     bg_sigma = 30, threshold = 25,
                                     min_area = 20L) {

  if (!setup_python_detector()) return(NULL)

  stopifnot(file.exists(image_path))

  target <- if (!is.null(expected_count) && expected_count > 0) {
    as.integer(expected_count)
  } else {
    0L
  }

  result <- tryCatch({
    run_full_pipeline(
      image_path = image_path,
      grid_rows = as.integer(grid_rows),
      grid_cols = as.integer(grid_cols),
      bg_sigma = as.double(bg_sigma),
      clip_sigma = 3.0,
      max_iter = 10L,
      threshold = as.double(threshold),
      min_area = as.integer(min_area),
      target_count = target
    )
  }, error = function(e) {
    log_message("  [ERROR] Python detection failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(result)) return(NULL)

  n <- as.integer(result$n_particles)
  log_message("  Python detector: ", n, " particles (threshold=",
              round(result$threshold_used, 1), ")")

  if (n == 0) return(NULL)

  # Convert Python dict-of-lists to R data.frame
  particles <- as.data.frame(result$particles, stringsAsFactors = FALSE)

  # Compute physical coordinates (µm)
  img_h <- as.integer(result$image_height)
  img_w <- as.integer(result$image_width)

  if (!is.null(scan_bounds)) {
    x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / img_w
    y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / img_h

    # Image convention: y=0 at top. Convert to Cartesian (y increases up)
    particles$x_um <- scan_bounds$x_min + particles$centroid_x * x_scale
    particles$y_um <- scan_bounds$y_max - particles$centroid_y * y_scale

    # Convert pixel sizes to µm (use average scale)
    um_per_px <- mean(c(x_scale, y_scale))
    particles$feret_max_um <- pmax(particles$bbox_width, particles$bbox_height) * um_per_px
    particles$feret_min_um <- pmin(particles$bbox_width, particles$bbox_height) * um_per_px
    particles$area_um2 <- particles$area_px * um_per_px^2
  } else {
    # No scan bounds — use pixel coordinates
    particles$x_um <- particles$centroid_x
    particles$y_um <- particles$centroid_y
    particles$feret_max_um <- pmax(particles$bbox_width, particles$bbox_height)
    particles$feret_min_um <- pmin(particles$bbox_width, particles$bbox_height)
    particles$area_um2 <- particles$area_px
  }

  attr(particles, "image_height") <- img_h
  attr(particles, "image_width") <- img_w
  attr(particles, "tile_height") <- as.integer(result$tile_height)
  attr(particles, "tile_width") <- as.integer(result$tile_width)
  attr(particles, "threshold_used") <- result$threshold_used

  particles
}
