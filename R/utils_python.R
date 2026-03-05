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

  # Verify required Python packages; auto-install via pip if missing
  # But first: only attempt install if Python is actually available.
  required  <- c("numpy", "scipy", "PIL")
  pip_names <- c(numpy = "numpy", scipy = "scipy", PIL = "Pillow")

  missing_pkgs <- required[!vapply(required, reticulate::py_module_available, logical(1))]
  if (length(missing_pkgs) > 0) {
    # Check if ANY Python is available before attempting install
    if (!reticulate::py_available()) {
      log_message("  [WARN] Python not installed on this system — Python detector disabled")
      log_message("  Falling back to R-based particle extraction")
      return(invisible(FALSE))
    }

    pip_missing <- unname(pip_names[missing_pkgs])
    log_message("  Python packages missing: ", paste(pip_missing, collapse = ", "),
                " — attempting auto-install via pip")
    tryCatch({
      reticulate::py_install(pip_missing, pip = TRUE)
      log_message("  Auto-install succeeded")
    }, error = function(e) {
      log_message("  [WARN] Auto-install failed: ", conditionMessage(e),
                  " — install manually with: pip install ",
                  paste(pip_missing, collapse = " "))
    })
    # Re-check after attempted install
    still_missing <- missing_pkgs[
      !vapply(missing_pkgs, reticulate::py_module_available, logical(1))
    ]
    if (length(still_missing) > 0) {
      log_message("  [WARN] Python packages still unavailable after install attempt: ",
                  paste(unname(pip_names[still_missing]), collapse = ", "),
                  " — Python detector disabled")
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
                                     min_area = 10L,
                                     circle_cx = -1.0, circle_cy = -1.0,
                                     circle_r = -1.0) {

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
      target_count = target,
      circle_cx = as.double(circle_cx),
      circle_cy = as.double(circle_cy),
      circle_r  = as.double(circle_r)
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


#' Detect the LDIR scan circle using the Python connected-component method
#'
#' Calls detect_scan_circle() from particle_detector.py, which is more robust
#' than the R algebraic edge-fit against bright particles near the image boundary.
#'
#' @param image_path Character. Path to the LDIR image file.
#' @return A list compatible with detect_ldir_scan_circle() output:
#'   cx_px, cy_px, radius_px, width, height, edge_gap_px, export_type, detected.
#'   Returns NULL if Python is unavailable.
detect_ldir_scan_circle_python <- function(image_path) {
  if (!setup_python_detector()) return(NULL)
  if (!file.exists(image_path)) return(NULL)

  result <- tryCatch({
    detect_scan_circle(image_path)
  }, error = function(e) {
    log_message("  [WARN] Python circle detection failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(result)) return(NULL)

  # Reject the Python fallback: fallback_center is image-centre + 95% radius,
  # indistinguishable from a real detection but unreliable.
  met <- as.character(result$method)
  if (identical(met, "fallback_center")) {
    log_message("  Python circle detection returned fallback_center — treating as failure",
                level = "WARN")
    return(NULL)
  }

  cx  <- as.numeric(result$cx)
  cy  <- as.numeric(result$cy)
  r   <- as.numeric(result$r)
  w   <- as.integer(result$width)
  h   <- as.integer(result$height)

  # Sanity: circle values must be inside image dimensions and radius positive
  if (!(cx > 0 && cx < w && cy > 0 && cy < h && r > 0)) {
    log_message("  Python circle detection: out-of-bounds result (cx=", round(cx, 1),
                ", cy=", round(cy, 1), ", r=", round(r, 1),
                ", w=", w, ", h=", h, ") — treating as failure",
                level = "WARN")
    return(NULL)
  }

  edge_gap <- min(cx, cy, w - cx, h - cy) - r
  export_type <- if (abs(edge_gap) <= 15) "scan_only" else "full_field"

  log_message("  Python circle detection (", met, "): center=(",
              round(cx, 1), ", ", round(cy, 1), "), radius=", round(r, 1),
              " px, edge_gap=", round(edge_gap, 1))

  if (edge_gap < -100) {
    log_message("  Python circle: edge_gap=", round(edge_gap, 1),
                " px (< -100) — detection likely wrong; results may be unreliable",
                level = "WARN")
  }

  list(
    cx_px       = cx,
    cy_px       = cy,
    radius_px   = r,
    width       = w,
    height      = h,
    edge_gap_px = edge_gap,
    export_type = export_type,
    detected    = TRUE,
    method      = met
  )
}
