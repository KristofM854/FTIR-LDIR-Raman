# =============================================================================
# utils_python.R — Reticulate bridge for Python particle detector
# =============================================================================

# Whether particle_detector.py has been loaded in this session.
.python_state <- new.env(parent = emptyenv())
.python_state$loaded <- FALSE

#' Python packages the LDIR image backend needs (import name = pip name).
#' Pinned versions: inst/python/requirements.txt.
.python_required <- c(numpy = "numpy", scipy = "scipy", PIL = "Pillow")

#' Bind reticulate to a real Python interpreter.
#'
#' Called once, lazily, before the first Python use -- not at source() time.
#' Left to itself reticulate provisions an ephemeral 'r-reticulate' virtualenv,
#' which on locked-down/offline machines tries to download `uv` and fails.
#'
#' Interpreter priority:
#'   1. RETICULATE_PYTHON environment variable (set it in ~/.Renviron if your
#'      interpreter is not found automatically).
#'   2. A standard Windows install that exists on this machine.
#'   3. python3 / python on PATH.
#'
#' @return The interpreter path, or "" when none was found.
configure_python <- function() {
  Sys.setenv(RETICULATE_USE_UV = "0")   # never auto-install uv
  py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (!nzchar(py)) {
    cands <- c(sprintf("C:/Program Files/Python3%d/python.exe", 14:11),
               unname(Sys.which("python3")), unname(Sys.which("python")))
    cands <- cands[nzchar(cands) & file.exists(cands)]
    if (length(cands) > 0) py <- cands[[1]]
  }
  if (nzchar(py)) {
    Sys.setenv(RETICULATE_PYTHON = py)
    try(reticulate::use_python(py, required = FALSE), silent = TRUE)
  }
  py
}

#' Load the Python particle detector, or stop with what is missing.
#'
#' LDIR image analysis requires Python 3 with numpy, scipy and Pillow (see
#' inst/python/requirements.txt). There is no R fallback and nothing is
#' installed automatically: a run must not change its own environment.
#' Idempotent.
#'
#' @return Invisible TRUE (errors otherwise).
setup_python_detector <- function() {
  if (isTRUE(.python_state$loaded)) return(invisible(TRUE))
  fail <- function(...) stop(
    "LDIR image analysis requires Python 3 with numpy, scipy and Pillow ",
    "(pip install -r inst/python/requirements.txt). ", ..., call. = FALSE)

  if (!requireNamespace("reticulate", quietly = TRUE))
    fail("The R package 'reticulate' is not installed.")
  py <- configure_python()
  if (!reticulate::py_available(initialize = TRUE))
    fail("No Python interpreter found", if (nzchar(py)) paste0(" (tried ", py, ")"),
         ". Set RETICULATE_PYTHON to your python executable.")

  missing <- names(.python_required)[!vapply(names(.python_required),
                                             reticulate::py_module_available, logical(1))]
  if (length(missing) > 0)
    fail("Missing Python packages: ", paste(.python_required[missing], collapse = ", "),
         " (interpreter: ", reticulate::py_config()$python, ").")

  # Repository root, shiny_app/ or tests/testthat/ as working directory.
  cands <- file.path(c(".", "..", file.path("..", "..")), "inst", "python",
                     "particle_detector.py")
  detector_path <- cands[file.exists(cands)][1]
  if (is.na(detector_path))
    fail("inst/python/particle_detector.py not found from ", getwd(), ".")

  reticulate::source_python(detector_path, envir = globalenv())
  .python_state$loaded <- TRUE
  log_message("  Python particle detector loaded (", reticulate::py_config()$python, ")")
  invisible(TRUE)
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
                                     circle_r = -1.0,
                                     overshoot_factor = 1.0,
                                     merge_fibers = TRUE,
                                     closing_radius = 2L) {

  setup_python_detector()

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
      circle_r  = as.double(circle_r),
      overshoot_factor = as.double(overshoot_factor),
      merge_fibers = as.logical(merge_fibers),
      closing_radius = as.integer(closing_radius)
    )
  }, error = function(e) {
    log_message("  [ERROR] Python detection failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(result)) return(NULL)

  n <- as.integer(result$n_particles)
  n_merges <- as.integer(result$n_fiber_merges %||% 0L)
  log_message("  Python detector: ", n, " particles (threshold=",
              round(result$threshold_used, 1),
              if (n_merges > 0) paste0(", fiber_merges=", n_merges) else "",
              ")")

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
  setup_python_detector()
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
