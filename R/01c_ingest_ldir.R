# =============================================================================
# 01c_ingest_ldir.R — LDIR data ingestion + image-based coordinate extraction
# =============================================================================
#
# The Agilent 8700 LDIR exports particle data (size, shape, material) in an
# Excel file with two sheets:
#   - "Info"      — sample name, instrument serial, etc.
#   - "Particles" — one row per particle with morphological + spectral data
#
# Critically, LDIR does NOT export particle X/Y centroid coordinates.
# This module:
#   1. Ingests the tabular data (particle sizes, materials, quality)
#   2. Extracts coordinates from the companion LDIR PNG image
#   3. Joins image-derived coordinates with Excel particles by size matching
# =============================================================================


#' Ingest LDIR particle data
#'
#' Reads the Agilent 8700 LDIR export and standardizes columns.
#'
#' @param filepath Path to the LDIR Excel file (.xlsx)
#' @param sheet    Sheet name to read (default "Particles")
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um (NA), area_um2, major_um, minor_um,
#'   feret_min_um, feret_max_um, material, quality, source_file
ingest_ldir <- function(filepath, sheet = "Particles") {
  log_message("Reading LDIR data from: ", filepath)

  raw <- readxl::read_excel(filepath, sheet = sheet)
  log_message("  Raw LDIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns: ", paste(names(raw), collapse = ", "))

  # Try to read sample info from "Info" sheet
  info_sheet <- tryCatch({
    readxl::read_excel(filepath, sheet = "Info")
  }, error = function(e) NULL)

  if (!is.null(info_sheet)) {
    sample_name <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Sample Name"]
    }, error = function(e) NA_character_)
    instrument <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Instrument Name"]
    }, error = function(e) NA_character_)
    log_message("  Sample: ", sample_name)
    log_message("  Instrument: ", instrument)
  }

  # Find columns by fuzzy matching
  id_col     <- find_column(raw, c("Id", "Identifier", "#"))
  width_col  <- find_column(raw, c("Width"))
  height_col <- find_column(raw, c("Height"))
  diam_col   <- find_column(raw, c("Diameter"))
  area_col   <- find_column(raw, c("Area"))
  perim_col  <- find_column(raw, c("Perimeter"))
  mat_col    <- find_column(raw, c("Identification", "Material"))
  qual_col   <- find_column(raw, c("Quality"))
  valid_col  <- find_column(raw, c("Is Valid"))
  aspect_col <- find_column(raw, c("Aspect Ratio"))

  # Particle IDs
  particle_ids <- if (!is.null(id_col)) {
    as.character(raw[[id_col]])
  } else {
    paste0("LDIR_", seq_len(nrow(raw)))
  }

  # Sizes
  width_um  <- safe_col_numeric(raw, width_col)
  height_um <- safe_col_numeric(raw, height_col)
  diam_um   <- safe_col_numeric(raw, diam_col)

  # For feret max/min: use width and height as proxies
  feret_max <- pmax(width_um, height_um, na.rm = TRUE)
  feret_min <- pmin(width_um, height_um, na.rm = TRUE)

  # Material
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Quality score
  quality <- safe_col_numeric(raw, qual_col)

  # Is Valid filter
  if (!is.null(valid_col)) {
    valid <- raw[[valid_col]]
    n_invalid <- sum(tolower(valid) != "true", na.rm = TRUE)
    if (n_invalid > 0) {
      log_message("  Flagged ", n_invalid, " invalid particles (keeping all for now)")
    }
  }

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = NA_real_,    # LDIR does not export coordinates
    y_um         = NA_real_,    # will be populated by image extraction
    area_um2     = safe_col_numeric(raw, area_col),
    major_um     = feret_max,
    minor_um     = feret_min,
    feret_min_um = feret_min,
    feret_max_um = feret_max,
    material     = material,
    quality      = quality,
    source_file  = basename(filepath),
    # Keep LDIR-specific columns
    diameter_um  = diam_um,
    aspect_ratio = safe_col_numeric(raw, aspect_col),
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " LDIR particles")
  log_message("  Coordinates: NOT available (LDIR does not export X/Y)")
  log_message("  Size range (feret max): ",
              round(min(df$feret_max_um, na.rm = TRUE), 1), " – ",
              round(max(df$feret_max_um, na.rm = TRUE), 1), " µm")
  log_message("  Materials: ", paste(names(head(sort(table(df$material), decreasing = TRUE), 10)),
                                     collapse = ", "))

  df
}


#' Detect the scan circle inside an LDIR PNG image
#'
#' Converts the image to a binary mask of the scan area, then fits a circle
#' using edge detection and least-squares circle fitting.
#'
#' @param image_path Path to LDIR PNG image file
#' @return List with cx_px, cy_px, radius_px, width, height, edge_gap_px,
#'   export_type ("scan_only" or "full_field")
detect_ldir_scan_circle <- function(image_path) {
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' required. Install with: install.packages('png')")
  }

  img <- png::readPNG(image_path)
  h <- nrow(img)
  w <- ncol(img)
  n_ch <- if (length(dim(img)) == 3) dim(img)[3] else 1

  # Convert to a foreground mask: any pixel that is not near-black
  if (n_ch >= 3) {
    brightness <- pmax(img[,,1], img[,,2], img[,,3])
  } else {
    brightness <- if (n_ch == 1) img else img[,,1]
  }

  # Threshold: pixels brighter than a low percentile are "scan area"
  thresh <- max(0.05, quantile(brightness, 0.10))
  mask <- brightness > thresh

  # Find bounding edge pixels of the mask (contour approximation)
  # For each row, find leftmost and rightmost foreground pixel
  # For each column, find topmost and bottommost foreground pixel
  edge_points_row <- integer(0)
  edge_points_col <- integer(0)

  # Sample rows and columns for speed
  row_sample <- seq(1, h, by = max(1, h %/% 200))
  for (r in row_sample) {
    fg_cols <- which(mask[r, ])
    if (length(fg_cols) >= 2) {
      edge_points_row <- c(edge_points_row, r, r)
      edge_points_col <- c(edge_points_col, min(fg_cols), max(fg_cols))
    }
  }
  col_sample <- seq(1, w, by = max(1, w %/% 200))
  for (cc in col_sample) {
    fg_rows <- which(mask[, cc])
    if (length(fg_rows) >= 2) {
      edge_points_row <- c(edge_points_row, min(fg_rows), max(fg_rows))
      edge_points_col <- c(edge_points_col, cc, cc)
    }
  }

  if (length(edge_points_row) < 10) {
    log_message("  Scan circle detection: too few edge points, falling back to image bounds",
                level = "WARN")
    return(list(
      cx_px = w / 2, cy_px = h / 2, radius_px = min(w, h) / 2,
      width = w, height = h,
      edge_gap_px = 0, export_type = "scan_only", detected = FALSE
    ))
  }

  # Least-squares circle fit (algebraic method)
  # Minimize sum((x-cx)^2 + (y-cy)^2 - r^2)^2
  # Linearized: x^2 + y^2 = 2*cx*x + 2*cy*y + (r^2 - cx^2 - cy^2)
  x_e <- as.numeric(edge_points_col)
  y_e <- as.numeric(edge_points_row)
  A <- cbind(x_e, y_e, 1)
  b_vec <- x_e^2 + y_e^2
  fit <- tryCatch(qr.solve(A, b_vec), error = function(e) NULL)

  if (is.null(fit)) {
    log_message("  Scan circle fit failed, falling back to image bounds", level = "WARN")
    return(list(
      cx_px = w / 2, cy_px = h / 2, radius_px = min(w, h) / 2,
      width = w, height = h,
      edge_gap_px = 0, export_type = "scan_only", detected = FALSE
    ))
  }

  cx_px <- fit[1] / 2
  cy_px <- fit[2] / 2
  radius_px <- sqrt(fit[3] + cx_px^2 + cy_px^2)

  # Classify export type
  edge_gap_px <- min(cx_px, cy_px, w - cx_px, h - cy_px) - radius_px
  export_type <- if (abs(edge_gap_px) <= 15) "scan_only" else "full_field"

  log_message("  Scan circle: center=(", round(cx_px, 1), ", ", round(cy_px, 1),
              "), radius=", round(radius_px, 1), " px")
  log_message("  Edge gap: ", round(edge_gap_px, 1), " px -> ", export_type)

  list(
    cx_px = cx_px, cy_px = cy_px, radius_px = radius_px,
    width = w, height = h,
    edge_gap_px = edge_gap_px, export_type = export_type, detected = TRUE
  )
}


#' Save scan circle debug diagnostic image
#'
#' Overlays detected circle + center crosshair on the original LDIR PNG.
#'
#' @param image_path Path to original LDIR PNG
#' @param circle_info Result from detect_ldir_scan_circle()
#' @param output_path Path to write the debug PNG
save_ldir_circle_debug <- function(image_path, circle_info, output_path) {
  tryCatch({
    img <- png::readPNG(image_path)
    h <- nrow(img)
    w <- ncol(img)

    cx <- circle_info$cx_px
    cy <- circle_info$cy_px
    r  <- circle_info$radius_px

    grDevices::png(output_path, width = w, height = h)
    par(mar = c(0, 0, 0, 0))

    # Plot image
    if (length(dim(img)) == 3) {
      plot(1, type = "n", xlim = c(1, w), ylim = c(h, 1),
           xlab = "", ylab = "", asp = 1, axes = FALSE)
      graphics::rasterImage(img, 1, h, w, 1)
    } else {
      plot(1, type = "n", xlim = c(1, w), ylim = c(h, 1),
           xlab = "", ylab = "", asp = 1, axes = FALSE)
    }

    # Draw detected circle
    theta <- seq(0, 2 * pi, length.out = 360)
    lines(cx + r * cos(theta), cy + r * sin(theta), col = "red", lwd = 2)

    # Draw center crosshair
    lines(c(cx - 30, cx + 30), c(cy, cy), col = "red", lwd = 2)
    lines(c(cx, cx), c(cy - 30, cy + 30), col = "red", lwd = 2)

    # Annotate
    text(cx, cy + 50, paste0("r=", round(r, 0), "px, gap=",
                              round(circle_info$edge_gap_px, 0), "px"),
         col = "yellow", cex = 1.5)

    grDevices::dev.off()
    log_message("  Saved scan circle debug: ", output_path)
  }, error = function(e) {
    log_message("  Could not save circle debug image: ", e$message, level = "WARN")
  })
}


#' Write scan circle export type info to a text file
#'
#' @param circle_info Result from detect_ldir_scan_circle()
#' @param output_path Path to write the text file
save_ldir_export_type <- function(circle_info, output_path) {
  lines <- c(
    paste0("image_width:  ", circle_info$width),
    paste0("image_height: ", circle_info$height),
    paste0("center_px:    (", round(circle_info$cx_px, 2), ", ",
           round(circle_info$cy_px, 2), ")"),
    paste0("radius_px:    ", round(circle_info$radius_px, 2)),
    paste0("edge_gap_px:  ", round(circle_info$edge_gap_px, 2)),
    paste0("export_type:  ", circle_info$export_type),
    paste0("detected:     ", circle_info$detected)
  )
  writeLines(lines, output_path)
}


#' Map particle pixel centroids to µm using scan-circle calibration
#'
#' Uses the detected scan circle center and radius to compute a
#' scale factor, rather than assuming full-image bounds = scan area.
#'
#' @param cx_particle_px Numeric vector, particle centroid x in pixels
#' @param cy_particle_px Numeric vector, particle centroid y in pixels
#' @param circle_info Result from detect_ldir_scan_circle()
#' @param scan_diameter_um Physical scan diameter in µm (default 13000)
#' @return Data frame with x_um, y_um columns
map_pixels_to_um_circle <- function(cx_particle_px, cy_particle_px,
                                     circle_info, scan_diameter_um = 13000) {
  scale_um_per_px <- (scan_diameter_um / 2) / circle_info$radius_px
  x_um <- (cx_particle_px - circle_info$cx_px) * scale_um_per_px
  # Pixel y increases downward; physical y increases upward
  y_um <- (circle_info$cy_px - cy_particle_px) * scale_um_per_px

  data.frame(x_um = x_um, y_um = y_um)
}


#' Extract LDIR particle coordinates from the companion PNG image
#'
#' The Agilent 8700 LDIR exports a particle map PNG where particles are
#' rendered as colored markers on a near-black background. This function
#' detects those markers and extracts centroids.
#'
#' Strategy:
#'   1. Primary: Python backend (scipy/numpy) with iterative sigma-clipping
#'      background correction per tile + global thresholding. Auto-tunes
#'      threshold to match expected_count. Achieves ~505 detections matching
#'      the LDIR software's particle count.
#'   2. Fallback (if Python unavailable): R-based saturation segmentation
#'      for color images, adaptive thresholding for grayscale.
#'
#' After extraction, pixel centroids are remapped to µm using scan-circle
#' calibration (not full-image bounds) for robustness to margins/padding.
#'
#' @param image_path Path to LDIR PNG image file
#' @param scan_bounds Physical scan bounds in µm (list with x_min, x_max, y_min, y_max)
#' @param expected_count Expected number of particles (from Excel data)
#' @param config Optional config list (for debug output and scan diameter)
#' @return Data frame with particle_id, x_um, y_um, area_um2, etc.
extract_ldir_image_coords <- function(image_path,
                                      scan_bounds = NULL,
                                      expected_count = NULL,
                                      config = NULL) {
  log_message("Extracting LDIR particle coordinates from image")

  # --- Step 1: Detect scan circle for calibrated mapping ---
  circle_info <- detect_ldir_scan_circle(image_path)

  scan_diam_um <- if (!is.null(config$ldir_scan_diameter_um)) {
    config$ldir_scan_diameter_um
  } else if (!is.null(scan_bounds)) {
    scan_bounds$x_max - scan_bounds$x_min
  } else {
    13000
  }

  # Save debug artifacts if debug mode enabled
  if (!is.null(config) && isTRUE(config$debug) && !is.null(config$debug_dir)) {
    save_ldir_circle_debug(image_path, circle_info,
                           file.path(config$debug_dir, "ldir_circle_debug.png"))
    save_ldir_export_type(circle_info,
                          file.path(config$debug_dir, "export_type.txt"))
  }

  # --- Step 2: Extract particle pixel centroids ---
  # Use a wrapper that returns particles with centroid_px_x, centroid_px_y columns
  # (pixel-space centroids before µm mapping)
  pixel_particles <- .extract_ldir_pixel_centroids(
    image_path, scan_bounds, expected_count
  )

  if (is.null(pixel_particles) || nrow(pixel_particles) == 0) {
    log_message("  No particles extracted from LDIR image")
    return(.empty_image_df())
  }

  # --- Step 3: Remap pixel centroids → µm using scan-circle calibration ---
  um_coords <- map_pixels_to_um_circle(
    pixel_particles$centroid_px_x,
    pixel_particles$centroid_px_y,
    circle_info, scan_diam_um
  )

  scale_um_per_px <- (scan_diam_um / 2) / circle_info$radius_px

  pixel_particles$x_um <- um_coords$x_um
  pixel_particles$y_um <- um_coords$y_um
  pixel_particles$coord_source <- "circle_calibrated"

  # Store calibration metadata as attributes
  attr(pixel_particles, "circle_cx_px") <- circle_info$cx_px
  attr(pixel_particles, "circle_cy_px") <- circle_info$cy_px
  attr(pixel_particles, "circle_radius_px") <- circle_info$radius_px
  attr(pixel_particles, "scale_um_per_px") <- scale_um_per_px

  log_message("  Circle-calibrated mapping: scale=",
              round(scale_um_per_px, 3), " µm/px, ",
              nrow(pixel_particles), " particles")

  pixel_particles
}


#' Internal: extract particle pixel centroids from LDIR image
#'
#' Returns particles with centroid_px_x, centroid_px_y in pixel coordinates
#' (not yet mapped to µm). Also includes size columns in µm using the
#' legacy scan_bounds for backward compatibility with size-based filtering.
#'
#' @param image_path Path to LDIR PNG image file
#' @param scan_bounds Physical scan bounds in µm
#' @param expected_count Expected number of particles
#' @return Data frame with centroid_px_x, centroid_px_y, and standard columns
.extract_ldir_pixel_centroids <- function(image_path, scan_bounds, expected_count) {
  # --- Try Python backend first (preferred) ---
  py_result <- tryCatch({
    detect_particles_python(
      image_path     = image_path,
      scan_bounds    = scan_bounds,
      expected_count = expected_count
    )
  }, error = function(e) {
    log_message("  Python detector error: ", conditionMessage(e))
    NULL
  })

  if (!is.null(py_result) && nrow(py_result) > 0) {
    log_message("  Extracted ", nrow(py_result), " particles from LDIR image (Python)")
    # Python result has centroid_x/centroid_y in pixel space (from result$particles)
    # The x_um/y_um were computed with old scan_bounds mapping — we keep pixel coords
    if ("centroid_x" %in% names(py_result)) {
      py_result$centroid_px_x <- py_result$centroid_x
      py_result$centroid_px_y <- py_result$centroid_y
    } else {
      # Reverse-compute pixel coords from x_um/y_um if centroid_x not available
      img_w <- attr(py_result, "image_width")
      img_h <- attr(py_result, "image_height")
      if (!is.null(scan_bounds) && !is.null(img_w)) {
        x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / img_w
        y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / img_h
        py_result$centroid_px_x <- (py_result$x_um - scan_bounds$x_min) / x_scale
        py_result$centroid_px_y <- (scan_bounds$y_max - py_result$y_um) / y_scale
      } else {
        py_result$centroid_px_x <- py_result$x_um
        py_result$centroid_px_y <- py_result$y_um
      }
    }
    return(py_result)
  }

  # --- Fallback: R-based extraction ---
  log_message("  Falling back to R-based extraction")

  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' required. Install with: install.packages('png')")
  }

  img <- png::readPNG(image_path)
  h_full <- nrow(img)
  w_full <- ncol(img)
  n_ch <- if (length(dim(img)) == 3) dim(img)[3] else 1
  log_message("  LDIR image: ", w_full, " x ", h_full, " px, ", n_ch, " channels")

  if (n_ch >= 3) {
    r_ch <- img[,,1]; g_ch <- img[,,2]; b_ch <- img[,,3]

    mx <- pmax(r_ch, g_ch, b_ch)
    mn <- pmin(r_ch, g_ch, b_ch)
    sat <- ifelse(mx > 0, (mx - mn) / mx, 0)
    mean_sat <- mean(sat)

    log_message("  Mean saturation: ", round(mean_sat, 3))

    if (mean_sat > 0.15) {
      log_message("  Using saturation-based extraction (colored LDIR image)")
      result <- .extract_ldir_saturation_px(img, h_full, w_full, scan_bounds,
                                             expected_count)
    } else {
      log_message("  Using adaptive-threshold extraction (grayscale image)")
      result <- .extract_ldir_adaptive_px(image_path, h_full, w_full,
                                           scan_bounds, expected_count)
    }
  } else {
    result <- .extract_ldir_adaptive_px(image_path, h_full, w_full,
                                         scan_bounds, expected_count)
  }

  result
}


#' Saturation-based extraction returning pixel centroids
#'
#' Same algorithm as .extract_ldir_saturation but returns centroid_px_x/y
#' in addition to x_um/y_um (which are computed for size filtering only).
.extract_ldir_saturation_px <- function(img, h, w, scan_bounds, expected_count) {
  r_ch <- img[,,1]; g_ch <- img[,,2]; b_ch <- img[,,3]

  mx <- pmax(r_ch, g_ch, b_ch)
  mn <- pmin(r_ch, g_ch, b_ch)
  sat <- ifelse(mx > 0, (mx - mn) / mx, 0)

  binary <- sat > 0.3 & mx > 0.08

  n_fg <- sum(binary)
  log_message("  Saturation threshold: ", n_fg, " foreground pixels (",
              round(n_fg / (h * w) * 100, 1), "%)")

  cc <- .two_pass_ccl(binary, h, w)
  lab_mat <- cc$labels
  n_components <- cc$n_components

  if (n_components == 0) return(.empty_image_df())

  min_pixels <- 15
  fg_idx    <- which(binary, arr.ind = TRUE)
  fg_labels <- lab_mat[binary]
  tab       <- tabulate(fg_labels, nbins = n_components)
  keep_ids  <- which(tab >= min_pixels)

  if (length(keep_ids) == 0) return(.empty_image_df())

  label_factor <- factor(fg_labels, levels = keep_ids)
  valid <- !is.na(label_factor)
  fg_rows <- fg_idx[valid, 1]
  fg_cols <- fg_idx[valid, 2]
  lf      <- droplevels(label_factor[valid])

  cy <- tapply(fg_rows, lf, mean)
  cx <- tapply(fg_cols, lf, mean)
  row_min <- tapply(fg_rows, lf, min)
  row_max <- tapply(fg_rows, lf, max)
  col_min <- tapply(fg_cols, lf, min)
  col_max <- tapply(fg_cols, lf, max)
  bb_h <- row_max - row_min + 1
  bb_w <- col_max - col_min + 1
  areas <- tab[keep_ids]

  # Pixel centroids (0-indexed convention: subtract 0.5)
  cx_full <- as.numeric(cx) - 0.5
  cy_full <- as.numeric(cy) - 0.5

  # Compute approximate µm sizes for filtering using scan_bounds
  if (!is.null(scan_bounds)) {
    x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / w
    y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / h
    size_scale <- max(x_scale, y_scale)
    area_scale <- x_scale * y_scale
  } else {
    size_scale <- 1
    area_scale <- 1
  }

  n_kept <- length(keep_ids)

  df <- data.frame(
    particle_id    = paste0("LDIR_IMG_", seq_len(n_kept)),
    centroid_px_x  = cx_full,
    centroid_px_y  = cy_full,
    x_um           = NA_real_,   # will be overwritten by circle calibration
    y_um           = NA_real_,
    area_um2       = as.numeric(areas) * area_scale,
    major_um       = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    minor_um       = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_min_um   = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_max_um   = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    material       = NA_character_,
    quality        = NA_real_,
    source_file    = "LDIR_image",
    stringsAsFactors = FALSE
  )

  # Size filter
  df <- df[df$feret_max_um >= 20, ]

  # Auto-trim if too many
  if (!is.null(expected_count) && expected_count > 0 && nrow(df) > expected_count * 1.5) {
    target_n <- round(expected_count * 1.3)
    if (nrow(df) > target_n) {
      size_cutoff <- sort(df$feret_max_um, decreasing = TRUE)[min(target_n, nrow(df))]
      df <- df[df$feret_max_um >= size_cutoff, ]
    }
  }

  log_message("  Final: ", nrow(df), " particles from saturation segmentation")
  df
}


#' Adaptive-threshold extraction returning pixel centroids (fallback)
.extract_ldir_adaptive_px <- function(image_path, h_full, w_full,
                                       scan_bounds, expected_count) {
  # Delegate to existing extract_particles_from_image, then recover pixel coords
  result <- extract_particles_from_image(
    image_path,
    scan_bounds     = scan_bounds,
    adaptive_radius = 50L,
    adaptive_offset = 0.05,
    min_pixels      = 15,
    min_size_um     = 20,
    downsample      = 2L,
    expected_count  = expected_count,
    instrument      = "LDIR"
  )

  if (nrow(result) > 0 && !is.null(scan_bounds)) {
    # Reverse the µm → pixel mapping to recover pixel centroids
    x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / w_full
    y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / h_full
    result$centroid_px_x <- (result$x_um - scan_bounds$x_min) / x_scale
    result$centroid_px_y <- (scan_bounds$y_max - result$y_um) / y_scale
  } else if (nrow(result) > 0) {
    result$centroid_px_x <- result$x_um
    result$centroid_px_y <- result$y_um
  }

  result
}


#' Join image-extracted coordinates with Excel particle table
#'
#' Since the LDIR Excel has IDs + sizes but no coordinates, and the image
#' has coordinates + sizes but no IDs, we use size-based optimal matching
#' (Hungarian algorithm) to pair them.
#'
#' @param excel_df Data frame from ingest_ldir() (has particle_id, sizes, material)
#' @param image_df Data frame from extract_ldir_image_coords() (has x_um, y_um, sizes)
#' @return excel_df with x_um and y_um populated, plus coord_match_cost column
join_ldir_coords <- function(excel_df, image_df) {
  log_message("Joining LDIR image coordinates with Excel particle table")

  n_excel <- nrow(excel_df)
  n_image <- nrow(image_df)
  log_message("  Excel: ", n_excel, " particles, Image: ", n_image, " particles")

  if (n_image == 0) {
    log_message("  No image particles to join — coordinates remain NA", level = "WARN")
    excel_df$coord_match_cost <- NA_real_
    excel_df$coord_source <- "none"
    return(excel_df)
  }

  # Build cost matrix using morphological similarity (vectorized)
  n_max <- max(n_excel, n_image)

  excel_area  <- excel_df$area_um2
  image_area  <- image_df$area_um2
  excel_feret <- excel_df$feret_max_um
  image_feret <- image_df$feret_max_um

  BIG <- 1e9
  cost <- matrix(BIG, nrow = n_max, ncol = n_max)

  # Vectorized log-ratio cost: outer() gives the full n_excel × n_image matrix
  # without any R-level loops.
  log_area <- outer(
    log(pmax(excel_area, 1e-6, na.rm = FALSE)),
    log(pmax(image_area,  1e-6, na.rm = FALSE)),
    FUN = function(a, b) abs(a - b)
  )
  log_feret <- outer(
    log(pmax(excel_feret, 1e-6, na.rm = FALSE)),
    log(pmax(image_feret,  1e-6, na.rm = FALSE)),
    FUN = function(a, b) abs(a - b)
  )

  # NA → 0 (missing size contributes no cost so we don't penalize)
  log_area[is.na(log_area)]   <- 0
  log_feret[is.na(log_feret)] <- 0

  cost[seq_len(n_excel), seq_len(n_image)] <- log_area + 0.5 * log_feret

  # Hungarian assignment
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- as.integer(clue::solve_LSAP(cost, maximum = FALSE))
  } else {
    # Greedy fallback: for each excel row, pick best remaining image
    assignment <- rep(NA_integer_, n_excel)
    taken <- logical(n_image)
    for (i in order(apply(cost[seq_len(n_excel), seq_len(n_image), drop = FALSE], 1, min))) {
      best_j <- which.min(cost[i, seq_len(n_image)] + ifelse(taken, BIG, 0))
      if (cost[i, best_j] < BIG && !taken[best_j]) {
        assignment[i] <- best_j
        taken[best_j] <- TRUE
      }
    }
  }

  # Apply coordinates (reject matches above cost threshold)
  max_cost <- 2.0
  excel_df$coord_match_cost <- NA_real_
  excel_df$coord_source <- "none"
  n_joined <- 0
  n_rejected <- 0

  for (i in seq_len(n_excel)) {
    j <- assignment[i]
    if (is.na(j) || j > n_image || cost[i, j] >= BIG) next

    if (cost[i, j] > max_cost) {
      n_rejected <- n_rejected + 1
      next
    }

    excel_df$x_um[i] <- image_df$x_um[j]
    excel_df$y_um[i] <- image_df$y_um[j]
    excel_df$coord_match_cost[i] <- cost[i, j]
    # Propagate coord_source from image extraction (e.g. "circle_calibrated")
    excel_df$coord_source[i] <- if ("coord_source" %in% names(image_df) &&
                                     !is.na(image_df$coord_source[j])) {
      image_df$coord_source[j]
    } else {
      "image"
    }
    n_joined <- n_joined + 1
  }

  n_missing <- n_excel - n_joined
  log_message("  Joined coordinates: ", n_joined, " of ", n_excel, " particles")
  if (n_rejected > 0) {
    log_message("  Rejected ", n_rejected, " joins with cost > ", max_cost,
                " (poor size match)")
  }
  if (n_missing > 0) {
    log_message("  Missing coordinates: ", n_missing, " particles (no image match)")
  }

  # Quality diagnostics
  costs <- excel_df$coord_match_cost[!is.na(excel_df$coord_match_cost)]
  if (length(costs) > 0) {
    log_message("  Join quality — cost median: ", round(median(costs), 3),
                ", mean: ", round(mean(costs), 3),
                ", max: ", round(max(costs), 3))
  }

  excel_df
}


#' Validate LDIR coordinate join by checking scan-order correlation
#'
#' If LDIR IDs follow a raster scan order, the ID sequence should correlate
#' with the spatial position of joined coordinates.
#'
#' @param df Data frame with particle_id and x_um, y_um (joined)
#' @return List with tau correlation coefficient and assessment
validate_ldir_scan_order <- function(df) {
  # Only use particles with valid coordinates
  valid <- !is.na(df$x_um) & !is.na(df$y_um)
  if (sum(valid) < 10) {
    return(list(tau = NA_real_, scan_order_consistent = NA,
                message = "Too few coordinated particles to test"))
  }

  df_valid <- df[valid, ]

  # Extract numeric part of ID for ordering
  id_num <- as.numeric(gsub("[^0-9]", "", df_valid$particle_id))
  if (all(is.na(id_num))) {
    return(list(tau = NA_real_, scan_order_consistent = NA,
                message = "Non-numeric IDs — cannot test scan order"))
  }

  # Test raster-scan correlation: sort by (y descending, x ascending)
  spatial_rank <- rank(-df_valid$y_um * 1e6 + df_valid$x_um)
  id_rank <- rank(id_num, na.last = "keep")

  tau <- cor(spatial_rank, id_rank, method = "kendall", use = "complete.obs")

  list(
    tau = tau,
    scan_order_consistent = !is.na(tau) && abs(tau) > 0.5,
    message = if (!is.na(tau)) {
      paste0("Kendall tau = ", round(tau, 3),
             if (abs(tau) > 0.7) " (strong)" else if (abs(tau) > 0.5) " (moderate)" else " (weak)")
    } else "Could not compute"
  )
}


#' Pre-filter LDIR particle data
#'
#' @param df Data frame from ingest_ldir()
#' @param min_quality Minimum quality score (0 = keep all)
#' @param min_size_um Minimum particle size in µm (applied to feret_max_um)
#' @param remove_invalid Logical, remove "Is Valid" = FALSE particles
#' @return Filtered data frame
prefilter_ldir <- function(df, min_quality = 0, min_size_um = 0,
                           remove_invalid = FALSE) {
  n_start <- nrow(df)
  log_message("Pre-filtering LDIR: ", n_start, " particles")

  # Quality filter
  if (min_quality > 0 && any(!is.na(df$quality))) {
    keep <- is.na(df$quality) | df$quality >= min_quality
    df <- df[keep, ]
    log_message("  Quality filter (>= ", min_quality, "): kept ", nrow(df),
                " of ", n_start, " particles")
  }

  # Size filter
  if (min_size_um > 0 && any(!is.na(df$feret_max_um))) {
    keep <- is.na(df$feret_max_um) | df$feret_max_um >= min_size_um
    df <- df[keep, ]
    log_message("  Size filter (>= ", min_size_um, " µm): kept ", nrow(df), " particles")
  }

  log_message("  LDIR after filtering: ", nrow(df), " particles")
  df
}
