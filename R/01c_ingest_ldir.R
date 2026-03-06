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

  # Sanitize column names to valid UTF-8 immediately
  names(raw) <- safe_colnames(names(raw))

  log_message("  Raw LDIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns (sanitized): ", paste(names(raw), collapse = ", "))

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
    valid <- safe_colnames(as.character(raw[[valid_col]]))
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


#' Detect the scan circle inside an LDIR image (any format)
#'
#' Converts the image to a binary mask of the scan area, then fits a circle
#' using edge detection and least-squares circle fitting.
#'
#' Accepts PNG, JPEG, TIFF, BMP, WEBP regardless of file extension.
#' Uses magick for robust multi-format reading (with png::readPNG fallback).
#'
#' @param image_path Path to LDIR image file
#' @return List with cx_px, cy_px, radius_px, width, height, edge_gap_px,
#'   export_type ("scan_only" or "full_field")
detect_ldir_scan_circle <- function(image_path, config = NULL) {
  detected_fmt <- guess_image_type(image_path)
  ext_type     <- toupper(tools::file_ext(image_path))
  if (detected_fmt != "unknown" && detected_fmt != ext_type) {
    log_message("  detect_ldir_scan_circle: extension=", ext_type,
                " but signature=", detected_fmt, " — using magick", level = "WARN")
  }

  img <- read_image_any(image_path, verbose = FALSE)
  if (is.null(img)) {
    stop("Could not read LDIR image: ", image_path,
         "\n  Detected format: ", detected_fmt,
         "\n  Make sure 'magick' package is installed.")
  }

  h <- nrow(img)
  w <- ncol(img)

  # --- Forced mosaic mode: skip circle detection entirely ---
  export_fmt <- if (!is.null(config$ldir_export_format)) config$ldir_export_format else "auto"
  if (identical(export_fmt, "mosaic")) {
    log_message("  Scan circle: ldir_export_format='mosaic' -- using full-image bounds")
    return(list(cx_px = w / 2, cy_px = h / 2, radius_px = min(w, h) / 2,
                width = w, height = h, edge_gap_px = 0,
                export_type = "mosaic_full_field",
                detected = TRUE, method = "mosaic_forced"))
  }

  # Helper: sanity-check a candidate circle and return TRUE if plausible
  .circle_sane <- function(cx, cy, r, w, h) {
    cx_ok <- abs(cx - w / 2) < 0.15 * w
    cy_ok <- abs(cy - h / 2) < 0.15 * h
    r_ok  <- r > 0.38 * min(w, h) && r < 0.54 * min(w, h)
    cx_ok && cy_ok && r_ok
  }

  # Quadrant-similarity check on float [0,1] brightness matrix
  .is_tiled_mosaic <- function(brightness) {
    bh <- nrow(brightness); bw <- ncol(brightness)
    hh <- bh %/% 2L; hw <- bw %/% 2L
    quad <- list(
      brightness[seq_len(hh),       seq_len(hw)],
      brightness[seq_len(hh),       seq(hw + 1L, bw)],
      brightness[seq(hh + 1L, bh),  seq_len(hw)],
      brightness[seq(hh + 1L, bh),  seq(hw + 1L, bw)]
    )
    stats_q <- lapply(quad, function(q)
      c(mean(q, na.rm = TRUE), stats::sd(q, na.rm = TRUE)))
    pairs <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
    n_near <- sum(vapply(pairs, function(ij) {
      abs(stats_q[[ij[1]]][1] - stats_q[[ij[2]]][1]) < 0.03 &&
      abs(stats_q[[ij[1]]][2] - stats_q[[ij[2]]][2]) < 0.03
    }, logical(1L)))
    n_near >= 2L
  }

  # Image-centre fallback used when every method fails sanity check.
  # If the image looks like a tiled mosaic (auto mode), use full-image bounds
  # with detected=TRUE so the pipeline continues.
  .fallback <- function(w, h, brightness = NULL) {
    if (!is.null(brightness) && !identical(export_fmt, "circular") &&
        .is_tiled_mosaic(brightness)) {
      log_message("  Scan circle: no circle found, image appears tiled -- ",
                  "using full-image bounds")
      return(list(cx_px = w / 2, cy_px = h / 2, radius_px = min(w, h) / 2,
                  width = w, height = h, edge_gap_px = 0,
                  export_type = "mosaic_full_field",
                  detected = TRUE, method = "mosaic_auto"))
    }
    r <- min(w, h) / 2 * 0.95
    edge_gap <- min(w, h) / 2 - r
    list(cx_px = w / 2, cy_px = h / 2, radius_px = r,
         width = w, height = h,
         edge_gap_px = edge_gap, export_type = "scan_only", detected = FALSE)
  }

  # --- Primary: Python connected-component method (most robust) ---
  py_circle <- tryCatch(
    detect_ldir_scan_circle_python(image_path),
    error = function(e) NULL
  )
  if (!is.null(py_circle)) {
    if (.circle_sane(py_circle$cx_px, py_circle$cy_px, py_circle$radius_px, w, h)) {
      log_message("  Scan circle (Python CC): center=(",
                  round(py_circle$cx_px, 1), ", ", round(py_circle$cy_px, 1),
                  "), radius=", round(py_circle$radius_px, 1), " px, gap=",
                  round(py_circle$edge_gap_px, 1), " -> ", py_circle$export_type)
      return(py_circle)
    }
    log_message("  Python circle failed sanity check (cx=",
                round(py_circle$cx_px, 1), ", cy=", round(py_circle$cy_px, 1),
                ", r=", round(py_circle$radius_px, 1), ") — trying R fallback",
                level = "WARN")
  }

  # --- Fallback: R algebraic edge-fit ---
  n_ch <- if (length(dim(img)) == 3) dim(img)[3] else 1
  if (n_ch >= 3) {
    brightness <- pmax(img[,,1], img[,,2], img[,,3])
  } else {
    brightness <- if (n_ch == 1) img else img[,,1]
  }

  thresh <- max(0.05, quantile(brightness, 0.10))
  mask <- brightness > thresh

  edge_points_row <- integer(0)
  edge_points_col <- integer(0)
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
    log_message("  Scan circle: too few edge points — using image-centre defaults",
                level = "WARN")
    return(.fallback(w, h, brightness))
  }

  x_e <- as.numeric(edge_points_col)
  y_e <- as.numeric(edge_points_row)
  A     <- cbind(x_e, y_e, 1)
  b_vec <- x_e^2 + y_e^2
  fit   <- tryCatch(qr.solve(A, b_vec), error = function(e) NULL)

  if (is.null(fit)) {
    log_message("  Scan circle: algebraic fit failed — using image-centre defaults",
                level = "WARN")
    return(.fallback(w, h, brightness))
  }

  cx_px    <- fit[1] / 2
  cy_px    <- fit[2] / 2
  radius_px <- sqrt(fit[3] + cx_px^2 + cy_px^2)

  # Sanity check on algebraic result
  if (!.circle_sane(cx_px, cy_px, radius_px, w, h)) {
    log_message("  Scan circle: algebraic result failed sanity check (cx=",
                round(cx_px, 1), ", cy=", round(cy_px, 1),
                ", r=", round(radius_px, 1), ") — using image-centre defaults",
                level = "WARN")
    return(.fallback(w, h, brightness))
  }

  edge_gap_px <- min(cx_px, cy_px, w - cx_px, h - cy_px) - radius_px
  export_type <- if (abs(edge_gap_px) <= 15) "scan_only" else "full_field"

  log_message("  Scan circle (R algebraic): center=(",
              round(cx_px, 1), ", ", round(cy_px, 1),
              "), radius=", round(radius_px, 1), " px, gap=",
              round(edge_gap_px, 1), " -> ", export_type)

  list(
    cx_px = cx_px, cy_px = cy_px, radius_px = radius_px,
    width = w, height = h,
    edge_gap_px = edge_gap_px, export_type = export_type,
    detected = TRUE, method = "r_algebraic"
  )
}


#' Save scan circle debug diagnostic image
#'
#' Overlays detected circle + center crosshair on the original LDIR image.
#' Accepts any image format supported by magick.
#'
#' @param image_path Path to original LDIR image (any format)
#' @param circle_info Result from detect_ldir_scan_circle()
#' @param output_path Path to write the debug PNG
save_ldir_circle_debug <- function(image_path, circle_info, output_path) {
  # Background image source: `image_path` argument (canonical PNG if canonicalization ran,
  # else original LDIR source file).  Rendered via read_image_canonical() -> [H,W,3] uint8 RGB
  # -> as.raster(img, max=255L) -> graphics::rasterImage().
  tryCatch({
    # If output_path is a directory, auto-create a filename
    if (dir.exists(output_path)) {
      output_path <- file.path(output_path, "ldir_circle_detection.png")
    }
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

    canon <- read_image_canonical(image_path, verbose = FALSE)
    if (is.null(canon)) {
      log_message("  save_ldir_circle_debug: could not read image", level = "WARN")
      return(invisible(NULL))
    }
    img <- canon$img_rgb
    h <- canon$height
    w <- canon$width

    assert_not_tiled_montage(img, tag = basename(image_path),
                              strict = isTRUE(getOption("ldir_strict_sanity")))

    log_message("  circle_debug img: dim=", paste(dim(img), collapse = "x"),
                " range=", paste(range(img), collapse = ".."))

    cx <- circle_info$cx_px
    cy <- circle_info$cy_px
    r  <- circle_info$radius_px

    grDevices::png(output_path, width = w, height = h)
    on.exit(grDevices::dev.off(), add = TRUE)
    par(mar = c(0, 0, 0, 0))

    plot(1, type = "n", xlim = c(1, w), ylim = c(h, 1),
         xlab = "", ylab = "", asp = 1, axes = FALSE)
    graphics::rasterImage(as.raster(img, max = 255L), 1, h, w, 1)

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

    log_message("  Saved scan circle debug: ", output_path)
  }, error = function(e) {
    log_message("  Could not save circle debug image: ", e$message, level = "WARN")
  })
}


#' Save pixel-space centroid overlay for extraction quality diagnosis
#'
#' Draws detected centroids (in pixel space) on top of the original LDIR image.
#' If centroids look correct here but µm coords are compressed/stretched,
#' the fault lies in the pixel→µm mapping (wrong radius_px or scan_diam_um).
#'
#' @param image_path   Path to original LDIR image
#' @param pixel_df     Data frame with centroid_px_x, centroid_px_y columns
#' @param circle_info  Result from detect_ldir_scan_circle()
#' @param output_path  Path for the debug PNG
save_ldir_centroids_px_debug <- function(image_path, pixel_df, circle_info, output_path) {
  # Background image source: `image_path` argument (canonical PNG if canonicalization ran,
  # else original LDIR source file).  Rendered via read_image_canonical() -> [H,W,3] uint8 RGB
  # -> as.raster(img, max=255L) -> graphics::rasterImage().
  tryCatch({
    canon <- read_image_canonical(image_path, verbose = FALSE)
    if (is.null(canon)) return(invisible(NULL))
    img <- canon$img_rgb
    h <- canon$height; w <- canon$width
    assert_not_tiled_montage(img, tag = basename(image_path),
                              strict = isTRUE(getOption("ldir_strict_sanity")))
    cx_pts <- pixel_df$centroid_px_x
    cy_pts <- pixel_df$centroid_px_y
    grDevices::png(output_path, width = w, height = h)
    par(mar = c(0, 0, 0, 0))
    plot(1, type = "n", xlim = c(1, w), ylim = c(h, 1), asp = 1, axes = FALSE,
         xlab = "", ylab = "")
    graphics::rasterImage(as.raster(img, max = 255L), 1, h, w, 1)
    if (!is.null(circle_info) && isTRUE(circle_info$detected)) {
      theta <- seq(0, 2 * pi, length.out = 360)
      lines(circle_info$cx_px + circle_info$radius_px * cos(theta),
            circle_info$cy_px + circle_info$radius_px * sin(theta),
            col = "red", lwd = 2)
    }
    points(cx_pts, cy_pts, pch = 3, cex = 0.6, col = "cyan", lwd = 1)
    title(main = paste0(nrow(pixel_df), " centroids (pixel space)"), line = -2,
          col.main = "yellow", cex.main = 1.2)
    grDevices::dev.off()
    log_message("  Saved centroid pixel overlay: ", output_path)
  }, error = function(e) {
    log_message("  Could not save centroid pixel debug: ", e$message, level = "WARN")
  })
}


#' Write LDIR calibration numbers to a text file for diagnosis
#'
#' Shows all the numbers that go into the pixel→µm mapping, plus a
#' fill_fraction diagnostic: if particles reach only ~90 % of the expected
#' half-extent, radius_px or scan_diam_um is likely wrong.
#'
#' @param circle_info     Result from detect_ldir_scan_circle()
#' @param scan_diam_um    Physical scan diameter used (µm)
#' @param scale_um_per_px Derived scale factor (µm/pixel)
#' @param n_particles     Number of particles after mapping
#' @param max_abs_um      Maximum absolute µm coordinate (diagnostic)
#' @param output_path     Path for the calibration text file
save_ldir_calibration <- function(circle_info, scan_diam_um, scale_um_per_px,
                                   n_particles, max_abs_um, output_path) {
  expected_half <- scan_diam_um / 2
  compression_pct <- if (!is.na(max_abs_um) && max_abs_um > 0)
    round((max_abs_um / expected_half) * 100, 1) else NA_real_
  lines <- c(
    paste0("# LDIR pixel-to-um calibration — generated ", Sys.time()),
    "",
    paste0("image_width_px:    ", circle_info$width),
    paste0("image_height_px:   ", circle_info$height),
    paste0("circle_cx_px:      ", round(circle_info$cx_px, 2)),
    paste0("circle_cy_px:      ", round(circle_info$cy_px, 2)),
    paste0("circle_radius_px:  ", round(circle_info$radius_px, 2)),
    paste0("edge_gap_px:       ", round(circle_info$edge_gap_px, 2)),
    paste0("export_type:       ", circle_info$export_type),
    paste0("circle_detected:   ", circle_info$detected),
    "",
    paste0("scan_diameter_um:  ", scan_diam_um),
    paste0("scale_um_per_px:   ", round(scale_um_per_px, 4)),
    "",
    paste0("n_particles:       ", n_particles),
    paste0("max_abs_um:        ", round(max_abs_um, 1)),
    paste0("fill_fraction_pct: ", compression_pct,
           "  (expected ~100 if scale correct, <90 => radius_px too large)")
  )
  writeLines(lines, output_path)
  log_message("  Saved calibration info: ", output_path)
}


#' Single source of truth for the LDIR µm-per-pixel scale factor
#'
#' Priority: (1) explicit config override, (2) pre-computed in circle_info,
#' (3) derived from scan diameter / detected radius.
#'
#' @param circle_info Result from detect_ldir_scan_circle()
#' @param config      Optional pipeline config list
#' @return Numeric: µm per pixel
ldir_um_per_px <- function(circle_info, config = NULL) {
  if (!is.null(config$ldir_um_per_px) && is.numeric(config$ldir_um_per_px) &&
      config$ldir_um_per_px > 0)
    return(config$ldir_um_per_px)
  if (!is.null(circle_info$scale_um_per_px) && circle_info$scale_um_per_px > 0)
    return(circle_info$scale_um_per_px)
  scan_r_um <- (config$ldir_scan_diameter_um %||% 13000) / 2
  scan_r_um / circle_info$radius_px
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


#' Extract LDIR particle coordinates from the companion image
#'
#' The Agilent 8700 LDIR exports a particle map image where particles are
#' rendered as colored markers on a near-black background.  This function
#' detects those markers and extracts centroids.
#'
#' Accepted formats: PNG, JPEG (.jpg/.jpeg), TIFF (.tif/.tiff), BMP, WEBP.
#' Format is detected from magic bytes, not the file extension, so a JPEG
#' file named .png is handled correctly.
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
#' @param image_path Path to LDIR image file (PNG/JPEG/TIFF/BMP/WEBP)
#' @param scan_bounds Physical scan bounds in µm (list with x_min, x_max, y_min, y_max)
#' @param expected_count Expected number of particles (from Excel data)
#' @param config Optional config list (for debug output and scan diameter)
#' @return Data frame with particle_id, x_um, y_um, area_um2, etc.
extract_ldir_image_coords <- function(image_path,
                                      scan_bounds = NULL,
                                      expected_count = NULL,
                                      config = NULL) {
  log_message("Extracting LDIR particle coordinates from image")

  # --- Step 1: Detect scan circle for calibrated µm mapping ---
  #
  # Circle detection can be disabled via config$ldir_use_circle_detection = FALSE.
  # When disabled (or when config is NULL and the flag is absent), the function
  # reads the image dimensions and uses full-image bounds as the "circle"
  # (equivalent to ldir_export_format = 'mosaic').  This is the correct behaviour
  # for LDIR exports that already contain the full scan field without a visible
  # circular crop.
  use_circle <- !isFALSE(config$ldir_use_circle_detection)  # default TRUE when absent

  if (use_circle) {
    circle_info <- detect_ldir_scan_circle(image_path, config = config)

    # --- Manual circle override + hard failure guard ---
    if (!isTRUE(circle_info$detected)) {
      mc <- if (!is.null(config)) config$ldir_circle_manual else NULL
      if (!is.null(mc) && is.numeric(mc$cx) && is.numeric(mc$cy) && is.numeric(mc$r)) {
        log_message("  Using manual circle override: cx=", mc$cx,
                    " cy=", mc$cy, " r=", mc$r)
        circle_info <- list(
          cx_px       = mc$cx, cy_px = mc$cy, radius_px = mc$r,
          width       = circle_info$width,  height    = circle_info$height,
          edge_gap_px = NA_real_,           export_type = "manual",
          detected    = TRUE,               method     = "manual"
        )
      } else {
        stop("LDIR circle detection failed for: ", image_path, "\n",
             "  Cannot compute reliable LDIR coordinates.\n",
             "  Options:\n",
             "    1. Disable circle detection: config$ldir_use_circle_detection <- FALSE\n",
             "    2. Provide manual coords:    config$ldir_circle_manual <- list(cx=, cy=, r=)")
      }
    }
  } else {
    # Circle detection disabled: read image dimensions and treat the full image
    # as the scan field (full-image bounds, no circular mask applied).
    img_hdr <- tryCatch(
      magick::image_info(magick::image_read(image_path)),
      error = function(e) NULL
    )
    if (is.null(img_hdr)) {
      # Fallback: read whole image to get dimensions
      tmp <- read_image_any(image_path, verbose = FALSE)
      img_w <- if (!is.null(tmp)) ncol(tmp) else 1L
      img_h <- if (!is.null(tmp)) nrow(tmp) else 1L
    } else {
      img_w <- img_hdr$width
      img_h <- img_hdr$height
    }
    log_message("  Circle detection disabled — using full-image bounds (",
                img_w, " x ", img_h, " px)")
    circle_info <- list(
      cx_px       = img_w / 2, cy_px = img_h / 2,
      radius_px   = min(img_w, img_h) / 2,
      width       = img_w, height = img_h,
      edge_gap_px = 0, export_type = "full_field",
      detected    = TRUE, method = "disabled"
    )
  }

  scan_diam_um <- if (!is.null(config$ldir_scan_diameter_um)) {
    config$ldir_scan_diameter_um
  } else if (!is.null(scan_bounds)) {
    scan_bounds$x_max - scan_bounds$x_min
  } else {
    13000
  }

  # Pixel-space diagnostic images are written unconditionally to the run's
  # debug/ subfolder (output_dir/debug/).  If debug=TRUE the same images are
  # also written to config$debug_dir for backward compatibility.
  debug_dir       <- if (!is.null(config) && !is.null(config$debug_dir)) config$debug_dir else NULL
  pixel_diag_dir  <- if (!is.null(config) && !is.null(config$output_dir)) {
    file.path(config$output_dir, "debug")
  } else {
    debug_dir   # fallback: use debug_dir if output_dir not set
  }
  if (!is.null(pixel_diag_dir)) {
    if (!dir.exists(pixel_diag_dir))
      dir.create(pixel_diag_dir, recursive = TRUE, showWarnings = FALSE)
    save_ldir_circle_debug(image_path, circle_info,
                           file.path(pixel_diag_dir, "ldir_circle_detection.png"))
  }
  # Also write to legacy debug_dir when debug=TRUE and it differs from pixel_diag_dir
  if (!is.null(debug_dir) && !identical(debug_dir, pixel_diag_dir)) {
    save_ldir_circle_debug(image_path, circle_info,
                           file.path(debug_dir, "ldir_circle_debug.png"))
  }

  # --- Step 2: Extract particle pixel centroids ---
  # Returns particles with centroid_px_x, centroid_px_y in pixel coordinates
  pixel_particles <- .extract_ldir_pixel_centroids(
    image_path, scan_bounds, expected_count, circle_info
  )

  if (is.null(pixel_particles) || nrow(pixel_particles) == 0) {
    log_message("  No particles extracted from LDIR image")
    return(list(particles = .empty_image_df(), circle_info = circle_info))
  }

  # --- Step 2b: Save pixel-space centroid overlay (always; key diagnostic) ---
  if (all(c("centroid_px_x", "centroid_px_y") %in% names(pixel_particles))) {
    if (!is.null(pixel_diag_dir)) {
      save_ldir_centroids_px_debug(
        image_path, pixel_particles, circle_info,
        file.path(pixel_diag_dir, "ldir_pixel_overlay.png")
      )
    }
    # Legacy debug_dir write
    if (!is.null(debug_dir) && !identical(debug_dir, pixel_diag_dir)) {
      save_ldir_centroids_px_debug(
        image_path, pixel_particles, circle_info,
        file.path(debug_dir, "ldir_centroids_px_debug.png")
      )
    }
  }

  # --- Step 2c: Remove particles outside the scan circle ---
  # The Python pipeline applies a circle mask, but the R fallback may not.
  # Guard with a generous 5 % tolerance to keep edge-touching particles.
  if (circle_info$radius_px > 0 &&
      all(c("centroid_px_x", "centroid_px_y") %in% names(pixel_particles))) {
    dist2 <- (pixel_particles$centroid_px_x - circle_info$cx_px)^2 +
              (pixel_particles$centroid_px_y - circle_info$cy_px)^2
    inside <- dist2 <= (circle_info$radius_px * 1.05)^2
    n_outside <- sum(!inside)
    if (n_outside > 0) {
      log_message("  Circle post-filter: removed ", n_outside,
                  " particles outside scan circle (kept ", sum(inside), ")")
      pixel_particles <- pixel_particles[inside, ]
    }
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

  # Augment circle_info with derived scale so callers have a single object
  circle_info$scale_um_per_px <- scale_um_per_px

  # Store calibration metadata as attributes (backward compat)
  attr(pixel_particles, "circle_cx_px") <- circle_info$cx_px
  attr(pixel_particles, "circle_cy_px") <- circle_info$cy_px
  attr(pixel_particles, "circle_radius_px") <- circle_info$radius_px
  attr(pixel_particles, "scale_um_per_px") <- scale_um_per_px

  # --- Step E: Quantitative compression check ---
  # max_abs_um should be ~= scan_diam_um/2. Consistently <90% means
  # radius_px is too large (scale_um_per_px too small) — likely the
  # circle detection is fitting a circle larger than the actual scan area.
  max_abs_um <- max(abs(c(pixel_particles$x_um, pixel_particles$y_um)),
                    na.rm = TRUE)
  fill_pct <- round(max_abs_um / (scan_diam_um / 2) * 100, 1)

  log_message("  Circle-calibrated mapping: scale=",
              round(scale_um_per_px, 3), " µm/px, ",
              nrow(pixel_particles), " particles")
  log_message("  Extent check: max_abs_um=", round(max_abs_um, 1),
              " µm  (", fill_pct, "% of ", scan_diam_um / 2, " µm half-extent)",
              if (fill_pct < 85) " ← LOW: circle radius may be too large" else "")

  # --- Save calibration text file (always, not just debug mode) ---
  output_dir <- if (!is.null(config) && !is.null(config$output_dir)) config$output_dir else NULL
  calib_dir  <- output_dir %||% debug_dir
  if (!is.null(calib_dir)) {
    save_ldir_calibration(
      circle_info     = circle_info,
      scan_diam_um    = scan_diam_um,
      scale_um_per_px = scale_um_per_px,
      n_particles     = nrow(pixel_particles),
      max_abs_um      = max_abs_um,
      output_path     = file.path(calib_dir, "ldir_calibration.txt")
    )
  }
  if (!is.null(debug_dir) && !is.null(output_dir) && debug_dir != output_dir) {
    # Also write to debug dir for backward compat when debug=TRUE
    save_ldir_calibration(
      circle_info     = circle_info,
      scan_diam_um    = scan_diam_um,
      scale_um_per_px = scale_um_per_px,
      n_particles     = nrow(pixel_particles),
      max_abs_um      = max_abs_um,
      output_path     = file.path(debug_dir, "ldir_calibration.txt")
    )
  }

  list(particles = pixel_particles, circle_info = circle_info)
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
#' @param circle_info Result from detect_ldir_scan_circle() (optional).
#'   When provided, the scan circle is applied as a mask in the Python
#'   backend and in the R saturation fallback.
#' @return Data frame with centroid_px_x, centroid_px_y, and standard columns
.extract_ldir_pixel_centroids <- function(image_path, scan_bounds, expected_count,
                                          circle_info = NULL) {
  # Extract circle parameters for Python (use -1 to signal "no mask")
  cx <- if (!is.null(circle_info) && isTRUE(circle_info$radius_px > 0))
          circle_info$cx_px else -1.0
  cy <- if (!is.null(circle_info) && isTRUE(circle_info$radius_px > 0))
          circle_info$cy_px else -1.0
  cr <- if (!is.null(circle_info) && isTRUE(circle_info$radius_px > 0))
          circle_info$radius_px else -1.0

  # --- Python backend (required — no R fallback) ---
  py_result <- tryCatch({
    detect_particles_python(
      image_path     = image_path,
      scan_bounds    = scan_bounds,
      expected_count = expected_count,
      circle_cx      = cx,
      circle_cy      = cy,
      circle_r       = cr
    )
  }, error = function(e) {
    stop(
      "LDIR image particle extraction requires the Python backend.\n",
      "  Python error: ", conditionMessage(e), "\n",
      "  Ensure Python is available and the ldir_detect script is installed.\n",
      "  LDIR image extraction cannot proceed without Python."
    )
  })

  if (is.null(py_result) || nrow(py_result) == 0) {
    stop(
      "Python particle detector returned no particles for: ", image_path, "\n",
      "  Check that the LDIR image is readable and the Python environment is correct."
    )
  }

  log_message("  Extracted ", nrow(py_result), " particles from LDIR image (Python)")

  # Normalise column names: ensure centroid_px_x / centroid_px_y exist
  if ("centroid_x" %in% names(py_result)) {
    py_result$centroid_px_x <- py_result$centroid_x
    py_result$centroid_px_y <- py_result$centroid_y
  } else {
    # Reverse-compute pixel coords from x_um/y_um when centroid_x not present
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

  py_result
}


#' Saturation-based extraction returning pixel centroids
#'
#' Same algorithm as .extract_ldir_saturation but returns centroid_px_x/y
#' in addition to x_um/y_um (which are computed for size filtering only).
.extract_ldir_saturation_px <- function(img, h, w, scan_bounds, expected_count,
                                         circle_info = NULL) {
  r_ch <- img[,,1]; g_ch <- img[,,2]; b_ch <- img[,,3]

  mx <- pmax(r_ch, g_ch, b_ch)
  mn <- pmin(r_ch, g_ch, b_ch)
  sat <- ifelse(mx > 0, (mx - mn) / mx, 0)

  binary <- sat > 0.3 & mx > 0.08

  # Mask out pixels outside the scan circle
  if (!is.null(circle_info) && isTRUE(circle_info$radius_px > 0)) {
    cx <- circle_info$cx_px
    cy <- circle_info$cy_px
    r  <- circle_info$radius_px
    rows_m <- matrix(seq_len(h), nrow = h, ncol = w)
    cols_m <- matrix(seq_len(w), nrow = h, ncol = w, byrow = TRUE)
    in_circle <- (rows_m - cy)^2 + (cols_m - cx)^2 <= (r * 1.05)^2
    binary <- binary & in_circle
  }

  n_fg <- sum(binary)
  log_message("  Saturation threshold: ", n_fg, " foreground pixels (",
              round(n_fg / (h * w) * 100, 1), "%)")

  cc <- .two_pass_ccl(binary, h, w)
  lab_mat <- cc$labels
  n_components <- cc$n_components

  if (n_components == 0) return(.empty_image_df())

  min_pixels <- 10
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
