# =============================================================================
# 01b_ingest_image.R — Image-based particle coordinate extraction
# =============================================================================
#
# When an instrument exports a chemical image (e.g., FTIR absorption heatmap)
# but NOT tabular particle coordinates, this module extracts particle
# centroids from the image and maps them to physical (um) coordinates.
#
# Approach:
#   1. Read the image as an RGB array
#   2. Convert to grayscale using mean(R,G,B) — color-agnostic
#   3. Downsample for performance (configurable factor)
#   4. Adaptive threshold: compare each pixel to local mean (box filter)
#   5. Two-pass connected component labeling (efficient, no recursion)
#   6. Compute centroid, area, bounding box per component
#   7. Scale pixel coordinates to um using known scan bounds
#   8. Filter by minimum physical size (um)
#   9. Sanity-check count against expected tabular count
# =============================================================================


#' Extract particle coordinates from an instrument chemical image
#'
#' @param image_path  Path to PNG image file
#' @param scan_bounds Named list with x_min, x_max, y_min, y_max in um
#'   defining the physical area the image covers. If NULL, uses pixel coords.
#' @param adaptive_radius  Radius (in downsampled pixels) for the local mean
#'   window. Must be larger than the largest expected particle in the
#'   downsampled image. Default 75 (= 151x151 window).
#' @param adaptive_offset  A pixel is foreground if it deviates from the local
#'   mean by more than this value (on 0-1 scale). Default 0.02.
#' @param min_pixels  Minimum component size (in downsampled pixels) to keep
#' @param min_size_um Minimum particle size in um (feret_max). Blobs smaller
#'   than this are removed after scaling to physical coordinates.
#' @param downsample  Integer downsample factor (e.g., 4 -> 1/4 resolution)
#' @param expected_count  Expected number of particles (e.g., from tabular
#'   data). If provided, triggers auto-trimming when the extracted count
#'   exceeds 1.5x the expected value.
#' @param instrument  Character label for source instrument
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um, area_um2, feret_max_um, feret_min_um, etc.
extract_particles_from_image <- function(image_path,
                                          scan_bounds = NULL,
                                          adaptive_radius = 75L,
                                          adaptive_offset = 0.02,
                                          min_pixels = 4,
                                          min_size_um = 20,
                                          downsample = 4L,
                                          expected_count = NULL,
                                          instrument = "IMAGE") {
  log_message("Extracting particles from image: ", basename(image_path))

  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' required. Install with: install.packages('png')")
  }

  # --- 1. Read image ---
  img <- png::readPNG(image_path)
  h_full <- nrow(img)
  w_full <- ncol(img)
  n_ch   <- if (length(dim(img)) == 3) dim(img)[3] else 1
  log_message("  Image dimensions: ", w_full, " x ", h_full, " px, ", n_ch, " channels")

  # --- 2. Convert to grayscale: mean(R,G,B), ignoring alpha (A) ---
  #     This is color-agnostic — no bias toward any channel, works on
  #     false-color maps, B&W images, and any instrument color scheme.
  if (n_ch >= 3) {
    gray <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  } else {
    gray <- if (n_ch == 1) img else img[, , 1]
  }

  # --- 3. Downsample for performance ---
  ds <- as.integer(max(1L, downsample))
  if (ds > 1) {
    row_idx <- seq(1, h_full, by = ds)
    col_idx <- seq(1, w_full, by = ds)
    gray_ds <- gray[row_idx, col_idx]
    h_ds    <- nrow(gray_ds)
    w_ds    <- ncol(gray_ds)
    log_message("  Downsampled ", ds, "x: ", w_ds, " x ", h_ds, " px")
  } else {
    gray_ds <- gray
    h_ds    <- h_full
    w_ds    <- w_full
  }

  # --- 4. Adaptive threshold ---
  #     Compare each pixel to the local mean in a (2*radius+1) window.
  #     A pixel is foreground if |pixel - local_mean| > offset.
  #     This is invariant to global color/intensity — only local contrast
  #     matters, making it robust across different instruments and color maps.
  local_mean <- .box_filter(gray_ds, radius = as.integer(adaptive_radius))
  binary <- abs(gray_ds - local_mean) > adaptive_offset
  n_fg   <- sum(binary)
  log_message("  Adaptive threshold (radius=", adaptive_radius,
              ", offset=", adaptive_offset, "): ",
              n_fg, " foreground pixels (",
              round(n_fg / (h_ds * w_ds) * 100, 2), "%)")

  if (n_fg == 0) {
    log_message("  No foreground pixels. Returning empty.", level = "WARN")
    return(.empty_image_df())
  }

  # --- 5. Two-pass connected component labeling ---
  cc <- .two_pass_ccl(binary, h_ds, w_ds)
  lab_mat <- cc$labels
  n_components <- cc$n_components
  log_message("  Found ", n_components, " raw components")

  if (n_components == 0) return(.empty_image_df())

  # --- 6. Compute properties per component ---
  fg_idx    <- which(binary, arr.ind = TRUE)  # row, col
  fg_labels <- lab_mat[binary]

  tab <- tabulate(fg_labels, nbins = n_components)
  keep_ids <- which(tab >= min_pixels)

  if (length(keep_ids) == 0) {
    log_message("  No components >= ", min_pixels, " px.", level = "WARN")
    return(.empty_image_df())
  }

  # Build a factor for tapply (only keep large components)
  label_factor <- factor(fg_labels, levels = keep_ids)
  valid <- !is.na(label_factor)
  fg_rows <- fg_idx[valid, 1]
  fg_cols <- fg_idx[valid, 2]
  lf      <- droplevels(label_factor[valid])

  # Centroids
  cy <- tapply(fg_rows, lf, mean)
  cx <- tapply(fg_cols, lf, mean)

  # Bounding boxes
  row_min <- tapply(fg_rows, lf, min)
  row_max <- tapply(fg_rows, lf, max)
  col_min <- tapply(fg_cols, lf, min)
  col_max <- tapply(fg_cols, lf, max)
  bb_h <- row_max - row_min + 1
  bb_w <- col_max - col_min + 1

  # Areas
  areas <- tab[keep_ids]

  # Mean intensity
  int_vals <- gray_ds[cbind(fg_rows, fg_cols)]
  mean_int <- tapply(int_vals, lf, mean)

  n_kept <- length(keep_ids)
  log_message("  Retained ", n_kept, " components (>= ", min_pixels, " px)")

  # --- 7. Scale pixel -> um (accounting for downsample) ---
  cx_full <- (as.numeric(cx) - 0.5) * ds
  cy_full <- (as.numeric(cy) - 0.5) * ds

  if (!is.null(scan_bounds)) {
    x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / w_full
    y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / h_full
    x_um <- scan_bounds$x_min + cx_full * x_scale
    y_um <- scan_bounds$y_min + cy_full * y_scale
    size_scale <- max(x_scale, y_scale) * ds
    area_scale <- x_scale * y_scale * ds * ds
    log_message("  Scale: ", round(x_scale, 2), " x ", round(y_scale, 2), " um/px")
  } else {
    x_um <- cx_full
    y_um <- cy_full
    size_scale <- ds
    area_scale <- ds * ds
  }

  df <- data.frame(
    particle_id  = paste0("IMG_", seq_len(n_kept)),
    x_um         = x_um,
    y_um         = y_um,
    area_um2     = as.numeric(areas) * area_scale,
    major_um     = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    minor_um     = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_min_um = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_max_um = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    material     = NA_character_,
    quality      = as.numeric(mean_int),
    source_file  = basename(image_path),
    stringsAsFactors = FALSE
  )

  log_message("  Extracted ", nrow(df), " particles (before size filter)")

  # --- 8. Filter by minimum physical size (D) ---
  if (min_size_um > 0) {
    n_before <- nrow(df)
    df <- df[df$feret_max_um >= min_size_um, ]
    log_message("  Size filter (>= ", min_size_um, " um): kept ", nrow(df),
                " of ", n_before, " particles")
  }

  # --- 9. Sanity check against expected count (E) ---
  if (!is.null(expected_count) && expected_count > 0 && nrow(df) > 0) {
    ratio <- nrow(df) / expected_count
    log_message("  Count check: ", nrow(df), " extracted vs ",
                expected_count, " expected (ratio ", round(ratio, 2), "x)")
    if (ratio > 1.5) {
      target_n <- round(expected_count * 1.3)
      if (nrow(df) > target_n) {
        size_cutoff <- sort(df$feret_max_um, decreasing = TRUE)[target_n]
        df <- df[df$feret_max_um >= size_cutoff, ]
        log_message("  Auto-trimmed to ", nrow(df),
                    " particles (kept largest, size >= ",
                    round(size_cutoff, 1), " um)")
      }
    } else if (ratio < 0.5) {
      log_message("  WARNING: extracted count much lower than expected.",
                  level = "WARN")
      log_message("  Consider lowering adaptive_offset or min_size_um.",
                  level = "WARN")
    }
  }

  if (!is.null(scan_bounds) && nrow(df) > 0) {
    log_message("  Coordinate range: X [", round(min(df$x_um)), ", ",
                round(max(df$x_um)), "], Y [", round(min(df$y_um)), ", ",
                round(max(df$y_um)), "]")
  }

  log_message("  Final: ", nrow(df), " particles from image")
  df
}


#' Fast box filter using integral image
#'
#' Computes the local mean over a (2*radius+1) x (2*radius+1) window for
#' every pixel. Uses an integral image (summed area table) for O(1) per-pixel
#' cost. Fully vectorized — no R loops.
#'
#' @param mat Numeric matrix
#' @param radius Integer radius of the box filter window
#' @return Numeric matrix of same size, containing local means
.box_filter <- function(mat, radius) {
  h <- nrow(mat)
  w <- ncol(mat)
  r <- as.integer(radius)

  # Build integral image: ii[i,j] = sum(mat[1:i, 1:j])
  ii <- apply(mat, 2, cumsum)        # cumsum down each column
  ii <- t(apply(ii, 1, cumsum))      # cumsum across each row

  # Zero-padded integral image (h+1 x w+1)
  # pad[i+1, j+1] = ii[i,j]; pad[1, :] = 0; pad[:, 1] = 0
  pad <- matrix(0, h + 1L, w + 1L)
  pad[2:(h + 1L), 2:(w + 1L)] <- ii

  # Window bounds per row/column (clamped to image edges)
  r_lo <- pmax(1L, seq_len(h) - r)
  r_hi <- pmin(h,  seq_len(h) + r)
  c_lo <- pmax(1L, seq_len(w) - r)
  c_hi <- pmin(w,  seq_len(w) + r)

  # Vectorized lookup: pad[rows, cols] gives h x w submatrix
  # sum(mat[r1:r2, c1:c2]) = pad[r2+1,c2+1] - pad[r1,c2+1] - pad[r2+1,c1] + pad[r1,c1]
  sums <- pad[r_hi + 1L, c_hi + 1L] -
          pad[r_lo,       c_hi + 1L] -
          pad[r_hi + 1L, c_lo      ] +
          pad[r_lo,       c_lo      ]

  # Number of pixels in each window (varies at edges)
  counts <- outer(r_hi - r_lo + 1L, c_hi - c_lo + 1L)

  sums / counts
}


#' Two-pass connected component labeling (4-connected)
#'
#' Efficient O(n) algorithm using union-find. Works on binary matrix.
#' @param binary Logical matrix (TRUE = foreground)
#' @param h Number of rows
#' @param w Number of columns
#' @return List with labels (integer matrix) and n_components (int)
.two_pass_ccl <- function(binary, h, w) {
  lab <- matrix(0L, nrow = h, ncol = w)

  # Union-find data structure (dynamic, grows as needed)
  uf_parent <- integer(0)
  next_label <- 1L

  # --- Pass 1: assign provisional labels ---
  for (r in seq_len(h)) {
    for (cc in seq_len(w)) {
      if (!binary[r, cc]) next

      up   <- if (r > 1 && binary[r - 1, cc]) lab[r - 1, cc] else 0L
      left <- if (cc > 1 && binary[r, cc - 1]) lab[r, cc - 1] else 0L

      if (up == 0L && left == 0L) {
        uf_parent <- c(uf_parent, next_label)
        lab[r, cc] <- next_label
        next_label <- next_label + 1L
      } else if (up != 0L && left == 0L) {
        lab[r, cc] <- up
      } else if (up == 0L && left != 0L) {
        lab[r, cc] <- left
      } else {
        lab[r, cc] <- min(up, left)
        if (up != left) {
          # Union
          ra <- up; while (uf_parent[ra] != ra) ra <- uf_parent[ra]
          rb <- left; while (uf_parent[rb] != rb) rb <- uf_parent[rb]
          if (ra != rb) uf_parent[max(ra, rb)] <- min(ra, rb)
        }
      }
    }
  }

  if (next_label == 1L) {
    return(list(labels = lab, n_components = 0L))
  }

  # --- Pass 2: resolve equivalences ---
  n_prov <- next_label - 1L
  roots <- integer(n_prov)
  for (i in seq_len(n_prov)) {
    x <- i
    while (uf_parent[x] != x) x <- uf_parent[x]
    roots[i] <- x
  }

  unique_roots <- unique(roots)
  remap <- integer(n_prov)
  remap[unique_roots] <- seq_along(unique_roots)
  final_map <- remap[roots]

  fg <- lab > 0L
  lab[fg] <- final_map[lab[fg]]

  list(labels = lab, n_components = length(unique_roots))
}


#' Return an empty data frame with standard image particle columns
.empty_image_df <- function() {
  data.frame(
    particle_id  = character(),
    x_um         = numeric(),
    y_um         = numeric(),
    area_um2     = numeric(),
    major_um     = numeric(),
    minor_um     = numeric(),
    feret_min_um = numeric(),
    feret_max_um = numeric(),
    material     = character(),
    quality      = numeric(),
    source_file  = character(),
    stringsAsFactors = FALSE
  )
}


#' Calibrate image-to-physical coordinate mapping
#'
#' Given known particle coordinates (tabular FTIR data) and image-extracted
#' positions, estimates the affine transform (pixel -> um) by matching the
#' largest particles between both sets.
#'
#' @param known_df  Data frame with x_um, y_um, feret_max_um (from tabular)
#' @param image_df  Data frame with x_um, y_um, feret_max_um (from image)
#' @param n_top     Number of largest particles to use for calibration
#' @return List with transform matrix, scale, rotation, residual_rms
calibrate_image_coords <- function(known_df, image_df, n_top = 20) {
  log_message("Calibrating image coordinates against known positions")

  n_top <- min(n_top, nrow(known_df), nrow(image_df))

  known_order <- order(known_df$feret_max_um, decreasing = TRUE)[seq_len(n_top)]
  image_order <- order(image_df$feret_max_um, decreasing = TRUE)[seq_len(n_top)]

  tf <- estimate_similarity_transform(
    src_x = image_df$x_um[image_order],
    src_y = image_df$y_um[image_order],
    dst_x = known_df$x_um[known_order],
    dst_y = known_df$y_um[known_order],
    allow_reflection = TRUE
  )

  log_message("  Calibration: scale=", round(tf$scale, 4),
              ", rot=", round(tf$rotation_deg, 2), "deg",
              ", reflected=", tf$reflected,
              ", rms=", round(tf$residual_rms, 1), " um")

  list(
    transform    = tf$matrix,
    scale        = tf$scale,
    rotation     = tf$rotation_deg,
    reflected    = tf$reflected,
    residual_rms = tf$residual_rms,
    n_matched    = n_top
  )
}
