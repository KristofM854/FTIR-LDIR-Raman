# =============================================================================
# global.R — Load pipeline output and prepare data for Shiny
# =============================================================================

library(shiny)
library(ggplot2)
library(png)

# ---------------------------------------------------------------------------
# Locate pipeline output.
# Supports two layouts:
#   1. Subdirectory runs:  output/2026-02-16_6/matched_particles.csv
#   2. Flat timestamped:   output/matched_particles_20260216_172156.csv
# Returns a list(dir, format) or NULL.
# ---------------------------------------------------------------------------
find_latest_run <- function(output_dir = file.path("..", "output")) {
  if (!dir.exists(output_dir)) return(NULL)

  # --- Try subdirectory format first ---
  runs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  # Keep only dirs that actually contain a matched_particles.csv
  runs <- runs[file.exists(file.path(runs, "matched_particles.csv"))]
  if (length(runs) > 0) {
    runs <- runs[order(file.mtime(runs), decreasing = TRUE)]
    return(list(dir = runs[1], format = "subdir"))
  }

  # --- Try flat timestamped format ---
  flat <- list.files(output_dir, pattern = "^matched_particles.*\\.csv$",
                     full.names = TRUE)
  if (length(flat) > 0) {
    flat <- flat[order(file.mtime(flat), decreasing = TRUE)]
    return(list(dir = output_dir, format = "flat",
                matched_file = flat[1]))
  }

  NULL
}

# ---------------------------------------------------------------------------
# Load pipeline CSVs from either subdirectory or flat format
# ---------------------------------------------------------------------------
load_run_data <- function(run_info) {
  data <- list(run_dir = run_info$dir)

  if (run_info$format == "subdir") {
    # Simple: files have fixed names in a subdirectory
    file_map <- list(
      matched         = "matched_particles.csv",
      unmatched_ftir  = "unmatched_ftir.csv",
      unmatched_raman = "unmatched_raman.csv",
      agreement       = "agreement_summary.csv"
    )
    for (nm in names(file_map)) {
      fp <- file.path(run_info$dir, file_map[[nm]])
      if (file.exists(fp)) data[[nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }
    tp <- file.path(run_info$dir, "transform_params.txt")
    if (file.exists(tp)) data$transform <- parse_transform_params(tp)

  } else {
    # Flat: files have timestamps in the name. Find each by prefix.
    find_flat <- function(prefix, ext = "csv") {
      pat <- paste0("^", prefix, ".*\\.", ext, "$")
      hits <- list.files(run_info$dir, pattern = pat, full.names = TRUE)
      if (length(hits) == 0) return(NULL)
      hits[order(file.mtime(hits), decreasing = TRUE)][1]
    }

    for (info in list(
      list(nm = "matched",         prefix = "matched_particles"),
      list(nm = "unmatched_ftir",  prefix = "unmatched_ftir"),
      list(nm = "unmatched_raman", prefix = "unmatched_raman"),
      list(nm = "agreement",       prefix = "agreement_summary")
    )) {
      fp <- find_flat(info$prefix)
      if (!is.null(fp)) data[[info$nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }
    tp <- find_flat("transform_params", "txt")
    if (!is.null(tp)) data$transform <- parse_transform_params(tp)
  }

  data
}

# ---------------------------------------------------------------------------
# Build run_data from user-uploaded CSVs (fallback when pipeline output is
# not available).  Minimum required: matched_particles.csv.
# ---------------------------------------------------------------------------
load_uploaded_data <- function(matched_path,
                                unmatched_ftir_path = NULL,
                                unmatched_raman_path = NULL,
                                transform_path = NULL) {
  data <- list(run_dir = dirname(matched_path))
  data$matched <- read.csv(matched_path, stringsAsFactors = FALSE)

  if (!is.null(unmatched_ftir_path) && file.exists(unmatched_ftir_path))
    data$unmatched_ftir <- read.csv(unmatched_ftir_path, stringsAsFactors = FALSE)
  if (!is.null(unmatched_raman_path) && file.exists(unmatched_raman_path))
    data$unmatched_raman <- read.csv(unmatched_raman_path, stringsAsFactors = FALSE)
  if (!is.null(transform_path) && file.exists(transform_path))
    data$transform <- parse_transform_params(transform_path)

  data
}

# ---------------------------------------------------------------------------
# Parse transform_params.txt into a structured list
# ---------------------------------------------------------------------------
parse_transform_params <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)

  get_val <- function(key) {
    pat <- paste0("^", key, ":\\s*")
    idx <- grep(pat, lines)
    if (length(idx) == 0) return(NA)
    trimws(sub(pat, "", lines[idx[1]]))
  }

  get_num <- function(key) as.numeric(get_val(key))

  get_row <- function(key) {
    v <- get_val(key)
    if (is.na(v)) return(NULL)
    as.numeric(strsplit(gsub("\\s+", "", v), ",")[[1]])
  }

  # 3x3 transform matrix (FTIR_norm -> Raman_norm)
  r1 <- get_row("matrix_row1")
  r2 <- get_row("matrix_row2")
  r3 <- get_row("matrix_row3")
  M <- NULL
  if (!is.null(r1) && !is.null(r2) && !is.null(r3)) {
    M <- rbind(r1, r2, r3)
  }

  # FTIR scan bounds (optional — present in newer pipeline output)
  scan_xmin <- get_num("ftir_scan_xmin")
  scan_xmax <- get_num("ftir_scan_xmax")
  scan_ymin <- get_num("ftir_scan_ymin")
  scan_ymax <- get_num("ftir_scan_ymax")
  scan_bounds <- NULL
  if (!is.na(scan_xmax)) {
    scan_bounds <- list(xmin = scan_xmin, xmax = scan_xmax,
                        ymin = scan_ymin, ymax = scan_ymax)
  }

  list(
    scale         = get_num("scale"),
    rotation_deg  = get_num("rotation_deg"),
    reflected     = tolower(get_val("reflected")) == "true",
    M             = M,
    ftir_centroid  = c(get_num("ftir_centroid_x"), get_num("ftir_centroid_y")),
    raman_centroid = c(get_num("raman_centroid_x"), get_num("raman_centroid_y")),
    ftir_scan_bounds = scan_bounds
  )
}

# ---------------------------------------------------------------------------
# Build the full 3x3 transform: FTIR original coords -> Raman coords
#   x_aligned = M * (x_orig - ftir_centroid) + raman_centroid
# ---------------------------------------------------------------------------
build_full_transform <- function(transform) {
  # Both FTIR aligned and Raman are displayed in the normalized (centered) frame.
  # Pipeline: x_norm = x_orig - ftir_centroid, then x_aligned = M * x_norm.
  # So full transform: subtract ftir centroid, then apply M.  No raman centroid added.
  T1 <- diag(3)
  T1[1, 3] <- -transform$ftir_centroid[1]
  T1[2, 3] <- -transform$ftir_centroid[2]

  # Full: M %*% T1  (maps FTIR original -> normalized aligned frame)
  transform$M %*% T1
}

# ---------------------------------------------------------------------------
# Transform 2D points with a 3x3 homogeneous matrix
# ---------------------------------------------------------------------------
transform_points <- function(x, y, M_full) {
  pts <- rbind(x, y, rep(1, length(x)))
  result <- M_full %*% pts
  list(x = result[1, ], y = result[2, ])
}

# ---------------------------------------------------------------------------
# Rotate a raster image array by an arbitrary angle (degrees)
#
# Returns: list(raster = rotated array, corners = transformed physical corners)
# The raster is suitable for annotation_raster; corners give the xmin/xmax/ymin/ymax.
# ---------------------------------------------------------------------------
rotate_raster <- function(img, angle_deg) {
  theta <- angle_deg * pi / 180
  cos_t <- cos(theta)
  sin_t <- sin(theta)

  h <- nrow(img)
  w <- ncol(img)
  nc <- if (length(dim(img)) == 3) dim(img)[3] else 1

  # Image center (rotation pivot)
  cx <- (w + 1) / 2
  cy <- (h + 1) / 2

  # Compute output dimensions from rotated bounding box
  # Corners of input image in centered coords
  corners_x <- c(1, w, 1, w) - cx
  corners_y <- c(1, 1, h, h) - cy
  rot_x <- corners_x * cos_t - corners_y * sin_t
  rot_y <- corners_x * sin_t + corners_y * cos_t

  out_w <- ceiling(max(rot_x) - min(rot_x))
  out_h <- ceiling(max(rot_y) - min(rot_y))
  out_cx <- (out_w + 1) / 2
  out_cy <- (out_h + 1) / 2

  # Build output coordinates
  out_col <- rep(seq_len(out_w), out_h)
  out_row <- rep(seq_len(out_h), each = out_w)

  # Inverse rotation: map output pixel -> input pixel
  dx <- out_col - out_cx
  dy <- out_row - out_cy
  src_x <- dx * cos_t + dy * sin_t + cx
  src_y <- -dx * sin_t + dy * cos_t + cy

  # Nearest-neighbour sampling (fast)
  src_col <- round(src_x)
  src_row <- round(src_y)
  valid <- src_col >= 1 & src_col <= w & src_row >= 1 & src_row <= h

  if (nc == 1) {
    out <- matrix(0, out_h, out_w)
    out[cbind(out_row[valid], out_col[valid])] <- img[cbind(src_row[valid], src_col[valid])]
  } else {
    out <- array(0, dim = c(out_h, out_w, nc))
    for (ch in seq_len(nc)) {
      ch_in <- img[, , ch]
      ch_out <- matrix(0, out_h, out_w)
      ch_out[cbind(out_row[valid], out_col[valid])] <- ch_in[cbind(src_row[valid], src_col[valid])]
      out[, , ch] <- ch_out
    }
    # Set alpha = 0 for background pixels (transparent)
    if (nc == 4) {
      alpha <- out[, , 4]
      alpha[!matrix(valid, out_h, out_w, byrow = FALSE)] <- 0
      # Reshape valid mask: out_row/out_col are in column-major order for rep()
      vm <- matrix(FALSE, out_h, out_w)
      vm[cbind(out_row[valid], out_col[valid])] <- TRUE
      out[, , 4] <- ifelse(vm, alpha, 0)
    }
  }

  list(raster = out, out_w = out_w, out_h = out_h,
       scale_x = w / out_w, scale_y = h / out_h)
}

# ---------------------------------------------------------------------------
# Transform an image for display via full affine warp.
#
# For each pixel in the output raster, the inverse of M_full maps back to
# the original image coordinates, which are then sampled.  This correctly
# handles the combination of translation (centroid subtraction) + rotation
# + scale in a single pass.
#
# img_raster : raster array (h x w x channels), already Y-flipped
# img_bounds : list(xmin, xmax, ymin, ymax) in original FTIR coords
# M_full     : 3x3 full transform matrix (original -> aligned)
# rotation_deg : (unused, kept for signature compat)
#
# Returns: list(raster, xmin, xmax, ymin, ymax) for annotation_raster
# ---------------------------------------------------------------------------
transform_image_for_display <- function(img_raster, img_bounds, M_full,
                                         rotation_deg = NULL) {
  h <- nrow(img_raster)
  w <- ncol(img_raster)
  nc <- if (length(dim(img_raster)) == 3) dim(img_raster)[3] else 1

  # Transform the 4 image corners to aligned coordinates
  cx <- c(img_bounds$xmin, img_bounds$xmax, img_bounds$xmin, img_bounds$xmax)
  cy <- c(img_bounds$ymin, img_bounds$ymin, img_bounds$ymax, img_bounds$ymax)
  tc <- transform_points(cx, cy, M_full)

  # Output extent: axis-aligned bounding box + small padding
  pad <- 50
  out_xmin <- min(tc$x) - pad
  out_xmax <- max(tc$x) + pad
  out_ymin <- min(tc$y) - pad
  out_ymax <- max(tc$y) + pad

  # Output raster resolution: match input pixels-per-µm
  px_per_um <- w / (img_bounds$xmax - img_bounds$xmin)
  out_w <- round((out_xmax - out_xmin) * px_per_um)
  out_h <- round((out_ymax - out_ymin) * px_per_um)
  # Cap resolution to avoid excessive memory
  max_dim <- 4000
  if (out_w > max_dim || out_h > max_dim) {
    scale <- max_dim / max(out_w, out_h)
    out_w <- round(out_w * scale)
    out_h <- round(out_h * scale)
  }

  # Inverse transform: aligned coords -> original coords
  M_inv <- solve(M_full)

  # For each output pixel, compute its aligned physical position, then
  # inverse-transform to original coords, then map to input pixel.
  out_col <- rep(seq_len(out_w), out_h)
  out_row <- rep(seq_len(out_h), each = out_w)

  # Output pixel -> aligned physical coords
  # annotation_raster: col 1 at xmin, col out_w at xmax
  #                    row 1 at ymax, row out_h at ymin  (Y-flipped image: row 1 = max Y)
  ax <- out_xmin + (out_col - 1) / max(out_w - 1, 1) * (out_xmax - out_xmin)
  ay <- out_ymax - (out_row - 1) / max(out_h - 1, 1) * (out_ymax - out_ymin)

  # Inverse transform to original coords
  orig <- M_inv %*% rbind(ax, ay, rep(1, length(ax)))
  ox <- orig[1, ]
  oy <- orig[2, ]

  # Map original coords to input pixel (Y-flipped image: row 1 = ymax, row h = ymin)
  src_col <- round(1 + (ox - img_bounds$xmin) / (img_bounds$xmax - img_bounds$xmin) * (w - 1))
  src_row <- round(1 + (img_bounds$ymax - oy) / (img_bounds$ymax - img_bounds$ymin) * (h - 1))

  valid <- src_col >= 1 & src_col <= w & src_row >= 1 & src_row <= h

  if (nc == 1) {
    out <- matrix(0, out_h, out_w)
    out[cbind(out_row[valid], out_col[valid])] <-
      img_raster[cbind(src_row[valid], src_col[valid])]
  } else {
    out <- array(0, dim = c(out_h, out_w, nc))
    for (ch in seq_len(nc)) {
      ch_in <- img_raster[, , ch]
      ch_out <- matrix(0, out_h, out_w)
      ch_out[cbind(out_row[valid], out_col[valid])] <-
        ch_in[cbind(src_row[valid], src_col[valid])]
      out[, , ch] <- ch_out
    }
    if (nc >= 4) {
      # Transparent background for unmapped pixels
      alpha <- out[, , 4]
      vm <- matrix(FALSE, out_h, out_w)
      vm[cbind(out_row[valid], out_col[valid])] <- TRUE
      out[, , 4] <- ifelse(vm, alpha, 0)
    }
  }

  list(raster = out, xmin = out_xmin, xmax = out_xmax,
       ymin = out_ymin, ymax = out_ymax)
}

# ---------------------------------------------------------------------------
# Estimate FTIR scan bounds from image dimensions and particle coordinates.
#
# The PerkinElmer Spotlight exports images at ~6 rendering pixels per 25µm
# grid cell.  From a 2993×2993 image: (2993+1)/6 ≈ 499 grid positions,
# giving a 499 * 25 = 12475 µm scan extent.  This function computes the
# bounds robustly from the image dimensions and grid step.
# ---------------------------------------------------------------------------
estimate_ftir_scan_bounds <- function(img_raster, particle_x_um = NULL,
                                       particle_y_um = NULL,
                                       grid_step_um = 25) {
  img_w <- ncol(img_raster)
  img_h <- nrow(img_raster)

  # Estimate grid positions from image pixel count
  # Empirical: (image_px + 1) / 6 gives the grid count
  render_px_per_cell <- 6
  grid_nx <- round((img_w + 1) / render_px_per_cell)
  grid_ny <- round((img_h + 1) / render_px_per_cell)

  # Scan extent: grid_count * grid_step
  x_extent <- grid_nx * grid_step_um
  y_extent <- grid_ny * grid_step_um

  # Sanity check against particle positions if available
  if (!is.null(particle_x_um)) {
    max_px <- max(particle_x_um, na.rm = TRUE)
    if (x_extent < max_px) x_extent <- ceiling(max_px / 500) * 500
  }
  if (!is.null(particle_y_um)) {
    max_py <- max(particle_y_um, na.rm = TRUE)
    if (y_extent < max_py) y_extent <- ceiling(max_py / 500) * 500
  }

  list(xmin = 0, xmax = x_extent, ymin = 0, ymax = y_extent)
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
      x = m$raman_x_norm, y = m$raman_y_norm,
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
      x = u$x_norm, y = u$y_norm,
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
# Coordinate bounds from particle data (safe for NULL / empty inputs)
# ---------------------------------------------------------------------------
compute_bounds <- function(...) {
  dfs <- list(...)
  all_x <- unlist(lapply(dfs, function(d) if (!is.null(d) && nrow(d) > 0) d$x))
  all_y <- unlist(lapply(dfs, function(d) if (!is.null(d) && nrow(d) > 0) d$y))
  all_x <- all_x[is.finite(all_x)]
  all_y <- all_y[is.finite(all_y)]
  if (length(all_x) == 0 || length(all_y) == 0) {
    return(list(x = c(-1000, 1000), y = c(-1000, 1000)))
  }
  pad <- 200
  list(
    x = c(min(all_x) - pad, max(all_x) + pad),
    y = c(min(all_y) - pad, max(all_y) + pad)
  )
}

# ---------------------------------------------------------------------------
# Load an image (PNG or JPEG) as a raster array for ggplot annotation.
# annotation_raster places row 1 at ymax (top of plot).  Standard images
# already have row 1 = top of the visual image, which corresponds to max-Y
# in Cartesian / stage coordinates.  So NO vertical flip is needed.
# ---------------------------------------------------------------------------
load_image_raster <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))

  raw <- NULL
  if (ext %in% c("jpg", "jpeg")) {
    raw <- tryCatch(jpeg::readJPEG(path), error = function(e) NULL)
    if (is.null(raw)) raw <- tryCatch(png::readPNG(path), error = function(e) NULL)
  } else {
    raw <- tryCatch(png::readPNG(path), error = function(e) NULL)
    if (is.null(raw) && requireNamespace("jpeg", quietly = TRUE))
      raw <- tryCatch(jpeg::readJPEG(path), error = function(e) NULL)
  }
  if (is.null(raw)) return(NULL)

  raw
}

# ---------------------------------------------------------------------------
# Find nearest particle to a hover coordinate.
# max_dist_frac: fraction of visible plot range used as snap radius.
# ---------------------------------------------------------------------------
nearest_particle <- function(df, hover_x, hover_y, max_dist_frac = 0.05) {
  if (is.null(df) || nrow(df) == 0 || is.null(hover_x) || is.null(hover_y)) return(NULL)
  dx <- df$x - hover_x
  dy <- df$y - hover_y
  dists <- sqrt(dx^2 + dy^2)
  idx <- which.min(dists)
  x_range <- diff(range(df$x, na.rm = TRUE))
  y_range <- diff(range(df$y, na.rm = TRUE))
  threshold <- max(x_range, y_range, 500) * max_dist_frac
  if (dists[idx] > threshold) return(NULL)
  df[idx, , drop = FALSE]
}
