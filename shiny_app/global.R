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
    matched         = "matched_particles.csv",
    unmatched_ftir  = "unmatched_ftir.csv",
    unmatched_raman = "unmatched_raman.csv",
    agreement       = "agreement_summary.csv"
  )

  for (nm in names(paths)) {
    fp <- file.path(run_dir, paths[[nm]])
    if (file.exists(fp)) data[[nm]] <- read.csv(fp, stringsAsFactors = FALSE)
  }

  tp <- file.path(run_dir, "transform_params.txt")
  if (file.exists(tp)) data$transform <- parse_transform_params(tp)

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

  list(
    scale         = get_num("scale"),
    rotation_deg  = get_num("rotation_deg"),
    reflected     = tolower(get_val("reflected")) == "true",
    M             = M,
    ftir_centroid  = c(get_num("ftir_centroid_x"), get_num("ftir_centroid_y")),
    raman_centroid = c(get_num("raman_centroid_x"), get_num("raman_centroid_y"))
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
# Transform an image for display: rotate + compute aligned physical bounds
#
# img_raster : PNG raster array (h x w x channels)
# img_bounds : list(xmin, xmax, ymin, ymax) in original FTIR coords
# M_full     : 3x3 full transform matrix (original -> aligned)
#
# Returns: list(raster, xmin, xmax, ymin, ymax) for annotation_raster
# ---------------------------------------------------------------------------
transform_image_for_display <- function(img_raster, img_bounds, M_full,
                                         rotation_deg) {
  # Transform the 4 image corners to aligned coordinates
  cx <- c(img_bounds$xmin, img_bounds$xmax, img_bounds$xmin, img_bounds$xmax)
  cy <- c(img_bounds$ymin, img_bounds$ymin, img_bounds$ymax, img_bounds$ymax)
  tc <- transform_points(cx, cy, M_full)

  # Axis-aligned bounding box in aligned coords
  aligned_xmin <- min(tc$x)
  aligned_xmax <- max(tc$x)
  aligned_ymin <- min(tc$y)
  aligned_ymax <- max(tc$y)

  # Rotate the raster by the transform angle
  # Note: ggplot annotation_raster has y-axis going up, but PNG row 1 = top.
  # annotation_raster handles this: ymin=bottom, ymax=top, raster row 1=top.
  # After rotation, the image is correctly oriented for the new coordinate frame.
  rot <- rotate_raster(img_raster, rotation_deg)

  list(
    raster = rot$raster,
    xmin   = aligned_xmin,
    xmax   = aligned_xmax,
    ymin   = aligned_ymin,
    ymax   = aligned_ymax
  )
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
# Tries PNG first, then JPEG as fallback if extension is ambiguous.
# ---------------------------------------------------------------------------
load_image_raster <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))

  # Try the format matching the extension first, then fallback
  if (ext %in% c("jpg", "jpeg")) {
    img <- tryCatch(jpeg::readJPEG(path), error = function(e) NULL)
    if (!is.null(img)) return(img)
    return(tryCatch(png::readPNG(path), error = function(e) NULL))
  }
  img <- tryCatch(png::readPNG(path), error = function(e) NULL)
  if (!is.null(img)) return(img)
  if (requireNamespace("jpeg", quietly = TRUE))
    return(tryCatch(jpeg::readJPEG(path), error = function(e) NULL))
  NULL
}

# ---------------------------------------------------------------------------
# Find nearest particle to a hover coordinate
# ---------------------------------------------------------------------------
nearest_particle <- function(df, hover_x, hover_y, max_dist_frac = 0.03) {
  if (is.null(df) || nrow(df) == 0 || is.null(hover_x) || is.null(hover_y)) return(NULL)
  dx <- df$x - hover_x
  dy <- df$y - hover_y
  dists <- sqrt(dx^2 + dy^2)
  idx <- which.min(dists)
  x_range <- diff(range(df$x, na.rm = TRUE))
  y_range <- diff(range(df$y, na.rm = TRUE))
  threshold <- max(x_range, y_range) * max_dist_frac
  if (dists[idx] > threshold) return(NULL)
  df[idx, , drop = FALSE]
}
