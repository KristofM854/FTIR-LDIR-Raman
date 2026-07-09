# =============================================================================
# utils.R — Shared utility functions for FTIR–Raman particle matching
# =============================================================================

# ---------------------------------------------------------------------------
# Output directory with timestamped subfolders
# ---------------------------------------------------------------------------

#' Create a timestamped run subfolder inside the base output directory
#'
#' Format: output/YYYY-MM-DD_1, output/YYYY-MM-DD_2, etc.
#' Automatically increments the run number for multiple runs on the same day.
#'
#' @param base_dir Base output directory (e.g., "output")
#' @return Path to the new run-specific subfolder (already created)
make_run_dir <- function(base_dir = "output") {
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  today <- format(Sys.Date(), "%Y-%m-%d")
  existing <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

  # Find existing run numbers for today
  pattern <- paste0("^", gsub("-", "\\\\-", today), "_(\\d+)$")
  matches <- regmatches(existing, regexec(pattern, existing))
  run_numbers <- as.integer(vapply(matches, function(m) {
    if (length(m) == 2) m[2] else NA_character_
  }, character(1)))
  run_numbers <- run_numbers[!is.na(run_numbers)]

  next_run <- if (length(run_numbers) == 0) 1L else max(run_numbers) + 1L
  run_dir <- file.path(base_dir, paste0(today, "_", next_run))
  dir.create(run_dir, recursive = TRUE)

  log_message("Run output directory: ", run_dir)
  run_dir
}

#' Build and create standardized stage-based output subdirectories for a run.
#'
#' @param run_dir Path to the run directory (e.g. "output/2026-03-06_1")
#' @return Named list of directory paths (all created on disk)
get_output_dirs <- function(run_dir) {
  dirs <- list(
    manifest    = file.path(run_dir, "00_manifest"),
    ingested    = file.path(run_dir, "01_ingested"),
    joined      = file.path(run_dir, "02_joined"),
    prefiltered = file.path(run_dir, "03_prefiltered"),
    alignment   = file.path(run_dir, "04_alignment"),
    matches     = file.path(run_dir, "05_matches"),
    agreement   = file.path(run_dir, "06_agreement"),
    diagnostics = file.path(run_dir, "07_diagnostics"),
    plots       = file.path(run_dir, "07_diagnostics", "plots")
  )
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  dirs
}

# ---------------------------------------------------------------------------
# Coordinate parsing
# ---------------------------------------------------------------------------

#' Parse FTIR coordinate string "[x;y]" into numeric x, y
#' @param coord_str Character vector of strings like "[1234;5678]" or "[1234.5;5678.9]"
#' @return Data frame with columns x_um and y_um
parse_ftir_coordinates <- function(coord_str) {
  # Remove brackets and split on semicolon

  cleaned <- gsub("\\[|\\]", "", coord_str)
  parts   <- strsplit(cleaned, ";")

  x_um <- vapply(parts, function(p) as.numeric(p[1]), numeric(1))
  y_um <- vapply(parts, function(p) as.numeric(p[2]), numeric(1))

  data.frame(x_um = x_um, y_um = y_um)
}

# ---------------------------------------------------------------------------
# 2D similarity transform utilities (3x3 homogeneous matrices)
# ---------------------------------------------------------------------------

#' Build a 3x3 homogeneous similarity transform matrix
#'
#' Without reflection:
#'   | a  -b  tx |      a = s*cos(theta), b = s*sin(theta)
#'   | b   a  ty |
#'   | 0   0   1 |
#'
#' With reflection (y-flip before rotation+scale):
#'   | a   b  tx |
#'   | b  -a  ty |
#'   | 0   0   1 |
#'
#' @param a Numeric, s*cos(theta)
#' @param b Numeric, s*sin(theta)
#' @param tx Numeric, x-translation
#' @param ty Numeric, y-translation
#' @param reflect Logical, whether the transform includes reflection
#' @return 3x3 matrix
build_transform_matrix <- function(a, b, tx, ty, reflect = FALSE) {
  if (!reflect) {
    matrix(c(a, b, 0,
             -b, a, 0,
             tx, ty, 1), nrow = 3, byrow = FALSE)
  } else {
    matrix(c(a, b, 0,
             b, -a, 0,
             tx, ty, 1), nrow = 3, byrow = FALSE)
  }
}

#' Estimate a 2D similarity transform from point correspondences
#'
#' Given source points P and destination points Q, find the similarity
#' transform T such that T(P) ≈ Q, minimizing sum of squared residuals.
#' Optional per-point weights allow down-weighting unreliable correspondences
#' (e.g. elongated particles whose centroid is uncertain).
#'
#' @param src_x Numeric vector, source x-coordinates
#' @param src_y Numeric vector, source y-coordinates
#' @param dst_x Numeric vector, destination x-coordinates
#' @param dst_y Numeric vector, destination y-coordinates
#' @param allow_reflection Logical, whether to also try reflection and pick best
#' @param weights Optional numeric vector of per-point weights (positive).
#'        Higher weight = more influence. NULL means equal weights.
#' @return List with: matrix (3x3), a, b, tx, ty, scale, rotation_deg,
#'         reflected, residual_rms
estimate_similarity_transform <- function(src_x, src_y, dst_x, dst_y,
                                          allow_reflection = TRUE,
                                          weights = NULL) {
  n <- length(src_x)
  stopifnot(n >= 2, n == length(src_y), n == length(dst_x), n == length(dst_y))

  # Response vector
  q_vec <- as.numeric(rbind(dst_x, dst_y))  # interleaved: qx1, qy1, qx2, qy2, ...

  # Build per-row weight vector (each point contributes 2 rows: x and y)
  if (!is.null(weights)) {
    stopifnot(length(weights) == n)
    w_sqrt <- sqrt(pmax(weights, 0))
    w_row <- rep(w_sqrt, each = 2)  # interleaved: w1, w1, w2, w2, ...
  } else {
    w_row <- NULL
  }

  # --- No-reflection model ---
  # qx_i = a*px_i - b*py_i + tx
  # qy_i = b*px_i + a*py_i + ty
  A_no <- matrix(0, nrow = 2 * n, ncol = 4)
  for (i in seq_len(n)) {
    row_x <- 2 * i - 1
    row_y <- 2 * i
    A_no[row_x, ] <- c(src_x[i], -src_y[i], 1, 0)
    A_no[row_y, ] <- c(src_y[i],  src_x[i], 0, 1)
  }

  # Apply weights via row scaling: W*A*x = W*q  (weighted least squares)
  if (!is.null(w_row)) {
    A_no_w <- A_no * w_row
    q_no_w <- q_vec * w_row
  } else {
    A_no_w <- A_no
    q_no_w <- q_vec
  }

  fit_no   <- qr.solve(A_no_w, q_no_w)
  res_no   <- q_vec - A_no %*% fit_no  # residuals on unweighted data
  rms_no   <- sqrt(mean(res_no^2))

  results <- list()
  results$no_reflect <- list(
    a = fit_no[1], b = fit_no[2], tx = fit_no[3], ty = fit_no[4],
    rms = rms_no, reflect = FALSE
  )

  if (allow_reflection) {
    # --- Reflection model (y-flip before rotation+scale) ---
    # qx_i =  a*px_i + b*py_i + tx
    # qy_i =  b*px_i - a*py_i + ty
    A_ref <- matrix(0, nrow = 2 * n, ncol = 4)
    for (i in seq_len(n)) {
      row_x <- 2 * i - 1
      row_y <- 2 * i
      A_ref[row_x, ] <- c( src_x[i],  src_y[i], 1, 0)
      A_ref[row_y, ] <- c(-src_y[i],  src_x[i], 0, 1)
    }

    if (!is.null(w_row)) {
      A_ref_w <- A_ref * w_row
      q_ref_w <- q_vec * w_row
    } else {
      A_ref_w <- A_ref
      q_ref_w <- q_vec
    }

    fit_ref <- qr.solve(A_ref_w, q_ref_w)
    res_ref <- q_vec - A_ref %*% fit_ref
    rms_ref <- sqrt(mean(res_ref^2))

    results$reflect <- list(
      a = fit_ref[1], b = fit_ref[2], tx = fit_ref[3], ty = fit_ref[4],
      rms = rms_ref, reflect = TRUE
    )
  }

  # Pick the best
  if (allow_reflection && results$reflect$rms < results$no_reflect$rms) {
    best <- results$reflect
  } else {
    best <- results$no_reflect
  }

  scale     <- sqrt(best$a^2 + best$b^2)
  rot_rad   <- atan2(best$b, best$a)
  rot_deg   <- rot_rad * 180 / pi

  M <- build_transform_matrix(best$a, best$b, best$tx, best$ty, best$reflect)

  list(
    matrix       = M,
    a            = best$a,
    b            = best$b,
    tx           = best$tx,
    ty           = best$ty,
    scale        = scale,
    rotation_deg = rot_deg,
    reflected    = best$reflect,
    residual_rms = best$rms
  )
}

#' Apply a 3x3 homogeneous transform to 2D points
#'
#' @param x Numeric vector of x-coordinates
#' @param y Numeric vector of y-coordinates
#' @param M 3x3 homogeneous transform matrix
#' @return Data frame with columns x_transformed, y_transformed
apply_transform_points <- function(x, y, M) {
  pts <- rbind(x, y, rep(1, length(x)))  # 3 x n
  result <- M %*% pts                     # 3 x n
  data.frame(
    x_transformed = result[1, ],
    y_transformed = result[2, ]
  )
}

#' Extract human-readable parameters from a 3x3 similarity transform matrix
#' @param M 3x3 homogeneous similarity transform matrix
#' @return List with scale, rotation_deg, tx, ty, reflected
extract_transform_params <- function(M) {
  a  <- M[1, 1]
  b  <- M[2, 1]
  tx <- M[1, 3]
  ty <- M[2, 3]

  # Check reflection: det of upper-left 2x2
  det_ul <- M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
  reflected <- det_ul < 0

  scale   <- sqrt(a^2 + b^2)
  rot_deg <- atan2(b, a) * 180 / pi

  list(
    scale        = scale,
    rotation_deg = rot_deg,
    tx           = tx,
    ty           = ty,
    reflected    = reflected
  )
}


#' Apply a multiple-of-90° in-plane rotation to centred (x, y) coordinate vectors
#'
#' Positive degrees = counter-clockwise; negative = clockwise (y-up convention).
#' Only exact multiples of 90 are accepted — no free-angle rotations.
#'
#' @param x,y Numeric vectors (centred, same length).
#' @param deg Integer. Must be in {0, 90, -90, 180}.
#' @return List with rotated \code{x} and \code{y} vectors.
rotate_coords_90 <- function(x, y, deg) {
  d <- as.integer(round(deg)) %% 360L
  if (d < 0L) d <- d + 360L
  switch(as.character(d),
    "0"   = list(x =  x, y =  y),
    "90"  = list(x = -y, y =  x),   # counter-clockwise
    "180" = list(x = -x, y = -y),
    "270" = list(x =  y, y = -x),   # equivalent to -90 / clockwise
    stop("rotate_coords_90: deg must be a multiple of 90, got: ", deg)
  )
}


# ---------------------------------------------------------------------------
# Image type detection and canonicalization (magick-based)
# ---------------------------------------------------------------------------

#' Detect image format from magic bytes
#'
#' Reads first 12 bytes of the file and returns the detected format
#' regardless of file extension.  Prevents errors when JPEG files carry
#' a .png extension (a known Agilent LDIR export quirk).
#'
#' @param path File path to inspect
#' @return Character: "PNG", "JPEG", "TIFF", "BMP", "WEBP", or "unknown"
guess_image_type <- function(path) {
  if (!file.exists(path)) return("unknown")
  hdr <- tryCatch(
    as.integer(readBin(path, "raw", n = 12)),
    error = function(e) integer(0)
  )
  if (length(hdr) < 4) return("unknown")

  # PNG: 0x89 P N G \r \n 0x1A \n
  if (hdr[1] == 137 && hdr[2] == 80 && hdr[3] == 78 && hdr[4] == 71)
    return("PNG")
  # JPEG: FF D8 FF
  if (hdr[1] == 255 && hdr[2] == 216 && hdr[3] == 255)
    return("JPEG")
  # TIFF little-endian (II) or big-endian (MM)
  if ((hdr[1] == 73 && hdr[2] == 73 && hdr[3] == 42 && hdr[4] == 0) ||
      (hdr[1] == 77 && hdr[2] == 77 && hdr[3] == 0  && hdr[4] == 42))
    return("TIFF")
  # BMP: BM
  if (hdr[1] == 66 && hdr[2] == 77) return("BMP")
  # WEBP: RIFF....WEBP
  if (length(hdr) >= 12 &&
      hdr[1] == 82 && hdr[2] == 73 && hdr[3] == 70 && hdr[4] == 70 &&
      hdr[9] == 87 && hdr[10] == 69 && hdr[11] == 66 && hdr[12] == 80)
    return("WEBP")

  "unknown"
}


#' Read any supported image via magick and return an R array (0–1 range)
#'
#' Works for PNG, JPEG, TIFF, BMP, WEBP regardless of file extension.
#' The returned array has the same layout as png::readPNG:
#'   dim = c(height, width, channels)  channels = 3 (RGB) or 4 (RGBA)
#'
#' If magick is unavailable the function falls back to png::readPNG /
#' jpeg::readJPEG based on the detected type (returns NULL on failure).
#'
#' @param path File path
#' @param verbose Emit log_message on detected format mismatch
#' @return Numeric array [0,1] or NULL on failure
read_image_any <- function(path, verbose = TRUE) {
  if (!file.exists(path)) {
    log_message("  read_image_any: file not found: ", path, level = "WARN")
    return(NULL)
  }

  detected <- guess_image_type(path)
  ext_type  <- toupper(tools::file_ext(path))

  if (verbose && detected != "unknown" && detected != ext_type) {
    log_message("  Image signature mismatch: extension says ", ext_type,
                " but magic bytes say ", detected, " — using magick",
                level = "WARN")
  }

  # Preferred path: magick (handles all formats transparently)
  if (requireNamespace("magick", quietly = TRUE)) {
    tryCatch({
      img_mg <- magick::image_read(path)
      if (length(img_mg) > 1) img_mg <- img_mg[1]
      info   <- magick::image_info(img_mg)
      log_message("  magick: ", info$format, " ", info$width, "x", info$height)

      # Use sRGB to preserve display-gamma colours (not linear RGB).
      # Keep alpha channel so images with transparent areas (e.g. FTIR RGBA PNGs)
      # remain transparent rather than being composited against black.
      img_rgb <- magick::image_convert(img_mg, colorspace = "sRGB")
      # Export as RGBA bitmap (4 channels; images without alpha get A=255)
      raw_data <- magick::image_data(img_rgb, channels = "rgba")
      vals <- if (is.character(raw_data)) strtoi(raw_data, base = 16L)
              else as.integer(raw_data)
      arr <- array(vals / 255, dim = dim(raw_data))
      arr <- aperm(arr, c(3, 2, 1))   # [4, W, H] -> [H, W, 4] (RGBA raster)
      return(arr)
    }, error = function(e) {
      log_message("  magick failed (", e$message, "), trying pkg fallback",
                  level = "WARN")
    })
  }

  # Fallback: dependency-free BMP reader (R/read_bmp.R) — instrument PCs
  # often export BMP, and magick may be unavailable on the analysis machine.
  if (detected == "BMP" && exists("read_bmp_raster")) {
    arr <- tryCatch(read_bmp_raster(path), error = function(e) NULL)
    if (!is.null(arr)) return(arr)
    log_message("  read_bmp_raster could not decode ", basename(path),
                " (compressed or unusual BMP variant)", level = "WARN")
  }

  # Fallback: pkg-specific readers
  if (detected == "JPEG" || ext_type %in% c("JPG", "JPEG")) {
    if (requireNamespace("jpeg", quietly = TRUE))
      return(tryCatch(jpeg::readJPEG(path), error = function(e) NULL))
  }
  if (requireNamespace("png", quietly = TRUE))
    return(tryCatch(png::readPNG(path),  error = function(e) NULL))

  NULL
}


#' Read an image and return a normalised [H, W, 3] uint8 RGB array
#'
#' Wraps read_image_any() and guarantees the output is an integer array with
#' dimensions [H, W, 3] (rows, cols, RGB channels) in the 0-255 range.
#' RGBA images are composited over a white background (alpha pre-multiplied).
#' Grayscale images are replicated to three channels.
#'
#' @param path    Image file path
#' @param verbose Passed to read_image_any()
#' @return List with img_rgb ([H,W,3] integer), width, height, format; or NULL on failure
read_image_canonical <- function(path, verbose = FALSE) {
  arr <- read_image_any(path, verbose = verbose)   # [H, W, C], values 0-1 float
  if (is.null(arr)) return(NULL)

  # Ensure 3-dim array even for a 2-D grayscale result
  if (length(dim(arr)) < 3) {
    arr <- array(arr, dim = c(nrow(arr), ncol(arr), 1L))
  }

  h    <- dim(arr)[1]
  w    <- dim(arr)[2]
  n_ch <- dim(arr)[3]

  if (n_ch == 1L) {
    base <- round(arr[,,1] * 255)
    rgb  <- array(0L, dim = c(h, w, 3L))
    rgb[,,1] <- base; rgb[,,2] <- base; rgb[,,3] <- base
  } else if (n_ch == 3L) {
    rgb <- round(arr[,,1:3, drop = FALSE] * 255)
  } else {
    # n_ch >= 4: composite over white background
    alpha <- arr[,,4]
    rgb   <- array(0, dim = c(h, w, 3L))
    for (k in 1:3) rgb[,,k] <- round((arr[,,k] * alpha + (1 - alpha)) * 255)
  }
  storage.mode(rgb) <- "integer"

  stopifnot(
    length(rgb) == h * w * 3L,
    dim(rgb)[1] == h, dim(rgb)[2] == w, dim(rgb)[3] == 3L
  )

  if (isTRUE(verbose)) {
    log_message("  read_image_canonical: dim=", paste(dim(rgb), collapse = "x"),
                " range=", paste(range(rgb), collapse = ".."))
  }

  list(img_rgb = rgb, width = w, height = h,
       format = guess_image_type(path))
}


#' Downsample a [H, W, 3] uint8 array to at most max_dim on its longest side
#'
#' Uses magick for quality-preserving resizing when available; falls back to
#' nearest-neighbour sub-sampling.
#'
#' @param img_rgb [H, W, 3] integer array (0-255)
#' @param max_dim Maximum size of the longest dimension
#' @return Downsampled [H2, W2, 3] integer array (or img_rgb unchanged if small enough)
make_preview_image <- function(img_rgb, max_dim = 1200L) {
  h <- nrow(img_rgb); w <- ncol(img_rgb)
  if (max(h, w) <= max_dim) return(img_rgb)
  scale <- max_dim / max(h, w)
  nh <- as.integer(round(h * scale))
  nw <- as.integer(round(w * scale))

  if (requireNamespace("magick", quietly = TRUE)) {
    tryCatch({
      # magick::image_read accepts [H,W,3] integer 0-255 directly
      mg  <- magick::image_read(img_rgb)
      mg  <- magick::image_scale(mg, paste0(nw, "x", nh, "!"))
      raw <- magick::image_data(
               magick::image_convert(mg, colorspace = "sRGB"),
               channels = "rgb")
      vals <- if (is.character(raw)) strtoi(raw, base = 16L) else as.integer(raw)
      out  <- array(vals, dim = c(3L, nw, nh))
      out  <- aperm(out, c(3L, 2L, 1L))        # [3,W,H] -> [H,W,3]
      storage.mode(out) <- "integer"
      return(out)
    }, error = function(e) NULL)
  }

  # Fallback: nearest-neighbour
  row_idx <- as.integer(round(seq(1, h, length.out = nh)))
  col_idx <- as.integer(round(seq(1, w, length.out = nw)))
  img_rgb[row_idx, col_idx, , drop = FALSE]
}


#' Warn (or stop) when an image appears to be a tiled / repeated montage
#'
#' Compares mean and SD of luminance in each quadrant.  If at least two pairs
#' of quadrants are near-identical the image is likely a repeated mosaic layout
#' rather than a single continuous scan field.
#'
#' @param img_rgb [H, W, 3] integer array (0-255), or a 2-D grayscale array
#' @param tag     Short name used in the warning message (e.g. basename)
#' @param strict  If TRUE, stop() instead of log a warning
#' @return Invisibly: number of near-identical quadrant pairs found (0-6)
assert_not_tiled_montage <- function(img_rgb, tag, strict = FALSE) {
  if (is.null(img_rgb) || length(dim(img_rgb)) < 2) return(invisible(0L))
  h <- nrow(img_rgb); w <- ncol(img_rgb)

  lum <- if (length(dim(img_rgb)) == 3) {
    (img_rgb[,,1] + img_rgb[,,2] + img_rgb[,,3]) / 3
  } else {
    img_rgb
  }

  quad <- list(
    lum[seq_len(h %/% 2L),        seq_len(w %/% 2L)],
    lum[seq_len(h %/% 2L),        seq(w %/% 2L + 1L, w)],
    lum[seq(h %/% 2L + 1L, h),    seq_len(w %/% 2L)],
    lum[seq(h %/% 2L + 1L, h),    seq(w %/% 2L + 1L, w)]
  )
  stats_q <- lapply(quad, function(q)
    c(mean(q, na.rm = TRUE), stats::sd(q, na.rm = TRUE)))

  pairs <- utils::combn(4L, 2L, simplify = FALSE)
  n_near <- sum(vapply(pairs, function(ij) {
    abs(stats_q[[ij[1]]][1] - stats_q[[ij[2]]][1]) < 0.03 * 255 &&
    abs(stats_q[[ij[1]]][2] - stats_q[[ij[2]]][2]) < 0.03 * 255
  }, logical(1L)))

  if (n_near >= 2L) {
    msg <- paste0("Image '", tag, "' appears tiled/repeated: ", n_near,
                  " quadrant pairs are near-identical. Check LDIR image source.")
    if (strict) stop(msg) else log_message("  WARN: ", msg, level = "WARN")
  }
  invisible(n_near)
}


#' Canonicalize an LDIR image to PNG and create a preview
#'
#' Canonicalize an instrument image: copy original, write lossless PNG, generate preview.
#'
#' Copies the original file, writes a canonical PNG (lossless), and a
#' downscaled preview for the Shiny viewer.  Returns a list of provenance
#' info suitable for inclusion in the run manifest.
#'
#' @param instrument Instrument name: "ftir", "raman", or "ldir"
#' @param src_path Source image path (any format)
#' @param inputs_dir Destination directory (output/<run>/inputs/)
#' @param max_preview_px Maximum dimension (width or height) of the preview
#' @return Named list: detected_format, orig_width, orig_height, md5,
#'   canonical_path, canonical_width, canonical_height,
#'   preview_path, preview_width, preview_height, preview_scale
canonicalize_instrument_image <- function(src_path, inputs_dir, instrument,
                                          max_preview_px = 2000L) {
  if (!dir.exists(inputs_dir)) dir.create(inputs_dir, recursive = TRUE)

  prefix   <- paste0(tolower(instrument), "_image_")
  detected <- guess_image_type(src_path)
  ext_orig <- tools::file_ext(src_path)
  if (nchar(ext_orig) == 0) ext_orig <- tolower(detected)
  inst <- tolower(instrument)

  result <- list(
    instrument       = instrument,
    detected_format  = detected,
    orig_path        = normalizePath(src_path, mustWork = FALSE),
    orig_basename    = basename(src_path),
    md5              = file_md5(src_path)
  )

  # Without magick, canonicalize with the native readers instead of skipping:
  # a skipped instrument gets no entry in manifest$image_assets and no PNGs in
  # inputs/, so the Shiny viewer never shows its image.  read_image_canonical()
  # handles PNG/JPEG via png/jpeg and BMP via read_bmp_raster (R/read_bmp.R).
  if (!requireNamespace("magick", quietly = TRUE)) {
    log_message("  magick not available; canonicalizing ", instrument,
                " image with native R readers", level = "WARN")
    if (!requireNamespace("png", quietly = TRUE)) {
      log_message("  png package unavailable; skipping ", instrument,
                  " image canonicalization", level = "WARN")
      return(result)
    }
    tryCatch({
      ci <- read_image_canonical(src_path)
      if (is.null(ci)) {
        log_message("  [", instrument, "] no native reader for this ",
                    detected, " image — skipping canonicalization ",
                    "(install magick for TIFF/WEBP/compressed-BMP support)",
                    level = "WARN")
        return(result)
      }
      result$orig_width  <- ci$width
      result$orig_height <- ci$height

      orig_dest <- file.path(inputs_dir,
                             paste0(inst, "_image_original.", tolower(ext_orig)))
      file.copy(src_path, orig_dest, overwrite = TRUE)
      result$orig_dest <- normalizePath(orig_dest, mustWork = FALSE)

      canon_path <- file.path(inputs_dir, paste0(inst, "_image_canonical.png"))
      if (identical(detected, "PNG")) {
        file.copy(src_path, canon_path, overwrite = TRUE)
      } else {
        png::writePNG(ci$img_rgb / 255, canon_path)
      }
      result$canonical_path   <- canon_path
      result$canonical_width  <- ci$width
      result$canonical_height <- ci$height
      log_message("  [", instrument, "] Canonical PNG (native): ", canon_path,
                  " (", ci$width, "x", ci$height, ")")

      prev <- make_preview_image(ci$img_rgb, max_dim = max_preview_px)
      prev_path <- file.path(inputs_dir, paste0(inst, "_image_preview.png"))
      png::writePNG(prev / 255, prev_path)
      result$preview_path   <- prev_path
      result$preview_width  <- ncol(prev)
      result$preview_height <- nrow(prev)
      result$preview_scale  <- ncol(prev) / ci$width
      log_message("  [", instrument, "] Preview PNG (native): ", prev_path,
                  " (", ncol(prev), "x", nrow(prev), ")")
    }, error = function(e) {
      log_message("  canonicalize_instrument_image [", instrument,
                  "] native fallback error: ", e$message, level = "WARN")
    })
    return(result)
  }

  tryCatch({
    img_mg <- magick::image_read(src_path)
    info   <- magick::image_info(img_mg)

    result$orig_width    <- info$width
    result$orig_height   <- info$height
    result$magick_format <- info$format

    # --- Copy original under inputs/ with detected extension ---
    orig_dest <- file.path(inputs_dir,
                           paste0(inst, "_image_original.", tolower(ext_orig)))
    file.copy(src_path, orig_dest, overwrite = TRUE)
    result$orig_dest <- normalizePath(orig_dest, mustWork = FALSE)

    # --- Write canonical PNG ---
    canon_path <- file.path(inputs_dir, paste0(inst, "_image_canonical.png"))
    if (identical(detected, "PNG")) {
      # Lossless: byte-identical copy avoids magick re-encode / colour-management drift
      file.copy(src_path, canon_path, overwrite = TRUE)
    } else {
      magick::image_write(img_mg, path = canon_path, format = "png")
    }
    canon_info <- magick::image_info(magick::image_read(canon_path))
    result$canonical_path   <- canon_path
    result$canonical_width  <- canon_info$width
    result$canonical_height <- canon_info$height
    log_message("  [", instrument, "] Canonical PNG: ", canon_path,
                " (", canon_info$width, "x", canon_info$height, ")")

    # --- Preview (downscaled) ---
    max_dim <- max(info$width, info$height)
    if (max_dim > max_preview_px) {
      scale_pct <- round(max_preview_px / max_dim * 100)
      prev_mg   <- magick::image_scale(img_mg, paste0(scale_pct, "%"))
    } else {
      prev_mg   <- img_mg
    }
    prev_path <- file.path(inputs_dir, paste0(inst, "_image_preview.png"))
    magick::image_write(prev_mg, path = prev_path, format = "png")
    prev_info <- magick::image_info(prev_mg)
    result$preview_path   <- prev_path
    result$preview_width  <- prev_info$width
    result$preview_height <- prev_info$height
    result$preview_scale  <- prev_info$width / info$width
    log_message("  [", instrument, "] Preview PNG: ", prev_path,
                " (", prev_info$width, "x", prev_info$height, ")")

  }, error = function(e) {
    log_message("  canonicalize_instrument_image [", instrument, "] error: ",
                e$message, level = "WARN")
  })

  result
}


canonicalize_ldir_image <- function(src_path, inputs_dir,
                                     max_preview_px = 2000L) {
  canonicalize_instrument_image(
    src_path = src_path,
    inputs_dir = inputs_dir,
    instrument = "ldir",
    max_preview_px = max_preview_px
  )
}


#' Compute MD5 hash of a file
#'
#' @param path File path
#' @return Hex MD5 string, or NA if file doesn't exist
file_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  tryCatch(
    as.character(tools::md5sum(path)),
    error = function(e) NA_character_
  )
}


# ---------------------------------------------------------------------------
# Run manifest (provenance)
# ---------------------------------------------------------------------------

#' Write a run manifest JSON file
#'
#' Resolve manifest.json path: prefer 00_manifest/ subfolder, fall back to flat
#' @param run_dir Run output directory
#' @return Path to manifest.json (may not exist yet)
resolve_manifest_path <- function(run_dir) {
  staged <- file.path(run_dir, "00_manifest", "manifest.json")
  if (file.exists(staged)) return(staged)
  legacy <- file.path(run_dir, "manifest.json")
  if (file.exists(legacy)) return(legacy)
  # Default to staged location for new writes
  staged
}

#' Records authoritative provenance for every pipeline run: run ID,
#' timestamp, git commit, R session, config snapshot, input file hashes,
#' and instrument image details.
#'
#' The manifest is written early (at run start) and updated by
#' update_manifest_stage() as the pipeline progresses.
#'
#' @param run_dir     Run output directory
#' @param run_id      Character run ID (e.g. "2026-02-27_1")
#' @param config      Config list from make_config()
#' @param input_paths Named list of input file paths (ftir, raman, ldir, ...)
#' @param images_info Named list of canonicalize_instrument_image() results
#'   keyed by instrument ("ftir", "raman", "ldir"). Supersedes ldir_image_info.
#' @param ldir_image_info Deprecated; use images_info$ldir instead
#' @param stage       Current pipeline stage label (default "started")
#' @return Invisible path to the manifest file
write_manifest <- function(run_dir, run_id, config,
                            input_paths = list(),
                            images_info = list(),
                            ldir_image_info = NULL,
                            image_infos = list(),
                            stage = "started") {
  # Backward compat: merge ldir_image_info into images_info if needed
  if (!is.null(ldir_image_info) && is.null(images_info$ldir)) {
    images_info$ldir <- ldir_image_info
  }
  manifest_path <- resolve_manifest_path(run_dir)
  manifest_dir <- dirname(manifest_path)
  if (!dir.exists(manifest_dir)) dir.create(manifest_dir, recursive = TRUE)

  # --- Git commit ---
  git_commit <- tryCatch(
    trimws(system("git rev-parse HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) NA_character_
  )
  if (length(git_commit) == 0 || startsWith(git_commit, "fatal")) {
    git_commit <- NA_character_
  }

  # --- Input file info ---
  image_names <- c("ftir_image", "raman_image", "ldir_image")
  detect_input_dims <- function(path) {
    if (!file.exists(path) || !requireNamespace("magick", quietly = TRUE)) {
      return(list(width = NA_integer_, height = NA_integer_))
    }
    tryCatch({
      info <- magick::image_info(magick::image_read(path))
      list(width = as.integer(info$width), height = as.integer(info$height))
    }, error = function(e) list(width = NA_integer_, height = NA_integer_))
  }

  inputs_info <- lapply(names(input_paths), function(nm) {
    p <- input_paths[[nm]]
    if (is.null(p) || !nzchar(p)) return(list(name = nm, path = NULL))
    abs_path <- tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE),
                         error = function(e) p)
    ext <- toupper(tools::file_ext(p))
    fmt <- if (nm %in% image_names) guess_image_type(p) else if (nzchar(ext)) ext else "unknown"
    dims <- if (nm %in% image_names) detect_input_dims(p) else list(width = NA_integer_, height = NA_integer_)
    list(
      name       = nm,
      path       = abs_path,
      basename   = basename(p),
      md5        = file_md5(p),
      size_bytes = if (file.exists(p)) file.info(p)$size else NA_integer_,
      format     = fmt,
      width      = dims$width,
      height     = dims$height
    )
  })
  names(inputs_info) <- names(input_paths)

  image_assets <- list()
  for (nm in names(image_infos)) {
    info_obj <- image_infos[[nm]]
    if (is.null(info_obj) || is.null(info_obj$orig_dest)) next
    image_assets[[nm]] <- list(
      original = list(
        path = normalizePath(info_obj$orig_dest, winslash = "/", mustWork = FALSE),
        md5 = file_md5(info_obj$orig_dest),
        format = if (!is.null(info_obj$detected_format)) info_obj$detected_format else "unknown",
        width = info_obj$orig_width,
        height = info_obj$orig_height
      ),
      canonical = list(
        path = if (!is.null(info_obj$canonical_path)) normalizePath(info_obj$canonical_path, winslash = "/", mustWork = FALSE) else NULL,
        md5 = if (!is.null(info_obj$canonical_path)) file_md5(info_obj$canonical_path) else NA_character_,
        format = "PNG",
        width = info_obj$canonical_width,
        height = info_obj$canonical_height
      ),
      preview = list(
        path = if (!is.null(info_obj$preview_path)) normalizePath(info_obj$preview_path, winslash = "/", mustWork = FALSE) else NULL,
        md5 = if (!is.null(info_obj$preview_path)) file_md5(info_obj$preview_path) else NA_character_,
        format = "PNG",
        width = info_obj$preview_width,
        height = info_obj$preview_height
      )
    )
  }

  # Build clean images dict for provenance panel (keys: "ftir", "raman", "ldir")
  images_clean <- list()
  for (nm in names(image_infos)) {
    clean_nm <- gsub("_image$", "", nm)
    if (!is.null(image_infos[[nm]])) images_clean[[clean_nm]] <- image_infos[[nm]]
  }

  # --- Config snapshot (key values only) ---
  cfg_keys <- c("ldir_scan_diameter_um", "ldir_flip_y_for_alignment",
                 "min_quality_ftir", "min_quality_raman", "min_size_um",
                 "ransac_iterations", "icp_max_iter", "match_radius_um",
                 # Raman image placement — the Shiny viewer reads these from
                 # the manifest snapshot to place the background image at its
                 # exact physical extent (raman_native_image_info Priority 1).
                 "raman_image_width_um", "raman_image_height_um",
                 "raman_image_center_x_um", "raman_image_center_y_um",
                 "raman_um_per_px")
  cfg_snap <- lapply(cfg_keys, function(k) config[[k]])
  names(cfg_snap) <- cfg_keys
  # Remove NULLs
  cfg_snap <- cfg_snap[!vapply(cfg_snap, is.null, logical(1))]

  manifest <- list(
    run_id          = run_id,
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stage           = stage,
    git_commit      = git_commit,
    r_version       = paste(R.version$major, R.version$minor, sep = "."),
    platform        = R.version$platform,
    user            = Sys.info()[["user"]],
    config_snapshot = cfg_snap,
    inputs          = inputs_info,
    image_assets    = image_assets,
    images          = if (length(images_clean) > 0) images_clean else NULL,
    ldir_image      = ldir_image_info
  )

  # Serialize to JSON (using jsonlite if available, else simple fallback)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    json_str <- jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE,
                                  null = "null", na = "null")
  } else {
    # Minimal fallback: key=value style text (not strict JSON but readable)
    json_str <- paste(
      "{",
      paste0('  "run_id": "', manifest$run_id, '",'),
      paste0('  "timestamp": "', manifest$timestamp, '",'),
      paste0('  "stage": "', manifest$stage, '",'),
      paste0('  "git_commit": "', ifelse(is.na(manifest$git_commit), "", manifest$git_commit), '"'),
      "}",
      sep = "\n"
    )
  }

  writeLines(json_str, manifest_path)
  log_message("  Manifest written: ", manifest_path,
              " (stage=", stage, ", git=",
              ifelse(is.na(git_commit), "N/A",
                     substr(git_commit, 1, 8)), ")")
  invisible(manifest_path)
}


#' Update the stage field in an existing manifest
#'
#' @param run_dir   Run output directory
#' @param stage     New stage label
#' @param error_msg Optional error message if pipeline failed
update_manifest_stage <- function(run_dir, stage, error_msg = NULL) {
  manifest_path <- resolve_manifest_path(run_dir)
  if (!file.exists(manifest_path)) return(invisible(NULL))

  tryCatch({
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      m <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
      m$stage <- stage
      if (!is.null(error_msg)) m$error <- error_msg
      m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
      writeLines(jsonlite::toJSON(m, pretty = TRUE, auto_unbox = TRUE,
                                   null = "null", na = "null"),
                 manifest_path)
    }
  }, error = function(e) {
    log_message("  update_manifest_stage failed: ", e$message, level = "WARN")
  })
  invisible(manifest_path)
}


#' Persist LDIR circle calibration in manifest
#'
#' Reads the run manifest, adds/replaces the `ldir_circle` key with the circle
#' calibration values needed by the Shiny viewer to compute correct image bounds,
#' then rewrites the manifest.
#'
#' @param run_dir Path to the run output directory containing manifest.json
#' @param circle_info Named list with cx_px, cy_px, radius_px, scale_um_per_px,
#'   width (image_width_px), height (image_height_px)
update_manifest_ldir_circle <- function(run_dir, circle_info) {
  manifest_path <- resolve_manifest_path(run_dir)
  if (!file.exists(manifest_path)) return(invisible(NULL))
  tryCatch({
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      m <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
      m$ldir_circle <- list(
        cx_px           = circle_info$cx_px,
        cy_px           = circle_info$cy_px,
        radius_px       = circle_info$radius_px,
        scale_um_per_px = circle_info$scale_um_per_px,
        image_width_px  = circle_info$width,
        image_height_px = circle_info$height,
        method          = circle_info$method %||% "unknown",
        detected        = isTRUE(circle_info$detected)
      )
      m$coord_frames <- list(
        ldir_native_px = list(
          origin  = "image_top_left",
          x_dir   = "right",
          y_dir   = "down",
          unit    = "pixels"
        ),
        ldir_centered_um = list(
          origin  = "scan_circle_center",
          x_dir   = "right",
          y_dir   = "up",
          unit    = "micrometres",
          formula = list(
            x_um = "(x_px - cx_px) * um_per_px",
            y_um = "(cy_px - y_px) * um_per_px"
          )
        )
      )
      writeLines(jsonlite::toJSON(m, pretty = TRUE, auto_unbox = TRUE,
                                   null = "null", na = "null"),
                 manifest_path)
    }
  }, error = function(e) {
    log_message("  update_manifest_ldir_circle failed: ", e$message, level = "WARN")
  })
  invisible(manifest_path)
}


# ---------------------------------------------------------------------------
# Encoding safety for column names
# ---------------------------------------------------------------------------

#' Sanitize column names to valid UTF-8
#'
#' Instrument exports on Windows often use latin1 encoding for special
#' characters like µ (\xb5), ² (\xb2), ³ (\xb3).  When R's locale
#' expects UTF-8 these appear as invalid multibyte strings, crashing
#' tolower(), grep(), and other string functions.
#'
#' Strategy: try system-default → UTF-8 first, then latin1 → UTF-8.
#' Pick whichever produces fewer NA/empty strings.  Any bytes that still
#' cannot convert are silently dropped (sub = "").
#'
#' @param x Character vector (typically column names)
#' @return Character vector with valid UTF-8 encoding
safe_colnames <- function(x) {
  # Attempt 1: system default → UTF-8
  x_sys <- tryCatch(
    iconv(x, from = "", to = "UTF-8", sub = ""),
    error = function(e) rep(NA_character_, length(x))
  )

  # Attempt 2: latin1 → UTF-8  (most common for Spotlight / Agilent exports)
  x_lat <- tryCatch(
    iconv(x, from = "latin1", to = "UTF-8", sub = ""),
    error = function(e) rep(NA_character_, length(x))
  )

  # Score: fewer NA/empty = better
  score <- function(v) sum(is.na(v) | nchar(v) == 0)
  result <- if (score(x_lat) < score(x_sys)) x_lat else x_sys

  # Last-resort: if still NA, keep original bytes (best effort)
  still_na <- is.na(result)
  if (any(still_na)) result[still_na] <- x[still_na]

  result
}


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log_message <- function(..., level = "INFO") {
  msg <- paste0("[", Sys.time(), "] [", level, "] ", paste0(..., collapse = ""))
  message(msg)
}

# ---------------------------------------------------------------------------
# Debug: A3 particle trace (Step 5)
# ---------------------------------------------------------------------------

#' Dump a single particle's coordinates at one pipeline stage to a CSV
#'
#' Writes one row per particle_id in \code{ids} to
#'   \code{<debug_dir>/<id>_stage_<stage>.csv}
#' If the particle is not found a "not_found" note is written instead.
#' Every write is verified with stopifnot(file.exists()).
#'
#' @param df Data frame at the current pipeline stage
#' @param ids Character vector of particle_id values to dump (e.g. c("A3","MP_11"))
#' @param stage Character label (used in filename and "stage" column)
#' @param debug_dir Path to debug directory
dump_particle <- function(df, ids, stage, debug_dir) {
  if (is.null(debug_dir) || !dir.exists(debug_dir)) return(invisible(NULL))

  all_cols <- c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
                "x_aligned", "y_aligned", "coord_source", "material",
                "feret_max_um", "quality")

  for (pid in ids) {
    row_idx <- which(df$particle_id == pid)
    out_path <- file.path(debug_dir,
                          paste0(gsub("[^A-Za-z0-9_-]", "_", pid),
                                 "_stage_", stage, ".csv"))

    if (length(row_idx) == 0) {
      row_df <- data.frame(
        particle_id = pid, stage = stage,
        note = "not_found_at_this_stage",
        stringsAsFactors = FALSE
      )
    } else {
      row_df <- df[row_idx[1], intersect(all_cols, names(df)), drop = FALSE]
      # Fill missing coord columns with NA
      for (col in all_cols) {
        if (!col %in% names(row_df)) row_df[[col]] <- NA
      }
      row_df$stage <- stage
    }

    tryCatch({
      write.csv(row_df, out_path, row.names = FALSE)
      stopifnot(file.exists(out_path))
    }, error = function(e) {
      log_message("  WARN: dump_particle write failed for '", pid,
                  "' stage '", stage, "': ", e$message, level = "WARN")
    })
  }

  log_message("  Dumped stage '", stage, "' for: ",
              paste(ids, collapse = ", "))
}


#' Trace a particle through pipeline stages
#'
#' Appends a snapshot row to debug_A3_trace.csv for each pipeline stage.
#' Traces ALL particles so the CSV can be filtered to any ID later.
#'
#' @param df Data frame with particle data at a given stage
#' @param stage Character label for this pipeline stage
#' @param config Config list (needs config$debug_dir)
trace_particle_snapshot <- function(df, stage, config) {
  if (!isTRUE(config$debug) || is.null(config$debug_dir)) return(invisible(NULL))

  trace_file <- file.path(config$debug_dir, "debug_A3_trace.csv")

  # Extract available coordinate columns
  cols_available <- intersect(
    c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
      "x_aligned", "y_aligned", "coord_source", "material",
      "feret_max_um"),
    names(df)
  )

  snapshot <- df[, cols_available, drop = FALSE]
  snapshot$stage <- stage

  # Fill missing columns with NA
  for (col in c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
                "x_aligned", "y_aligned", "coord_source", "material",
                "feret_max_um")) {
    if (!col %in% names(snapshot)) snapshot[[col]] <- NA
  }

  tryCatch({
    if (file.exists(trace_file)) {
      write.table(snapshot, trace_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    } else {
      write.csv(snapshot, trace_file, row.names = FALSE)
    }
    stopifnot(file.exists(trace_file))
  }, error = function(e) {
    log_message("  WARN: trace_particle_snapshot write failed: ", e$message,
                level = "WARN")
  })

  log_message("  Traced ", nrow(snapshot), " particles at stage: ", stage)
}


# ---------------------------------------------------------------------------
# Debug: Branch A/B comparison + residual vectors (Steps 2 & 6)
# ---------------------------------------------------------------------------

#' Run Branch A/B Y-flip comparison and save debug artifacts
#'
#' Generates overlay plots for both Y-flip settings and computes quality
#' metrics (ICP RMS, median match distance) for each.
#'
#' @param ldir_aligned LDIR data frame with x_aligned, y_aligned
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param ldir_raman_match Match result for LDIR-Raman
#' @param ldir_icp ICP result for LDIR alignment
#' @param config Configuration list
debug_ldir_branches <- function(ldir_aligned, raman_df, ldir_raman_match,
                                 ldir_icp, config) {
  if (!isTRUE(config$debug) || is.null(config$debug_dir)) return(invisible(NULL))

  debug_dir <- config$debug_dir

  tryCatch({
    # Current branch overlay (the active flip setting)
    branch_label <- if (isTRUE(config$ldir_flip_y_for_alignment)) "A" else "B"

    # Step 5: Write overlay_plot_data_ldir.csv so we can verify plotted columns
    overlay_csv <- file.path(debug_dir, "overlay_plot_data_ldir.csv")
    plot_cols <- intersect(c("particle_id", "x_aligned", "y_aligned",
                              "x_norm", "y_norm", "x_um", "y_um",
                              "material", "feret_max_um", "coord_source"),
                           names(ldir_aligned))
    write.csv(ldir_aligned[seq_len(min(200, nrow(ldir_aligned))), plot_cols,
                            drop = FALSE],
              overlay_csv, row.names = FALSE)
    stopifnot(file.exists(overlay_csv))
    log_message("  Debug: wrote overlay_plot_data_ldir.csv (",
                min(200, nrow(ldir_aligned)), " rows, cols: ",
                paste(plot_cols, collapse = ","), ")")

    # A3 in the final overlay (what is actually plotted)
    trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
    for (pid in trace_ids) {
      row_idx <- which(ldir_aligned$particle_id == pid)
      if (length(row_idx) > 0) {
        r <- ldir_aligned[row_idx[1], plot_cols, drop = FALSE]
        log_message("  Debug overlay coords for '", pid, "': ",
                    "x_aligned=", round(r$x_aligned, 2),
                    " y_aligned=", round(r$y_aligned, 2))
      } else {
        log_message("  Debug: '", pid, "' NOT FOUND in ldir_aligned", level = "WARN")
      }
    }

    # Save current branch overlay
    icp_rms_val <- if (length(ldir_icp$rms_history) > 0)
      round(tail(ldir_icp$rms_history, 1), 1) else NA_real_

    p_current <- plot_overlay(
      ldir_aligned, raman_df, ldir_raman_match,
      ftir_color = "darkgreen", raman_color = "steelblue",
      src_label = "ldir"
    ) + ggplot2::labs(
      title = paste0("LDIR-Raman Overlay — Branch ", branch_label,
                      " (flip_y=", config$ldir_flip_y_for_alignment, ")"),
      subtitle = paste0("Green=LDIR(x_aligned/y_aligned), Blue=Raman(x_norm/y_norm) | ",
                        "ICP RMS=", icp_rms_val, " \u00b5m")
    )
    overlay_png <- file.path(debug_dir,
                              paste0("overlay_branch", branch_label, ".png"))
    ggplot2::ggsave(overlay_png, p_current, width = 10, height = 8, dpi = 150)
    stopifnot(file.exists(overlay_png))

    # Compute metrics for current branch
    matched <- ldir_raman_match$matched
    metrics <- data.frame(
      branch = branch_label,
      flip_y = config$ldir_flip_y_for_alignment,
      icp_rms = icp_rms_val,
      n_matched = nrow(matched),
      median_match_dist = if (nrow(matched) > 0)
        round(median(matched$match_distance), 2) else NA_real_,
      stringsAsFactors = FALSE
    )

    # Write branch summary
    summary_file <- file.path(debug_dir, "branch_summary.txt")
    summary_lines <- c(
      paste0("Branch ", branch_label, " (ldir_flip_y_for_alignment = ",
             config$ldir_flip_y_for_alignment, "):"),
      paste0("  ICP RMS:             ", metrics$icp_rms, " \u00b5m"),
      paste0("  Matched pairs:       ", metrics$n_matched),
      paste0("  Median match dist:   ", metrics$median_match_dist, " \u00b5m"),
      ""
    )
    writeLines(summary_lines, summary_file)
    stopifnot(file.exists(summary_file))
    log_message("  Debug: saved Branch ", branch_label, " overlay + metrics")

    # --- Residual vectors (Step 6) ---
    if (nrow(matched) > 0) {
      src_x_col <- "ldir_x_aligned"
      src_y_col <- "ldir_y_aligned"
      ref_x_col <- "raman_x_norm"
      ref_y_col <- "raman_y_norm"

      if (all(c(src_x_col, src_y_col, ref_x_col, ref_y_col) %in% names(matched))) {
        resid_df <- data.frame(
          x = matched[[src_x_col]],
          y = matched[[src_y_col]],
          dx = matched[[ref_x_col]] - matched[[src_x_col]],
          dy = matched[[ref_y_col]] - matched[[src_y_col]],
          dist = matched$match_distance
        )
        resid_df <- resid_df[complete.cases(resid_df), ]

        if (nrow(resid_df) > 0) {
          p_resid <- ggplot2::ggplot(resid_df) +
            ggplot2::geom_segment(
              ggplot2::aes(x = x, y = y,
                           xend = x + dx * 5, yend = y + dy * 5,
                           colour = dist),
              arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")),
              linewidth = 0.5, alpha = 0.7
            ) +
            ggplot2::scale_colour_viridis_c(name = "Distance (µm)") +
            ggplot2::coord_equal() +
            ggplot2::labs(
              title = "LDIR Residual Vector Field",
              subtitle = "Arrows magnified 5x | Systematic patterns = local distortion",
              x = "X (µm)", y = "Y (µm)"
            ) +
            ggplot2::theme_minimal()

          ggplot2::ggsave(file.path(debug_dir, "ldir_residual_vectors.png"),
                          p_resid, width = 10, height = 8, dpi = 150)
          log_message("  Debug: saved ldir_residual_vectors.png")
        }
      }
    }

  }, error = function(e) {
    log_message("  Debug branch comparison failed: ", e$message, level = "WARN")
  })
}
