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

#' Compose two 3x3 transforms: apply T1 first, then T2
#' @param T1 3x3 matrix (applied first)
#' @param T2 3x3 matrix (applied second)
#' @return 3x3 matrix T2 %*% T1
compose_transforms <- function(T1, T2) {
  T2 %*% T1
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

#' Create a translation-only 3x3 matrix
make_translation_matrix <- function(tx, ty) {
  M <- diag(3)
  M[1, 3] <- tx
  M[2, 3] <- ty
  M
}

#' Create a rotation-only 3x3 matrix (rotation about origin)
#' @param angle_deg Rotation angle in degrees
make_rotation_matrix <- function(angle_deg) {
  theta <- angle_deg * pi / 180
  ct <- cos(theta)
  st <- sin(theta)
  matrix(c(ct, st, 0,
           -st, ct, 0,
           0, 0, 1), nrow = 3, byrow = FALSE)
}

#' Create a y-axis reflection matrix
make_mirror_y_matrix <- function() {
  matrix(c(1, 0, 0,
           0, -1, 0,
           0, 0, 1), nrow = 3, byrow = FALSE)
}

# ---------------------------------------------------------------------------
# Nearest-neighbor helpers
# ---------------------------------------------------------------------------

#' Find nearest neighbors using RANN
#' @param query_x, query_y Coordinates of query points
#' @param ref_x, ref_y Coordinates of reference points
#' @param k Number of neighbors to return
#' @return List with nn_idx (indices into ref) and nn_dist (distances)
find_nearest_neighbors <- function(query_x, query_y, ref_x, ref_y, k = 1) {
  query_mat <- cbind(query_x, query_y)
  ref_mat   <- cbind(ref_x, ref_y)

  nn <- RANN::nn2(data = ref_mat, query = query_mat, k = k)

  list(
    nn_idx  = nn$nn.idx,
    nn_dist = nn$nn.dists
  )
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log_message <- function(..., level = "INFO") {
  msg <- paste0("[", Sys.time(), "] [", level, "] ", paste0(..., collapse = ""))
  message(msg)
}
