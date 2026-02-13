# =============================================================================
# 04_ransac.R — RANSAC-based alignment estimation
# =============================================================================

#' Estimate the geometric transform from FTIR to Raman coordinate frame
#'
#' Two-phase approach:
#'   1. Coarse grid search over rotation angles (± mirror) to find the
#'      approximate alignment by counting spatially overlapping particles.
#'   2. RANSAC refinement: from tentative nearest-neighbor correspondences,
#'      robustly estimate the optimal similarity transform.
#'
#' @param ftir_df FTIR data frame with x_norm, y_norm (centered coordinates)
#' @param raman_df Raman data frame with x_norm, y_norm (centered coordinates)
#' @param config Configuration list (see 00_config.R)
#' @return List with:
#'   transform    — 3x3 similarity transform matrix (FTIR_norm → Raman_norm)
#'   params       — human-readable transform parameters
#'   n_inliers    — number of inlier correspondences
#'   inlier_pairs — data frame of inlier FTIR–Raman index pairs
#'   diagnostics  — list of diagnostic info (coarse scores, etc.)
ransac_align <- function(ftir_df, raman_df, config) {
  log_message("Starting RANSAC alignment")

  step_deg       <- config$ransac_coarse_step_deg
  n_ransac       <- config$ransac_n_iterations
  min_samples    <- config$ransac_min_samples
  inlier_dist    <- config$ransac_inlier_dist_um
  allow_mirror   <- config$ransac_allow_mirror

  ftir_x  <- ftir_df$x_norm
  ftir_y  <- ftir_df$y_norm
  raman_x <- raman_df$x_norm
  raman_y <- raman_df$y_norm
  n_ftir  <- length(ftir_x)
  n_raman <- length(raman_x)

  if (n_ftir < 3 || n_raman < 3) {
    stop("Need at least 3 particles in each dataset for alignment. ",
         "FTIR: ", n_ftir, ", Raman: ", n_raman)
  }

  # =========================================================================
  # Phase 1: Coarse grid search
  # =========================================================================
  log_message("  Phase 1: Coarse rotation grid search (step = ", step_deg, "°)")

  angles <- seq(0, 360 - step_deg, by = step_deg)
  mirror_opts <- if (allow_mirror) c(FALSE, TRUE) else FALSE

  best_score <- 0
  best_angle <- 0
  best_mirror <- FALSE
  coarse_scores <- data.frame(angle = numeric(), mirror = logical(),
                              n_inliers = integer())

  for (mirror in mirror_opts) {
    for (angle in angles) {
      # Build rotation (+ optional mirror) matrix
      theta <- angle * pi / 180
      ct <- cos(theta)
      st <- sin(theta)

      if (!mirror) {
        tx <- ct * ftir_x - st * ftir_y
        ty <- st * ftir_x + ct * ftir_y
      } else {
        # y-flip then rotate: reflect y, then rotate
        tx <- ct * ftir_x + st * ftir_y
        ty <- st * ftir_x - ct * ftir_y
      }

      # Count inliers (FTIR points with a Raman neighbor within threshold)
      nn <- RANN::nn2(cbind(raman_x, raman_y), cbind(tx, ty), k = 1)
      n_inliers <- sum(nn$nn.dists[, 1] <= inlier_dist)

      coarse_scores <- rbind(coarse_scores,
                             data.frame(angle = angle, mirror = mirror,
                                        n_inliers = n_inliers))

      if (n_inliers > best_score) {
        best_score  <- n_inliers
        best_angle  <- angle
        best_mirror <- mirror
      }
    }
  }

  log_message("  Coarse search best: angle = ", best_angle,
              "°, mirror = ", best_mirror,
              ", inliers = ", best_score, " / ", n_ftir)

  if (best_score < min_samples) {
    warning("Coarse alignment found very few inliers (", best_score,
            "). Results may be unreliable. Check that the datasets are from ",
            "the same physical sample.")
  }

  # =========================================================================
  # Phase 2: RANSAC refinement
  # =========================================================================
  log_message("  Phase 2: RANSAC refinement (", n_ransac, " iterations)")

  # Apply coarse alignment to get tentative correspondences
  theta_best <- best_angle * pi / 180
  ct <- cos(theta_best)
  st <- sin(theta_best)
  if (!best_mirror) {
    coarse_x <- ct * ftir_x - st * ftir_y
    coarse_y <- st * ftir_x + ct * ftir_y
  } else {
    coarse_x <- ct * ftir_x + st * ftir_y
    coarse_y <- st * ftir_x - ct * ftir_y
  }

  # Find tentative nearest-neighbor correspondences
  nn_coarse <- RANN::nn2(cbind(raman_x, raman_y), cbind(coarse_x, coarse_y), k = 1)
  tent_mask <- nn_coarse$nn.dists[, 1] <= inlier_dist * 2  # generous threshold for candidates
  tent_ftir_idx  <- which(tent_mask)
  tent_raman_idx <- nn_coarse$nn.idx[tent_mask, 1]
  n_tentative <- length(tent_ftir_idx)

  log_message("  Tentative correspondences: ", n_tentative)

  if (n_tentative < min_samples) {
    # Fall back to coarse alignment only
    log_message("  Too few tentative correspondences for RANSAC. Using coarse alignment.", level = "WARN")
    M_coarse <- build_coarse_transform(best_angle, best_mirror)
    params <- extract_transform_params(M_coarse)
    return(list(
      transform    = M_coarse,
      params       = params,
      n_inliers    = best_score,
      inlier_pairs = data.frame(ftir_idx = tent_ftir_idx,
                                raman_idx = tent_raman_idx),
      diagnostics  = list(coarse_scores = coarse_scores,
                          method = "coarse_only")
    ))
  }

  # RANSAC loop
  best_ransac_inliers <- 0
  best_ransac_M       <- NULL
  best_ransac_pairs   <- NULL

  for (iter in seq_len(n_ransac)) {
    # Sample min_samples random correspondences from tentative set
    sample_idx <- sample(n_tentative, min(min_samples, n_tentative))
    s_ftir_idx  <- tent_ftir_idx[sample_idx]
    s_raman_idx <- tent_raman_idx[sample_idx]

    # Estimate similarity transform from these pairs
    tryCatch({
      tf <- estimate_similarity_transform(
        src_x = ftir_x[s_ftir_idx],
        src_y = ftir_y[s_ftir_idx],
        dst_x = raman_x[s_raman_idx],
        dst_y = raman_y[s_raman_idx],
        allow_reflection = allow_mirror
      )

      # Apply to all FTIR points and count inliers
      transformed <- apply_transform_points(ftir_x, ftir_y, tf$matrix)
      nn_check <- RANN::nn2(cbind(raman_x, raman_y),
                            cbind(transformed$x_transformed, transformed$y_transformed),
                            k = 1)
      inlier_mask <- nn_check$nn.dists[, 1] <= inlier_dist
      n_in <- sum(inlier_mask)

      if (n_in > best_ransac_inliers) {
        best_ransac_inliers <- n_in
        best_ransac_M       <- tf$matrix
        best_ransac_pairs   <- data.frame(
          ftir_idx  = which(inlier_mask),
          raman_idx = nn_check$nn.idx[inlier_mask, 1]
        )
      }
    }, error = function(e) {
      # Skip degenerate samples (e.g., collinear points)
    })
  }

  # Refit transform using ALL inliers of the best RANSAC model
  if (!is.null(best_ransac_pairs) && nrow(best_ransac_pairs) >= min_samples) {
    final_tf <- estimate_similarity_transform(
      src_x = ftir_x[best_ransac_pairs$ftir_idx],
      src_y = ftir_y[best_ransac_pairs$ftir_idx],
      dst_x = raman_x[best_ransac_pairs$raman_idx],
      dst_y = raman_y[best_ransac_pairs$raman_idx],
      allow_reflection = allow_mirror
    )

    # Re-evaluate inliers with the refined transform
    transformed <- apply_transform_points(ftir_x, ftir_y, final_tf$matrix)
    nn_final <- RANN::nn2(cbind(raman_x, raman_y),
                          cbind(transformed$x_transformed, transformed$y_transformed),
                          k = 1)
    inlier_mask <- nn_final$nn.dists[, 1] <= inlier_dist
    final_pairs <- data.frame(
      ftir_idx  = which(inlier_mask),
      raman_idx = nn_final$nn.idx[inlier_mask, 1]
    )

    log_message("  RANSAC result: ", sum(inlier_mask), " inliers, ",
                "scale = ", round(final_tf$scale, 4),
                ", rotation = ", round(final_tf$rotation_deg, 2), "°",
                ", reflected = ", final_tf$reflected)

    params <- extract_transform_params(final_tf$matrix)

    return(list(
      transform    = final_tf$matrix,
      params       = params,
      n_inliers    = sum(inlier_mask),
      inlier_pairs = final_pairs,
      diagnostics  = list(coarse_scores = coarse_scores,
                          ransac_best_inliers = best_ransac_inliers,
                          method = "ransac_refined")
    ))
  }

  # Fallback: use coarse alignment
  log_message("  RANSAC did not improve over coarse. Using coarse alignment.", level = "WARN")
  M_coarse <- build_coarse_transform(best_angle, best_mirror)
  params <- extract_transform_params(M_coarse)

  list(
    transform    = M_coarse,
    params       = params,
    n_inliers    = best_score,
    inlier_pairs = data.frame(ftir_idx = tent_ftir_idx,
                              raman_idx = tent_raman_idx),
    diagnostics  = list(coarse_scores = coarse_scores,
                        method = "coarse_fallback")
  )
}


#' Build a 3x3 transform matrix from coarse rotation + optional mirror
build_coarse_transform <- function(angle_deg, mirror) {
  theta <- angle_deg * pi / 180
  a <- cos(theta)
  b <- sin(theta)
  build_transform_matrix(a, b, tx = 0, ty = 0, reflect = mirror)
}
