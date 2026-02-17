# =============================================================================
# 06_icp_refine.R — Iterative Closest Point (ICP) refinement
# =============================================================================

#' Refine alignment using ICP
#'
#' Starting from the RANSAC transform, iteratively:
#'   1. Find nearest Raman neighbor for each transformed FTIR particle
#'   2. Filter pairs by maximum distance
#'   3. Re-estimate similarity transform from the filtered pairs
#'   4. Apply the updated transform
#'   5. Stop when RMS improvement is below threshold or max iterations reached
#'
#' @param ftir_df FTIR data frame with x_norm, y_norm
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param initial_transform 3x3 matrix from RANSAC
#' @param config Configuration list
#' @return List with:
#'   transform     — refined 3x3 transform matrix
#'   params        — human-readable parameters
#'   residuals     — vector of per-pair distances at final iteration
#'   rms_history   — RMS error at each iteration
#'   converged     — logical
#'   n_iterations  — number of iterations performed
icp_refine <- function(ftir_df, raman_df, initial_transform, config) {
  log_message("Starting ICP refinement")

  max_iter       <- config$icp_max_iterations
  conv_thresh    <- config$icp_convergence_thresh
  max_pair_dist  <- config$icp_max_pair_dist_um
  allow_mirror   <- config$ransac_allow_mirror
  use_reciprocal <- isTRUE(config$icp_reciprocal)
  trim_pct       <- if (!is.null(config$icp_trim_pct)) config$icp_trim_pct else 0

  use_elong_weight <- isTRUE(config$icp_elongation_downweight)
  elong_alpha      <- if (!is.null(config$icp_elongation_alpha)) config$icp_elongation_alpha else 0.5

  ftir_x  <- ftir_df$x_norm
  ftir_y  <- ftir_df$y_norm
  raman_x <- raman_df$x_norm
  raman_y <- raman_df$y_norm
  n_ftir  <- length(ftir_x)
  n_raman <- length(raman_x)

  ftir_mat  <- cbind(ftir_x, ftir_y)
  raman_mat <- cbind(raman_x, raman_y)

  # Compute per-particle elongation weights for FTIR
  # Round particles (aspect ~1) get weight ~1; fibers (aspect ~5) get lower weight
  ftir_weights <- rep(1.0, n_ftir)
  if (use_elong_weight && "major_um" %in% names(ftir_df) && "minor_um" %in% names(ftir_df)) {
    major <- ftir_df$major_um
    minor <- ftir_df$minor_um
    aspect <- ifelse(!is.na(major) & !is.na(minor) & minor > 0,
                     major / minor, 1.0)
    ftir_weights <- 1.0 / (1.0 + elong_alpha * pmax(aspect - 1, 0))
    log_message("  Elongation weighting: ON (alpha = ", elong_alpha,
                ", weight range: ", round(min(ftir_weights), 3), " - ",
                round(max(ftir_weights), 3), ")")
  }

  if (use_reciprocal) log_message("  Reciprocal nearest-neighbor filtering: ON")
  if (trim_pct > 0)   log_message("  Trimming worst ", trim_pct * 100, "% of pairs each iteration")

  current_M   <- initial_transform
  rms_history <- numeric()
  prev_rms    <- Inf
  converged   <- FALSE

  for (iter in seq_len(max_iter)) {
    # Apply current transform
    transformed <- apply_transform_points(ftir_x, ftir_y, current_M)
    tx <- transformed$x_transformed
    ty <- transformed$y_transformed
    transformed_mat <- cbind(tx, ty)

    # Forward: for each FTIR particle, find nearest Raman
    nn_fwd <- RANN::nn2(raman_mat, transformed_mat, k = 1)
    dists     <- nn_fwd$nn.dists[, 1]
    raman_idx <- nn_fwd$nn.idx[, 1]

    # Filter by max distance
    keep <- dists <= max_pair_dist

    # Reciprocal filter: only keep pairs where Raman→FTIR also agrees
    if (use_reciprocal && sum(keep) > 3) {
      nn_rev <- RANN::nn2(transformed_mat, raman_mat, k = 1)
      # For each FTIR[i] → Raman[j] pair, check that Raman[j] → FTIR[i]
      reciprocal <- logical(n_ftir)
      for (fi in which(keep)) {
        ri <- raman_idx[fi]
        reciprocal[fi] <- (nn_rev$nn.idx[ri, 1] == fi)
      }
      n_before_recip <- sum(keep)
      keep <- keep & reciprocal
      if (iter == 1) {
        log_message("  Reciprocal filter: ", n_before_recip, " → ", sum(keep), " pairs")
      }
    }

    if (sum(keep) < 3) {
      log_message("  ICP iteration ", iter, ": too few pairs (", sum(keep),
                  "). Stopping.", level = "WARN")
      break
    }

    # Trim the worst N% of remaining pairs (by distance)
    if (trim_pct > 0 && sum(keep) > 5) {
      keep_idx <- which(keep)
      keep_dists <- dists[keep_idx]
      n_keep <- length(keep_idx)
      n_trim <- max(0, floor(n_keep * trim_pct))
      if (n_trim > 0 && (n_keep - n_trim) >= 3) {
        cutoff <- sort(keep_dists)[n_keep - n_trim]
        keep[keep_idx[keep_dists > cutoff]] <- FALSE
      }
    }

    # Current RMS
    current_rms <- sqrt(mean(dists[keep]^2))
    rms_history <- c(rms_history, current_rms)

    log_message("  ICP iter ", iter, ": pairs = ", sum(keep),
                ", RMS = ", round(current_rms, 3), " µm")

    # Check convergence
    improvement <- prev_rms - current_rms
    if (abs(improvement) < conv_thresh) {
      log_message("  ICP converged (improvement = ", round(improvement, 4), " µm)")
      converged <- TRUE
      break
    }

    # Re-estimate transform from filtered pairs (with elongation weights)
    pair_weights <- if (use_elong_weight) ftir_weights[keep] else NULL
    new_tf <- estimate_similarity_transform(
      src_x = ftir_x[keep],
      src_y = ftir_y[keep],
      dst_x = raman_x[raman_idx[keep]],
      dst_y = raman_y[raman_idx[keep]],
      allow_reflection = allow_mirror,
      weights = pair_weights
    )

    current_M <- new_tf$matrix
    prev_rms  <- current_rms
  }

  # Final evaluation
  transformed_final <- apply_transform_points(ftir_x, ftir_y, current_M)
  nn_final <- RANN::nn2(cbind(raman_x, raman_y),
                        cbind(transformed_final$x_transformed,
                              transformed_final$y_transformed),
                        k = 1)

  params <- extract_transform_params(current_M)
  final_rms <- if (length(rms_history) > 0) tail(rms_history, 1) else NA_real_

  log_message("  ICP complete: ", length(rms_history), " iterations, ",
              "final RMS = ", round(final_rms, 3), " µm, ",
              "converged = ", converged)

  list(
    transform    = current_M,
    params       = params,
    residuals    = nn_final$nn.dists[, 1],
    rms_history  = rms_history,
    converged    = converged,
    n_iterations = length(rms_history)
  )
}
