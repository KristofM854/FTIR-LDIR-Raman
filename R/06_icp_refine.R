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

  ftir_x  <- ftir_df$x_norm
  ftir_y  <- ftir_df$y_norm
  raman_x <- raman_df$x_norm
  raman_y <- raman_df$y_norm

  current_M   <- initial_transform
  rms_history <- numeric()
  prev_rms    <- Inf
  converged   <- FALSE

  for (iter in seq_len(max_iter)) {
    # Apply current transform
    transformed <- apply_transform_points(ftir_x, ftir_y, current_M)
    tx <- transformed$x_transformed
    ty <- transformed$y_transformed

    # Find nearest Raman neighbor for each transformed FTIR particle
    nn <- RANN::nn2(cbind(raman_x, raman_y), cbind(tx, ty), k = 1)
    dists     <- nn$nn.dists[, 1]
    raman_idx <- nn$nn.idx[, 1]

    # Filter by max distance
    keep <- dists <= max_pair_dist
    if (sum(keep) < 3) {
      log_message("  ICP iteration ", iter, ": too few pairs (", sum(keep),
                  "). Stopping.", level = "WARN")
      break
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

    # Re-estimate transform from filtered pairs
    new_tf <- estimate_similarity_transform(
      src_x = ftir_x[keep],
      src_y = ftir_y[keep],
      dst_x = raman_x[raman_idx[keep]],
      dst_y = raman_y[raman_idx[keep]],
      allow_reflection = allow_mirror
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
