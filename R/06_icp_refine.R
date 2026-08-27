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
icp_refine <- function(ftir_df, raman_df, initial_transform, config,
                       anchor_pairs = NULL) {
  log_message("Starting ICP refinement")

  max_iter       <- config$icp_max_iterations
  conv_thresh    <- config$icp_convergence_thresh
  max_pair_dist  <- config$icp_max_pair_dist_um
  allow_mirror   <- config$ransac_allow_mirror
  use_reciprocal <- isTRUE(config$icp_reciprocal)
  trim_pct       <- if (!is.null(config$icp_trim_pct)) config$icp_trim_pct else 0

  use_elong_weight <- isTRUE(config$icp_elongation_downweight)
  elong_alpha      <- if (!is.null(config$icp_elongation_alpha)) config$icp_elongation_alpha else 0.5

  # --- Scale guard ----------------------------------------------------------
  # ICP re-estimates an unconstrained similarity from nearest-neighbour pairs
  # every iteration, and that objective REWARDS collapse: shrinking the source
  # cloud packs points into dense regions of the target, so the nearest
  # neighbour distance falls as the fit gets more wrong. Measured on a real
  # LDIR<->Raman run: the correct pose (scale 1.06) scores RMS 301 um while a
  # collapsed one (scale 0.23) scores 67 um and a fully collapsed one (scale
  # 0.06) scores 38 um. Left unconstrained ICP walks straight down that hill.
  #
  # Both clouds are in physical micrometres, so the true scale is near 1.
  # ransac_align() already enforces [0.8, 1.25]; ICP had no such guard, which
  # is how a run ended up with scale 0.2336. Clamp each re-estimate to the same
  # band and keep the rotation/translation from that step.
  scale_min <- if (!is.null(config$icp_scale_min)) config$icp_scale_min else 0.8
  scale_max <- if (!is.null(config$icp_scale_max)) config$icp_scale_max else 1.25
  n_clamped <- 0L

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
    # Anchor pairs (if any) are pinned with very high weight (100×) so ICP
    # cannot move the landmark correspondences away from the Procrustes solution.
    pair_weights <- if (use_elong_weight) ftir_weights[keep] else rep(1.0, sum(keep))

    est_src_x <- ftir_x[keep]
    est_src_y <- ftir_y[keep]
    est_dst_x <- raman_x[raman_idx[keep]]
    est_dst_y <- raman_y[raman_idx[keep]]
    est_w     <- pair_weights

    if (!is.null(anchor_pairs) && nrow(anchor_pairs) > 0) {
      anc_src_x <- ftir_x[anchor_pairs$src_idx]
      anc_src_y <- ftir_y[anchor_pairs$src_idx]
      anc_dst_x <- raman_x[anchor_pairs$tgt_idx]
      anc_dst_y <- raman_y[anchor_pairs$tgt_idx]
      anc_w     <- rep(100.0, nrow(anchor_pairs))

      est_src_x <- c(est_src_x, anc_src_x)
      est_src_y <- c(est_src_y, anc_src_y)
      est_dst_x <- c(est_dst_x, anc_dst_x)
      est_dst_y <- c(est_dst_y, anc_dst_y)
      est_w     <- c(est_w, anc_w)
    }

    new_tf <- estimate_similarity_transform(
      src_x = est_src_x,
      src_y = est_src_y,
      dst_x = est_dst_x,
      dst_y = est_dst_y,
      allow_reflection = allow_mirror,
      weights = est_w
    )

    # Clamp the scale back into the plausible band (see the scale guard note
    # above). Rescaling about the source centroid keeps the rotation and the
    # correspondence geometry this step found, and only removes the collapse
    # component -- rebuilding the matrix from scratch would discard the step.
    new_M <- new_tf$matrix
    s_new <- tryCatch(extract_transform_params(new_M)$scale,
                      error = function(e) NA_real_)
    if (is.finite(s_new) && (s_new < scale_min || s_new > scale_max)) {
      s_clamped <- min(max(s_new, scale_min), scale_max)
      k <- s_clamped / s_new
      cx <- mean(est_src_x); cy <- mean(est_src_y)
      # Map the source centroid, scale the linear part by k, then re-anchor so
      # the centroid still lands where this iteration put it.
      p_before <- apply_transform_points(cx, cy, new_M)
      new_M[1:2, 1:2] <- new_M[1:2, 1:2] * k
      p_after <- apply_transform_points(cx, cy, new_M)
      new_M[1, 3] <- new_M[1, 3] + (p_before$x_transformed - p_after$x_transformed)
      new_M[2, 3] <- new_M[2, 3] + (p_before$y_transformed - p_after$y_transformed)
      n_clamped <- n_clamped + 1L
      log_message("  ICP iter ", iter, ": scale ", round(s_new, 4),
                  " outside [", scale_min, ", ", scale_max, "] - clamped to ",
                  round(s_clamped, 4), level = "WARN")
    }

    current_M <- new_M
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

  # A run that needed clamping was actively trying to collapse -- almost always
  # a sign the initial pose was wrong, so say so rather than silently returning
  # a transform that merely looks plausible.
  if (n_clamped > 0)
    log_message("  ICP: scale clamped on ", n_clamped, " of ",
                length(rms_history), " iterations. The starting pose is ",
                "probably wrong - check the alignment inlier count.",
                level = "WARN")

  list(
    transform    = current_M,
    params       = params,
    residuals    = nn_final$nn.dists[, 1],
    rms_history  = rms_history,
    converged    = converged,
    n_iterations = length(rms_history),
    n_scale_clamped = n_clamped
  )
}
