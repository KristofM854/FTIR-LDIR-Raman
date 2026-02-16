# =============================================================================
# 03b_landmark_align.R — Landmark-first alignment using large particles & fibers
# =============================================================================
#
# Strategy: large particles and fibers are easily identifiable in both FTIR
# and Raman images. By aligning on these "landmarks" first, we can often
# determine the spatial transform quickly and with high confidence, skipping
# the expensive full RANSAC grid search entirely.
#
# A particle qualifies as a landmark if:
#   - feret_max_um >= landmark_min_size_um  (large particle), OR
#   - aspect ratio (major/minor) >= landmark_fiber_aspect_ratio AND
#     feret_max_um >= landmark_fiber_min_size_um  (fiber)
# =============================================================================

#' Select landmark particles from a data frame
#'
#' @param df Data frame with feret_max_um, major_um, minor_um columns
#' @param config Configuration list
#' @return Logical vector (TRUE = landmark)
select_landmarks <- function(df, config) {
  n <- nrow(df)
  is_landmark <- logical(n)

  size   <- df$feret_max_um
  major  <- df$major_um
  minor  <- df$minor_um

  # Criterion 1: large particles
  if (any(!is.na(size))) {
    is_landmark <- is_landmark | (!is.na(size) & size >= config$landmark_min_size_um)
  }

  # Criterion 2: fibers (high aspect ratio, above minimum size)
  aspect <- ifelse(!is.na(major) & !is.na(minor) & minor > 0,
                   major / minor, NA_real_)
  if (any(!is.na(aspect))) {
    is_fiber <- !is.na(aspect) &
      aspect >= config$landmark_fiber_aspect_ratio &
      !is.na(size) &
      size >= config$landmark_fiber_min_size_um
    is_landmark <- is_landmark | is_fiber
  }

  is_landmark
}


#' Run landmark-based alignment (Tier 1)
#'
#' Extracts landmarks from both datasets, runs RANSAC on them, and evaluates
#' confidence. Returns the transform and a confidence assessment.
#'
#' @param ftir_df FTIR data frame with x_norm, y_norm (centered)
#' @param raman_df Raman data frame with x_norm, y_norm (centered)
#' @param config Configuration list
#' @return List with:
#'   success       — logical, TRUE if enough landmarks were found and aligned
#'   confident     — logical, TRUE if alignment passes confidence thresholds
#'   transform     — 3x3 transform matrix (NULL if !success)
#'   params        — human-readable transform params
#'   n_ftir_landmarks  — number of FTIR landmarks
#'   n_raman_landmarks — number of Raman landmarks
#'   n_inliers     — number of inlier landmark pairs
#'   inlier_ratio  — fraction of landmarks that are inliers
#'   mean_residual — mean distance of inlier pairs (µm)
#'   inlier_pairs  — data frame of matched landmark pairs
landmark_align <- function(ftir_df, raman_df, config) {
  log_message(strrep("-", 50))
  log_message("Tier 1: Landmark-based alignment")

  # --- Select landmarks ---
  ftir_lm_mask  <- select_landmarks(ftir_df, config)
  raman_lm_mask <- select_landmarks(raman_df, config)

  n_ftir_lm  <- sum(ftir_lm_mask)
  n_raman_lm <- sum(raman_lm_mask)

  log_message("  FTIR landmarks:  ", n_ftir_lm, " / ", nrow(ftir_df))
  log_message("  Raman landmarks: ", n_raman_lm, " / ", nrow(raman_df))

  min_count <- config$landmark_min_count

  if (n_ftir_lm < min_count || n_raman_lm < min_count) {
    log_message("  Too few landmarks (need >= ", min_count, " per dataset). Skipping Tier 1.")
    return(list(success = FALSE, confident = FALSE, transform = NULL,
                params = NULL, n_ftir_landmarks = n_ftir_lm,
                n_raman_landmarks = n_raman_lm, n_inliers = 0,
                inlier_ratio = 0, mean_residual = NA_real_,
                inlier_pairs = data.frame()))
  }

  ftir_lm  <- ftir_df[ftir_lm_mask, ]
  raman_lm <- raman_df[raman_lm_mask, ]

  # Log size stats for landmarks
  if (any(!is.na(ftir_lm$feret_max_um))) {
    log_message("  FTIR landmark sizes: ",
                round(min(ftir_lm$feret_max_um, na.rm = TRUE)), " – ",
                round(max(ftir_lm$feret_max_um, na.rm = TRUE)), " µm")
  }
  if (any(!is.na(raman_lm$feret_max_um))) {
    log_message("  Raman landmark sizes: ",
                round(min(raman_lm$feret_max_um, na.rm = TRUE)), " – ",
                round(max(raman_lm$feret_max_um, na.rm = TRUE)), " µm")
  }

  # --- Run RANSAC on landmarks only ---
  # Use tighter inlier distance for landmarks (they should match closely)
  landmark_config <- config
  landmark_config$ransac_inlier_dist_um <- min(config$ransac_inlier_dist_um,
                                                config$landmark_confidence_max_residual_um * 2)

  ransac_result <- ransac_align(ftir_lm, raman_lm, landmark_config)

  # --- Evaluate confidence ---
  # Apply the transform and measure how well landmarks match
  transformed <- apply_transform_points(ftir_lm$x_norm, ftir_lm$y_norm,
                                        ransac_result$transform)
  raman_mat <- cbind(raman_lm$x_norm, raman_lm$y_norm)
  nn <- RANN::nn2(raman_mat,
                  cbind(transformed$x_transformed, transformed$y_transformed),
                  k = 1)

  inlier_dist <- config$landmark_confidence_max_residual_um
  inlier_mask <- nn$nn.dists[, 1] <= inlier_dist
  n_inliers   <- sum(inlier_mask)
  n_landmarks <- n_ftir_lm  # denominator = FTIR landmarks (we're transforming FTIR → Raman)
  inlier_ratio <- n_inliers / n_landmarks
  mean_residual <- if (n_inliers > 0) mean(nn$nn.dists[inlier_mask, 1]) else NA_real_

  log_message("  Landmark inliers: ", n_inliers, " / ", n_landmarks,
              " (", round(inlier_ratio * 100, 1), "%)")
  if (!is.na(mean_residual)) {
    log_message("  Mean landmark residual: ", round(mean_residual, 2), " µm")
  }

  # Build inlier pairs
  inlier_pairs <- data.frame(
    ftir_idx  = which(inlier_mask),
    raman_idx = nn$nn.idx[inlier_mask, 1],
    distance  = nn$nn.dists[inlier_mask, 1]
  )

  # Confidence check
  confident <- (inlier_ratio >= config$landmark_confidence_min_inlier_ratio) &&
    (!is.na(mean_residual) && mean_residual <= config$landmark_confidence_max_residual_um) &&
    (n_inliers >= config$landmark_min_count)

  if (confident) {
    log_message("  >>> Landmark alignment CONFIDENT — Tier 1 succeeded!")
  } else {
    reasons <- character()
    if (inlier_ratio < config$landmark_confidence_min_inlier_ratio)
      reasons <- c(reasons, paste0("inlier ratio too low (",
                                   round(inlier_ratio * 100, 1), "%)"))
    if (!is.na(mean_residual) && mean_residual > config$landmark_confidence_max_residual_um)
      reasons <- c(reasons, paste0("mean residual too high (",
                                   round(mean_residual, 1), " µm)"))
    if (n_inliers < config$landmark_min_count)
      reasons <- c(reasons, paste0("too few inliers (", n_inliers, ")"))
    log_message("  Landmark alignment not confident: ", paste(reasons, collapse = "; "))
  }

  params <- extract_transform_params(ransac_result$transform)

  list(
    success           = TRUE,
    confident         = confident,
    transform         = ransac_result$transform,
    params            = params,
    n_ftir_landmarks  = n_ftir_lm,
    n_raman_landmarks = n_raman_lm,
    n_inliers         = n_inliers,
    inlier_ratio      = inlier_ratio,
    mean_residual     = mean_residual,
    inlier_pairs      = inlier_pairs
  )
}
