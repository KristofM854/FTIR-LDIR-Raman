# =============================================================================
# 04_ransac.R — RANSAC-based alignment estimation
# =============================================================================

#' Estimate the geometric transform from FTIR to Raman coordinate frame
#'
#' Two-phase approach:
#'   1. Coarse grid search over rotation angles (+/- mirror). For each candidate,
#'      estimate the residual translation from nearest-neighbor pairs, then count
#'      inliers. This handles centroid mismatch between different particle
#'      populations.
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
ransac_align <- function(ftir_df, raman_df, config, src_label = NULL) {
  # src_label names the instrument pair in the log and in the low-inlier
  # warning. Without it the warning read only "Coarse alignment found very few
  # inliers", with R attributing it to `ransac_align(ftir_lm, raman_lm, ...)`
  # -- which says "ftir" for BOTH FTIR instruments and gives no way to tell
  # which alignment was weak.
  .lbl <- if (is.null(src_label) || !nzchar(src_label)) "" else
    paste0(" (", src_label, " -> Raman)")
  log_message("Starting RANSAC alignment", .lbl)

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

  # Deterministic sampling (B3): seed a local RNG stream and restore the
  # caller's global RNG state on exit, so repeated runs on identical data yield
  # the same transform (and match counts) without perturbing randomness
  # elsewhere in the pipeline. on.exit evaluates in this frame, so .old_seed is
  # in scope when it fires.
  .seed <- if (is.null(config$align_seed)) 1L else config$align_seed
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
              rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(.seed)

  raman_mat <- cbind(raman_x, raman_y)

  # =========================================================================
  # Phase 0: Candidate scale factors
  # =========================================================================
  # Instruments do not always deliver true-µm coordinates — e.g. an LDIR
  # export covering only the deposit region gets inflated ~2.6x by the 13 mm
  # scan-circle assumption. A robust span ratio between the clouds seeds
  # additional scale candidates so Phase 1 can find poses far from scale 1;
  # Phase 2 then re-estimates the exact scale from correspondences. For
  # same-scale datasets the ratio is ~1 and this collapses to the old search.
  rspan <- function(v) {
    q <- stats::quantile(v, c(0.05, 0.95), na.rm = TRUE)
    max(q[2] - q[1], 1e-9)
  }
  span_ratio <- (rspan(raman_x) + rspan(raman_y)) / (rspan(ftir_x) + rspan(ftir_y))
  scale_cands <- c(1, span_ratio * c(0.85, 1, 1.15))
  scale_cands <- sort(unique(round(scale_cands[scale_cands > 0.05], 3)))
  # Merge candidates within 10% of each other
  keep <- c(TRUE, diff(scale_cands) / head(scale_cands, -1) > 0.1)
  scale_cands <- scale_cands[keep]
  log_message("  Phase 0: scale candidates = ",
              paste(scale_cands, collapse = ", "),
              " (robust span ratio = ", round(span_ratio, 3), ")")

  # =========================================================================
  # Phase 1: Coarse grid search with translation estimation
  # =========================================================================
  log_message("  Phase 1: Coarse rotation grid search (step = ", step_deg, "deg)")

  angles <- seq(0, 360 - step_deg, by = step_deg)
  mirror_opts <- if (allow_mirror) c(FALSE, TRUE) else FALSE

  best_score  <- 0
  best_angle  <- 0
  best_mirror <- FALSE
  best_scale  <- 1
  best_tx     <- 0
  best_ty     <- 0
  coarse_scores <- data.frame(angle = numeric(), mirror = logical(),
                              scale = numeric(), n_inliers = integer())

  # Pick a subset of FTIR indices to use as translation anchors
  # (use all if small, sample if large)
  anchor_indices <- if (n_ftir <= 30) seq_len(n_ftir) else sample(n_ftir, 30)

  for (s_cand in scale_cands) {
  for (mirror in mirror_opts) {
    for (angle in angles) {
      theta <- angle * pi / 180
      ct <- cos(theta)
      st <- sin(theta)

      # Apply scale + rotation (+ optional mirror) to FTIR points
      if (!mirror) {
        rx <- s_cand * (ct * ftir_x - st * ftir_y)
        ry <- s_cand * (st * ftir_x + ct * ftir_y)
      } else {
        rx <- s_cand * (ct * ftir_x + st * ftir_y)
        ry <- s_cand * (st * ftir_x - ct * ftir_y)
      }

      # Anchor-based translation estimation:
      # For each FTIR anchor particle, use its nearest Raman neighbor to
      # define a translation hypothesis, then count inliers. This handles
      # arbitrary centroid offsets between datasets.
      nn <- RANN::nn2(raman_mat, cbind(rx, ry), k = 1)

      angle_best_n <- 0
      angle_best_tx <- 0
      angle_best_ty <- 0

      for (ai in anchor_indices) {
        # Translation that maps rotated FTIR[ai] onto its nearest Raman neighbor
        ri <- nn$nn.idx[ai, 1]
        est_tx <- raman_x[ri] - rx[ai]
        est_ty <- raman_y[ri] - ry[ai]

        # Quick inlier count with this translation
        shifted_x <- rx + est_tx
        shifted_y <- ry + est_ty
        nn2 <- RANN::nn2(raman_mat, cbind(shifted_x, shifted_y), k = 1)
        n_inliers <- sum(nn2$nn.dists[, 1] <= inlier_dist)

        if (n_inliers > angle_best_n) {
          angle_best_n  <- n_inliers
          angle_best_tx <- est_tx
          angle_best_ty <- est_ty
        }
      }

      coarse_scores <- rbind(coarse_scores,
                             data.frame(angle = angle, mirror = mirror,
                                        scale = s_cand,
                                        n_inliers = angle_best_n))

      if (angle_best_n > best_score) {
        best_score  <- angle_best_n
        best_angle  <- angle
        best_mirror <- mirror
        best_scale  <- s_cand
        best_tx     <- angle_best_tx
        best_ty     <- angle_best_ty
      }
    }
  }
  }

  log_message("  Coarse search best: angle = ", best_angle,
              " deg, mirror = ", best_mirror,
              ", scale = ", best_scale,
              ", inliers = ", best_score, " / ", n_ftir,
              ", translation = (", round(best_tx, 1), ", ", round(best_ty, 1), ")")

  # Log top 5 candidates for diagnostics
  top5 <- coarse_scores[order(-coarse_scores$n_inliers), ][1:min(5, nrow(coarse_scores)), ]
  for (i in seq_len(nrow(top5))) {
    log_message("    #", i, ": angle=", top5$angle[i],
                " deg, mirror=", top5$mirror[i],
                ", inliers=", top5$n_inliers[i])
  }

  if (best_score < min_samples) {
    warning("Coarse alignment", .lbl, " found very few inliers (", best_score,
            "). Results may be unreliable. Check that the datasets are from ",
            "the same physical sample.", call. = FALSE)
  }

  # =========================================================================
  # Phase 2: RANSAC refinement
  # =========================================================================
  log_message("  Phase 2: RANSAC refinement (", n_ransac, " iterations)")

  # Apply coarse alignment (scale + rotation + mirror + translation) to get
  # tentative correspondences
  theta_best <- best_angle * pi / 180
  ct <- cos(theta_best)
  st <- sin(theta_best)
  if (!best_mirror) {
    coarse_x <- best_scale * (ct * ftir_x - st * ftir_y) + best_tx
    coarse_y <- best_scale * (st * ftir_x + ct * ftir_y) + best_ty
  } else {
    coarse_x <- best_scale * (ct * ftir_x + st * ftir_y) + best_tx
    coarse_y <- best_scale * (st * ftir_x - ct * ftir_y) + best_ty
  }

  # Find tentative nearest-neighbor correspondences with generous threshold
  nn_coarse <- RANN::nn2(raman_mat, cbind(coarse_x, coarse_y), k = 1)
  tent_mask <- nn_coarse$nn.dists[, 1] <= inlier_dist * 2
  tent_ftir_idx  <- which(tent_mask)
  tent_raman_idx <- nn_coarse$nn.idx[tent_mask, 1]
  n_tentative <- length(tent_ftir_idx)

  log_message("  Tentative correspondences: ", n_tentative)

  if (n_tentative < min_samples) {
    log_message("  Too few tentative correspondences for RANSAC. Using coarse alignment.", level = "WARN")
    M_coarse <- build_coarse_transform_with_translation(best_angle, best_mirror, best_tx, best_ty, scale = best_scale)
    params <- extract_transform_params(M_coarse)
    return(list(
      transform    = M_coarse,
      params       = params,
      n_inliers    = best_score,
      inlier_pairs = data.frame(ftir_idx = tent_ftir_idx,
                                raman_idx = tent_raman_idx),
      diagnostics  = list(coarse_scores = coarse_scores, method = "coarse_only")
    ))
  }

  # RANSAC loop
  best_ransac_inliers <- 0
  best_ransac_M       <- NULL
  best_ransac_pairs   <- NULL

  for (iter in seq_len(n_ransac)) {
    sample_idx <- sample(n_tentative, min(min_samples, n_tentative))
    s_ftir_idx  <- tent_ftir_idx[sample_idx]
    s_raman_idx <- tent_raman_idx[sample_idx]

    tryCatch({
      tf <- estimate_similarity_transform(
        src_x = ftir_x[s_ftir_idx],
        src_y = ftir_y[s_ftir_idx],
        dst_x = raman_x[s_raman_idx],
        dst_y = raman_y[s_raman_idx],
        allow_reflection = allow_mirror
      )

      # Reject transforms whose scale strays far from the coarse estimate.
      # NOT a fixed [0.9, 1.1] window: coordinate scales can legitimately be
      # far from 1 when an instrument's export doesn't cover the assumed
      # physical extent (e.g. LDIR deposit-region export -> true scale ~0.38).
      if (tf$scale < best_scale * 0.8 || tf$scale > best_scale * 1.25) next

      transformed <- apply_transform_points(ftir_x, ftir_y, tf$matrix)
      nn_check <- RANN::nn2(raman_mat,
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
      # Skip degenerate samples
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
    nn_final <- RANN::nn2(raman_mat,
                          cbind(transformed$x_transformed, transformed$y_transformed),
                          k = 1)
    inlier_mask <- nn_final$nn.dists[, 1] <= inlier_dist
    final_pairs <- data.frame(
      ftir_idx  = which(inlier_mask),
      raman_idx = nn_final$nn.idx[inlier_mask, 1]
    )

    log_message("  RANSAC result: ", sum(inlier_mask), " inliers, ",
                "scale = ", round(final_tf$scale, 4),
                ", rotation = ", round(final_tf$rotation_deg, 2), " deg",
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
  M_coarse <- build_coarse_transform_with_translation(best_angle, best_mirror, best_tx, best_ty, scale = best_scale)
  params <- extract_transform_params(M_coarse)

  list(
    transform    = M_coarse,
    params       = params,
    n_inliers    = best_score,
    inlier_pairs = data.frame(ftir_idx = tent_ftir_idx,
                              raman_idx = tent_raman_idx),
    diagnostics  = list(coarse_scores = coarse_scores, method = "coarse_fallback")
  )
}


# =============================================================================
# Descriptor-based RANSAC alignment (optional Tier 2 replacement)
# =============================================================================

#' Compute per-particle descriptors for descriptor-based RANSAC
#'
#' Adds two columns to the data frame:
#'   log_size        — log1p(feret_max_um), robust size descriptor
#'   local_density_um — mean distance to k nearest neighbours (k=5 by default)
#'
#' @param df    Data frame with feret_max_um and coordinate columns x_col/y_col
#' @param k     Number of nearest neighbours for density (default 5)
#' @param x_col Column name for x-coordinates (default "x_norm")
#' @param y_col Column name for y-coordinates (default "y_norm")
#' @return df with log_size and local_density_um added
compute_particle_descriptors <- function(df, k = 5,
                                          x_col = "x_norm",
                                          y_col = "y_norm") {
  df$log_size <- log1p(pmax(0, df$feret_max_um))

  coords  <- cbind(df[[x_col]], df[[y_col]])
  n_valid <- nrow(coords)
  if (n_valid < 2) {
    df$local_density_um <- NA_real_
    return(df)
  }
  k_actual <- min(k + 1L, n_valid)  # +1 because nn2 includes self
  nn       <- RANN::nn2(data = coords, query = coords, k = k_actual)
  # exclude self (first column distance = 0)
  dists_ex_self <- nn$nn.dists[, -1, drop = FALSE]
  df$local_density_um <- if (ncol(dists_ex_self) > 0) {
    rowMeans(dists_ex_self)
  } else {
    NA_real_
  }
  df
}


#' Descriptor-based RANSAC similarity alignment (optional Tier 2)
#'
#' Generates candidate correspondences using log-size similarity, then runs
#' RANSAC to find the best similarity transform, refines with all inliers,
#' and applies scale / rotation guardrails to avoid degenerate solutions.
#'
#' Returns the same structure as ransac_align() for drop-in compatibility.
#' Returns NULL if RANSAC cannot find a transform meeting the constraints.
#'
#' @param ldir_df   LDIR data frame with x_norm, y_norm, feret_max_um
#' @param raman_df  Raman data frame with x_norm, y_norm, feret_max_um
#' @param config    Pipeline config (for logging)
#' @param size_tol  Maximum |log_size_ldir - log_size_raman| for candidates
#' @param n_ransac  Number of RANSAC iterations
#' @param n_sample  Points per RANSAC sample (minimum 3)
#' @param inlier_threshold_um  Distance threshold for inlier counting (µm)
#' @param scale_min  Minimum acceptable scale factor
#' @param scale_max  Maximum acceptable scale factor
#' @param rot_limit_deg  Maximum |rotation| in degrees
descriptor_ransac_align <- function(ldir_df, raman_df, config,
                                     size_tol            = 0.5,
                                     n_ransac            = 500L,
                                     n_sample            = 3L,
                                     inlier_threshold_um = 500,
                                     scale_min           = 0.8,
                                     scale_max           = 1.25,
                                     rot_limit_deg       = 45) {
  log_message("  Descriptor RANSAC: computing particle descriptors…")

  ldir_d  <- compute_particle_descriptors(ldir_df,  x_col = "x_norm", y_col = "y_norm")
  raman_d <- compute_particle_descriptors(raman_df, x_col = "x_norm", y_col = "y_norm")

  # --- Candidate pair lists: for each LDIR particle, Raman candidates by size ---
  candidates <- lapply(seq_len(nrow(ldir_d)), function(i) {
    which(abs(raman_d$log_size - ldir_d$log_size[i]) < size_tol)
  })
  eligible <- which(lengths(candidates) >= 1L)

  if (length(eligible) < n_sample) {
    log_message("  Descriptor RANSAC: too few eligible pairs (",
                length(eligible), ") — aborting", level = "WARN")
    return(NULL)
  }

  log_message("  Descriptor RANSAC: ", nrow(ldir_d), " LDIR pts, ",
              nrow(raman_d), " Raman pts, ",
              length(eligible), " eligible LDIR pts")

  # Deterministic sampling (B3): local seed + global RNG restore on exit.
  .seed <- if (is.null(config$align_seed)) 1L else config$align_seed
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
              rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(.seed)

  best_inliers <- 0L
  best_tf      <- NULL

  for (iter in seq_len(n_ransac)) {
    # Sample n_sample distinct LDIR particles with candidates
    sel <- sample(eligible, min(n_sample, length(eligible)), replace = FALSE)

    # For each selected LDIR point draw one random Raman candidate
    raman_sel <- vapply(sel, function(i) {
      cands <- candidates[[i]]
      cands[sample.int(length(cands), 1L)]
    }, integer(1L))

    tf <- tryCatch(
      estimate_similarity_transform(
        ldir_d$x_norm[sel],  ldir_d$y_norm[sel],
        raman_d$x_norm[raman_sel], raman_d$y_norm[raman_sel],
        allow_reflection = FALSE   # avoid reflection ambiguity in RANSAC sampling
      ),
      # Expected control flow, not an error to surface: a degenerate random
      # triplet (collinear / coincident points) makes the fit unsolvable — we
      # simply skip this iteration. Runs thousands of times, so no logging.
      error = function(e) NULL
    )
    if (is.null(tf)) next
    if (tf$scale < scale_min || tf$scale > scale_max) next
    if (abs(tf$rotation_deg) > rot_limit_deg) next

    # Count inliers across all LDIR particles
    tf_pts <- apply_transform_points(ldir_d$x_norm, ldir_d$y_norm, tf$matrix)
    nn <- RANN::nn2(
      cbind(raman_d$x_norm, raman_d$y_norm),
      cbind(tf_pts$x_transformed, tf_pts$y_transformed),
      k = 1
    )
    n_in <- sum(nn$nn.dists[, 1] < inlier_threshold_um)
    if (n_in > best_inliers) {
      best_inliers <- n_in
      best_tf      <- tf
    }
  }

  if (is.null(best_tf)) {
    log_message("  Descriptor RANSAC: no valid transform found after ",
                n_ransac, " iterations", level = "WARN")
    return(NULL)
  }

  log_message("  Descriptor RANSAC best: scale=", round(best_tf$scale, 4),
              ", rot=", round(best_tf$rotation_deg, 2),
              "°, inliers=", best_inliers, "/", nrow(ldir_d))

  # --- Refine with all inliers from the best model ---
  tf_pts_best <- apply_transform_points(ldir_d$x_norm, ldir_d$y_norm, best_tf$matrix)
  nn_final    <- RANN::nn2(
    cbind(raman_d$x_norm, raman_d$y_norm),
    cbind(tf_pts_best$x_transformed, tf_pts_best$y_transformed),
    k = 1
  )
  inlier_mask <- nn_final$nn.dists[, 1] < inlier_threshold_um
  raman_inlier_idx <- nn_final$nn.idx[inlier_mask, 1]

  if (sum(inlier_mask) >= 3L) {
    refined_tf <- tryCatch(
      estimate_similarity_transform(
        ldir_d$x_norm[inlier_mask],  ldir_d$y_norm[inlier_mask],
        raman_d$x_norm[raman_inlier_idx], raman_d$y_norm[raman_inlier_idx]
      ),
      error = function(e) best_tf
    )
  } else {
    refined_tf <- best_tf
  }

  log_message("  Descriptor RANSAC refined: scale=", round(refined_tf$scale, 4),
              ", rot=", round(refined_tf$rotation_deg, 2),
              "°, rms=", round(refined_tf$residual_rms, 1), " \u00b5m")

  params <- list(
    scale        = refined_tf$scale,
    rotation_deg = refined_tf$rotation_deg,
    reflected    = isTRUE(refined_tf$reflected),
    tx           = refined_tf$tx,
    ty           = refined_tf$ty
  )
  list(
    transform    = refined_tf$matrix,
    params       = params,
    n_inliers    = sum(inlier_mask),
    inlier_frac  = mean(inlier_mask),
    rms          = refined_tf$residual_rms,
    inlier_pairs = data.frame(
      ldir_idx  = which(inlier_mask),
      raman_idx = raman_inlier_idx,
      stringsAsFactors = FALSE
    ),
    diagnostics  = list(method = "descriptor_ransac",
                        n_ransac_iterations = n_ransac,
                        size_tol = size_tol)
  )
}


#' Build a 3x3 transform matrix from coarse scale + rotation + optional mirror + translation
build_coarse_transform_with_translation <- function(angle_deg, mirror, tx, ty,
                                                    scale = 1) {
  theta <- angle_deg * pi / 180
  a <- scale * cos(theta)
  b <- scale * sin(theta)
  build_transform_matrix(a, b, tx = tx, ty = ty, reflect = mirror)
}


# ============================================================================
# Global registration — robust alignment for sparse / heavily-transformed
# point clouds (e.g. LDIR -> Raman with large scale + rotation differences).
# ============================================================================
#
# The coarse RANSAC in ransac_align() anchors translation on single nearest-
# neighbour guesses, which is fragile when few particles overlap and the
# scale/rotation are far from 1/0: it can lock onto a local optimum that fits
# only a handful of points (observed on real LDIR data: 5 inliers where 23
# are achievable).  This routine instead sweeps rotation x scale globally and,
# for each candidate, recovers the translation by VOTING over all pairwise
# offset vectors, scoring by ONE-TO-ONE inlier count (each particle used at
# most once) so a degenerate collapse cannot win.  The winning pose's
# correspondences are then refit with estimate_similarity_transform().
#
# Same return contract as ransac_align().
global_register_align <- function(src_df, ref_df, config,
                                  allow_mirror = TRUE) {
  sx <- src_df$x_norm; sy <- src_df$y_norm
  rx <- ref_df$x_norm; ry <- ref_df$y_norm
  n_src <- length(sx); n_ref <- length(rx)
  if (n_src < 3 || n_ref < 3)
    stop("global_register_align needs >=3 points per cloud (src ", n_src,
         ", ref ", n_ref, ")")

  tol      <- if (is.null(config$ransac_inlier_dist_um)) 200 else config$ransac_inlier_dist_um
  step_deg <- if (is.null(config$ransac_coarse_step_deg)) 2 else config$ransac_coarse_step_deg
  if (step_deg < 2) step_deg <- 2   # 1-deg grid is needless here and slow

  scx <- mean(sx); scy <- mean(sy); rcx <- mean(rx); rcy <- mean(ry)
  sxc <- sx - scx; syc <- sy - scy
  rxc <- rx - rcx; ryc <- ry - rcy

  span_ratio <- (align_rspan(rxc) + align_rspan(ryc)) /
                (align_rspan(sxc) + align_rspan(syc))

  # Scale candidates. Both clouds are in physical micrometres, so the true
  # scale is near 1 and near span_ratio (the two clouds cover the same field).
  # The old sweep ran seq(0.2, 1.3, 0.05), which offered the search a range of
  # collapsed poses that shrink the source into a dense part of the target:
  # those score well on nearest-neighbour distance and then hand ICP a wrong
  # starting pose. Keep the same band ransac_align() enforces, widened only as
  # far as span_ratio actually implies. Override with config$align_scale_min /
  # _max for data where the two frames genuinely differ in scale.
  s_min <- if (!is.null(config$align_scale_min)) config$align_scale_min else 0.8
  s_max <- if (!is.null(config$align_scale_max)) config$align_scale_max else 1.25
  if (is.finite(span_ratio) && span_ratio > 0) {
    s_min <- min(s_min, span_ratio * 0.75)
    s_max <- max(s_max, span_ratio * 1.25)
  }
  scales <- sort(unique(round(c(1, seq(s_min, s_max, by = 0.05),
                                span_ratio * seq(0.75, 1.25, 0.05)), 3)))
  scales <- scales[scales >= s_min & scales <= s_max & scales > 0.02]
  mirrors <- if (allow_mirror) c(FALSE, TRUE) else FALSE

  # Subsample for the O(n_src*n_ref) search; re-score winner on full clouds.
  # Deterministic sampling (B3): seed a local RNG stream and restore the
  # caller's global RNG state on exit (on.exit fires in this frame), so this is
  # reproducible without leaving the global RNG reset for downstream code.
  .seed <- if (is.null(config$align_seed)) 1L else config$align_seed
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
              rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(.seed)
  sub <- function(n, m = 150L) if (n > m) sort(sample.int(n, m)) else seq_len(n)
  si <- sub(n_src); ri <- sub(n_ref)

  # Pose search: align_score_pose() / align_pose_xy() live in R/align_helpers.R
  # (shared with tools/diagnose_*.R). `tol` is threaded through explicitly.
  best <- list(n = -1)
  angles <- seq(0, 360 - step_deg, by = step_deg)
  for (mir in mirrors) for (s in scales) for (deg in angles) {
    r <- align_score_pose(deg, s, mir, sxc[si], syc[si], rxc[ri], ryc[ri], tol)
    if (r$n > best$n) best <- c(r, list(deg = deg, s = s, mir = mir))
  }

  # Correspondences at the winning pose, on the FULL clouds
  p <- align_pose_xy(best$deg, best$s, best$mir, sxc, syc)
  px <- p$x + best$tx; py <- p$y + best$ty
  d  <- sqrt(outer(rxc, px, "-")^2 + outer(ryc, py, "-")^2)
  ok <- which(d <= tol, arr.ind = TRUE)
  pairs <- data.frame()
  if (nrow(ok) > 0) {
    ok <- ok[order(d[ok]), , drop = FALSE]
    ur <- logical(n_ref); uc <- logical(n_src); rows <- list()
    for (r in seq_len(nrow(ok))) {
      i <- ok[r, 1]; j <- ok[r, 2]
      if (!ur[i] && !uc[j]) {
        ur[i] <- TRUE; uc[j] <- TRUE
        rows[[length(rows) + 1]] <- data.frame(src_idx = j, ref_idx = i)
      }
    }
    pairs <- do.call(rbind, rows)
  }

  # Refit a clean similarity from the correspondences (in original x_norm space)
  if (nrow(pairs) >= 2) {
    tf <- estimate_similarity_transform(
      src_x = sx[pairs$src_idx], src_y = sy[pairs$src_idx],
      dst_x = rx[pairs$ref_idx], dst_y = ry[pairs$ref_idx],
      allow_reflection = allow_mirror)
    trans <- apply_transform_points(sx, sy, tf$matrix)
    nn <- RANN::nn2(cbind(rx, ry),
                    cbind(trans$x_transformed, trans$y_transformed), k = 1)
    inl <- nn$nn.dists[, 1] <= tol
    final_pairs <- data.frame(src_idx = which(inl),
                              ref_idx = nn$nn.idx[inl, 1])
    log_message("  Global registration: rot=", round(tf$rotation_deg, 2),
                " deg, scale=", round(tf$scale, 4),
                ", reflected=", tf$reflected,
                ", inliers=", sum(inl), " / ", n_src,
                " (coarse best ", best$n, ")")
    return(list(transform = tf$matrix,
                params = extract_transform_params(tf$matrix),
                n_inliers = sum(inl),
                inlier_pairs = final_pairs,
                diagnostics = list(method = "global_register",
                                   coarse_inliers = best$n,
                                   coarse_scale = best$s,
                                   coarse_deg = best$deg,
                                   coarse_mirror = best$mir)))
  }

  # No usable correspondences — hand back the coarse pose transform
  M <- build_coarse_transform_with_translation(best$deg, best$mir,
         best$tx + rcx - best$s * 0, best$ty + rcy - best$s * 0,
         scale = best$s)
  log_message("  Global registration: only ", best$n,
              " coarse inliers, no refit — using coarse pose", level = "WARN")
  list(transform = M, params = extract_transform_params(M),
       n_inliers = best$n, inlier_pairs = data.frame(),
       diagnostics = list(method = "global_register_coarse"))
}
