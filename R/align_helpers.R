# =============================================================================
# align_helpers.R — shared primitives for coarse pose search & inlier counting
# =============================================================================
# Extracted from global_register_align() (04_ransac.R). The alignment core and
# the tools/diagnose_*.R scripts previously each carried their own copy of the
# rotation/scale pose transform, the one-to-one inlier count, the robust span,
# and the translation-voting score. They now share this one implementation.
#
# The pose/inlier primitives are pure (no global state, no RNG) so they are safe
# to reuse from any context. select_material_anchors() is the one exception: it
# logs, and therefore needs log_message() from utils.R.

#' Robust span of a coordinate vector: 5th–95th percentile range, floored at
#' 1e-9 to avoid divide-by-zero on degenerate (near-constant) inputs.
#' @param v Numeric vector.
#' @return Positive numeric span.
align_rspan <- function(v) {
  q <- stats::quantile(v, c(0.05, 0.95), na.rm = TRUE)
  max(q[2] - q[1], 1e-9)
}

#' Apply rotation (degrees) + uniform scale (+ optional mirror) to coordinates.
#' Mirror flips the sign convention on the cross terms (reflected pose).
#' @param deg Rotation in degrees.
#' @param s Uniform scale factor.
#' @param mir Logical; TRUE applies a reflection.
#' @param X,Y Numeric coordinate vectors.
#' @return list(x, y) of transformed coordinates.
align_pose_xy <- function(deg, s, mir, X, Y) {
  th <- deg * pi / 180; ct <- cos(th); st <- sin(th)
  if (!mir) list(x = s * (ct * X - st * Y), y = s * (st * X + ct * Y))
  else      list(x = s * (ct * X + st * Y), y = s * (st * X - ct * Y))
}

#' Count greedy one-to-one correspondences within `tol` between two clouds.
#' Points (px,py) are matched against (qx,qy); nearest pairs first, each point
#' on either side used at most once. Returns the integer number of mutual
#' matches (a lower bound that resists the "collapse everything into the densest
#' region" degeneracy that plain nearest-neighbour counting rewards).
#' @param px,py Numeric coordinate vectors of the first cloud.
#' @param qx,qy Numeric coordinate vectors of the second cloud.
#' @param tol Distance tolerance (same units as the coordinates).
#' @return Integer count of one-to-one inliers.
align_one_to_one <- function(px, py, qx, qy, tol) {
  d <- sqrt(outer(qx, px, "-")^2 + outer(qy, py, "-")^2)
  ok <- which(d <= tol, arr.ind = TRUE)
  if (nrow(ok) == 0) return(0L)
  ok <- ok[order(d[ok]), , drop = FALSE]
  ur <- logical(length(qx)); uc <- logical(length(px)); n <- 0L
  for (r in seq_len(nrow(ok))) {
    i <- ok[r, 1]; j <- ok[r, 2]
    if (!ur[i] && !uc[j]) { ur[i] <- TRUE; uc[j] <- TRUE; n <- n + 1L }
  }
  n
}

#' Score a (rotation, scale, mirror) pose by translation voting + one-to-one
#' inliers. Bins the pairwise src→ref displacement at the given tolerance, tries
#' the top few translation candidates, and returns the best.
#' @param deg,s,mir Pose parameters (see align_pose_xy).
#' @param sX,sY Source coordinate vectors.
#' @param rX,rY Reference coordinate vectors.
#' @param tol Distance tolerance for binning and inlier counting.
#' @return list(n, tx, ty): inlier count and winning translation.
align_score_pose <- function(deg, s, mir, sX, sY, rX, rY, tol) {
  p <- align_pose_xy(deg, s, mir, sX, sY)
  dx <- outer(rX, p$x, "-"); dy <- outer(rY, p$y, "-")
  key <- paste(round(dx / tol), round(dy / tol))
  tb  <- sort(table(key), decreasing = TRUE)
  best <- list(n = 0L, tx = 0, ty = 0)
  for (k in names(tb)[seq_len(min(5, length(tb)))]) {
    sel <- key == k; tx <- mean(dx[sel]); ty <- mean(dy[sel])
    n <- align_one_to_one(p$x + tx, p$y + ty, rX, rY, tol)
    if (n > best$n) best <- list(n = n, tx = tx, ty = ty)
  }
  best
}


#' Restrict a particle set to material anchors, falling back to the full set
#'
#' Material-based anchoring was a lab-sample convenience: on spiked samples you
#' know PET/PP are present in both instruments, and restricting to them raises
#' the fraction of true correspondences. On a field sample — or any sample
#' dominated by a polymer that is not on the list — the intersection can be
#' small or empty, and an empty anchor set is far worse than no filter at all.
#'
#' This mirrors the fallback the LDIR path already uses: apply the material
#' mask only when it leaves enough particles to actually anchor on, otherwise
#' keep the full cloud and say so.
#'
#' @param df        Particle data frame with a `material` column.
#' @param patterns  Character vector of case-insensitive regex patterns, or
#'                  NULL / empty to skip material filtering entirely.
#' @param min_count Minimum anchors required before the filter is honoured.
#' @param label     Dataset name used in log messages.
#' @return The filtered data frame, or `df` unchanged when the filter would
#'         leave fewer than `min_count` particles.
select_material_anchors <- function(df, patterns, min_count = 4, label = "") {
  if (is.null(patterns) || length(patterns) == 0) {
    log_message("  ", label, " anchors: no material filter configured — ",
                "using all ", nrow(df), " particles")
    return(df)
  }

  mask <- grepl(paste(patterns, collapse = "|"), df$material, ignore.case = TRUE)
  mask[is.na(mask)] <- FALSE
  n_match <- sum(mask)

  if (n_match < min_count) {
    log_message("  ", label, " anchors: only ", n_match, " particle(s) match ",
                "the configured anchor materials (", paste(patterns, collapse = ", "),
                ") — need >= ", min_count, ". Using all ", nrow(df),
                " particles instead; alignment is geometric and does not ",
                "require a specific polymer.", level = "WARN")
    return(df)
  }

  out <- df[mask, ]
  log_message("  ", label, " anchors: ", n_match, " of ", nrow(df),
              " particles (materials: ",
              paste(sort(unique(out$material)), collapse = ", "), ")")
  out
}

#' Count one-to-one inliers for a candidate transform
#'
#' Applies `M` to the source cloud and counts how many source particles pair
#' uniquely with a reference particle within `tol`. This is the acceptance
#' criterion that survives a collapse, unlike RMS.
#'
#' RMS actively REWARDS collapse: shrinking the source cloud packs points into
#' dense regions of the reference, so mean nearest-neighbour distance falls as
#' the fit gets more wrong. On a real LDIR<->Raman run the correct pose scored
#' RMS 301 um and 31 inliers, while a collapsed one scored RMS 67 um and 18
#' inliers. Selecting on RMS picks the wrong transform; selecting on inliers
#' picks the right one.
#'
#' @param M 3x3 transform matrix applied to the source cloud
#' @param src_x,src_y Source coordinates
#' @param ref_x,ref_y Reference coordinates
#' @param tol Pairing tolerance in the coordinate unit (um)
#' @return Integer inlier count, or NA when the inputs are unusable
align_transform_inliers <- function(M, src_x, src_y, ref_x, ref_y, tol) {
  if (is.null(M) || !is.matrix(M) || length(src_x) == 0 || length(ref_x) == 0)
    return(NA_integer_)
  p <- tryCatch(apply_transform_points(src_x, src_y, M), error = function(e) NULL)
  if (is.null(p)) return(NA_integer_)
  ok <- is.finite(p$x_transformed) & is.finite(p$y_transformed)
  if (!any(ok)) return(NA_integer_)
  align_one_to_one(p$x_transformed[ok], p$y_transformed[ok],
                   ref_x, ref_y, tol)
}


#' Chance of a random partner within r: 1 - exp(-rho * pi * r^2), with rho the
#' reference density over the area of its convex hull (Poisson, uniform).
#'
#' @param ref_x,ref_y Reference particle positions (um)
#' @param r Distance (um); a vector is allowed
#' @return list(p = probability per particle, density = particles per um^2,
#'   area_um2 = hull area, n = reference particles)
chance_partner_rate <- function(ref_x, ref_y, r) {
  ok <- is.finite(ref_x) & is.finite(ref_y)
  x <- ref_x[ok]; y <- ref_y[ok]; n <- length(x)
  if (n < 3) return(list(p = rep(NA_real_, length(r)), density = NA_real_,
                         area_um2 = NA_real_, n = n))
  h <- grDevices::chull(x, y)
  hx <- x[h]; hy <- y[h]
  area <- abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy)) / 2
  rho <- if (area > 0) n / area else NA_real_
  list(p = 1 - exp(-rho * pi * r^2), density = rho, area_um2 = area, n = n)
}


#' Is a coarse-search pose more than a search artefact? Two pairs fix a
#' similarity pose, so only pairs beyond two are evidence; each source anchor
#' finds a chance reference anchor within r with p (chance_partner_rate()),
#' and a search over correspondence pairs tries about (n_src x n_ref)^2 x 2
#' (mirror) poses (Bonferroni). Few anchors in a dense reference give p near
#' 1: the data cannot fix the pose.
#'
#' @param n_pairs One-to-one anchor pairs of the chosen pose
#' @param n_src Source anchors
#' @param ref_x,ref_y Reference anchors (um)
#' @param r Pair distance (um)
#' @return list(p, expected = chance pairs per pose, n_pairs)
pose_significance <- function(n_pairs, n_src, ref_x, ref_y, r) {
  cr <- chance_partner_rate(ref_x, ref_y, r)
  lam <- n_src * cr$p
  extra <- n_pairs - 2
  p <- if (!is.finite(lam) || n_pairs < 3) 1 else
    min(1, (n_src * cr$n)^2 * 2 * stats::ppois(extra - 1, lam, lower.tail = FALSE))
  list(p = p, expected = lam, n_pairs = n_pairs)
}


#' Which particles pass their instrument's agreement quality cutoff.
#'
#' The cutoff of the agreement scoring (.agreement_quality_cutoff(): FTIR AAU,
#' Raman HQI, LDIR quality). Without a cutoff, or without quality scores,
#' every particle passes.
#' @param df Particle table with `quality`
#' @param prefix Instrument prefix ("ftir", "ftir_bruker", "raman", "ldir", ...)
#' @return logical vector
quality_passes <- function(df, prefix, config) {
  if (is.null(df) || nrow(df) == 0) return(logical(0))
  cut <- .agreement_quality_cutoff(prefix, config)
  q <- suppressWarnings(as.numeric(df$quality))
  if (is.null(df$quality) || cut <= 0 || all(is.na(q))) return(rep(TRUE, nrow(df)))
  !is.na(q) & q >= cut
}


#' Tier 2 alignment anchors.
#'
#' Explicitly configured anchor materials (align_<role>_materials) win.
#' Otherwise, with align_anchor_plastics, the plastics (synthetic polymer
#' families) that pass the instrument's quality cutoff; with fewer than
#' align_min_anchor_count of those, every quality-passing particle; with
#' fewer of those, every particle.
#' @param prefix Instrument prefix for the quality cutoff
#' @return The anchor rows of df
select_alignment_anchors <- function(df, patterns, prefix, config, label = "") {
  min_n <- config$align_min_anchor_count %||% 4
  if ((!is.null(patterns) && length(patterns) > 0) ||
      !isTRUE(cfg_value(config, "align_anchor_plastics")))
    return(select_material_anchors(df, patterns, min_n, label))
  q_ok <- quality_passes(df, prefix, config)
  plastic <- classify_category_vec(classify_family_vec(df$material)) == "Synthetic"
  if (sum(q_ok & plastic) >= min_n) {
    log_message("  ", label, " anchors: ", sum(q_ok & plastic), " of ", nrow(df),
                " particles (plastics passing the quality cutoff)")
    return(df[q_ok & plastic, , drop = FALSE])
  }
  if (sum(q_ok) >= min_n) {
    log_message("  ", label, " anchors: only ", sum(q_ok & plastic), " plastic(s) pass ",
                "the quality cutoff; using all ", sum(q_ok), " particles that pass it",
                level = "WARN")
    return(df[q_ok, , drop = FALSE])
  }
  log_message("  ", label, " anchors: fewer than ", min_n, " particles pass the quality ",
              "cutoff; using all ", nrow(df), " particles", level = "WARN")
  df
}
