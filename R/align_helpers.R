# =============================================================================
# align_helpers.R — shared primitives for coarse pose search & inlier counting
# =============================================================================
# Extracted from global_register_align() (04_ransac.R). The alignment core and
# the tools/diagnose_*.R scripts previously each carried their own copy of the
# rotation/scale pose transform, the one-to-one inlier count, the robust span,
# and the translation-voting score. They now share this one implementation.
#
# All functions are pure (no global state, no RNG) so they are safe to reuse
# from any context.

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
