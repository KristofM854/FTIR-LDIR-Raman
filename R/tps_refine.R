# =============================================================================
# tps_refine.R — thin-plate-spline local refinement of LDIR->Raman alignment
# =============================================================================
# A single global similarity transform cannot overlay EVERY particle when the
# two instruments' coordinate systems differ by a small non-rigid distortion
# (stage nonlinearity, optical distortion, magnification drift across the
# field). Most particles land within the match gate; a few peripheral ones
# carry a local residual of a few hundred µm and stay unmatched even though
# they clearly correspond.
#
# This fits a regularized thin-plate spline to the DISPLACEMENT field
# (Raman - aligned-LDIR) sampled at the confident matched pairs, then applies
# it to all LDIR particles to correct each one's LOCAL offset using its
# confident neighbours — catching the stragglers without globally loosening
# the gate. Safeguards against distorting non-overlapping debris:
#   * the spline models the small residual, not absolute position;
#   * positions are normalized so the smoothing parameter is scale-free;
#   * the applied displacement is capped at the control-point residual scale,
#     so a particle far from any control (debris) cannot be flung onto a
#     spurious partner.
# =============================================================================

# Fit one TPS component mapping normalized control positions -> scalar values.
# lambda regularizes (0 = exact interpolation; larger = smoother/approximate).
.tps_fit_component <- function(nx, ny, val, lambda) {
  n <- length(nx)
  r2 <- outer(nx, nx, "-")^2 + outer(ny, ny, "-")^2
  K  <- ifelse(r2 > 0, 0.5 * r2 * log(r2), 0)     # U(r) = r^2 log r
  Pm <- cbind(1, nx, ny)
  A  <- rbind(cbind(K + lambda * diag(n), Pm),
              cbind(t(Pm), matrix(0, 3, 3)))
  b  <- c(val, 0, 0, 0)
  sol <- tryCatch(solve(A, b), error = function(e) {
    # fall back to a least-squares solve if the system is singular
    qr.solve(A, b)
  })
  list(w = sol[seq_len(n)], a = sol[n + 1:3])
}

#' Fit a regularized TPS displacement field from matched control pairs
#'
#' @param src_x,src_y Aligned LDIR coordinates at the matched pairs (µm)
#' @param dst_x,dst_y Raman coordinates at the matched pairs (µm)
#' @param lambda Smoothing (normalized units); larger = smoother
#' @return a fitted-warp object for tps_apply(), or NULL if under-determined
tps_fit_warp <- function(src_x, src_y, dst_x, dst_y, lambda = 0.5) {
  ok <- is.finite(src_x) & is.finite(src_y) & is.finite(dst_x) & is.finite(dst_y)
  sx <- src_x[ok]; sy <- src_y[ok]; dx <- dst_x[ok]; dy <- dst_y[ok]
  n <- length(sx)
  if (n < 4) return(NULL)                          # need >=4 controls
  # Reject near-collinear controls (TPS affine part degenerate)
  if (min(diff(range(sx)), diff(range(sy))) < 1e-6) return(NULL)

  cx <- mean(sx); cy <- mean(sy)
  s  <- sqrt(mean((sx - cx)^2 + (sy - cy)^2))
  if (!is.finite(s) || s <= 0) return(NULL)
  nx <- (sx - cx) / s; ny <- (sy - cy) / s          # normalized control coords

  # Model the residual displacement (Raman - aligned LDIR), not absolute pos
  vdx <- dx - sx; vdy <- dy - sy
  fx <- .tps_fit_component(nx, ny, vdx, lambda)
  fy <- .tps_fit_component(nx, ny, vdy, lambda)

  cap <- 1.5 * max(sqrt(vdx^2 + vdy^2), na.rm = TRUE)  # displacement cap (µm)
  structure(list(nx = nx, ny = ny, cx = cx, cy = cy, s = s,
                 fx = fx, fy = fy, cap = cap, n = n,
                 ctrl_res = sqrt(mean(vdx^2 + vdy^2))),
            class = "tps_warp")
}

.tps_eval <- function(fit, nx, ny) {
  r2 <- outer(nx, fit$nx_ctrl, "-")^2 + outer(ny, fit$ny_ctrl, "-")^2
  U  <- ifelse(r2 > 0, 0.5 * r2 * log(r2), 0)
  as.vector(fit$comp$a[1] + fit$comp$a[2] * nx + fit$comp$a[3] * ny +
            U %*% fit$comp$w)
}

#' Apply a fitted TPS warp to LDIR coordinates
#' @return list(x, y) of warped coordinates (same length as input)
tps_apply <- function(warp, x, y) {
  if (is.null(warp)) return(list(x = x, y = y))
  nx <- (x - warp$cx) / warp$s; ny <- (y - warp$cy) / warp$s
  ev <- function(comp) .tps_eval(
    list(nx_ctrl = warp$nx, ny_ctrl = warp$ny, comp = comp), nx, ny)
  ddx <- ev(warp$fx); ddy <- ev(warp$fy)
  # Cap displacement magnitude so far-from-control particles can't be flung
  mag <- sqrt(ddx^2 + ddy^2)
  scale <- ifelse(mag > warp$cap & mag > 0, warp$cap / mag, 1)
  list(x = x + ddx * scale, y = y + ddy * scale)
}
