# =============================================================================
# ftir_image_placement.R -- place an FTIR chemical image under its particles
# =============================================================================
#
# ONE function decides where the FTIR "Average Abs" image (or a Bruker FTIR
# image) sits in the native FTIR coordinate frame: place_ftir_image(). The
# Shiny FTIR / FTIR (Bruker) tabs, the PDF/HTML report and the Multi-Run
# backdrop all call it, and tools/diag_ftir_overlay.R measures it.
#
# Why this exists. The viewer used to stretch the PNG onto the bounding box of
# the particle coordinates (min/max of x_orig / y_orig). The particles never
# reach the edges of the scan, so the image was squeezed into a box that is
# smaller than the area it actually covers: the misfit is a per-axis SCALE
# error that grows with distance from a fixed point, not a constant shift.
# The bounding box also depended on which particles were in the set, so the
# report (which filters first) placed the image differently from the viewer.
#
# The PNG export carries no geometry metadata (no origin, no pixel size), so
# the extent has to come from one of:
#
#   P1  config / run metadata: ftir_image_width_um + ftir_image_height_um
#       (+ optional ftir_image_center_x_um / _y_um, default = scan origin at
#       0,0). Exact and resize-invariant when known. Describes the WHOLE PNG.
#   P2  registration: detect the particles in the image itself (threshold +
#       connected components), match them to the CSV particles and fit the
#       pixel -> um mapping per axis (scale + offset). Margins, colour bars and
#       an unknown origin drop out because the fit measures the mapping
#       directly. Fitted on ALL particles, so it does not move with filters.
#   P3  legacy fallback: aspect-preserving fit to the particle extent. Only
#       approximate (the particles do not reach the scan edges); flagged.
#
# Pixel convention: an image of W x H pixels spans the edge coordinates
# [0, W] x [0, H]; pixel column c (1-based) covers [c-1, c] and its centre is
# c - 0.5. annotation_raster() uses the same convention (column 1 fills
# [xmin, xmin + dx]), so an extent computed from edge coordinates places the
# raster with no half-pixel error. Row 1 is the TOP of the image; when the
# CSV Y axis points up (the normal case) y_um decreases with row.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# Image -> grayscale intensity
# ---------------------------------------------------------------------------

#' Grayscale intensity in [0, 1] from an H x W (x channels) raster.
#' Colour-agnostic mean of R, G, B (alpha ignored) -- the same reduction
#' R/01b_ingest_image.R uses.
.ftir_gray <- function(raw) {
  d <- dim(raw)
  if (length(d) == 2) return(raw * 1)
  if (d[3] >= 3) return((raw[, , 1] + raw[, , 2] + raw[, , 3]) / 3)
  raw[, , 1] * 1
}


# ---------------------------------------------------------------------------
# Connected components on a binary image, by run-length encoding
# ---------------------------------------------------------------------------
# Labels 8-connected foreground runs row by row and unions runs that touch the
# previous row. Cost scales with the number of runs, not the number of pixels,
# so a sparse 3000 x 3000 absorption map is labelled in about a second without
# compiled code (01b's per-pixel two-pass labeller takes minutes at that size).
#
# Returns a data frame with one row per component:
#   cx_px, cy_px  centroid in EDGE pixel coordinates (see header)
#   area_px       pixel count
#   col_min, col_max, row_min, row_max   1-based bounding box
.ftir_label_runs <- function(fg) {
  h <- nrow(fg)
  run_row <- list(); run_s <- list(); run_e <- list()
  k <- 0L
  for (r in seq_len(h)) {
    v <- fg[r, ]
    if (!any(v)) next
    dv <- diff(c(FALSE, v, FALSE))
    s <- which(dv == 1L); e <- which(dv == -1L) - 1L
    k <- k + 1L
    run_row[[k]] <- rep.int(r, length(s)); run_s[[k]] <- s; run_e[[k]] <- e
  }
  if (k == 0L) return(NULL)
  row <- unlist(run_row); s <- unlist(run_s); e <- unlist(run_e)
  n <- length(row)

  # Union-find over runs.
  parent <- seq_len(n)
  find <- function(i) {
    root <- i
    while (parent[root] != root) root <- parent[root]
    while (parent[i] != root) { nx <- parent[i]; parent[i] <<- root; i <- nx }
    root
  }
  # Runs are ordered by row then by start, so each row is a contiguous block.
  row_first <- match(seq_len(h), row)
  row_last  <- n + 1L - match(seq_len(h), rev(row))
  for (r in seq_len(h)[-1]) {
    if (is.na(row_first[r]) || is.na(row_first[r - 1L])) next
    cur  <- row_first[r]:row_last[r]
    prev <- row_first[r - 1L]:row_last[r - 1L]
    pe <- e[prev]; ps <- s[prev]
    for (i in cur) {
      # 8-connectivity: previous-row runs overlapping [s-1, e+1].
      lo <- findInterval(s[i] - 1L, pe, left.open = TRUE) + 1L
      hi <- findInterval(e[i] + 1L, ps)
      if (lo > hi) next
      ri <- find(i)
      for (j in prev[lo:hi]) {
        rj <- find(j)
        if (ri != rj) { parent[max(ri, rj)] <- min(ri, rj); ri <- min(ri, rj) }
      }
    }
  }
  lab <- vapply(seq_len(n), find, integer(1))
  lab <- match(lab, unique(lab))

  len  <- e - s + 1L
  area <- as.numeric(rowsum(len, lab))
  # Sum of pixel-centre columns over a run [s, e]: len * ((s + e) / 2 - 0.5)
  sx <- as.numeric(rowsum(len * ((s + e) / 2 - 0.5), lab))
  sy <- as.numeric(rowsum(len * (row - 0.5), lab))
  data.frame(
    cx_px   = sx / area,
    cy_px   = sy / area,
    area_px = area,
    col_min = as.numeric(tapply(s, lab, min)),
    col_max = as.numeric(tapply(e, lab, max)),
    row_min = as.numeric(tapply(row, lab, min)),
    row_max = as.numeric(tapply(row, lab, max))
  )
}


#' Detect particle blobs in an FTIR chemical image
#'
#' Simple and deliberately dumb: global threshold on grayscale intensity,
#' connected components, area centroids. The background level is the median
#' pixel and the noise level its MAD; a pixel is foreground when it deviates
#' from the background by more than `k_mad` MADs (and by at least `min_contrast`
#' of the intensity range). Polarity ("bright" particles on a dark map, or the
#' reverse) is picked automatically as the side with the longer tail.
#'
#' Components larger than `max_frac` of the shorter image side are dropped:
#' those are frames, colour bars and legends, not particles.
#'
#' @param raw H x W (x channels) raster in [0, 1] (png::readPNG layout)
#' @return data frame from .ftir_label_runs() (+ attr "polarity"), or NULL
detect_image_blobs <- function(raw, polarity = c("auto", "bright", "dark"),
                               k_mad = 6, min_contrast = 0.08,
                               min_px = 4, max_frac = 0.2) {
  polarity <- match.arg(polarity)
  g <- .ftir_gray(raw)
  if (is.null(dim(g)) || any(dim(g) < 2)) return(NULL)
  med <- stats::median(g)
  q <- stats::quantile(g, c(0.0005, 0.9995), names = FALSE)
  if (polarity == "auto")
    polarity <- if ((q[2] - med) >= (med - q[1])) "bright" else "dark"
  dev <- if (polarity == "bright") g - med else med - g
  spread <- stats::mad(g)
  tail <- if (polarity == "bright") q[2] - med else med - q[1]
  thr <- max(k_mad * spread, min_contrast * tail, 1e-6)
  fg <- dev > thr
  if (!any(fg)) return(NULL)
  cc <- .ftir_label_runs(fg)
  if (is.null(cc)) return(NULL)
  lim <- max_frac * min(dim(g))
  keep <- cc$area_px >= min_px &
          (cc$col_max - cc$col_min + 1) <= lim &
          (cc$row_max - cc$row_min + 1) <= lim
  cc <- cc[keep, , drop = FALSE]
  if (nrow(cc) == 0) return(NULL)
  rownames(cc) <- NULL
  attr(cc, "polarity")  <- polarity
  attr(cc, "threshold") <- thr
  cc
}


# ---------------------------------------------------------------------------
# Nearest neighbours and mutual matching
# ---------------------------------------------------------------------------

.nn1 <- function(ref, query) {
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(ref, query, k = 1)
    return(list(idx = as.integer(nn$nn.idx[, 1]), dist = as.numeric(nn$nn.dists[, 1])))
  }
  idx <- integer(nrow(query)); dist <- numeric(nrow(query))
  for (i in seq_len(nrow(query))) {
    d2 <- (ref[, 1] - query[i, 1])^2 + (ref[, 2] - query[i, 2])^2
    idx[i] <- which.min(d2); dist[i] <- sqrt(d2[idx[i]])
  }
  list(idx = idx, dist = dist)
}

#' Mutual nearest-neighbour pairs between point sets a and b within `gate`.
#' Returns data frame(ia, ib, dist).
mutual_nn_pairs <- function(a, b, gate = Inf) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (nrow(a) == 0 || nrow(b) == 0) return(data.frame(ia = integer(), ib = integer(), dist = numeric()))
  ab <- .nn1(b, a); ba <- .nn1(a, b)
  ia <- seq_len(nrow(a))
  ok <- ba$idx[ab$idx] == ia & ab$dist <= gate
  data.frame(ia = ia[ok], ib = ab$idx[ok], dist = ab$dist[ok])
}


# ---------------------------------------------------------------------------
# Registration: pixel -> um, per axis
# ---------------------------------------------------------------------------
#   x_um = ax + bx * cx_px
#   y_um = ay + by * cy_px      (by < 0 when CSV Y points up, the usual case)

.apply_px_map <- function(m, cx, cy) cbind(m$ax + m$bx * cx, m$ay + m$by * cy)

.fit_px_map <- function(cx, cy, x, y) {
  fx <- stats::lm.fit(cbind(1, cx), x)$coefficients
  fy <- stats::lm.fit(cbind(1, cy), y)$coefficients
  list(ax = fx[[1]], bx = fx[[2]], ay = fy[[1]], by = fy[[2]])
}

#' Register detected blobs to the CSV particle coordinates
#'
#' Initial guesses: (a) the blobs' bounding box maps onto the particles'
#' bounding box (both are bounded by the same outermost particles, whatever
#' margins the PNG carries); (b) translation voting between the largest blobs
#' and the largest particles over a grid of scales, which also works when the
#' particle set covers only part of the image. Each is tried with the image Y
#' axis pointing down and up; the fit matching most particles wins. Refinement:
#' iterate mutual-nearest-neighbour matching with a shrinking gate and a
#' trimmed per-axis least-squares refit. The coarse stage matches only the
#' larger half of both sets -- big particles are unambiguous -- then the fit
#' is refined on everything.
#'
#' @param blobs  output of detect_image_blobs()
#' @param x, y   CSV particle coordinates (um), ALL particles of the run
#' @param size_um optional particle size (Feret max, um) for the large-first stage
#' @param img_w, img_h raster dimensions (pixels)
#' @return list(ok, ax, bx, ay, by, n_matched, rms_um, rms_px, pairs, ...)
register_blobs_to_particles <- function(blobs, x, y, size_um = NULL,
                                        img_w, img_h,
                                        min_matches = 6L, max_rms_px = 3,
                                        max_aniso = 0.2, n_iter = 25L) {
  fail <- function(reason) list(ok = FALSE, reason = reason)
  okp <- is.finite(x) & is.finite(y)
  x <- x[okp]; y <- y[okp]
  size_um <- if (is.null(size_um)) rep(NA_real_, length(okp))[okp] else size_um[okp]
  if (is.null(blobs) || nrow(blobs) < min_matches) return(fail("too few image blobs"))
  if (length(x) < min_matches) return(fail("too few particles"))

  P  <- cbind(x, y)
  big_p <- if (any(is.finite(size_um)))
    which(size_um >= stats::median(size_um, na.rm = TRUE) | !is.finite(size_um))
  else seq_along(x)
  big_b <- which(blobs$area_px >= stats::median(blobs$area_px))

  run_from <- function(m0) {
    m <- m0
    span <- max(diff(range(x)), diff(range(y)))
    gate <- 0.05 * span
    pairs <- NULL
    for (it in seq_len(n_iter)) {
      coarse <- it <= 5L
      bi <- if (coarse) big_b else seq_len(nrow(blobs))
      pi <- if (coarse) big_p else seq_along(x)
      B <- .apply_px_map(m, blobs$cx_px[bi], blobs$cy_px[bi])
      pr <- mutual_nn_pairs(B, P[pi, , drop = FALSE], gate)
      if (nrow(pr) < 3) return(NULL)
      pr$ib_blob <- bi[pr$ia]; pr$ip <- pi[pr$ib]
      # Trim gross outliers (wrong pairings) before refitting.
      if (nrow(pr) >= 8) {
        cut <- stats::median(pr$dist) + 3 * stats::mad(pr$dist)
        pr <- pr[pr$dist <= max(cut, 1e-9), , drop = FALSE]
      }
      m_new <- .fit_px_map(blobs$cx_px[pr$ib_blob], blobs$cy_px[pr$ib_blob],
                           x[pr$ip], y[pr$ip])
      if (any(!is.finite(unlist(m_new)))) return(NULL)
      m <- m_new
      res <- sqrt(stats::median(pr$dist^2))
      gate <- max(min(gate, 4 * res), 3 * abs(m$bx), 1e-6)
      pairs <- pr
    }
    # Final residuals under the final map, on all mutual pairs within gate.
    B  <- .apply_px_map(m, blobs$cx_px, blobs$cy_px)
    pr <- mutual_nn_pairs(B, P, gate)
    if (nrow(pr) < 3) return(NULL)
    list(m = m, pairs = data.frame(blob = pr$ia, particle = which(okp)[pr$ib],
                                   dist = pr$dist),
         rms_um = sqrt(mean(pr$dist^2)), n = nrow(pr))
  }

  bx_rng <- range(blobs$cx_px); by_rng <- range(blobs$cy_px)
  if (diff(bx_rng) <= 0 || diff(by_rng) <= 0) return(fail("degenerate blob extent"))
  inits <- list()
  for (lo_hi in list(c(0, 1), c(0.02, 0.98))) {
    qbx <- stats::quantile(blobs$cx_px, lo_hi, names = FALSE)
    qby <- stats::quantile(blobs$cy_px, lo_hi, names = FALSE)
    qx  <- stats::quantile(x, lo_hi, names = FALSE)
    qy  <- stats::quantile(y, lo_hi, names = FALSE)
    bx  <- diff(qx) / diff(qbx)
    for (ysign in c(-1, 1)) {
      by <- ysign * diff(qy) / diff(qby)
      ay <- if (ysign < 0) qy[2] - by * qby[1] else qy[1] - by * qby[1]
      inits[[length(inits) + 1]] <- list(ax = qx[1] - bx * qbx[1], bx = bx,
                                         ay = ay, by = by)
    }
  }
  # Voting inits, robust to a particle set that covers only part of the image
  # (the bounding boxes then disagree): for a grid of isotropic scales around
  # the bounding-box estimate, every pairing of the K largest blobs with the K
  # largest particles votes for a translation; the densest vote wins.
  K <- min(20L, nrow(blobs), length(x))
  kb <- order(-blobs$area_px)[seq_len(K)]
  kp <- if (any(is.finite(size_um))) order(-size_um)[seq_len(K)] else seq_len(K)
  s0 <- sqrt(abs(inits[[1]]$bx * inits[[1]]$by))
  span <- max(diff(range(x)), diff(range(y)))
  votes <- list()
  for (s in s0 * 2^seq(-1.5, 1.5, by = 0.05)) for (ysign in c(-1, 1)) {
    ox <- outer(x[kp], s * blobs$cx_px[kb], "-")
    oy <- outer(y[kp], ysign * s * blobs$cy_px[kb], "-")
    bin <- max(0.01 * span, 2 * s)
    key <- paste(round(ox / bin), round(oy / bin))
    tb <- table(key)
    top <- names(tb)[which.max(tb)]
    sel <- key == top
    votes[[length(votes) + 1]] <- list(n = max(tb), m = list(
      ax = mean(ox[sel]), bx = s, ay = mean(oy[sel]), by = ysign * s))
  }
  nv <- vapply(votes, `[[`, numeric(1), "n")
  for (v in votes[order(-nv)][seq_len(min(4L, length(votes)))]) inits[[length(inits) + 1]] <- v$m

  fits <- Filter(Negate(is.null), lapply(inits, function(m0)
    tryCatch(run_from(m0), error = function(e) NULL)))
  if (length(fits) == 0) return(fail("matching did not converge"))
  score <- vapply(fits, function(f) f$n - f$rms_um / (abs(f$m$bx) * 1e3), numeric(1))
  best <- fits[[which.max(score)]]
  m <- best$m
  px_um  <- sqrt(abs(m$bx * m$by))
  rms_px <- best$rms_um / px_um
  aniso  <- abs(abs(m$bx) / abs(m$by) - 1)
  out <- c(m, list(n_matched = best$n, rms_um = best$rms_um, rms_px = rms_px,
                   anisotropy = aniso, img_w = img_w, img_h = img_h,
                   n_blobs = nrow(blobs), n_particles = length(x),
                   pairs = best$pairs))
  need <- max(min_matches, ceiling(0.25 * min(nrow(blobs), length(x))))
  out$ok <- best$n >= need && rms_px <= max_rms_px && aniso <= max_aniso &&
            m$bx > 0
  if (!out$ok)
    out$reason <- sprintf("rejected: n=%d (need %d), rms=%.2f px (max %.1f), anisotropy=%.1f%%",
                          best$n, need, rms_px, max_rms_px, 100 * aniso)
  out
}


#' Extent (edge coordinates) of the full raster under a pixel -> um map.
#' Returns list(xmin, xmax, ymin, ymax, flip_rows); flip_rows = TRUE when the
#' image Y axis points the same way as the um Y axis, i.e. raster row 1 must
#' be drawn at the BOTTOM.
px_map_extent <- function(m, img_w, img_h) {
  xs <- m$ax + m$bx * c(0, img_w)
  ys <- m$ay + m$by * c(0, img_h)
  list(xmin = min(xs), xmax = max(xs), ymin = min(ys), ymax = max(ys),
       flip_rows = m$by > 0)
}


# ---------------------------------------------------------------------------
# Config / metadata extent (P1)
# ---------------------------------------------------------------------------

.cfg_num <- function(cfg, key) {
  if (is.null(cfg)) return(NA_real_)
  v <- if (is.data.frame(cfg)) { if (key %in% names(cfg)) cfg[[key]][1] else NULL }
       else cfg[[key]]
  v <- suppressWarnings(as.numeric(unlist(v)))
  if (length(v) == 0 || !is.finite(v[1])) NA_real_ else v[1]
}

#' Physical extent of the whole FTIR PNG from config/run metadata, or NULL.
#' prefix "ftir_image" (PerkinElmer) or "ftir_bruker_image" (Bruker).
ftir_extent_from_config <- function(cfg, prefix = "ftir_image") {
  w <- .cfg_num(cfg, paste0(prefix, "_width_um"))
  h <- .cfg_num(cfg, paste0(prefix, "_height_um"))
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) return(NULL)
  cx <- .cfg_num(cfg, paste0(prefix, "_center_x_um"))
  cy <- .cfg_num(cfg, paste0(prefix, "_center_y_um"))
  if (!is.finite(cx)) cx <- w / 2      # default: scan origin at (0, 0)
  if (!is.finite(cy)) cy <- h / 2
  list(xmin = cx - w / 2, xmax = cx + w / 2, ymin = cy - h / 2, ymax = cy + h / 2)
}


# ---------------------------------------------------------------------------
# P3: aspect-preserving particle-extent fit (approximate)
# ---------------------------------------------------------------------------
.particle_extent_fit <- function(raw, x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NULL)
  if (is.null(raw)) return(list(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)))
  sx <- diff(range(x)); sy <- diff(range(y))
  if (sx <= 0 && sy <= 0) return(NULL)
  asp <- ncol(raw) / nrow(raw)
  if (sy <= 0 || sx / sy > asp) sy <- sx / asp else sx <- sy * asp
  cx <- mean(range(x)); cy <- mean(range(y))
  list(xmin = cx - sx / 2, xmax = cx + sx / 2, ymin = cy - sy / 2, ymax = cy + sy / 2)
}


# ---------------------------------------------------------------------------
# The single entry point
# ---------------------------------------------------------------------------

#' Raster pixel dimensions the extent refers to. downsample_raster() (Shiny
#' upload path) block-averages by an integer factor k and drops the remainder
#' rows/columns; it records the original size as attributes. Returns the
#' original size and the fraction of it the raster still covers.
.raster_coverage <- function(raw) {
  w <- ncol(raw); h <- nrow(raw)
  ow <- attr(raw, "orig_width_px") %||% w
  oh <- attr(raw, "orig_height_px") %||% h
  k <- max(1, round(ow / w))
  list(orig_w = ow, orig_h = oh,
       fx = min(1, (w * k) / ow), fy = min(1, (h * k) / oh))
}

#' Place an FTIR image in the native FTIR (CSV) coordinate frame
#'
#' @param raw   H x W (x channels) raster
#' @param x, y  CSV particle coordinates (um). Pass the FULL particle set of the
#'              run -- the result must not depend on display filters. Used by
#'              P2 (registration, when no stored one is given) and P3.
#' @param size_um optional Feret max per particle (helps P2 match big ones first)
#' @param cfg   config list / one-row meta data frame (P1 fields)
#' @param registration stored registration (read_ftir_image_placement()) to
#'              reuse; recomputed when NULL or when it refers to another image size
#' @param prefix "ftir_image" or "ftir_bruker_image" (P1 config field prefix)
#' @param register FALSE skips P2 (e.g. no image analysis wanted)
#' @return list(raster, xmin, xmax, ymin, ymax, method, registration) or NULL
place_ftir_image <- function(raw, x, y, size_um = NULL, cfg = NULL,
                             registration = NULL, prefix = "ftir_image",
                             register = TRUE) {
  if (is.null(raw)) return(NULL)
  cov <- .raster_coverage(raw)
  finish <- function(ext, method, reg = NULL, flip_rows = FALSE) {
    r <- raw
    if (isTRUE(flip_rows)) {
      r <- if (length(dim(raw)) == 2) raw[nrow(raw):1, , drop = FALSE]
           else raw[nrow(raw):1, , , drop = FALSE]
    }
    # A downsampled raster that dropped remainder pixels covers slightly
    # less than the original image: trim the extent to match (left/top
    # edges are kept, right/bottom edges move in).
    if (cov$fx < 1) ext$xmax <- ext$xmin + (ext$xmax - ext$xmin) * cov$fx
    if (cov$fy < 1) {
      span <- (ext$ymax - ext$ymin) * cov$fy
      if (isTRUE(flip_rows)) ext$ymax <- ext$ymin + span else ext$ymin <- ext$ymax - span
    }
    list(raster = r, xmin = ext$xmin, xmax = ext$xmax,
         ymin = ext$ymin, ymax = ext$ymax, method = method,
         registration = reg)
  }

  # P1: explicit physical extent.
  ext <- ftir_extent_from_config(cfg, prefix)
  if (!is.null(ext)) return(finish(ext, "config_extent"))

  # P2: registration (stored, else computed now on the full particle set).
  reg <- registration
  if (!is.null(reg) && (!isTRUE(reg$ok) ||
      !isTRUE(all.equal(c(reg$img_w, reg$img_h), c(cov$orig_w, cov$orig_h)))))
    reg <- NULL
  if (is.null(reg) && isTRUE(register) && length(x) > 0) {
    reg <- tryCatch(register_ftir_image(raw, x, y, size_um),
                    error = function(e) NULL)
    if (!is.null(reg) && !isTRUE(reg$ok)) reg <- NULL
  }
  if (!is.null(reg)) {
    ext <- px_map_extent(reg, reg$img_w, reg$img_h)
    return(finish(ext[c("xmin", "xmax", "ymin", "ymax")], "registration",
                  reg, flip_rows = ext$flip_rows))
  }

  # P3: approximate fallback.
  ext <- .particle_extent_fit(raw, x, y)
  if (is.null(ext)) return(NULL)
  finish(ext, "particle_extent_approx")
}

#' Detect + register in one call, on the raster's ORIGINAL pixel grid.
#' For a downsampled raster the blob coordinates are rescaled to the original
#' pixel grid so the stored map stays valid for the full-resolution file.
register_ftir_image <- function(raw, x, y, size_um = NULL, ...) {
  blobs <- detect_image_blobs(raw)
  cov <- .raster_coverage(raw)
  if (!is.null(blobs)) {
    kx <- cov$orig_w / ncol(raw) * cov$fx; ky <- cov$orig_h / nrow(raw) * cov$fy
    blobs$cx_px <- blobs$cx_px * kx; blobs$cy_px <- blobs$cy_px * ky
  }
  register_blobs_to_particles(blobs, x, y, size_um,
                              img_w = cov$orig_w, img_h = cov$orig_h, ...)
}


# ---------------------------------------------------------------------------
# Persistence: computed once by the pipeline, reused by every plot
# ---------------------------------------------------------------------------

FTIR_PLACEMENT_FILES <- c(ftir = "ftir_image_placement.json",
                          ftir_bruker = "ftir_bruker_image_placement.json")

write_ftir_image_placement <- function(placement, path) {
  if (is.null(placement) || !requireNamespace("jsonlite", quietly = TRUE))
    return(invisible(NULL))
  reg <- placement$registration
  rec <- list(
    method = placement$method,
    extent_um = list(xmin = placement$xmin, xmax = placement$xmax,
                     ymin = placement$ymin, ymax = placement$ymax),
    registration = if (!is.null(reg)) reg[c("ok", "ax", "bx", "ay", "by",
      "n_matched", "rms_um", "rms_px", "anisotropy", "img_w", "img_h",
      "n_blobs", "n_particles")] else NULL)
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE,
                              digits = NA, null = "null"), path)
  invisible(path)
}

#' Stored registration (list with ok/ax/bx/ay/by/img_w/img_h...) or NULL.
read_ftir_image_placement <- function(path) {
  if (is.null(path) || !file.exists(path) ||
      !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  rec <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE),
                  error = function(e) NULL)
  reg <- rec$registration
  if (is.null(reg) || !isTRUE(reg$ok)) return(NULL)
  reg
}
