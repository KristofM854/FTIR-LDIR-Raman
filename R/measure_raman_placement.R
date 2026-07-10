# =============================================================================
# measure_raman_placement.R — measure the Raman image's physical extent
# =============================================================================
# Shared core used by BOTH the pipeline (main.R, automatic calibration) and
# the CLI diagnostic (tools/diagnose_raman_placement.R) so the two can never
# drift.  Given a Raman micrograph and the Raman particle stage coordinates,
# it finds the scale + translation that best places the image under the
# particles, scoring each candidate by the fraction of ALL particles landing
# on a bright pixel (particles are bright blobs on a dark membrane).
#
# The WITec panel Width/Height seed the scale search; the returned extent may
# differ from them when the export is cropped/zoomed (scale != 1).
# =============================================================================

# 5x5-style separable running-max dilation (radius in px)
.mrp_run_max <- function(mat, radius, along_rows) {
  out <- mat
  for (k in setdiff(seq(-radius, radius), 0)) {
    n <- if (along_rows) nrow(mat) else ncol(mat)
    idx <- pmin(pmax(seq_len(n) + k, 1), n)
    out <- if (along_rows) pmax(out, mat[idx, , drop = FALSE])
           else                 pmax(out, mat[, idx, drop = FALSE])
  }
  out
}
.mrp_dilate <- function(mat, radius)
  .mrp_run_max(.mrp_run_max(mat, radius, TRUE), radius, FALSE)

#' Measure the best-fit physical extent of a Raman image under its particles
#'
#' @param lum    Numeric matrix [H, W] of image luminance in [0,1]
#' @param x,y    Particle stage coordinates (µm); non-finite entries ignored
#' @param W,H    WITec panel Width/Height (µm) — the scale-search seed
#' @param min_frac Minimum bright-fraction to accept the fit
#' @return list(width_um, height_um, center_x_um, center_y_um, scale,
#'   frac_bright, mirrored) or NULL when no confident, non-mirrored fit found.
measure_raman_placement_core <- function(lum, x, y, W, H, min_frac = 0.6) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n_total <- length(x)
  if (n_total < 4 || is.null(dim(lum)) || W <= 0 || H <= 0) return(NULL)

  Hpx <- nrow(lum); Wpx <- ncol(lum)
  lum_max    <- .mrp_dilate(lum, 2)
  r_coarse   <- max(4L, as.integer(ceiling(min(Hpx, Wpx) / 60)))
  lum_coarse <- .mrp_dilate(lum, r_coarse)

  extent_center <- function(cx, cy, w, h)
    list(xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2, ymax = cy + h/2)

  # frac_bright over ALL particles (outside the extent counts as a miss),
  # so a shrunken placement covering a few lucky particles cannot win.
  score <- function(ext, flip_h = FALSE, flip_v = FALSE, mat = lum_max) {
    fx <- (x - ext$xmin) / (ext$xmax - ext$xmin)
    fy <- (ext$ymax - y) / (ext$ymax - ext$ymin)
    if (flip_h) fx <- 1 - fx
    if (flip_v) fy <- 1 - fy
    inside <- fx >= 0 & fx <= 1 & fy >= 0 & fy <= 1
    if (!any(inside)) return(c(frac = 0, mean = 0))
    ci <- pmin(pmax(ceiling(fx[inside] * Wpx), 1), Wpx)
    ri <- pmin(pmax(ceiling(fy[inside] * Hpx), 1), Hpx)
    v <- mat[cbind(ri, ci)]
    c(frac = sum(v > 0.5) / n_total, mean = mean(v))
  }
  better <- function(a, b)
    a["frac"] > b$frac || (a["frac"] == b$frac && a["mean"] > b$mean)

  pcx <- mean(range(x)); pcy <- mean(range(y))

  # Coarse: scale x mirror x translation on the coarsely-dilated luminance
  coarse_hits <- list()
  for (fh in c(FALSE, TRUE)) for (fv in c(FALSE, TRUE)) {
    for (s in seq(0.3, 2.4, by = 0.1)) {
      w_s <- W * s; h_s <- H * s
      step <- r_coarse * (w_s / Wpx); span <- max(w_s, h_s)
      hit <- list(frac = -1, mean = -1)
      for (dx in seq(-span/2, span/2, by = step))
        for (dy in seq(-span/2, span/2, by = step)) {
          sc <- score(extent_center(pcx + dx, pcy + dy, w_s, h_s), fh, fv,
                      mat = lum_coarse)
          if (better(sc, hit))
            hit <- list(frac = sc[["frac"]], mean = sc[["mean"]],
                        flip_h = fh, flip_v = fv, scale = s,
                        cx = pcx + dx, cy = pcy + dy)
        }
      coarse_hits[[length(coarse_hits) + 1]] <- hit
    }
  }
  ord <- order(-vapply(coarse_hits, `[[`, 0, "frac"),
               -vapply(coarse_hits, `[[`, 0, "mean"))
  top <- coarse_hits[ord[seq_len(min(6, length(ord)))]]

  # Fine: refine each leading coarse hit on the 5x5 luminance
  best <- list(frac = -1, mean = -1)
  for (h in top) {
    for (s in seq(h$scale - 0.06, h$scale + 0.06, by = 0.02)) {
      w_s <- W * s; h_s <- H * s
      for (dx in seq(-3*r_coarse, 3*r_coarse, by = max(1, r_coarse/4)) * (w_s/Wpx))
        for (dy in seq(-3*r_coarse, 3*r_coarse, by = max(1, r_coarse/4)) * (h_s/Hpx)) {
          sc <- score(extent_center(h$cx + dx, h$cy + dy, w_s, h_s),
                      h$flip_h, h$flip_v)
          if (better(sc, best))
            best <- list(frac = sc[["frac"]], mean = sc[["mean"]],
                         flip_h = h$flip_h, flip_v = h$flip_v, scale = s,
                         cx = h$cx + dx, cy = h$cy + dy)
        }
    }
  }
  # Sub-cell polish
  for (dx in seq(-2, 2, by = 0.5) * (W * best$scale / Wpx))
    for (dy in seq(-2, 2, by = 0.5) * (H * best$scale / Hpx)) {
      sc <- score(extent_center(best$cx + dx, best$cy + dy,
                                W * best$scale, H * best$scale),
                  best$flip_h, best$flip_v)
      if (better(sc, best)) {
        best$frac <- sc[["frac"]]; best$mean <- sc[["mean"]]
        best$cx <- best$cx + dx; best$cy <- best$cy + dy
      }
    }

  if (best$frac < min_frac) return(NULL)
  list(width_um    = W * best$scale,
       height_um   = H * best$scale,
       center_x_um = best$cx,
       center_y_um = best$cy,
       scale       = best$scale,
       frac_bright = best$frac,
       mirrored    = isTRUE(best$flip_h) || isTRUE(best$flip_v),
       flip_h      = isTRUE(best$flip_h),
       flip_v      = isTRUE(best$flip_v))
}

#' Read an image to a luminance matrix (PNG fast path, else read_image_any)
#' @return [H, W] numeric matrix in [0,1], or NULL on failure
read_image_luminance <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  arr <- NULL
  ext <- tolower(tools::file_ext(path))
  if (ext == "png" && requireNamespace("png", quietly = TRUE))
    arr <- tryCatch(png::readPNG(path), error = function(e) NULL)
  if (is.null(arr) && exists("read_image_any"))
    arr <- tryCatch(read_image_any(path, verbose = FALSE), error = function(e) NULL)
  if (is.null(arr)) return(NULL)
  if (length(dim(arr)) == 3)
    0.2126 * arr[,,1] + 0.7152 * arr[,,2] + 0.0722 * arr[,,3]
  else arr
}

#' Measure Raman image placement from a file path (pipeline entry point)
#'
#' @param image_path Path to the Raman micrograph (canonical PNG preferred)
#' @param x,y Raman particle stage coordinates (µm)
#' @param W,H WITec panel Width/Height (µm)
#' @param min_frac Minimum bright-fraction to accept
#' @return same as measure_raman_placement_core(), or NULL
measure_raman_image_placement <- function(image_path, x, y, W, H,
                                          min_frac = 0.6) {
  if (is.null(W) || is.null(H) || !is.numeric(W) || !is.numeric(H) ||
      W <= 0 || H <= 0) return(NULL)
  lum <- read_image_luminance(image_path)
  if (is.null(lum)) return(NULL)
  measure_raman_placement_core(lum, x, y, W, H, min_frac = min_frac)
}

#' Merge key/value pairs into a run manifest's config_snapshot in place
#'
#' Reads the manifest JSON, updates config_snapshot with `updates`, and writes
#' it back (keeping a .bak on first change). Used to persist an auto-measured
#' Raman extent so the viewer places the image correctly.
#' @return TRUE on success, FALSE otherwise (never throws)
update_manifest_config_snapshot <- function(run_dir, updates) {
  tryCatch({
    mp <- file.path(run_dir, "00_manifest", "manifest.json")
    if (!file.exists(mp)) mp <- file.path(run_dir, "manifest.json")
    if (!file.exists(mp)) return(FALSE)
    if (!requireNamespace("jsonlite", quietly = TRUE)) return(FALSE)
    man <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
    if (is.null(man$config_snapshot)) man$config_snapshot <- list()
    for (k in names(updates)) man$config_snapshot[[k]] <- updates[[k]]
    if (!file.exists(paste0(mp, ".bak")))
      file.copy(mp, paste0(mp, ".bak"), overwrite = FALSE)
    jsonlite::write_json(man, mp, auto_unbox = TRUE, pretty = TRUE,
                         null = "null", digits = 8)
    TRUE
  }, error = function(e) FALSE)
}
