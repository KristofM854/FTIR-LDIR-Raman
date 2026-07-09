# =============================================================================
# diagnose_raman_placement.R — measure, don't guess, the Raman image placement
# =============================================================================
# Scores candidate mappings between Raman particle stage coordinates and the
# background image by sampling image brightness at each particle position
# (particles are bright blobs on a dark membrane, so the correct mapping
# maximizes the fraction of particles landing on bright pixels).
#
# Tested candidates:
#   - WITec extent, Center Y negated (what the viewer currently does)
#   - WITec extent, Center Y as reported
#   - image-relative placement (particle coords measured from the image
#     corner instead of the stage origin), Y up and Y down
#   - each of the above with the raster mirrored horizontally / vertically
#   - a free scale + translation registration search, which reveals whether
#     the exported image's footprint differs from the WITec Width/Height
#     (e.g. export region != scan region) and by how much
#
# Usage:
#   Rscript tools/diagnose_raman_placement.R output/<run_dir>
#
# Output: a ranked table on stdout and a rendered comparison image at
#   <run_dir>/debug/raman_placement_diagnostic.png
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript tools/diagnose_raman_placement.R <run_dir>")
run_dir <- args[1]
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

# --- Load manifest config (WITec values) -------------------------------------
m_path <- file.path(run_dir, "00_manifest", "manifest.json")
if (!file.exists(m_path)) m_path <- file.path(run_dir, "manifest.json")
cfg <- list()
if (file.exists(m_path)) {
  man <- jsonlite::fromJSON(m_path, simplifyVector = FALSE)
  cfg <- man$config_snapshot
}
W  <- cfg$raman_image_width_um
H  <- cfg$raman_image_height_um
CX <- cfg$raman_image_center_x_um
CY <- cfg$raman_image_center_y_um
if (is.null(W) || is.null(H) || is.null(CX) || is.null(CY))
  stop("Manifest has no raman_image_* values in config_snapshot — ",
       "re-run the pipeline with them set in 00_config.R, or paste them ",
       "into ", m_path)
cat(sprintf("WITec panel values: W=%.1f H=%.1f Center=(%.1f, %.1f)\n", W, H, CX, CY))

# --- Locate the image the viewer displays (same priority as the app) ---------
img_path <- NULL
for (cand in c(file.path(run_dir, "inputs", "raman_image_uploaded.png"),
               file.path(run_dir, "inputs", "raman_image_canonical.png"),
               file.path(run_dir, "inputs", "raman_image_preview.png"))) {
  if (file.exists(cand)) { img_path <- cand; break }
}
if (is.null(img_path)) stop("No raman image found under ", file.path(run_dir, "inputs"))
img <- png::readPNG(img_path)
lum <- if (length(dim(img)) == 3)
  0.2126 * img[,,1] + 0.7152 * img[,,2] + 0.0722 * img[,,3] else img
Hpx <- nrow(lum); Wpx <- ncol(lum)
cat(sprintf("Image: %s (%d x %d px)\n", basename(img_path), Wpx, Hpx))
cat(sprintf("Aspect check: image px W/H = %.4f vs WITec um W/H = %.4f %s\n",
            Wpx / Hpx, W / H,
            if (abs(Wpx / Hpx - W / H) > 0.02)
              "<-- MISMATCH: exported image footprint differs from WITec extent!"
            else "(consistent)"))

# --- Load Raman particle stage coordinates ------------------------------------
p_path <- file.path(run_dir, "01_ingested", "raman_ingested.csv")
if (!file.exists(p_path)) stop("Missing ", p_path)
pts <- read.csv(p_path)
x <- pts$x_um[is.finite(pts$x_um) & is.finite(pts$y_um)]
y <- pts$y_um[is.finite(pts$x_um) & is.finite(pts$y_um)]
cat(sprintf("Particles: %d, X [%.0f, %.0f], Y [%.0f, %.0f]\n\n",
            length(x), min(x), max(x), min(y), max(y)))

# --- Scoring: 5x5-max luminance at each particle position ---------------------
# Precompute a 5x5 max-dilated luminance matrix once (separable running max)
# so each placement score is pure vectorized indexing — the registration
# search below evaluates ~100k placements.
run_max <- function(mat, radius, along_rows) {
  out <- mat
  ks <- setdiff(seq(-radius, radius), 0)
  for (k in ks) {
    n <- if (along_rows) nrow(mat) else ncol(mat)
    idx <- pmin(pmax(seq_len(n) + k, 1), n)
    out <- if (along_rows) pmax(out, mat[idx, , drop = FALSE])
           else                 pmax(out, mat[, idx, drop = FALSE])
  }
  out
}
dilate <- function(mat, radius) run_max(run_max(mat, radius, TRUE), radius, FALSE)
lum_max <- dilate(lum, 2)                                   # fine score (5x5)
r_coarse   <- max(4L, as.integer(ceiling(min(Hpx, Wpx) / 60)))
lum_coarse <- dilate(lum, r_coarse)                         # coarse-search score

# extent: list(xmin,xmax,ymin,ymax); flip_h/flip_v mirror the raster in place.
score_placement <- function(ext, flip_h = FALSE, flip_v = FALSE, mat = lum_max) {
  fx <- (x - ext$xmin) / (ext$xmax - ext$xmin)
  fy <- (ext$ymax - y) / (ext$ymax - ext$ymin)   # row direction (top -> bottom)
  if (flip_h) fx <- 1 - fx
  if (flip_v) fy <- 1 - fy
  inside <- fx >= 0 & fx <= 1 & fy >= 0 & fy <= 1
  if (!any(inside)) return(c(frac_bright = 0, mean_lum = 0, n_inside = 0))
  ci <- pmin(pmax(ceiling(fx[inside] * Wpx), 1), Wpx)
  ri <- pmin(pmax(ceiling(fy[inside] * Hpx), 1), Hpx)
  v <- mat[cbind(ri, ci)]
  c(frac_bright = mean(v > 0.5), mean_lum = mean(v), n_inside = sum(inside))
}

set.seed(1)
base_v <- lum_max[cbind(sample.int(Hpx, 3000, TRUE), sample.int(Wpx, 3000, TRUE))]
cat(sprintf("Random baseline: frac_bright = %.3f  mean_lum = %.3f\n\n",
            mean(base_v > 0.5), mean(base_v)))

# --- Named candidates ----------------------------------------------------------
extent_center <- function(cx, cy, w = W, h = H)
  list(xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2, ymax = cy + h/2)
cands <- list(
  witec_negY     = extent_center(CX, -CY),
  witec_rawY     = extent_center(CX,  CY),
  corner_topleft = list(xmin = CX, xmax = CX + W, ymin = CY - H, ymax = CY),
  image_relative = list(xmin = 0, xmax = W, ymin = 0, ymax = H)
)
rows <- list()
for (nm in names(cands)) for (fh in c(FALSE, TRUE)) for (fv in c(FALSE, TRUE)) {
  s <- score_placement(cands[[nm]], fh, fv)
  rows[[length(rows) + 1]] <- data.frame(
    candidate = nm, flip_h = fh, flip_v = fv, scale = 1,
    frac_bright = s["frac_bright"], mean_lum = s["mean_lum"],
    n_inside = s["n_inside"])
}
tab <- do.call(rbind, rows)
tab <- tab[order(-tab$frac_bright, -tab$mean_lum), ]
rownames(tab) <- NULL
cat("=== Named candidates (scale = 1) ===\n")
print(tab, digits = 3)

# --- Free registration: scale + translation search -----------------------------
# Coarse-to-fine: the coarse pass scores against a heavily dilated luminance
# (tolerance ~ r_coarse px) with a grid step matched to that tolerance, so
# sharp optima cannot fall between grid points; the top coarse hits are then
# refined locally against the 5x5-dilated luminance.
cat("\n=== Free scale + translation search ===\n")
pcx <- mean(range(x)); pcy <- mean(range(y))

better <- function(a, b)   # is score-vector a better than list b?
  a["frac_bright"] > b$frac_bright ||
  (a["frac_bright"] == b$frac_bright && a["mean_lum"] > b$mean_lum)

coarse_hits <- list()
for (fh in c(FALSE, TRUE)) for (fv in c(FALSE, TRUE)) {
  for (s in seq(0.3, 2.4, by = 0.1)) {
    w_s <- W * s; h_s <- H * s
    step <- r_coarse * (w_s / Wpx)            # µm per coarse-tolerance unit
    span <- max(w_s, h_s)
    hit <- list(frac_bright = -1, mean_lum = -1)
    for (dx in seq(-span/2, span/2, by = step))
      for (dy in seq(-span/2, span/2, by = step)) {
        sc <- score_placement(extent_center(pcx + dx, pcy + dy, w_s, h_s),
                              fh, fv, mat = lum_coarse)
        if (better(sc, hit))
          hit <- list(frac_bright = sc[["frac_bright"]], mean_lum = sc[["mean_lum"]],
                      flip_h = fh, flip_v = fv, scale = s,
                      cx = pcx + dx, cy = pcy + dy)
      }
    coarse_hits[[length(coarse_hits) + 1]] <- hit
  }
}
ord <- order(-vapply(coarse_hits, `[[`, 0, "frac_bright"),
             -vapply(coarse_hits, `[[`, 0, "mean_lum"))
top <- coarse_hits[ord[seq_len(min(6, length(ord)))]]

best <- list(frac_bright = -1, mean_lum = -1)
for (h in top) {
  for (s in seq(h$scale - 0.06, h$scale + 0.06, by = 0.02)) {
    w_s <- W * s; h_s <- H * s
    fine_step <- max(1, r_coarse / 4) * (w_s / Wpx)
    for (dx in seq(-3 * r_coarse, 3 * r_coarse, by = max(1, r_coarse / 4)) * (w_s / Wpx))
      for (dy in seq(-3 * r_coarse, 3 * r_coarse, by = max(1, r_coarse / 4)) * (h_s / Hpx)) {
        sc <- score_placement(extent_center(h$cx + dx, h$cy + dy, w_s, h_s),
                              h$flip_h, h$flip_v)
        if (better(sc, best))
          best <- list(frac_bright = sc[["frac_bright"]], mean_lum = sc[["mean_lum"]],
                       n_inside = sc[["n_inside"]], flip_h = h$flip_h,
                       flip_v = h$flip_v, scale = s,
                       cx = h$cx + dx, cy = h$cy + dy)
      }
  }
}
# Final sub-pixel polish around the winner
for (dx in seq(-2, 2, by = 0.5) * (W * best$scale / Wpx))
  for (dy in seq(-2, 2, by = 0.5) * (H * best$scale / Hpx)) {
    sc <- score_placement(extent_center(best$cx + dx, best$cy + dy,
                                        W * best$scale, H * best$scale),
                          best$flip_h, best$flip_v)
    if (better(sc, best)) {
      best$frac_bright <- sc[["frac_bright"]]; best$mean_lum <- sc[["mean_lum"]]
      best$n_inside <- sc[["n_inside"]]
      best$cx <- best$cx + dx; best$cy <- best$cy + dy
    }
  }

cat(sprintf(paste0(
  "Best registration: frac_bright = %.2f (baseline %.2f)\n",
  "  raster mirrored: horizontal=%s vertical=%s\n",
  "  scale = %.2f x WITec extent (=> image spans %.0f x %.0f um)\n",
  "  image center = (%.0f, %.0f)\n",
  "  vs WITec-negY center (%.0f, %.0f): offset dx=%.0f dy=%.0f um\n"),
  best$frac_bright, mean(base_v > 0.5),
  best$flip_h, best$flip_v,
  best$scale, W * best$scale, H * best$scale,
  best$cx, best$cy, CX, -CY, best$cx - CX, best$cy - (-CY)))

# --- Render comparison image ----------------------------------------------------
dbg <- file.path(run_dir, "debug"); dir.create(dbg, showWarnings = FALSE)
out_png <- file.path(dbg, "raman_placement_diagnostic.png")
draw_panel <- function(ext, fh, fv, title) {
  ras <- img
  if (fv) ras <- if (length(dim(ras)) == 3) ras[rev(seq_len(Hpx)), , , drop = FALSE]
                 else ras[rev(seq_len(Hpx)), , drop = FALSE]
  if (fh) ras <- if (length(dim(ras)) == 3) ras[, rev(seq_len(Wpx)), , drop = FALSE]
                 else ras[, rev(seq_len(Wpx)), drop = FALSE]
  xr <- range(c(ext$xmin, ext$xmax, x)); yr <- range(c(ext$ymin, ext$ymax, y))
  plot(NA, xlim = xr, ylim = yr, asp = 1, xlab = "X (um)", ylab = "Y (um)",
       main = title, cex.main = 0.9)
  rasterImage(ras, ext$xmin, ext$ymin, ext$xmax, ext$ymax)
  points(x, y, col = "#FF3333", pch = 1, lwd = 1.5, cex = 1.1)
}
png(out_png, width = 2000, height = 1050, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
draw_panel(cands$witec_negY, FALSE, FALSE,
           "Current viewer placement (WITec, Center Y negated)")
draw_panel(extent_center(best$cx, best$cy, W * best$scale, H * best$scale),
           best$flip_h, best$flip_v,
           sprintf("Best found (scale %.2f, mirrorH=%s, mirrorV=%s)",
                   best$scale, best$flip_h, best$flip_v))
dev.off()
cat("\nComparison image written to: ", out_png, "\n")
cat("Open it: red circles should sit on bright blobs in the right panel.\n")
