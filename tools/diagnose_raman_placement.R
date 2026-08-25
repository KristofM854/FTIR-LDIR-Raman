# =============================================================================
# diagnose_raman_placement.R -- measure, don't guess, the Raman image placement
# =============================================================================
# Scores candidate mappings between Raman particle stage coordinates and the
# background image by sampling image brightness at each particle position.
# Particles are bright blobs on a dark membrane, so the correct mapping
# maximizes the fraction of ALL particles landing on bright pixels (particles
# outside the image count as misses -- this keeps the registration search from
# collapsing onto a tiny image that covers two lucky particles).
#
# Analyzes EVERY Raman image present in <run_dir>/inputs -- the viewer-uploaded
# image (which overrides everything in the app) and the pipeline canonical --
# because they may be different files with different footprints.
#
# Tested candidates per image:
#   - WITec extent, Center Y negated (what the viewer currently does)
#   - WITec extent, Center Y as reported
#   - WITec values as top-left corner
#   - image-relative placement (extent x[0,W], y[0,H])
#   - each of the above with the raster mirrored horizontally / vertically
#   - a free scale + translation registration search (coarse-to-fine),
#     which reveals whether the image's real footprint differs from the
#     WITec Width/Height (e.g. cropped/zoomed export) and by how much
#
# Usage:
#   Rscript tools/diagnose_raman_placement.R output/<run_dir> [--apply]
#
# --apply writes the measured best-fit extent (of the image the viewer
# actually displays) into the run manifest's config_snapshot, replacing the
# stale WITec values for THIS run only, so the viewer places the image
# correctly on the next app start.  A manifest.json.bak backup is kept.
# Refused when the fit is weak or requires mirroring (which the viewer's
# placement config cannot express).
#
# Output: ranked tables on stdout and one comparison PNG per analyzed image at
#   <run_dir>/debug/raman_placement_diagnostic_<which>.png
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
apply_fit <- "--apply" %in% args
args <- setdiff(args, "--apply")
if (length(args) < 1)
  stop("Usage: Rscript tools/diagnose_raman_placement.R <run_dir> [--apply]")
run_dir <- args[1]
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

# Shared measurement core (also used by the pipeline's auto-calibration)
.self <- tryCatch(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                           value = TRUE)[1]), error = function(e) NA)
.core <- file.path(if (!is.na(.self)) dirname(dirname(normalizePath(.self))) else ".",
                   "R", "measure_raman_placement.R")
if (!file.exists(.core)) .core <- "R/measure_raman_placement.R"
source(.core)

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

# --- Load Raman particle stage coordinates ------------------------------------
p_path <- file.path(run_dir, "01_ingested", "raman_ingested.csv")
if (!file.exists(p_path)) stop("Missing ", p_path)
pts <- read.csv(p_path)
x <- pts$x_um[is.finite(pts$x_um) & is.finite(pts$y_um)]
y <- pts$y_um[is.finite(pts$x_um) & is.finite(pts$y_um)]
n_total <- length(x)
cat(sprintf("Particles: %d, X [%.0f, %.0f], Y [%.0f, %.0f]\n",
            n_total, min(x), max(x), min(y), max(y)))

# --- Collect all Raman images in the run --------------------------------------
img_files <- c(uploaded  = file.path(run_dir, "inputs", "raman_image_uploaded.png"),
               canonical = file.path(run_dir, "inputs", "raman_image_canonical.png"),
               preview   = file.path(run_dir, "inputs", "raman_image_preview.png"))
img_files <- img_files[file.exists(img_files)]
if ("canonical" %in% names(img_files))          # preview = downscaled canonical
  img_files <- img_files[names(img_files) != "preview"]
if (length(img_files) == 0)
  stop("No raman image found under ", file.path(run_dir, "inputs"))
if ("uploaded" %in% names(img_files)) {
  cat("\nNOTE: raman_image_uploaded.png exists — the viewer displays THIS file\n",
      "and ignores the pipeline canonical. If the upload was a cropped or\n",
      "zoomed export, no metadata-based placement can match it. Delete the\n",
      "file to make the viewer use the canonical image again.\n", sep = "")
}

dbg <- file.path(run_dir, "debug"); dir.create(dbg, showWarnings = FALSE)

run_max <- function(mat, radius, along_rows) {
  out <- mat
  for (k in setdiff(seq(-radius, radius), 0)) {
    n <- if (along_rows) nrow(mat) else ncol(mat)
    idx <- pmin(pmax(seq_len(n) + k, 1), n)
    out <- if (along_rows) pmax(out, mat[idx, , drop = FALSE])
           else                 pmax(out, mat[, idx, drop = FALSE])
  }
  out
}
dilate <- function(mat, radius) run_max(run_max(mat, radius, TRUE), radius, FALSE)
extent_center <- function(cx, cy, w = W, h = H)
  list(xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2, ymax = cy + h/2)

analyze_image <- function(img_path, which_img) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("=== Analyzing: %s ===\n", which_img))
  img <- png::readPNG(img_path)
  lum <- if (length(dim(img)) == 3)
    0.2126 * img[,,1] + 0.7152 * img[,,2] + 0.0722 * img[,,3] else img
  Hpx <- nrow(lum); Wpx <- ncol(lum)
  cat(sprintf("Image: %s (%d x %d px)\n", basename(img_path), Wpx, Hpx))
  cat(sprintf("Aspect check: image px W/H = %.4f vs WITec um W/H = %.4f %s\n",
              Wpx / Hpx, W / H,
              if (abs(Wpx / Hpx - W / H) > 0.02)
                "<-- MISMATCH: this file's footprint differs from the WITec extent!"
              else "(consistent)"))

  lum_max    <- dilate(lum, 2)                                 # fine (5x5)
  r_coarse   <- max(4L, as.integer(ceiling(min(Hpx, Wpx) / 60)))
  lum_coarse <- dilate(lum, r_coarse)                          # coarse search

  # frac_bright counts over ALL particles; ones outside the extent are misses.
  score_placement <- function(ext, flip_h = FALSE, flip_v = FALSE, mat = lum_max) {
    fx <- (x - ext$xmin) / (ext$xmax - ext$xmin)
    fy <- (ext$ymax - y) / (ext$ymax - ext$ymin)   # row direction (top->bottom)
    if (flip_h) fx <- 1 - fx
    if (flip_v) fy <- 1 - fy
    inside <- fx >= 0 & fx <= 1 & fy >= 0 & fy <= 1
    if (!any(inside)) return(c(frac_bright = 0, mean_lum = 0, n_inside = 0))
    ci <- pmin(pmax(ceiling(fx[inside] * Wpx), 1), Wpx)
    ri <- pmin(pmax(ceiling(fy[inside] * Hpx), 1), Hpx)
    v <- mat[cbind(ri, ci)]
    c(frac_bright = sum(v > 0.5) / n_total,      # denominator: ALL particles
      mean_lum = mean(v), n_inside = sum(inside))
  }

  set.seed(1)
  base_v <- lum_max[cbind(sample.int(Hpx, 3000, TRUE), sample.int(Wpx, 3000, TRUE))]
  cat(sprintf("Random baseline: frac_bright = %.3f  mean_lum = %.3f\n\n",
              mean(base_v > 0.5), mean(base_v)))

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

  # --- Free registration: coarse-to-fine scale + translation search ---
  # Delegates to the shared core (R/measure_raman_placement.R) so the tool
  # and the pipeline's auto-calibration measure identically.  Both gates are
  # disabled here so the tool reports weak fits instead of returning NULL --
  # seeing the weak numbers is the point of running the diagnostic.
  cat("\n=== Free scale + translation search ===\n")
  core <- measure_raman_placement_core(lum, x, y, W, H,
                                       min_frac = -1, min_lift = -1)
  best <- list(frac_bright = core$frac_bright, flip_h = core$flip_h,
               flip_v = core$flip_v, scale = core$scale,
               cx = core$center_x_um, cy = core$center_y_um)
  disp <- score_placement(extent_center(best$cx, best$cy,
                                        W * best$scale, H * best$scale),
                          best$flip_h, best$flip_v)
  best$mean_lum <- disp[["mean_lum"]]; best$n_inside <- disp[["n_inside"]]

  cat(sprintf(paste0(
    "Best registration: frac_bright = %.2f of ALL %d particles ",
    "(%d inside the image; baseline %.2f)\n",
    "  raster mirrored: horizontal=%s vertical=%s\n",
    "  scale = %.2f x WITec extent (=> image spans %.0f x %.0f um)\n",
    "  image center = (%.0f, %.0f)\n",
    "  vs WITec-negY center (%.0f, %.0f): offset dx=%.0f dy=%.0f um\n"),
    best$frac_bright, n_total, best$n_inside, mean(base_v > 0.5),
    best$flip_h, best$flip_v,
    best$scale, W * best$scale, H * best$scale,
    best$cx, best$cy, CX, -CY, best$cx - CX, best$cy - (-CY)))

  # --- Render comparison ---
  out_png <- file.path(dbg, paste0("raman_placement_diagnostic_", which_img, ".png"))
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
             sprintf("[%s] Current viewer placement (WITec, Center Y negated)", which_img))
  draw_panel(extent_center(best$cx, best$cy, W * best$scale, H * best$scale),
             best$flip_h, best$flip_v,
             sprintf("[%s] Best found (frac %.2f, scale %.2f, mirrorH=%s, mirrorV=%s)",
                     which_img, best$frac_bright, best$scale, best$flip_h, best$flip_v))
  dev.off()
  cat("Comparison image written to: ", out_png, "\n")
  invisible(best)
}

results <- list()
for (which_img in names(img_files))
  results[[which_img]] <- analyze_image(img_files[[which_img]], which_img)

cat("\nDone. Open the PNG(s) under ", dbg,
    " — red circles should sit on bright blobs in the right panel.\n", sep = "")

# --- Optionally write the measured extent into the run manifest ---------------
if (apply_fit) {
  viewer_img <- names(img_files)[1]   # same priority order as the app
  best <- results[[viewer_img]]
  cat(sprintf("\n--apply: using the '%s' image's best fit\n", viewer_img))
  if (best$frac_bright < 0.3) {
    stop("Refusing to apply: best fit puts only ",
         round(best$frac_bright * 100), "% of particles on bright pixels — ",
         "too weak to trust.")
  }
  if (isTRUE(best$flip_h) || isTRUE(best$flip_v)) {
    stop("Refusing to apply: best fit requires a mirrored raster, which the ",
         "viewer's placement config cannot express. Re-export the image ",
         "without mirroring, or report this case.")
  }
  file.copy(m_path, paste0(m_path, ".bak"), overwrite = TRUE)
  man$config_snapshot$raman_image_width_um    <- W * best$scale
  man$config_snapshot$raman_image_height_um   <- H * best$scale
  man$config_snapshot$raman_image_center_x_um <- best$cx
  # Stage-frame center: the viewer scores both Y conventions and will pick
  # this value as-reported (negating it would throw the particles outside).
  man$config_snapshot$raman_image_center_y_um <- best$cy
  jsonlite::write_json(man, m_path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", digits = 8)
  cat(sprintf(paste0(
    "Manifest updated (%s; backup at %s):\n",
    "  raman_image_width_um    = %.1f\n",
    "  raman_image_height_um   = %.1f\n",
    "  raman_image_center_x_um = %.1f\n",
    "  raman_image_center_y_um = %.1f\n",
    "Restart the Shiny app (or re-select the run) to see the corrected ",
    "placement.\n"),
    m_path, paste0(m_path, ".bak"),
    W * best$scale, H * best$scale, best$cx, best$cy))
}
