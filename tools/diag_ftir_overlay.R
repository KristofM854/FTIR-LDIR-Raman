#!/usr/bin/env Rscript
# =============================================================================
# diag_ftir_overlay.R -- measure how well the FTIR image sits under its particles
# =============================================================================
#
# Usage (from the repo root):
#   Rscript tools/diag_ftir_overlay.R --csv "<FTIR particle CSV>" --png "<Average Abs PNG>"
#   Rscript tools/diag_ftir_overlay.R --synthetic reported   # built-in scene
#   Rscript tools/diag_ftir_overlay.R --synthetic hard       # margins + colour bar
#   Optional: --out output/diag   --top 20
#
# What it does
#   1. Loads the CSV with ingest_ftir() + prefilter_ftir() and the PNG with the
#      viewer's load_image_raster() -- the same functions the app uses.
#   2. Reports PNG size, the extent each placement produces, the particle range
#      (all particles and the top-N filtered set) and the implied um/px.
#   3. Detects blobs in the PNG independently of the CSV (global threshold +
#      connected components; R/ftir_image_placement.R::detect_image_blobs()).
#   4. Maps blob centroids into the plotted frame with a given placement.
#   5. Matches blobs to CSV particles by mutual nearest neighbour, largest first.
#   6. Fits CSV -> blob with translation / per-axis scale / full affine models.
#   7. Writes a residual vector plot and a before/after overlay to --out.
#   8. Repeats with the particle set reduced to the top-N largest particles.
#
# Placements compared
#   legacy_tab     the FTIR tab before this fix: raster stretched onto the
#                  bounding box of ALL run particles (app.R ftir_native_image_info)
#   legacy_report  the pipeline report before this fix: aspect-preserving fit
#                  to the bounding box of the FILTERED particles
#                  (global.R place_image_particle_extent via report_instrument_image)
#   fixed          place_ftir_image() -- the shared function every plot now uses
#
# On a synthetic scene the true extent is known, so the script also reports the
# placement error against ground truth (independent of the blob detector), and
# a split-half check: register on half the particles, score the other half.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default = NULL) {
  i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1]
}
csv_path  <- .opt("csv")
png_path  <- .opt("png")
synthetic <- .opt("synthetic")
out_dir   <- .opt("out", file.path("output", "diag"))
top_n     <- as.integer(.opt("top", "20"))
if (is.null(synthetic) && (is.null(csv_path) || is.null(png_path)))
  synthetic <- "reported"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- App code ----------------------------------------------------------------
for (f in c("R/utils.R", "R/00_config.R", "R/01_ingest.R", "R/02_prefilter.R",
            "R/ftir_image_placement.R", "tools/ftir_synthetic.R"))
  sys.source(f, envir = globalenv())
.viewer <- local({
  env <- new.env(parent = globalenv())
  txt <- paste(readLines("shiny_app/global.R", warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  env
})
log_message <- function(...) invisible(NULL)   # keep the report readable

truth <- NULL
if (!is.null(synthetic)) {
  spec  <- synthetic_ftir_scenes()[[synthetic]]
  if (is.null(spec)) stop("unknown --synthetic scene: ", synthetic)
  scene <- do.call(make_synthetic_ftir, spec)
  files <- write_synthetic_ftir(scene, file.path(out_dir, paste0("synthetic_", synthetic)),
                                stem = paste0("synthetic ", synthetic))
  csv_path <- files$csv; png_path <- files$png; truth <- scene$truth
}
label <- if (!is.null(synthetic)) paste0("synthetic:", synthetic) else basename(csv_path)

ftir <- prefilter_ftir(ingest_ftir(csv_path), min_quality = 0, min_size_um = 0)
raw  <- .viewer$load_image_raster(png_path)
if (is.null(raw)) stop("could not read ", png_path)
W <- ncol(raw); H <- nrow(raw)
x <- ftir$x_um; y <- ftir$y_um; sz <- ftir$feret_max_um
top <- order(-sz)[seq_len(min(top_n, nrow(ftir)))]

cat(sprintf("\n=== FTIR overlay diagnostic: %s ===\n", label))
cat(sprintf("PNG: %d x %d px | particles: %d (top-%d set: %d)\n", W, H, nrow(ftir), top_n, length(top)))
rng <- function(v) sprintf("[%.1f, %.1f] (span %.1f)", min(v), max(v), diff(range(v)))
cat("Particle X range, all:   ", rng(x), "\n")
cat("Particle Y range, all:   ", rng(y), "\n")
cat("Particle X range, top-N: ", rng(x[top]), "\n")
cat("Particle Y range, top-N: ", rng(y[top]), "\n")

# --- Placements ---------------------------------------------------------------
legacy_tab <- function(xx, yy) list(raster = raw, xmin = min(xx), xmax = max(xx),
                                    ymin = min(yy), ymax = max(yy), method = "legacy_tab")
legacy_report <- function(xx, yy) c(list(raster = raw, method = "legacy_report"),
  .viewer$compute_image_bounds(raw, xx, yy, padding_um = 0))
fixed <- function(xx, yy, ss) place_ftir_image(raw, xx, yy, ss)

# --- Blob detection (independent of the CSV) -----------------------------------
blobs <- detect_image_blobs(raw)
cat(sprintf("Blobs detected: %d (polarity %s, threshold %.3f)\n", NROW(blobs),
            attr(blobs, "polarity"), attr(blobs, "threshold")))

# Blob centroids (edge px) -> plotted frame under a placement. The placement's
# raster may be row-flipped; centroids were measured on the unflipped image.
blob_xy <- function(pl) {
  fx <- (pl$xmax - pl$xmin) / W; fy <- (pl$ymax - pl$ymin) / H
  flipped <- isTRUE(pl$registration$by > 0)
  cbind(pl$xmin + blobs$cx_px * fx,
        if (flipped) pl$ymin + blobs$cy_px * fy else pl$ymax - blobs$cy_px * fy)
}

# Mutual-NN matching, largest first: stage 1 pairs the larger half of each
# set and estimates a rough shift; stage 2 pairs everything after removing it.
match_blobs <- function(B, idx) {
  P <- cbind(x[idx], y[idx])
  span <- max(diff(range(P[, 1])), diff(range(P[, 2])))
  bigb <- which(blobs$area_px >= stats::median(blobs$area_px))
  bigp <- which(sz[idx] >= stats::median(sz[idx]))
  s1 <- mutual_nn_pairs(B[bigb, , drop = FALSE], P[bigp, , drop = FALSE], 0.1 * span)
  shift <- if (nrow(s1) >= 3)
    apply(B[bigb[s1$ia], , drop = FALSE] - P[bigp[s1$ib], , drop = FALSE], 2, stats::median)
  else c(0, 0)
  Bs <- sweep(B, 2, shift)
  s2 <- mutual_nn_pairs(Bs, P, 0.05 * span)
  data.frame(blob = s2$ia, particle = idx[s2$ib])
}

fit_models <- function(B, pr) {
  cx <- x[pr$particle]; cy <- y[pr$particle]
  bx <- B[pr$blob, 1];  by <- B[pr$blob, 2]
  rms <- function(rx, ry) sqrt(mean(rx^2 + ry^2))
  dx <- mean(bx - cx); dy <- mean(by - cy)
  t_rms <- rms(bx - cx - dx, by - cy - dy)
  fx <- stats::lm(bx ~ cx); fy <- stats::lm(by ~ cy)
  s_rms <- rms(stats::residuals(fx), stats::residuals(fy))
  ax <- stats::lm(bx ~ cx + cy); ay <- stats::lm(by ~ cx + cy)
  a_rms <- rms(stats::residuals(ax), stats::residuals(ay))
  sx <- unname(coef(fx)[2]); sy <- unname(coef(fy)[2])
  list(n = nrow(pr), dx = dx, dy = dy, rms_translation = t_rms,
       sx = sx, sy = sy, dx_s = unname(coef(fx)[1]), dy_s = unname(coef(fy)[1]),
       fixed_x = unname(coef(fx)[1]) / (1 - sx), fixed_y = unname(coef(fy)[1]) / (1 - sy),
       rms_scale = s_rms, affine = rbind(coef(ax), coef(ay)), rms_affine = a_rms,
       raw_rms = rms(bx - cx, by - cy))
}

truth_error <- function(pl) {
  if (is.null(truth)) return(NULL)
  # Where the placement puts true particle positions, vs where they are.
  # Each particle sits at pixel (u, v) of the image; the placement draws that
  # pixel at xmin + u * (xmax - xmin) / W.
  u <- (x - truth$xmin) / (truth$xmax - truth$xmin) * W
  v <- (truth$ymax - y) / (truth$ymax - truth$ymin) * H
  flipped <- isTRUE(pl$registration$by > 0)
  px <- pl$xmin + u * (pl$xmax - pl$xmin) / W
  py <- if (flipped) pl$ymin + v * (pl$ymax - pl$ymin) / H else pl$ymax - v * (pl$ymax - pl$ymin) / H
  e <- sqrt((px - x)^2 + (py - y)^2)
  list(edge_err_um = c(xmin = pl$xmin - truth$xmin, xmax = pl$xmax - truth$xmax,
                       ymin = pl$ymin - truth$ymin, ymax = pl$ymax - truth$ymax),
       rms_um = sqrt(mean(e^2)), max_um = max(e),
       rms_px = sqrt(mean(e^2)) / sqrt(truth$px_x * truth$px_y))
}

px_um <- NA_real_
results <- list()
# Blob <-> particle correspondence is a property of the image and the CSV, not
# of the placement under test, so it is established ONCE (by the scale-aware
# registration, which matches the larger particles first) and every placement
# is scored on the same pairs. A placement that is off by more than the
# particle spacing would otherwise be scored on wrong pairings. Without a
# registration (it failed on this dataset) each placement is matched on its own.
corr <- NULL
report <- function(tag, pl, idx) {
  B <- blob_xy(pl)
  pr <- if (!is.null(corr)) corr[corr$particle %in% idx, c("blob", "particle")]
        else match_blobs(B, idx)
  fm <- fit_models(B, pr)
  te <- truth_error(pl)
  upx <- (pl$xmax - pl$xmin) / W; upy <- (pl$ymax - pl$ymin) / H
  cat(sprintf("\n--- %s ---\n", tag))
  cat(sprintf("  extent: X [%.1f, %.1f]  Y [%.1f, %.1f]  -> %.4f x %.4f um/px  (method %s)\n",
              pl$xmin, pl$xmax, pl$ymin, pl$ymax, upx, upy, pl$method))
  cat(sprintf("  particle span / image span: X %.4f  Y %.4f\n",
              diff(range(x[idx])) / (pl$xmax - pl$xmin), diff(range(y[idx])) / (pl$ymax - pl$ymin)))
  pxs <- sqrt(upx * upy)
  cat(sprintf("  matched: %d of %d particles / %d blobs | raw RMS %.1f um (%.2f px)\n",
              fm$n, length(idx), nrow(blobs), fm$raw_rms, fm$raw_rms / pxs))
  cat(sprintf("  translation : dx %+.1f um  dy %+.1f um  | RMS %.1f um (%.2f px)\n",
              fm$dx, fm$dy, fm$rms_translation, fm$rms_translation / pxs))
  cat(sprintf("  axis scale  : sx %.5f (%+.3f%%)  sy %.5f (%+.3f%%)  dx %+.1f dy %+.1f | RMS %.1f um (%.2f px)\n",
              fm$sx, 100 * (fm$sx - 1), fm$sy, 100 * (fm$sy - 1), fm$dx_s, fm$dy_s,
              fm$rms_scale, fm$rms_scale / pxs))
  if (abs(fm$sx - 1) > 1e-3 || abs(fm$sy - 1) > 1e-3)
    cat(sprintf("                zero-error point: X %.0f um, Y %.0f um\n", fm$fixed_x, fm$fixed_y))
  cat(sprintf("  affine      : RMS %.1f um (%.2f px)\n", fm$rms_affine, fm$rms_affine / pxs))
  if (!is.null(te))
    cat(sprintf("  vs TRUTH    : edge error xmin %+.1f xmax %+.1f ymin %+.1f ymax %+.1f um | particle RMS %.2f um (%.3f px), max %.2f um\n",
                te$edge_err_um[1], te$edge_err_um[2], te$edge_err_um[3], te$edge_err_um[4],
                te$rms_um, te$rms_px, te$max_um))
  results[[tag]] <<- c(list(tag = tag, method = pl$method, xmin = pl$xmin, xmax = pl$xmax,
                            ymin = pl$ymin, ymax = pl$ymax), fm[setdiff(names(fm), "affine")],
                       if (!is.null(te)) list(truth_rms_um = te$rms_um, truth_rms_px = te$rms_px))
  invisible(list(pl = pl, B = B, pr = pr, fm = fm))
}

all_idx <- seq_along(x)
pl_tab  <- legacy_tab(x, y)
pl_fix  <- fixed(x, y, sz)
if (identical(pl_fix$method, "registration")) {
  corr <- pl_fix$registration$pairs
  cat(sprintf("Correspondence: %d blob/particle pairs from the registration\n", nrow(corr)))
} else {
  cat("Registration failed (", pl_fix$registration$reason %||% "no usable blobs",
      "); matching each placement separately\n")
}
r_tab   <- report("legacy_tab   | placement from ALL particles", pl_tab, all_idx)
r_rep   <- report("legacy_report| placement from ALL particles", legacy_report(x, y), all_idx)
r_fix   <- report("fixed        | placement from ALL particles", pl_fix, all_idx)

# Step 8: placement recomputed from the filtered (top-N) set.
report(sprintf("legacy_tab   | placement from top-%d", top_n), legacy_tab(x[top], y[top]), all_idx)
report(sprintf("legacy_report| placement from top-%d", top_n), legacy_report(x[top], y[top]), all_idx)
pl_fix_top <- fixed(x[top], y[top], sz[top])
report(sprintf("fixed        | placement from top-%d", top_n), pl_fix_top, all_idx)
cat(sprintf("\nExtent shift when the particle set changes (all -> top-%d), um:\n", top_n))
shift <- function(a, b) sprintf("xmin %+.1f xmax %+.1f ymin %+.1f ymax %+.1f",
                                b$xmin - a$xmin, b$xmax - a$xmax, b$ymin - a$ymin, b$ymax - a$ymax)
cat("  legacy_tab   :", shift(pl_tab, legacy_tab(x[top], y[top])), "\n")
cat("  legacy_report:", shift(legacy_report(x, y), legacy_report(x[top], y[top])), "\n")
cat("  fixed        :", shift(pl_fix, pl_fix_top), "\n")
cat("  (the viewer and report now compute the placement once from ALL particles,\n",
    "  so a display filter cannot change it; the top-N row shows robustness to a sparse sample)\n")

# Split-half: register on even-indexed particles, score the odd ones.
even <- all_idx[all_idx %% 2 == 0]; odd <- all_idx[all_idx %% 2 == 1]
pl_half <- fixed(x[even], y[even], sz[even])
report("fixed        | registered on even half, scored on odd half", pl_half, odd)

# Particle-reference check (hypothesis F): residual vs particle size.
if (nrow(r_fix$pr) >= 5) {
  d <- r_fix$B[r_fix$pr$blob, , drop = FALSE] - cbind(x[r_fix$pr$particle], y[r_fix$pr$particle])
  s <- sz[r_fix$pr$particle]
  cat(sprintf("\nResidual vs particle size (fixed placement): cor(dx, size) %.3f, cor(dy, size) %.3f\n",
              suppressWarnings(stats::cor(d[, 1], s)), suppressWarnings(stats::cor(d[, 2], s))))
}

# --- Plots ---------------------------------------------------------------------
circles <- function(cx, cy, r, n = 36) {
  t <- seq(0, 2 * pi, length.out = n + 1)[-1]
  data.frame(id = rep(seq_along(cx), each = n),
             x = rep(cx, each = n) + rep(r, each = n) * cos(t),
             y = rep(cy, each = n) + rep(r, each = n) * sin(t))
}
overlay <- function(pl, title) {
  ci <- circles(x, y, sz / 2)
  ggplot() +
    annotation_raster(pl$raster, pl$xmin, pl$xmax, pl$ymin, pl$ymax) +
    geom_polygon(data = ci, aes(x, y, group = id), fill = NA, colour = "cyan", linewidth = 0.3) +
    geom_point(aes(x = x, y = y), colour = "cyan", size = 0.4) +
    coord_fixed(xlim = range(c(pl$xmin, pl$xmax, x)), ylim = range(c(pl$ymin, pl$ymax, y))) +
    labs(title = title, subtitle = sprintf("circles: Feret Max diameter in um | extent X [%.0f, %.0f] Y [%.0f, %.0f]",
                                           pl$xmin, pl$xmax, pl$ymin, pl$ymax),
         x = "X (um)", y = "Y (um)") +
    theme_minimal(base_size = 9)
}
stem <- gsub("[^A-Za-z0-9]+", "_", label)
ggsave(file.path(out_dir, paste0("overlay_before_", stem, ".png")),
       overlay(pl_tab, paste("BEFORE (legacy FTIR tab):", label)), width = 7, height = 7, dpi = 150, bg = "white")
ggsave(file.path(out_dir, paste0("overlay_after_", stem, ".png")),
       overlay(pl_fix, paste("AFTER (place_ftir_image):", label)), width = 7, height = 7, dpi = 150, bg = "white")

vec <- function(r, title, mag) {
  d <- data.frame(x = x[r$pr$particle], y = y[r$pr$particle],
                  ex = r$B[r$pr$blob, 1], ey = r$B[r$pr$blob, 2])
  ggplot(d) +
    geom_segment(aes(x = x, y = y, xend = x + mag * (ex - x), yend = y + mag * (ey - y)),
                 arrow = grid::arrow(length = grid::unit(0.08, "cm")), linewidth = 0.3) +
    geom_point(aes(x, y), size = 0.5) + coord_fixed() +
    labs(title = title, subtitle = sprintf("arrow: CSV particle -> image blob, magnified %gx", mag),
         x = "X (um)", y = "Y (um)") + theme_minimal(base_size = 9)
}
ggsave(file.path(out_dir, paste0("residuals_before_", stem, ".png")),
       vec(r_tab, paste("Residuals BEFORE:", label), 5), width = 7, height = 7, dpi = 150, bg = "white")
ggsave(file.path(out_dir, paste0("residuals_after_", stem, ".png")),
       vec(r_fix, paste("Residuals AFTER:", label), 100), width = 7, height = 7, dpi = 150, bg = "white")

res_df <- do.call(rbind, lapply(results, function(r) as.data.frame(r[!vapply(r, is.null, TRUE)])))
utils::write.csv(res_df, file.path(out_dir, paste0("diag_", stem, ".csv")), row.names = FALSE)
cat("\nWrote plots and diag_", stem, ".csv to ", out_dir, "\n", sep = "")
