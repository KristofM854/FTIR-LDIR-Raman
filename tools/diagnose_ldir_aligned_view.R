# =============================================================================
# diagnose_ldir_aligned_view.R — why aligned LDIR doesn't sit on the Raman image
# =============================================================================
# The LDIR "aligned (Raman space)" view plots LDIR particles at their aligned
# coordinates over the Raman micrograph (overlay_image_info). If the dots
# don't sit on the image blobs, the cause is one of three things, and this
# tool measures which:
#   1. Transform:   do matched aligned-LDIR coincide with Raman particles?
#                   (small residual = transform is correct)
#   2. Unmatched:   does unmatched_ldir_vs_raman.csv carry x_aligned? If not,
#                   the viewer plots unmatched LDIR at RAW coords in aligned
#                   mode — a gross visual mismatch for most particles.
#   3. Image frame: does the Raman image extent (WITec metadata, shifted into
#                   the normalized/aligned frame by the Raman centroid, as
#                   overlay_image_info does) actually contain the aligned-LDIR
#                   and Raman-particle clouds?
#
# Usage: Rscript tools/diagnose_ldir_aligned_view.R output/<run_dir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript tools/diagnose_ldir_aligned_view.R <run_dir>")
run_dir <- args[1]
mpath <- file.path(run_dir, "05_matches", "matched_ldir_raman.csv")
if (!file.exists(mpath)) stop("Missing ", mpath)
m <- read.csv(mpath)

# --- 1. Transform correctness -----------------------------------------------
res <- sqrt((m$ldir_x_aligned - m$raman_x_norm)^2 +
            (m$ldir_y_aligned - m$raman_y_norm)^2)
cat(sprintf("[1] Transform: %d matched pairs, aligned-LDIR vs Raman-norm residual median %.1f um, 90th %.1f um\n",
            nrow(m), median(res, na.rm = TRUE), quantile(res, .9, na.rm = TRUE)))
cat(sprintf("    %s\n", if (median(res, na.rm = TRUE) < 50)
  "-> transform is CORRECT: matched LDIR coincide with Raman particles."
  else "-> transform residual is large; alignment itself is off."))

# Raman centroid (overlay_image_info shifts the stage image by this)
cenx <- mean(m$raman_x_um - m$raman_x_norm, na.rm = TRUE)
ceny <- mean(m$raman_y_um - m$raman_y_norm, na.rm = TRUE)
cat(sprintf("    Raman centroid (stage->norm shift): (%.0f, %.0f) um\n", cenx, ceny))

# --- 2. Unmatched LDIR carry aligned coords? --------------------------------
upath <- file.path(run_dir, "05_matches", "unmatched_ldir_vs_raman.csv")
if (file.exists(upath)) {
  u <- read.csv(upath)
  has_al <- all(c("x_aligned", "y_aligned") %in% names(u))
  cat(sprintf("\n[2] Unmatched LDIR: %d particles; x_aligned present: %s\n",
              nrow(u), has_al))
  if (!has_al) {
    cat("    -> BUG: viewer falls back to RAW x_um for these in aligned mode,\n")
    cat("       so most unmatched LDIR are drawn at the wrong place.\n")
  } else {
    axr <- range(u$x_aligned, na.rm = TRUE); ayr <- range(u$y_aligned, na.rm = TRUE)
    cat(sprintf("    aligned range x[%.0f,%.0f] y[%.0f,%.0f]; raw x_um range x[%.0f,%.0f]\n",
                axr[1], axr[2], ayr[1], ayr[2],
                min(u$x_um, na.rm=TRUE), max(u$x_um, na.rm=TRUE)))
  }
} else {
  cat("\n[2] No unmatched_ldir_vs_raman.csv (all LDIR matched, or force-complete).\n")
}

# --- 3. Raman image extent in the aligned/normalized frame ------------------
cfg <- tryCatch(jsonlite::fromJSON(
  file.path(run_dir, "00_manifest", "manifest.json"),
  simplifyVector = FALSE)$config_snapshot, error = function(e) NULL)
W <- cfg$raman_image_width_um; H <- cfg$raman_image_height_um
CX <- cfg$raman_image_center_x_um; CY <- cfg$raman_image_center_y_um
cat("\n[3] Raman image placement in the aligned (normalized) frame:\n")
if (is.null(W) || is.null(CX)) {
  cat("    No raman_image_* in manifest — overlay uses the particle-bbox\n")
  cat("    fallback, which will NOT line up with the true image footprint.\n")
} else {
  # Same Y-convention resolution the viewer uses: pick the center-Y sign that
  # puts the Raman particles inside the stage extent.
  px <- c(m$raman_x_um); py <- c(m$raman_y_um)
  score <- function(cy) mean(px >= CX-W/2 & px <= CX+W/2 & py >= cy-H/2 & py <= cy+H/2)
  cy_use <- if (score(-CY) >= score(CY)) -CY else CY
  # stage extent -> normalized (subtract centroid), as overlay_image_info does
  ext_n <- list(xmin = (CX-W/2)-cenx, xmax = (CX+W/2)-cenx,
                ymin = (cy_use-H/2)-ceny, ymax = (cy_use+H/2)-ceny)
  cat(sprintf("    image extent (normalized): x[%.0f,%.0f] y[%.0f,%.0f]\n",
              ext_n$xmin, ext_n$xmax, ext_n$ymin, ext_n$ymax))
  # aligned clouds
  alx <- c(m$ldir_x_aligned, if (exists("u") && all(c("x_aligned") %in% names(u))) u$x_aligned)
  aly <- c(m$ldir_y_aligned, if (exists("u") && all(c("y_aligned") %in% names(u))) u$y_aligned)
  cat(sprintf("    aligned-LDIR range:        x[%.0f,%.0f] y[%.0f,%.0f]\n",
              min(alx,na.rm=TRUE), max(alx,na.rm=TRUE), min(aly,na.rm=TRUE), max(aly,na.rm=TRUE)))
  cat(sprintf("    Raman-norm particle range: x[%.0f,%.0f] y[%.0f,%.0f]\n",
              min(m$raman_x_norm,na.rm=TRUE), max(m$raman_x_norm,na.rm=TRUE),
              min(m$raman_y_norm,na.rm=TRUE), max(m$raman_y_norm,na.rm=TRUE)))
  frac_in <- mean(m$raman_x_norm >= ext_n$xmin & m$raman_x_norm <= ext_n$xmax &
                  m$raman_y_norm >= ext_n$ymin & m$raman_y_norm <= ext_n$ymax, na.rm = TRUE)
  cat(sprintf("    -> %.0f%% of Raman particles fall inside the normalized image extent.\n",
              100*frac_in))
  cat(sprintf("    %s\n", if (frac_in > 0.8)
    "Image frame is consistent: aligned LDIR (= Raman particles) should sit on the image."
    else "MISMATCH: the Raman image is NOT placed over its own particles in the aligned frame -> image-placement bug."))
}

# --- 4. Is the residual reducible? (systematic vs scatter vs smooth warp) ----
# The matched residual can be: a constant offset (fixable by translation),
# a smooth distortion field (fixable by TPS), or random scatter from coarse
# LDIR centroids (irreducible — 19 matches is then the honest ceiling).
cat("\n[4] Residual structure (matched pairs):\n")
dx <- m$raman_x_norm - m$ldir_x_aligned
dy <- m$raman_y_norm - m$ldir_y_aligned
ok <- is.finite(dx) & is.finite(dy)
dx <- dx[ok]; dy <- dy[ok]
sys_mag <- sqrt(mean(dx)^2 + mean(dy)^2)
scatter <- sqrt(mean((dx - mean(dx))^2 + (dy - mean(dy))^2))
cat(sprintf("    systematic offset |mean| = %.0f um; scatter (after removing it) = %.0f um\n",
            sys_mag, scatter))
if (sys_mag > scatter)
  cat("    -> dominated by a CONSTANT offset: a translation tweak recovers most of it.\n")

# Leave-one-out TPS: fit the warp on N-1 controls, predict the held-out one.
# If the held-out residual drops well below the raw residual, the field is
# SMOOTH (TPS will help). If not, it is random scatter (TPS cannot help).
if (exists("tps_fit_warp", mode = "function") ||
    file.exists("R/tps_refine.R")) {
  if (!exists("tps_fit_warp", mode = "function")) source("R/tps_refine.R")
  ax <- m$ldir_x_aligned[ok]; ay <- m$ldir_y_aligned[ok]
  rx <- m$raman_x_norm[ok];   ry <- m$raman_y_norm[ok]
  n <- length(ax)
  if (n > 60) {
    cat(sprintf("    (%d matched pairs — dense run, residual already tiny; skipping LOO-TPS)\n", n))
  } else if (n >= 8) {
    loo <- vapply(seq_len(n), function(i) {
      w <- tryCatch(tps_fit_warp(ax[-i], ay[-i], rx[-i], ry[-i], lambda = 0.5),
                    error = function(e) NULL)
      if (is.null(w)) return(NA_real_)
      p <- tps_apply(w, ax[i], ay[i])
      sqrt((p$x - rx[i])^2 + (p$y - ry[i])^2)
    }, numeric(1))
    raw_i <- sqrt((ax - rx)^2 + (ay - ry)^2)
    cat(sprintf("    leave-one-out TPS: median residual %.0f um -> %.0f um (%.0f%% reduction)\n",
                median(raw_i, na.rm = TRUE), median(loo, na.rm = TRUE),
                100 * (1 - median(loo, na.rm = TRUE) / median(raw_i, na.rm = TRUE))))
    if (median(loo, na.rm = TRUE) < 0.6 * median(raw_i, na.rm = TRUE))
      cat("    -> SMOOTH distortion: TPS should recover stragglers. If the pipeline\n       reported 'no gain', the TPS adopt criterion (count-only) is too strict.\n")
    else
      cat("    -> mostly RANDOM scatter (coarse LDIR centroids): the ~130 um is\n       largely irreducible; 19 matches is close to the honest ceiling.\n")
  }
}
