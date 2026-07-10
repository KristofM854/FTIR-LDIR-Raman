# =============================================================================
# diagnose_ldir_rotation.R — measure the true LDIR -> Raman orientation
# =============================================================================
# Registers the LDIR particle point cloud (raw circle-calibrated µm from
# 02_joined/ldir_joined_raw.csv) onto the Raman point cloud
# (01_ingested/raman_ingested.csv) over a full grid of rotation (1° steps),
# mirror, and scale, using translation-voting: for every candidate
# orientation the best translation is found by binning all pairwise
# LDIR->Raman offset vectors and taking the densest bin, then counting exact
# inliers.  No dependence on the pipeline's RANSAC — this is an independent
# measurement of the true relative orientation of the two coordinate frames.
#
# It then reads the transform the pipeline actually chose
# (04_alignment/transform_params_ldir_raman.txt) and says whether the
# pipeline landed on the measured consensus or on something spurious.
#
# NOTE: the pipeline applies config$ldir_rotate_deg_for_alignment (-90 by
# default) BEFORE its RANSAC, and the RANSAC itself searches all rotations.
# The rotation reported here is the TOTAL rotation from raw LDIR coordinates
# to the Raman frame — i.e. pre-rotation + RANSAC residual combined.
#
# Usage:
#   Rscript tools/diagnose_ldir_rotation.R output/<run_dir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript tools/diagnose_ldir_rotation.R <run_dir>")
run_dir <- args[1]
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

l_path <- file.path(run_dir, "02_joined", "ldir_joined_raw.csv")
r_path <- file.path(run_dir, "01_ingested", "raman_ingested.csv")
if (!file.exists(l_path)) stop("Missing ", l_path)
if (!file.exists(r_path)) stop("Missing ", r_path)

ldir  <- read.csv(l_path)
raman <- read.csv(r_path)
lok <- is.finite(ldir$x_um)  & is.finite(ldir$y_um)
rok <- is.finite(raman$x_um) & is.finite(raman$y_um)
lx <- ldir$x_um[lok];  ly <- ldir$y_um[lok]
rx <- raman$x_um[rok]; ry <- raman$y_um[rok]
cat(sprintf("LDIR:  %d particles, X [%.0f, %.0f], Y [%.0f, %.0f]\n",
            length(lx), min(lx), max(lx), min(ly), max(ly)))
cat(sprintf("Raman: %d particles, X [%.0f, %.0f], Y [%.0f, %.0f]\n\n",
            length(rx), min(rx), max(rx), min(ry), max(ry)))
if (length(lx) < 4 || length(rx) < 4) stop("Too few particles to register.")

# Center both clouds (translation is recovered by voting, so only the
# relative offset matters)
lx <- lx - mean(lx); ly <- ly - mean(ly)
rx <- rx - mean(rx); ry <- ry - mean(ry)

INLIER_UM <- 200   # same tolerance class as ransac_inlier_dist_um

# Subsample large clouds for the orientation search (the spatial pattern
# survives subsampling); the winners are re-scored on the full sets below.
set.seed(42)
subsample_idx <- function(v, n_max = 150)
  if (length(v) > n_max) sort(sample(seq_along(v), n_max)) else seq_along(v)
li <- subsample_idx(lx); ri <- subsample_idx(rx)
slx <- lx[li]; sly <- ly[li]; srx <- rx[ri]; sry <- ry[ri]

rotate_cloud <- function(px, py, angle_deg, mirror, s) {
  th <- angle_deg * pi / 180
  ct <- cos(th); st <- sin(th)
  if (!mirror) list(x = s * (ct * px - st * py), y = s * (st * px + ct * py))
  else         list(x = s * (ct * px + st * py), y = s * (st * px - ct * py))
}

# Inlier count of transformed LDIR (X,Y) against a Raman cloud (ax,ay)
count_inliers <- function(X, Y, ax, ay) {
  d2 <- outer(ax, X, "-")^2 + outer(ay, Y, "-")^2   # [n_raman, n_ldir]
  dmin <- d2[cbind(max.col(-t(d2)), seq_along(X))]  # per-LDIR min distance^2
  sum(sqrt(dmin) <= INLIER_UM)
}

score_orientation <- function(angle_deg, mirror, s) {
  p <- rotate_cloud(slx, sly, angle_deg, mirror, s)
  # translation voting: bin all pairwise raman - ldir offset vectors
  dx <- outer(srx, p$x, "-"); dy <- outer(sry, p$y, "-")
  key <- (round(dx / INLIER_UM) + 4096) * 8192 + (round(dy / INLIER_UM) + 4096)
  uk  <- unique(as.vector(key))
  cnt <- tabulate(match(as.vector(key), uk))
  ord <- order(-cnt)[seq_len(min(5, length(uk)))]
  best_n <- 0; best_tx <- 0; best_ty <- 0
  for (k in uk[ord]) {
    sel <- key == k
    tx <- mean(dx[sel]); ty <- mean(dy[sel])
    n <- count_inliers(p$x + tx, p$y + ty, srx, sry)
    if (n > best_n) { best_n <- n; best_tx <- tx; best_ty <- ty }
  }
  list(n = best_n, tx = best_tx, ty = best_ty)
}

cat("Searching rotation (1° grid) x mirror x scale ... (~30-90 s)\n")
results <- list()
for (mirror in c(FALSE, TRUE)) {
  for (s in seq(0.85, 1.15, by = 0.05)) {
    for (a in seq(0, 359, by = 1)) {
      sc <- score_orientation(a, mirror, s)
      results[[length(results) + 1]] <- data.frame(
        rotation_deg = a, mirror = mirror, scale = s,
        n_inliers = sc$n, tx = sc$tx, ty = sc$ty)
    }
  }
}
res <- do.call(rbind, results)

# Re-score the leading orientations on the FULL point sets
res <- res[order(-res$n_inliers), ]
for (i in seq_len(min(25, nrow(res)))) {
  p <- rotate_cloud(lx, ly, res$rotation_deg[i], res$mirror[i], res$scale[i])
  res$n_inliers[i] <- count_inliers(p$x + res$tx[i], p$y + res$ty[i], rx, ry)
}
res$rotation_signed <- ifelse(res$rotation_deg > 180,
                              res$rotation_deg - 360, res$rotation_deg)
res <- res[order(-res$n_inliers), ]
rownames(res) <- NULL

cat("\n=== Top 10 orientation hypotheses ===\n")
print(head(res[, c("rotation_signed", "mirror", "scale", "n_inliers", "tx", "ty")], 10),
      digits = 4)

best <- res[1, ]
n_l <- length(lx)
cat(sprintf(paste0(
  "\nBest: rotation = %.0f deg, mirror = %s, scale = %.2f -> ",
  "%d of %d LDIR particles land within %d um of a Raman particle\n",
  "Nearest 90-degree convention for config$ldir_rotate_deg_for_alignment: %d\n"),
  best$rotation_signed, best$mirror, best$scale, best$n_inliers, n_l,
  INLIER_UM,
  (round(best$rotation_signed / 90) * 90) %% 360 -
    ifelse((round(best$rotation_signed / 90) * 90) %% 360 > 180, 360, 0)))

# --- Compare with what the pipeline chose --------------------------------------
tp <- file.path(run_dir, "04_alignment", "transform_params_ldir_raman.txt")
if (file.exists(tp)) {
  cat("\n=== Pipeline's chosen transform (", tp, ") ===\n", sep = "")
  writeLines(readLines(tp))
  ln <- readLines(tp)
  g <- function(p) as.numeric(sub(".*:\\s*", "", grep(p, ln, value = TRUE)[1]))
  pipe_rot <- g("^rotation_deg")
  pipe_ref <- grepl("TRUE", grep("^reflected", ln, value = TRUE)[1])
  # The pipeline's rotation_deg is the RESIDUAL after the configured
  # pre-rotation (normalize_coords_ldir). The pre-rotation is only recorded
  # when the run had debug=TRUE (debug/ldir_norm_params.json); otherwise we
  # check which of the four allowed conventions reconciles the pipeline with
  # the measured consensus.
  pre <- NA
  np <- file.path(run_dir, "debug", "ldir_norm_params.json")
  if (file.exists(np)) {
    nj <- tryCatch(jsonlite::fromJSON(np), error = function(e) NULL)
    if (!is.null(nj$rotate_deg_applied)) pre <- nj$rotate_deg_applied
  }
  if (!is.na(pre)) {
    total_pipe <- ((pipe_rot + pre + 180) %% 360) - 180
    d_rot <- abs(((total_pipe - best$rotation_signed + 180) %% 360) - 180)
    cat(sprintf(paste0(
      "\nPipeline total rotation (pre-rotation %s deg + RANSAC/ICP %.1f deg) ",
      "= %.1f deg, reflected = %s\n"), pre, pipe_rot, total_pipe, pipe_ref))
    if (d_rot <= 15 && identical(pipe_ref, best$mirror)) {
      cat("=> Pipeline orientation AGREES with the measured consensus.\n",
          "   If matching is still poor, the problem is downstream",
          " (match radius, join, or filters), not the rotation.\n", sep = "")
    } else {
      cat(sprintf(paste0(
        "=> MISMATCH: pipeline landed %.0f deg / mirror=%s away from the ",
        "measured consensus (%.0f deg, mirror=%s).\n",
        "   The aligner picked a spurious optimum for this dataset.\n"),
        d_rot, pipe_ref, best$rotation_signed, best$mirror))
    }
  } else {
    cat(sprintf(paste0(
      "\nPipeline RANSAC/ICP residual rotation = %.1f deg, reflected = %s.\n",
      "The pre-rotation (config$ldir_rotate_deg_for_alignment) was not ",
      "recorded in this run\n(debug/ldir_norm_params.json missing — run with ",
      "config$debug = TRUE to record it).\nChecking which convention would ",
      "reconcile the pipeline with the measured consensus (%.0f deg):\n"),
      pipe_rot, pipe_ref, best$rotation_signed))
    for (cand in c(0L, 90L, -90L, 180L)) {
      total <- ((pipe_rot + cand + 180) %% 360) - 180
      d_rot <- abs(((total - best$rotation_signed + 180) %% 360) - 180)
      cat(sprintf("  pre-rotation %4d deg -> total %6.1f deg : %s\n",
                  cand, total,
                  if (d_rot <= 15 && identical(pipe_ref, best$mirror))
                    "MATCHES the consensus" else "off"))
    }
    cat("If the run used the config default (-90) and that row says 'off',\n",
        "the aligner picked a spurious optimum for this dataset.\n", sep = "")
  }
} else {
  cat("\n(No transform_params_ldir_raman.txt in this run — pipeline comparison skipped.)\n")
}
