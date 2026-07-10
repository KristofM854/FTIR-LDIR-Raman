# =============================================================================
# diagnose_matching.R — quantify LDIR<->Raman match headroom for a run
# =============================================================================
# Question: are genuinely co-located LDIR and Raman particles being left
# unmatched only because the distance gate (match_dist_threshold_um) is too
# tight for LDIR's coarser spatial precision?
#
# Method (self-contained per run, no re-run needed):
#   1. From 05_matches/matched_ldir_raman.csv, fit the affine map that the
#      pipeline used to bring LDIR into the Raman reference frame (raw
#      ldir_x_um/y_um -> ldir_x_aligned/y_aligned) and the Raman recentring
#      (raman_x_um/y_um -> raman_x_norm/y_norm).
#   2. Apply those maps to ALL particles (LDIR from 02_joined/
#      ldir_joined_raw.csv, Raman from 01_ingested/raman_ingested.csv), so
#      every particle sits in the common aligned frame.
#   3. For each currently-UNMATCHED LDIR particle, find its nearest
#      unmatched Raman particle. Report, per candidate distance threshold,
#      how many additional pairs would form — and how many of those are
#      "safe" (mutual nearest neighbour + size-consistent + clear margin to
#      the next-nearest), i.e. very likely true matches rather than
#      coincidences.
#
# Usage:
#   Rscript tools/diagnose_matching.R output/<run_dir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript tools/diagnose_matching.R <run_dir>")
run_dir <- args[1]
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

mpath <- file.path(run_dir, "05_matches", "matched_ldir_raman.csv")
lpath <- file.path(run_dir, "02_joined", "ldir_joined_raw.csv")
rpath <- file.path(run_dir, "01_ingested", "raman_ingested.csv")
for (p in c(mpath, lpath, rpath))
  if (!file.exists(p)) stop("Missing required file: ", p)

matched <- read.csv(mpath)
ldir    <- read.csv(lpath)
raman   <- read.csv(rpath)
if (nrow(matched) < 3)
  stop("Only ", nrow(matched), " matched LDIR-Raman pairs — need >=3 to fit ",
       "the alignment transform. Nothing to diagnose.")

# --- Fit affine map raw -> aligned frame from the matched pairs --------------
fit_affine <- function(src_x, src_y, dst_x, dst_y) {
  ok <- is.finite(src_x) & is.finite(src_y) & is.finite(dst_x) & is.finite(dst_y)
  X <- cbind(1, src_x[ok], src_y[ok])
  bx <- qr.solve(X, dst_x[ok]); by <- qr.solve(X, dst_y[ok])
  function(x, y) list(x = bx[1] + bx[2]*x + bx[3]*y,
                      y = by[1] + by[2]*x + by[3]*y)
}
ldir_to_aligned <- fit_affine(matched$ldir_x_um, matched$ldir_y_um,
                              matched$ldir_x_aligned, matched$ldir_y_aligned)
raman_to_frame  <- fit_affine(matched$raman_x_um, matched$raman_y_um,
                              matched$raman_x_norm, matched$raman_y_norm)

# Residual of the fit on the matched pairs (sanity: should be small)
la <- ldir_to_aligned(matched$ldir_x_um, matched$ldir_y_um)
fit_rms <- sqrt(mean((la$x - matched$ldir_x_aligned)^2 +
                     (la$y - matched$ldir_y_aligned)^2, na.rm = TRUE))
cat(sprintf("Matched pairs: %d | affine refit RMS: %.1f um\n",
            nrow(matched), fit_rms))

# --- Project ALL particles into the common frame ----------------------------
lok <- is.finite(ldir$x_um) & is.finite(ldir$y_um)
rok <- is.finite(raman$x_um) & is.finite(raman$y_um)
LA <- ldir_to_aligned(ldir$x_um[lok], ldir$y_um[lok])
RF <- raman_to_frame(raman$x_um[rok], raman$y_um[rok])
ldir_id  <- ldir$particle_id[lok]
raman_id <- raman$particle_id[rok]
lf_max <- ldir$feret_max_um[lok]
rf_max <- raman$feret_max_um[rok]

matched_ldir  <- unique(as.character(matched$ldir_particle_id))
matched_raman <- unique(as.character(matched$raman_particle_id))
u_l <- !(as.character(ldir_id)  %in% matched_ldir)
u_r <- !(as.character(raman_id) %in% matched_raman)
cat(sprintf("LDIR: %d with coords (%d matched, %d unmatched) | Raman: %d (%d matched, %d unmatched)\n",
            sum(lok), sum(!u_l), sum(u_l), sum(rok), sum(!u_r), sum(u_r)))
cat(sprintf("Current match distances: median %.0f um, 90th pct %.0f um, max %.0f um\n\n",
            median(matched$match_distance, na.rm = TRUE),
            quantile(matched$match_distance, 0.9, na.rm = TRUE),
            max(matched$match_distance, na.rm = TRUE)))

# ===========================================================================
# Independent re-registration: is the pipeline's alignment the limit?
# ===========================================================================
# The headroom section below uses the pipeline's transform (fit from the
# matched pairs), so a poor pipeline alignment caps the headroom it reveals.
# Here we register the RAW clouds from scratch (full rotation x scale x
# translation, one-to-one inliers) to find the best ACHIEVABLE alignment,
# independent of the pipeline, and count matches under it. If this pairs far
# more particles than the pipeline matched, the ALIGNMENT is the bottleneck,
# not the match gate.
cat("\n=== Best achievable alignment (independent re-registration) ===\n")
LX <- ldir$x_um[lok]; LY <- ldir$y_um[lok]
RX <- raman$x_um[rok]; RY <- raman$y_um[rok]
LXc <- LX - mean(LX); LYc <- LY - mean(LY)
RXc <- RX - mean(RX); RYc <- RY - mean(RY)
rspan <- function(v) { q <- quantile(v, c(.05, .95), na.rm = TRUE); max(q[2]-q[1], 1e-9) }
span_ratio <- (rspan(RXc) + rspan(RYc)) / (rspan(LXc) + rspan(LYc))
TOL <- 250

# The rotation x scale search is O(n_ldir * n_raman) per candidate; subsample
# large clouds for the search (the spatial pattern survives), then the reported
# inlier count is recomputed on the FULL clouds under the winning transform.
set.seed(7)
sidx <- function(v, m = 150) if (length(v) > m) sort(sample(length(v), m)) else seq_along(v)
si <- sidx(LXc); ri2 <- sidx(RXc)
sLXc <- LXc[si]; sLYc <- LYc[si]; sRXc <- RXc[ri2]; sRYc <- RYc[ri2]

one_to_one <- function(px, py, qx, qy, tol) {  # greedy mutual, count only
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
# translation via voting on the subsampled clouds; returns best translation
score_tr <- function(deg, s, lx, ly, rx2, ry2) {
  th <- deg * pi / 180
  X <- s * (cos(th) * lx - sin(th) * ly)
  Y <- s * (sin(th) * lx + cos(th) * ly)
  dx <- outer(rx2, X, "-"); dy <- outer(ry2, Y, "-")
  key <- paste(round(dx / TOL), round(dy / TOL))
  tb <- sort(table(key), decreasing = TRUE)
  best <- list(n = 0L, tx = 0, ty = 0)
  for (k in names(tb)[seq_len(min(5, length(tb)))]) {
    sel <- key == k; tx <- mean(dx[sel]); ty <- mean(dy[sel])
    n <- one_to_one(X + tx, Y + ty, rx2, ry2, TOL)
    if (n > best$n) best <- list(n = n, tx = tx, ty = ty)
  }
  best
}
apply_tr <- function(deg, s, tx, ty, lx, ly) {
  th <- deg * pi / 180
  list(x = s*(cos(th)*lx - sin(th)*ly) + tx, y = s*(sin(th)*lx + cos(th)*ly) + ty)
}
scales <- sort(unique(round(c(seq(0.2, 1.3, 0.05),
                             span_ratio * seq(0.8, 1.2, 0.05)), 3)))
reg_best <- list(n = -1)
for (s in scales) for (deg in seq(0, 358, by = 2)) {
  r <- score_tr(deg, s, sLXc, sLYc, sRXc, sRYc)
  if (r$n > reg_best$n) reg_best <- c(r, list(deg = deg, s = s))
}
# Recompute inlier count on the FULL clouds under the winning transform
full <- apply_tr(reg_best$deg, reg_best$s, reg_best$tx, reg_best$ty, LXc, LYc)
reg_full_n <- one_to_one(full$x, full$y, RXc, RYc, TOL)
reg_best$n <- reg_full_n
cat(sprintf("Best independent alignment: rotation %d deg, scale %.2f -> %d of %d LDIR within %d um\n",
            ((reg_best$deg + 180) %% 360) - 180, reg_best$s, reg_best$n,
            length(LX), TOL))
cat(sprintf("Pipeline matched %d. ", nrow(matched)))
if (reg_best$n >= nrow(matched) + 3) {
  cat("=> ALIGNMENT is the bottleneck: a better-fit transform pairs far more\n")
  cat("   particles than the pipeline found. Fix sparse LDIR->Raman alignment\n")
  cat("   (more anchors / global registration), not just the match gate.\n")
} else {
  cat("=> Even the best independent alignment finds few pairs, so the two\n")
  cat("   instruments largely detect different particles (limited true\n")
  cat("   overlap) — not an alignment or gate problem.\n")
}

if (sum(u_l) == 0 || sum(u_r) == 0) {
  cat("\nNo unmatched particles on one side — no per-particle headroom to analyse.\n")
  quit(save = "no")
}

# --- For each unmatched LDIR, nearest + 2nd-nearest unmatched Raman ----------
ux <- LA$x[u_l]; uy <- LA$y[u_l]; ulf <- lf_max[u_l]; uli <- ldir_id[u_l]
vx <- RF$x[u_r]; vy <- RF$y[u_r]; vrf <- rf_max[u_r]; vri <- raman_id[u_r]
near <- t(vapply(seq_along(ux), function(i) {
  d <- sqrt((vx - ux[i])^2 + (vy - uy[i])^2)
  o <- order(d)
  c(j = o[1], d1 = d[o[1]], d2 = if (length(d) > 1) d[o[2]] else Inf)
}, numeric(3)))
# mutual nearest neighbour: is unmatched-LDIR i also the nearest to Raman j?
rev_near <- vapply(seq_along(vx), function(j) {
  d <- sqrt((ux - vx[j])^2 + (uy - vy[j])^2); which.min(d)
}, integer(1))
mutual <- vapply(seq_along(ux), function(i) rev_near[near[i, "j"]] == i, logical(1))
size_ok <- {
  a <- ulf; b <- vrf[near[, "j"]]
  ok <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  r <- ifelse(ok, pmax(a, b) / pmin(a, b), Inf)
  r <= 3         # within 3x in Feret max
}

cat("=== Additional LDIR->Raman pairs recoverable by a looser gate ===\n")
cat(sprintf("%-10s %-12s %-14s %-18s\n", "gate_um", "would_match",
            "mutual_NN", "safe(mutual+size+margin)"))
for (thr in c(100, 150, 200, 250, 300, 400)) {
  would <- near[, "d1"] <= thr
  mut   <- would & mutual
  # margin: nearest clearly closer than 2nd nearest (avoids ambiguous grabs)
  margin_ok <- near[, "d1"] <= 0.6 * near[, "d2"]
  safe  <- mut & size_ok & margin_ok
  cat(sprintf("%-10d %-12d %-14d %-18d\n", thr, sum(would), sum(mut), sum(safe)))
}

cat("\nInterpretation: 'safe' pairs are unmatched particles that are each\n")
cat("other's nearest neighbour, size-consistent, and unambiguous — almost\n")
cat("certainly true matches the current gate is rejecting. If 'safe' grows\n")
cat("well beyond the current match count as the gate opens to ~200-300 um,\n")
cat("the gate is the bottleneck; if it stays flat, alignment or genuine\n")
cat("non-overlap is the limit.\n")

# List the safe recoverable pairs at 250 um for spot-checking
thr <- 250
would <- near[, "d1"] <= thr
margin_ok <- near[, "d1"] <= 0.6 * near[, "d2"]
safe <- would & mutual & size_ok & margin_ok
if (any(safe)) {
  cat(sprintf("\nRecoverable 'safe' pairs at %d um (LDIR -> Raman, dist um):\n", thr))
  idx <- which(safe)
  idx <- idx[order(near[idx, "d1"])]
  for (i in idx)
    cat(sprintf("  %-10s -> %-10s  %.0f um (feret %.0f / %.0f)\n",
                uli[i], vri[near[i, "j"]], near[i, "d1"],
                ulf[i], vrf[near[i, "j"]]))
}
