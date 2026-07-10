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

if (sum(u_l) == 0 || sum(u_r) == 0) {
  cat("No unmatched particles on one side — no headroom to analyse.\n"); quit(save = "no")
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
