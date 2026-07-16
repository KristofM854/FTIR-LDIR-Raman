# =============================================================================
# reproducibility.R — intra-instrument replicate comparison
# =============================================================================
# Compares N repeat measurements ("runs") of the SAME filter on the SAME
# instrument to quantify reproducibility (precision) and, when the true
# reference plastic is known, accuracy (trueness).
#
# The engine is instrument-agnostic: it operates on the standardized particle
# frame every ingester emits (particle_id, x_um, y_um, area_um2, major_um,
# minor_um, feret_max_um, material, quality). Because the runs share one
# instrument coordinate frame, registration is a light rigid ICP (rotation +
# translation, scale locked to 1) that absorbs a slight filter re-seat.
#
# Matching is SPATIAL ONLY — material is never a gate — so a particle called
# PET in run 1 and PP in run 3 is still linked as one physical particle and the
# disagreement is surfaced rather than hidden. Consensus particles are the
# connected components of the pairwise run-to-run matches (unbiased: a particle
# seen only in runs 2 & 3 still counts).
#
# Depends on: RANN (nearest neighbours), clue (optimal assignment; greedy
# fallback), and classify_family_vec() from R/08b_material_map.R for polymer
# family comparison. No side effects at source time.

# ---- 2D rigid registration (Kabsch / ICP, scale = 1) ------------------------

# Best-fit rigid transform mapping rows of P onto rows of Q (both n x 2).
.repro_kabsch2d <- function(P, Q) {
  cp <- colMeans(P); cq <- colMeans(Q)
  Pc <- sweep(P, 2, cp); Qc <- sweep(Q, 2, cq)
  H  <- t(Pc) %*% Qc
  sv <- svd(H)
  d  <- sign(det(sv$v %*% t(sv$u)))          # reflection guard
  R  <- sv$v %*% diag(c(1, d)) %*% t(sv$u)
  list(R = R, t = as.numeric(cq - R %*% cp))
}

.repro_apply_rt <- function(P, R, t) {
  out <- P %*% t(R)
  out[, 1] <- out[, 1] + t[1]
  out[, 2] <- out[, 2] + t[2]
  out
}

#' Rigidly register `src` onto `ref` (same instrument, slight offset).
#'
#' Iterative closest point with a full re-fit each step and a correspondence
#' gate to reject particles present in only one run. Adds x_aligned/y_aligned
#' to `src` (ref's frame); records the transform as an attribute.
repro_rigid_register <- function(src, ref, gate = 300, max_iter = 50, tol = 1e-3) {
  P      <- cbind(src$x_um, src$y_um)
  ref_xy <- cbind(ref$x_um, ref$y_um)
  if (nrow(P) == 0 || nrow(ref_xy) == 0) {
    src$x_aligned <- src$x_um; src$y_aligned <- src$y_um
    attr(src, "transform") <- list(R = diag(2), t = c(0, 0))
    return(src)
  }
  R  <- diag(2)
  t  <- colMeans(ref_xy) - colMeans(P)       # centroid init
  cur <- .repro_apply_rt(P, R, t)
  prev_err <- Inf
  for (it in seq_len(max_iter)) {
    nn   <- RANN::nn2(ref_xy, cur, k = 1)
    keep <- nn$nn.dists[, 1] <= gate
    if (sum(keep) < 2) break
    kf  <- .repro_kabsch2d(P[keep, , drop = FALSE],
                           ref_xy[nn$nn.idx[keep, 1], , drop = FALSE])
    R <- kf$R; t <- kf$t
    cur <- .repro_apply_rt(P, R, t)
    err <- mean(nn$nn.dists[keep, 1])
    if (is.finite(prev_err) && abs(prev_err - err) < tol) break
    prev_err <- err
  }
  src$x_aligned <- cur[, 1]
  src$y_aligned <- cur[, 2]
  attr(src, "transform") <- list(R = R, t = t)
  src
}

#' Register runs 2..N onto run 1. Run 1 keeps its native coordinates.
repro_align_runs <- function(runs, gate = 300) {
  stopifnot(length(runs) >= 1)
  runs[[1]]$x_aligned <- runs[[1]]$x_um
  runs[[1]]$y_aligned <- runs[[1]]$y_um
  for (i in seq_along(runs)[-1]) {
    runs[[i]] <- repro_rigid_register(runs[[i]], runs[[1]], gate = gate)
  }
  runs
}

# ---- spatial run-to-run matching --------------------------------------------

#' Optimal 1-to-1 spatial match between two registered runs.
#'
#' Cost = aligned-coordinate distance + a mild size penalty; pairs beyond `gate`
#' are never linked. Material is deliberately ignored. Returns a data frame of
#' a_idx / b_idx / dist (row indices into runs a and b).
repro_match_pair <- function(a, b, gate = 75, size_weight = 0.2) {
  na <- nrow(a); nb <- nrow(b)
  if (na == 0 || nb == 0) return(data.frame(a_idx = integer(), b_idx = integer(),
                                             dist = numeric()))
  ax <- a$x_aligned; ay <- a$y_aligned
  bx <- b$x_aligned; by <- b$y_aligned
  k  <- min(nb, 10L)
  nn <- RANN::nn2(cbind(bx, by), cbind(ax, ay), k = k)

  BIG <- 1e12
  n   <- max(na, nb)
  cost <- matrix(BIG, n, n)
  for (i in seq_len(na)) {
    for (kk in seq_len(k)) {
      j <- nn$nn.idx[i, kk]; d <- nn$nn.dists[i, kk]
      if (d > gate) next
      pen <- 0
      fa <- a$feret_max_um[i]; fb <- b$feret_max_um[j]
      if (size_weight > 0 && !is.na(fa) && !is.na(fb) && fa > 0 && fb > 0)
        pen <- size_weight * abs(log(fa / fb)) * d
      cost[i, j] <- min(cost[i, j], d + pen)
    }
  }
  if (requireNamespace("clue", quietly = TRUE)) {
    asg <- as.integer(clue::solve_LSAP(cost))
  } else {
    asg <- apply(cost, 1, which.min)         # greedy fallback
  }
  out <- list()
  for (i in seq_len(na)) {
    j <- asg[i]
    if (j > nb || cost[i, j] >= BIG) next
    out[[length(out) + 1]] <- data.frame(
      a_idx = i, b_idx = j,
      dist  = sqrt((ax[i] - bx[j])^2 + (ay[i] - by[j])^2))
  }
  if (length(out) == 0) data.frame(a_idx = integer(), b_idx = integer(),
                                   dist = numeric())
  else do.call(rbind, out)
}

# ---- union-find over all pairwise matches -> consensus particles ------------

.repro_uf_new  <- function(n) seq_len(n)
.repro_uf_find <- function(parent, i) {
  i <- as.integer(i)
  while (parent[i] != i) { parent[i] <- parent[parent[i]]; i <- parent[i] }
  as.integer(i)
}

#' Link runs into consensus particles.
#'
#' Matches every unordered pair of runs, unions matched nodes, and returns a
#' membership list: consensus id -> per-run row index (NA when absent).
repro_link_runs <- function(runs, gate = 75, size_weight = 0.2) {
  n_runs  <- length(runs)
  offsets <- as.integer(cumsum(c(0L, vapply(runs, nrow, integer(1)))))
  total   <- offsets[n_runs + 1]
  parent  <- .repro_uf_new(total)

  union <- function(a, b) {
    ra <- .repro_uf_find(parent, as.integer(a))
    rb <- .repro_uf_find(parent, as.integer(b))
    if (ra != rb) parent[rb] <<- ra
  }
  for (i in seq_len(n_runs - 1)) {
    for (j in (i + 1):n_runs) {
      mm <- repro_match_pair(runs[[i]], runs[[j]], gate = gate, size_weight = size_weight)
      if (nrow(mm) == 0) next
      for (r in seq_len(nrow(mm)))
        union(offsets[i] + mm$a_idx[r], offsets[j] + mm$b_idx[r])
    }
  }
  # Group global nodes by root.
  roots <- vapply(seq_len(total), function(i) .repro_uf_find(parent, i), integer(1))
  comps <- split(seq_len(total), roots)

  # Membership matrix: one row per consensus particle, one column per run.
  memb <- matrix(NA_integer_, nrow = length(comps), ncol = n_runs)
  node_run <- rep(seq_len(n_runs), times = diff(offsets))
  node_row <- unlist(lapply(runs, function(d) seq_len(nrow(d))))
  ci <- 0L
  for (comp in comps) {
    ci <- ci + 1L
    for (nd in comp) {
      r <- node_run[nd]
      # Keep the first occurrence if a run somehow contributes twice.
      if (is.na(memb[ci, r])) memb[ci, r] <- node_row[nd]
    }
  }
  memb
}

# ---- consensus attributes + aggregate metrics -------------------------------

.repro_family <- function(mat) {
  if (exists("classify_family_vec")) classify_family_vec(mat) else as.character(mat)
}

#' Per-consensus-particle table: presence, material and size per run, plus
#' concordance / CV / positional-jitter summaries.
repro_consensus_table <- function(runs, memb) {
  n_runs <- length(runs)
  rn <- seq_len(n_runs)
  get <- function(ci, r, col) {
    idx <- memb[ci, r]
    if (is.na(idx)) NA else runs[[r]][[col]][idx]
  }
  rows <- lapply(seq_len(nrow(memb)), function(ci) {
    present <- rn[!is.na(memb[ci, ])]
    mats  <- vapply(rn, function(r) as.character(get(ci, r, "material")), character(1))
    feret <- vapply(rn, function(r) as.numeric(get(ci, r, "feret_max_um")), numeric(1))
    xs    <- vapply(rn, function(r) as.numeric(get(ci, r, "x_aligned")), numeric(1))
    ys    <- vapply(rn, function(r) as.numeric(get(ci, r, "y_aligned")), numeric(1))
    fams  <- .repro_family(mats)
    fams_p <- fams[present]
    feret_p <- feret[present]
    concordant <- length(unique(fams_p[!is.na(fams_p)])) <= 1
    feret_cv <- if (sum(is.finite(feret_p)) >= 2)
                  stats::sd(feret_p, na.rm = TRUE) / mean(feret_p, na.rm = TRUE) else NA_real_
    # Positional jitter: max pairwise distance among present runs (aligned).
    jitter <- NA_real_
    if (length(present) >= 2) {
      pxy <- cbind(xs[present], ys[present])
      dd  <- as.matrix(stats::dist(pxy))
      jitter <- max(dd)
    }
    row <- list(
      consensus_id = ci, n_runs_detected = length(present),
      material_consensus = if (all(is.na(fams_p))) NA_character_
                           else names(sort(table(fams_p), decreasing = TRUE))[1],
      material_concordant = concordant,
      feret_cv = feret_cv, pos_jitter_um = jitter)
    for (r in rn) {
      row[[paste0("run", r, "_detected")]] <- !is.na(memb[ci, r])
      row[[paste0("run", r, "_material")]] <- mats[r]
      row[[paste0("run", r, "_feret_um")]] <- feret[r]
    }
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Aggregate reproducibility (and, if reference_family given, accuracy) metrics.
repro_summary <- function(runs, consensus, reference_family = NULL) {
  n_runs  <- length(runs)
  counts  <- vapply(runs, nrow, integer(1))
  det     <- consensus$n_runs_detected
  multi   <- consensus[det >= 2, , drop = FALSE]

  concord_rate <- if (nrow(multi) > 0) mean(multi$material_concordant) else NA_real_

  accuracy <- NA_real_
  if (!is.null(reference_family)) {
    all_fams <- unlist(lapply(runs, function(d) .repro_family(d$material)))
    all_fams <- all_fams[!is.na(all_fams)]
    accuracy <- if (length(all_fams) > 0) mean(all_fams == reference_family) else NA_real_
  }

  list(
    n_runs                 = n_runs,
    counts_per_run         = counts,
    count_mean             = mean(counts),
    count_cv               = if (mean(counts) > 0) stats::sd(counts) / mean(counts) else NA_real_,
    n_consensus            = nrow(consensus),
    detected_in_all        = sum(det == n_runs),
    detected_in_all_frac   = if (nrow(consensus) > 0) mean(det == n_runs) else NA_real_,
    detection_breakdown    = table(factor(det, levels = seq_len(n_runs))),
    material_concordance   = concord_rate,
    accuracy_vs_reference  = accuracy,
    median_feret_cv        = stats::median(consensus$feret_cv, na.rm = TRUE),
    median_pos_jitter_um   = stats::median(consensus$pos_jitter_um, na.rm = TRUE)
  )
}

#' One-call driver: align -> link -> tabulate -> summarize.
run_reproducibility <- function(runs, gate = 75, align_gate = 300,
                                 size_weight = 0.2, reference_family = NULL) {
  stopifnot(length(runs) >= 2)
  runs      <- repro_align_runs(runs, gate = align_gate)
  memb      <- repro_link_runs(runs, gate = gate, size_weight = size_weight)
  consensus <- repro_consensus_table(runs, memb)
  summary   <- repro_summary(runs, consensus, reference_family = reference_family)
  list(runs = runs, membership = memb, consensus = consensus, summary = summary)
}

# ---- overview plots ---------------------------------------------------------

#' Write the four overview plots (detection consistency, per-run material
#' counts, Feret CV, positional jitter) as PNGs into out_dir. No-op with a
#' warning when ggplot2 is unavailable. Returns the file paths (invisibly).
repro_plots <- function(res, out_dir, title_prefix = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available — skipping reproducibility plots")
    return(invisible(character(0)))
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  g  <- ggplot2::ggplot
  aes <- ggplot2::aes
  written <- character(0)
  save <- function(name, plot, w = 7, h = 5) {
    fp <- file.path(out_dir, name)
    ggplot2::ggsave(fp, plot, width = w, height = h, dpi = 120)
    written <<- c(written, fp)
  }
  n_runs    <- res$summary$n_runs
  consensus <- res$consensus
  ttl <- function(s) if (nzchar(title_prefix)) paste0(title_prefix, " — ", s) else s

  # 1. Detection consistency: how many particles were seen in k of N runs.
  det_df <- as.data.frame(table(factor(consensus$n_runs_detected,
                                        levels = seq_len(n_runs))))
  names(det_df) <- c("runs_detected", "n")
  save("repro_detection_consistency.png",
       g(det_df, aes(x = runs_detected, y = n)) +
         ggplot2::geom_col(fill = "#1f77b4", width = 0.7) +
         ggplot2::geom_text(aes(label = n), vjust = -0.3) +
         ggplot2::labs(title = ttl("Detection consistency"),
                       x = paste0("Detected in k of ", n_runs, " runs"),
                       y = "Particles") +
         ggplot2::theme_minimal(base_size = 14))

  # 2. Per-run particle count by polymer family.
  fam_rows <- do.call(rbind, lapply(seq_len(n_runs), function(r) {
    fam <- .repro_family(res$runs[[r]]$material)
    d <- as.data.frame(table(family = fam))
    if (nrow(d) == 0) return(NULL)
    data.frame(run = paste0("run", r), family = d$family, n = d$Freq)
  }))
  if (!is.null(fam_rows) && nrow(fam_rows) > 0) {
    save("repro_material_counts.png",
         g(fam_rows, aes(x = family, y = n, fill = run)) +
           ggplot2::geom_col(position = "dodge") +
           ggplot2::labs(title = ttl("Particle count per run by material"),
                         x = "Polymer family", y = "Particles", fill = NULL) +
           ggplot2::theme_minimal(base_size = 14) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
  }

  # 3. Feret CV distribution across matched particles.
  cv <- consensus$feret_cv[is.finite(consensus$feret_cv)]
  if (length(cv) > 0) {
    save("repro_feret_cv.png",
         g(data.frame(cv = cv), aes(x = cv)) +
           ggplot2::geom_histogram(bins = 20, fill = "#2ca02c", colour = "white") +
           ggplot2::labs(title = ttl("Size reproducibility (Feret CV)"),
                         x = "Coefficient of variation", y = "Particles") +
           ggplot2::theme_minimal(base_size = 14))
  }

  # 4. Positional jitter across matched particles.
  jit <- consensus$pos_jitter_um[is.finite(consensus$pos_jitter_um)]
  if (length(jit) > 0) {
    save("repro_positional_jitter.png",
         g(data.frame(j = jit), aes(x = j)) +
           ggplot2::geom_histogram(bins = 20, fill = "#d62728", colour = "white") +
           ggplot2::labs(title = ttl("Localization precision (positional jitter)"),
                         x = "Max inter-run distance (µm)", y = "Particles") +
           ggplot2::theme_minimal(base_size = 14))
  }
  invisible(written)
}
