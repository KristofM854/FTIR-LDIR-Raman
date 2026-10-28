# =============================================================================
# reproducibility.R — intra-instrument replicate comparison
# =============================================================================
# Compares N repeat measurements ("runs") of the SAME filter on the SAME
# instrument to quantify reproducibility (precision) and, when the true
# reference plastic is known, accuracy (trueness). A remount may turn the
# filter by any angle: the starting pose comes from the pipeline's landmark
# search; matching, gate and pair confidence are the pipeline's.
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
# groups the pairwise run-to-run matches form, merged most confident first and
# never two particles of one run (unbiased: a particle seen only in runs 2 & 3
# still counts).
#
# Depends on the pipeline modules (R/load.R: landmark_align, match_particles,
# .position_error, chance_partner_rate, make_config) and classify_family_vec()
# from R/08b_material_map.R. No side effects at source time.

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

#' Pipeline settings for replicate comparison: one instrument, so no mirror
#' and scale fixed at 1. NULL = make_config().
repro_config <- function(config = NULL) {
  cfg <- config %||% make_config()
  cfg$ransac_allow_mirror <- FALSE
  cfg$align_scale_min <- 1
  cfg$align_scale_max <- 1
  cfg
}

# Particles per run used for the starting pose (the largest)
REPRO_POSE_LANDMARKS <- 200

#' Starting pose of `src` on `ref`: the pipeline's landmark search over all
#' rotations (large particles; scale 1, no mirror), accepted only when its
#' pairs are far above chance (landmark_align()). A remounted filter can be
#' turned by any angle; the rigid ICP alone only corrects a few degrees.
#' @return list(R, t, rotation_deg, method = "landmark" / "none", p_chance)
.repro_initial_pose <- function(src, ref, config) {
  none <- list(R = diag(2), t = colMeans(cbind(ref$x_um, ref$y_um)) -
                 colMeans(cbind(src$x_um, src$y_um)),
               rotation_deg = 0, method = "none", p_chance = NA_real_)
  if (is.null(config) || !exists("landmark_align", mode = "function") ||
      nrow(src) < 3 || nrow(ref) < 3) return(none)
  cs <- colMeans(cbind(src$x_um, src$y_um)); cr <- colMeans(cbind(ref$x_um, ref$y_um))
  # The same instrument sizes a particle the same way in every run, so its
  # largest particles are the same ones: the pose search runs on those only
  # (a dense field has thousands >= 100 um), the ICP then on all particles.
  top <- function(d) {
    k <- order(-d$feret_max_um)[seq_len(min(REPRO_POSE_LANDMARKS, nrow(d)))]
    d[k[!is.na(d$feret_max_um[k])], , drop = FALSE]
  }
  s <- top(src); s$x_norm <- s$x_um - cs[1]; s$y_norm <- s$y_um - cs[2]
  r <- top(ref); r$x_norm <- r$x_um - cr[1]; r$y_norm <- r$y_um - cr[2]
  lcfg <- config
  lcfg$landmark_min_size_um <- 0; lcfg$landmark_fiber_min_size_um <- 0
  lm <- tryCatch(suppressMessages(landmark_align(s, r, lcfg, src_label = "replicate",
                                                 ref_label = "run 1")),
                 error = function(e) NULL)
  if (is.null(lm) || !isTRUE(lm$confident) || is.null(lm$transform)) {
    none$p_chance <- lm$p_chance %||% NA_real_
    return(none)
  }
  M <- lm$transform
  R <- M[1:2, 1:2]
  R <- R / sqrt(abs(det(R)))                    # scale 1
  list(R = R, t = as.numeric(cr + M[1:2, 3] - R %*% cs),
       rotation_deg = atan2(R[2, 1], R[1, 1]) * 180 / pi,
       method = "landmark", p_chance = lm$p_chance)
}

#' Rigidly register `src` onto `ref` (same instrument).
#'
#' Iterative closest point with a full re-fit each step and a correspondence
#' gate to reject particles present in only one run, started from `init`
#' (.repro_initial_pose()) or, without it, the centroid shift. Adds
#' x_aligned/y_aligned to `src` (ref's frame); records the transform and the
#' starting pose as attributes.
repro_rigid_register <- function(src, ref, gate = 300, max_iter = 50, tol = 1e-3,
                                 init = NULL) {
  P      <- cbind(src$x_um, src$y_um)
  ref_xy <- cbind(ref$x_um, ref$y_um)
  if (nrow(P) == 0 || nrow(ref_xy) == 0) {
    src$x_aligned <- src$x_um; src$y_aligned <- src$y_um
    attr(src, "transform") <- list(R = diag(2), t = c(0, 0))
    return(src)
  }
  R  <- init$R %||% diag(2)
  t  <- init$t %||% (colMeans(ref_xy) - colMeans(P))       # centroid init
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
  attr(src, "transform") <- list(R = R, t = t, rotation_deg = atan2(R[2, 1], R[1, 1]) * 180 / pi)
  attr(src, "initial_pose") <- init
  src
}

#' Register runs 2..N onto run 1 (run 1 keeps its native coordinates), each
#' from the landmark starting pose when `config` is given.
repro_align_runs <- function(runs, gate = 300, config = NULL) {
  stopifnot(length(runs) >= 1)
  runs[[1]]$x_aligned <- runs[[1]]$x_um
  runs[[1]]$y_aligned <- runs[[1]]$y_um
  for (i in seq_along(runs)[-1]) {
    init <- if (!is.null(config)) .repro_initial_pose(runs[[i]], runs[[1]], config) else NULL
    runs[[i]] <- repro_rigid_register(runs[[i]], runs[[1]], gate = gate, init = init)
    tf <- attr(runs[[i]], "transform")
    if (!is.null(init))
      log_message(sprintf("  run %d -> run 1: rotation %.1f deg (%s)", i, tf$rotation_deg,
                          if (identical(init$method, "landmark"))
                            sprintf("landmark start, p = %.1g", init$p_chance)
                          else "no landmark pose above chance: started at 0 deg"),
                  level = if (identical(init$method, "landmark")) "INFO" else "WARN")
  }
  runs
}

# ---- spatial run-to-run matching --------------------------------------------

#' Optimal 1-to-1 spatial match between two registered runs.
#'
#' The pipeline's matcher (match_particles(): optimal assignment per group of
#' nearby particles, size-ratio guard, distance + size / shape cost). The
#' gate is 3 x the positional error fitted to this pair of runs
#' (.position_error()), within 10 um .. `gate`; `gate` itself when there is
#' no fit. Material is deliberately ignored. Returns a_idx / b_idx / dist /
#' p_true (row indices into runs a and b), with attributes gate_um, sigma_um,
#' chance_share.
repro_match_pair <- function(a, b, gate = 75, config = NULL) {
  empty <- data.frame(a_idx = integer(), b_idx = integer(), dist = numeric(), p_true = numeric())
  if (nrow(a) == 0 || nrow(b) == 0) return(empty)
  cfg <- repro_config(config)
  b$x_norm <- b$x_aligned; b$y_norm <- b$y_aligned
  rho <- chance_partner_rate(b$x_norm, b$y_norm, 1)$density
  pe <- .position_error(a$x_aligned, a$y_aligned, b$x_norm, b$y_norm, rho, gate)
  g <- if (is.finite(pe$sigma)) min(max(3 * pe$sigma, 10), gate) else gate
  cfg$match_dist_threshold_um <- g
  cfg$match_gate_from_sigma <- FALSE
  res <- suppressMessages(match_particles(a, b, cfg, src_label = "a", ref_label = "b"))
  m <- res$matched
  out <- if (is.null(m) || nrow(m) == 0) empty else
    data.frame(a_idx = m$a_idx, b_idx = m$b_idx, dist = m$match_distance,
               p_true = m$p_true %||% NA_real_)
  attr(out, "gate_um") <- g
  attr(out, "sigma_um") <- pe$sigma
  attr(out, "chance_share") <- res$match_stats$chance_share %||% NA_real_
  out
}

# ---- consensus particles -----------------------------------------------------

#' Link runs into consensus particles.
#'
#' Matches every unordered pair of runs and merges the links into groups,
#' most confident first (p_true, then distance). A link that would put two
#' particles of the same run into one group is skipped and counted
#' (attr "n_link_conflicts"), so each particle belongs to exactly one
#' consensus particle and none is dropped. Returns the membership matrix:
#' consensus id x run -> row index (NA when absent); attr "pairs" lists the
#' gate, positional error and chance share of every run pair.
repro_link_runs <- function(runs, gate = 75, config = NULL) {
  n_runs  <- length(runs)
  offsets <- as.integer(cumsum(c(0L, vapply(runs, nrow, integer(1)))))
  total   <- offsets[n_runs + 1]
  node_run <- rep(seq_len(n_runs), times = diff(offsets))
  node_row <- unlist(lapply(runs, function(d) seq_len(nrow(d))))

  links <- list(); pairs <- list()
  for (i in seq_len(n_runs - 1)) for (j in (i + 1):n_runs) {
    mm <- repro_match_pair(runs[[i]], runs[[j]], gate = gate, config = config)
    pairs[[length(pairs) + 1]] <- data.frame(run_a = i, run_b = j, n_links = nrow(mm),
      gate_um = attr(mm, "gate_um") %||% NA_real_, sigma_um = attr(mm, "sigma_um") %||% NA_real_,
      chance_share = attr(mm, "chance_share") %||% NA_real_)
    if (nrow(mm) > 0)
      links[[length(links) + 1]] <- data.frame(u = offsets[i] + mm$a_idx, v = offsets[j] + mm$b_idx,
                                               p = mm$p_true, d = mm$dist)
  }
  links <- if (length(links)) do.call(rbind, links) else data.frame(u = integer(), v = integer(),
                                                                    p = numeric(), d = numeric())
  links <- links[order(-ifelse(is.na(links$p), -1, links$p), links$d), , drop = FALSE]

  # Union-find with the set of runs in each group
  parent <- seq_len(total)
  runs_in <- lapply(seq_len(total), function(k) node_run[k])
  find <- function(k) { while (parent[k] != k) { parent[k] <<- parent[parent[k]]; k <- parent[k] }; k }
  n_conflict <- 0L
  for (r in seq_len(nrow(links))) {
    a <- find(links$u[r]); b <- find(links$v[r])
    if (a == b) next
    if (length(intersect(runs_in[[a]], runs_in[[b]]))) { n_conflict <- n_conflict + 1L; next }
    parent[b] <- a
    runs_in[[a]] <- c(runs_in[[a]], runs_in[[b]])
  }
  roots <- vapply(seq_len(total), find, integer(1))
  comp <- match(roots, unique(roots))
  memb <- matrix(NA_integer_, nrow = max(c(0L, comp)), ncol = n_runs)
  memb[cbind(comp, node_run)] <- node_row
  attr(memb, "n_link_conflicts") <- n_conflict
  attr(memb, "n_links") <- nrow(links)
  attr(memb, "pairs") <- if (length(pairs)) do.call(rbind, pairs) else NULL
  memb
}

# ---- consensus attributes + aggregate metrics -------------------------------

.repro_family <- function(mat) {
  if (exists("classify_family_vec")) classify_family_vec(mat) else as.character(mat)
}

#' Per-consensus-particle table: presence, material and size per run, plus
#' concordance / CV / positional-jitter summaries (vectorised over particles).
repro_consensus_table <- function(runs, memb) {
  n_runs <- length(runs); nc <- nrow(memb)
  col <- function(name, cast) {
    m <- matrix(cast(NA), nc, n_runs)
    for (r in seq_len(n_runs)) {
      idx <- memb[, r]; ok <- !is.na(idx)
      v <- runs[[r]][[name]]
      if (!is.null(v)) m[ok, r] <- cast(v[idx[ok]])
    }
    m
  }
  mats <- col("material", as.character); feret <- col("feret_max_um", as.numeric)
  xs <- col("x_aligned", as.numeric); ys <- col("y_aligned", as.numeric)
  fams <- matrix(.repro_family(as.vector(mats)), nc, n_runs)
  fams[is.na(mats)] <- NA
  present <- !is.na(memb)
  n_det <- rowSums(present)
  concordant <- apply(fams, 1, function(f) length(unique(f[!is.na(f)])) <= 1)
  consensus_mat <- apply(fams, 1, function(f) {
    f <- f[!is.na(f)]; if (!length(f)) NA_character_ else names(sort(table(f), decreasing = TRUE))[1] })
  fm <- rowMeans(feret, na.rm = TRUE)
  fsd <- apply(feret, 1, function(v) if (sum(is.finite(v)) >= 2) stats::sd(v, na.rm = TRUE) else NA_real_)
  # Positional jitter: max pairwise distance among the runs that saw it
  jitter <- rep(NA_real_, nc)
  for (r1 in seq_len(n_runs - 1)) for (r2 in (r1 + 1):n_runs) {
    d <- sqrt((xs[, r1] - xs[, r2])^2 + (ys[, r1] - ys[, r2])^2)
    jitter <- pmax(jitter, d, na.rm = TRUE)
  }
  out <- data.frame(consensus_id = seq_len(nc), n_runs_detected = n_det,
                    material_consensus = consensus_mat, material_concordant = concordant,
                    feret_cv = fsd / fm, pos_jitter_um = jitter, stringsAsFactors = FALSE)
  for (r in seq_len(n_runs)) {
    out[[paste0("run", r, "_detected")]] <- present[, r]
    out[[paste0("run", r, "_material")]] <- mats[, r]
    out[[paste0("run", r, "_feret_um")]] <- feret[, r]
  }
  out
}

#' Aggregate reproducibility (and, if reference_family given, accuracy) metrics.
#'
#' `reference_family` may be a single polymer-family string (monotype filter) or
#' a character vector of expected families (mixed-polymer reference standard).
#' Accuracy = fraction of all per-run calls whose family is in that set.
repro_summary <- function(runs, consensus, reference_family = NULL, memb = NULL) {
  n_runs  <- length(runs)
  counts  <- vapply(runs, nrow, integer(1))
  det     <- consensus$n_runs_detected
  multi   <- consensus[det >= 2, , drop = FALSE]

  concord_rate <- if (nrow(multi) > 0) mean(multi$material_concordant) else NA_real_

  accuracy <- NA_real_
  if (!is.null(reference_family) && length(reference_family) > 0) {
    all_fams <- unlist(lapply(runs, function(d) .repro_family(d$material)))
    all_fams <- all_fams[!is.na(all_fams)]
    accuracy <- if (length(all_fams) > 0) mean(all_fams %in% reference_family) else NA_real_
  }
  pairs <- if (!is.null(memb)) attr(memb, "pairs") else NULL

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
    median_pos_jitter_um   = stats::median(consensus$pos_jitter_um, na.rm = TRUE),
    n_link_conflicts       = if (!is.null(memb)) attr(memb, "n_link_conflicts") %||% 0L else NA_integer_,
    max_chance_share       = if (!is.null(pairs) && any(is.finite(pairs$chance_share)))
                               max(pairs$chance_share, na.rm = TRUE) else NA_real_,
    run_pairs              = pairs
  )
}

#' One-call driver: align -> link -> tabulate -> summarize.
#' @param config Pipeline config (NULL = make_config()); replicate settings
#'   (no mirror, scale 1) are applied on top (repro_config()).
run_reproducibility <- function(runs, gate = 75, align_gate = 300,
                                 reference_family = NULL, config = NULL) {
  stopifnot(length(runs) >= 2)
  cfg       <- repro_config(config)
  runs      <- repro_align_runs(runs, gate = align_gate, config = cfg)
  memb      <- repro_link_runs(runs, gate = gate, config = cfg)
  consensus <- repro_consensus_table(runs, memb)
  summary   <- repro_summary(runs, consensus, reference_family = reference_family, memb = memb)
  if (summary$n_link_conflicts > 0)
    log_message("  ", summary$n_link_conflicts, " run-to-run link(s) skipped: they would have ",
                "joined two particles of the same run", level = "WARN")
  list(runs = runs, membership = memb, consensus = consensus, summary = summary)
}

# ---- viewer-ready long table ------------------------------------------------

#' One row per (consensus particle x run it was detected in), in the aligned
#' frame. Drives the Multi-Run viewer: colour by `run`, group by `consensus_id`
#' for lines/hover, encode reproducibility from `n_runs_detected` /
#' `material_concordant`.
#'
#' `coord_match_cost` (LDIR only; NA for instruments with native coordinates)
#' is the Hungarian-assignment cost from join_ldir_coords() — how confidently
#' that particle's Excel row (material/quality) was matched to an
#' image-detected blob purely by size similarity (no position is involved in
#' that join; see multi-run placement design notes, `git show 1983824:docs/multirun_image_placement_plan.md`). A high cost here
#' flags particles whose MATERIAL LABEL could plausibly have been assigned to
#' the wrong blob, which the viewer surfaces as a false "material
#' disagreement" (red ring) between runs even when the physical particle
#' genuinely reproduced.
repro_long_table <- function(res) {
  runs <- res$runs; memb <- res$membership; cons <- res$consensus
  rows <- lapply(seq_along(runs), function(r) {
    ci <- which(!is.na(memb[, r]))
    if (!length(ci)) return(NULL)
    d <- runs[[r]]; idx <- memb[ci, r]
    opt <- function(name) if (name %in% names(d)) as.numeric(d[[name]][idx]) else NA_real_
    data.frame(
      consensus_id        = ci,
      run                 = r,
      particle_id         = as.character(d$particle_id[idx]),
      x_aligned           = d$x_aligned[idx],
      y_aligned           = d$y_aligned[idx],
      material            = as.character(d$material[idx]),
      material_family     = .repro_family(d$material[idx]),
      # Carried so the Multi-Run viewer can gate on identification quality
      # (AAU / HQI, scale depends on instrument). NA for a run whose ingester
      # did not supply it; the viewer hides its slider in that case.
      quality             = opt("quality"),
      feret_max_um        = d$feret_max_um[idx],
      area_um2            = opt("area_um2"),
      major_um            = opt("major_um"),
      image_area_um2      = opt("image_area_um2"),
      image_feret_um      = opt("image_feret_um"),
      n_runs_detected     = cons$n_runs_detected[ci],
      material_concordant = cons$material_concordant[ci],
      coord_match_cost    = opt("coord_match_cost"),
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$consensus_id, out$run), , drop = FALSE]
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
         ggplot2::geom_col(fill = "#4D4D4D", width = 0.7) +   # not an instrument colour
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
           ggplot2::geom_histogram(bins = 20, fill = "#4D4D4D", colour = "white") +
           ggplot2::labs(title = ttl("Size reproducibility (Feret CV)"),
                         x = "Coefficient of variation", y = "Particles") +
           ggplot2::theme_minimal(base_size = 14))
  }

  # 4. Positional jitter across matched particles.
  jit <- consensus$pos_jitter_um[is.finite(consensus$pos_jitter_um)]
  if (length(jit) > 0) {
    save("repro_positional_jitter.png",
         g(data.frame(j = jit), aes(x = j)) +
           ggplot2::geom_histogram(bins = 20, fill = "#4D4D4D", colour = "white") +
           ggplot2::labs(title = ttl("Localization precision (positional jitter)"),
                         x = "Max inter-run distance (µm)", y = "Particles") +
           ggplot2::theme_minimal(base_size = 14))
  }
  invisible(written)
}
