# =============================================================================
# validation.R -- how far the results can be trusted, for docs/validation.md
# =============================================================================
#
#   source("R/load.R"); load_pipeline(); source("tools/validation.R")
#   result <- run_pipeline(config, report = FALSE)
#   v <- validate_run(result, out_dir = "validation")
#
# or, for the bundled example:  Rscript examples/validate_example_A.R
#
# Four checks, each on the pipeline's own aligned particles:
#
#   null model     How many matches would arise by chance at the gate? The
#                  source cloud is displaced so that true correspondences are
#                  destroyed but density and overlap are kept, and re-matched:
#                    rotation  random rotation about the reference centroid,
#                              excluding angles that move the median particle
#                              by less than two gates
#                    shift     toroidal shift (>= 2 gates) within the reference
#                              bounding box
#                    Poisson   n_src * (1 - exp(-density * pi * gate^2))
#   gate           matches and agreement as the distance gate varies (the
#                  alignment is kept)
#   alignment      the alignment re-run with one parameter changed: transform,
#                  matches, and the share of the baseline pairs that survive
#   LDIR join      whether the LDIR table rows (material, size) were attached
#                  to the right image positions: the agreement of each LDIR
#                  particle with its Raman partner, against the same with the
#                  LDIR rows randomly re-assigned to the positions
#
# Nothing here changes the pipeline or its outputs.
# =============================================================================


# ---- Pairs in the shared frame ----------------------------------------------

#' Every instrument pair of a run, as (aligned source, reference) particle sets.
#'
#' @param result run_pipeline() result
#' @return Named list of list(src, ref, src_label, ref_label, matched)
validation_pairs <- function(result) {
  ref <- result$ref$for_match
  pairs <- list()
  for (nm in names(result$pairs)) {
    p <- result$pairs[[nm]]
    if (is.data.frame(p$src) && nrow(p$src) > 0 && !is.null(p$match))
      pairs[[nm]] <- list(src = p$src, ref = ref, src_label = p$labels$src_prefix,
                          ref_label = p$labels$ref_prefix, matched = p$match$matched)
  }
  pairs
}

# The distance gate match_particles() uses for a pair.
.pair_gate <- function(config, src_label, ref_label) {
  if ((src_label == "ldir" || ref_label == "ldir") && !is.null(config$match_dist_threshold_ldir_um))
    config$match_dist_threshold_ldir_um
  else config$match_dist_threshold_um
}

.n_matched <- function(src, ref, config, src_label, ref_label) {
  m <- suppressMessages(match_particles(src, ref, config,
                                        src_label = src_label, ref_label = ref_label))
  m$match_stats$n_matched
}

# Quiet: the pipeline logs every step; a sweep would print thousands of lines.
.quietly <- function(expr) suppressMessages(expr)


# ---- Null model -------------------------------------------------------------

#' Expected number of chance matches for one pair.
#'
#' @param p One entry of validation_pairs()
#' @param config Pipeline config
#' @param n_iter Displacements per null
#' @param seed RNG seed (the function restores the caller's RNG state)
#' @return One-row data.frame
null_model_matches <- function(p, config, n_iter = 200, seed = 2026) {
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  src <- p$src[is.finite(p$src$x_aligned) & is.finite(p$src$y_aligned), ]
  ref <- p$ref[is.finite(p$ref$x_norm) & is.finite(p$ref$y_norm), ]
  gate <- .pair_gate(config, p$src_label, p$ref_label)
  match_n <- function(s) .n_matched(s, ref, config, p$src_label, p$ref_label)
  observed <- match_n(src)

  # Rotations small enough to leave a typical particle within two gates of
  # where it was would keep the true correspondences: excluded.
  cx <- mean(ref$x_norm); cy <- mean(ref$y_norm)
  r_med <- stats::median(sqrt((src$x_aligned - cx)^2 + (src$y_aligned - cy)^2))
  th_min <- min(pi / 2, 2 * asin(min(1, gate / max(r_med, 1e-9))))   # chord = 2 gates
  rot <- vapply(seq_len(n_iter), function(i) {
    th <- stats::runif(1, th_min, 2 * pi - th_min)
    x <- src$x_aligned - cx; y <- src$y_aligned - cy; s <- src
    s$x_aligned <- cos(th) * x - sin(th) * y + cx
    s$y_aligned <- sin(th) * x + cos(th) * y + cy
    match_n(s)
  }, numeric(1))

  rx <- range(ref$x_norm); ry <- range(ref$y_norm)
  wx <- diff(rx); wy <- diff(ry); m <- 2 * gate
  shift <- if (wx > 2 * m && wy > 2 * m) vapply(seq_len(n_iter), function(i) {
    s <- src
    s$x_aligned <- (src$x_aligned - rx[1] + stats::runif(1, m, wx - m)) %% wx + rx[1]
    s$y_aligned <- (src$y_aligned - ry[1] + stats::runif(1, m, wy - m)) %% wy + ry[1]
    match_n(s)
  }, numeric(1)) else NA_real_

  density <- nrow(ref) / (wx * wy)
  poisson <- nrow(src) * (1 - exp(-density * pi * gate^2))
  q95 <- function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.95, na.rm = TRUE))
  data.frame(
    n_src = nrow(src), n_ref = nrow(ref), gate_um = gate, observed = observed,
    chance_rotation_mean = mean(rot), chance_rotation_q95 = q95(rot),
    chance_shift_mean = mean(shift), chance_shift_q95 = q95(shift),
    chance_poisson = poisson,
    chance_share_pct = round(100 * max(mean(rot), mean(shift, na.rm = TRUE), poisson,
                                       na.rm = TRUE) / max(observed, 1), 1),
    n_iter = n_iter)
}


# ---- Gate sensitivity -------------------------------------------------------

#' Matches, distances and agreement as the distance gate varies.
#'
#' @param p One entry of validation_pairs()
#' @param gates Gates to try (um); NULL = 0.5x to 2x the pair's gate
#' @return data.frame, one row per gate
sensitivity_gate <- function(p, config, gates = NULL) {
  g0 <- .pair_gate(config, p$src_label, p$ref_label)
  if (is.null(gates)) gates <- sort(unique(round(g0 * c(0.5, 0.75, 1, 1.25, 1.5, 2))))
  key <- c(paste0(p$src_label, "_particle_id"), paste0(p$ref_label, "_particle_id"))
  base <- paste(p$matched[[key[1]]], p$matched[[key[2]]])
  do.call(rbind, lapply(gates, function(g) {
    cfg <- config
    if (p$src_label == "ldir" || p$ref_label == "ldir") cfg$match_dist_threshold_ldir_um <- g
    else cfg$match_dist_threshold_um <- g
    m <- suppressMessages(match_particles(p$src, p$ref, cfg,
                                          src_label = p$src_label, ref_label = p$ref_label))
    a <- suppressMessages(analyze_agreement(m, cfg,
                                            instrument_a = toupper(p$src_label),
                                            instrument_b = toupper(p$ref_label)))
    tr <- a$tiered_rates
    pairs <- paste(m$matched[[key[1]]], m$matched[[key[2]]])
    data.frame(gate_um = g, n_matched = m$match_stats$n_matched,
               median_distance_um = m$match_stats$median_distance,
               baseline_pairs_kept_pct = round(100 * mean(base %in% pairs), 1),
               n_scored = tr$n_total %||% 0, family_pct = tr$family_or_better_pct %||% NA,
               exact_pct = tr$exact_pct %||% NA)
  }))
}


# ---- Alignment sensitivity --------------------------------------------------

# One parameter at a time, from the run's config. landmark_skip_full_ransac =
# FALSE runs Tier 2 (RANSAC vs global registration) even when the landmark
# tier is confident: do both tiers give the same answer?
default_alignment_sweep <- list(
  landmark_skip_full_ransac = c(TRUE, FALSE),
  landmark_min_size_um      = c(50, 100, 200),
  align_raman_min_size_um   = c(10, 20, 40),
  ldir_tps_refine           = c(FALSE, TRUE)
)
# The Tier 2 parameters matter only when Tier 2 runs: swept with it forced.
tier2_alignment_sweep <- list(
  ransac_inlier_dist_um = c(100, 150, 200, 300)
)

#' Re-run the alignment of every pair with one parameter changed at a time.
#'
#' FTIR-type pairs re-run align_and_match(); LDIR re-runs its stage from the
#' joined coordinates (the image extraction and join are kept, so this varies
#' the alignment only).
#'
#' @param result run_pipeline() result; result$config is the baseline
#' @param sweep Named list: config key -> values
#' @param base Config values applied before every sweep value (e.g.
#'   list(landmark_skip_full_ransac = FALSE))
#' @return data.frame, one row per (pair, parameter, value)
sensitivity_alignment <- function(result, sweep = default_alignment_sweep, base = list()) {
  config <- result$config
  base_label <- if (length(base)) paste(names(base), unlist(base), sep = "=", collapse = ", ") else ""
  baseline <- validation_pairs(result)
  pair_keys <- function(m, sl, rl) paste(m[[paste0(sl, "_particle_id")]], m[[paste0(rl, "_particle_id")]])
  rows <- list()
  for (key in names(sweep)) for (val in sweep[[key]]) {
    cfg <- utils::modifyList(config, base); cfg[[key]] <- val
    cfg$output_dir <- tempfile("validation_"); dir.create(cfg$output_dir)
    out_dirs <- get_output_dirs(cfg$output_dir)
    ref_role <- result$reference %||% "raman"
    ref <- .quietly(prepare_reference(result$raw[[ref_role]], cfg, ref_role))
    runs <- list()
    ftir_run <- NULL
    for (role in c("ftir", "ftir_bruker")) {
      if (identical(role, ref_role) || is.null(result$raw[[role]])) next
      nm <- pair_labels(role, ref_role)$stem
      if (!nm %in% names(baseline)) next
      runs[[nm]] <- .quietly(align_and_match(result$raw[[role]], ref, cfg, ROLE_NAME[[role]], role))
      if (role == "ftir") ftir_run <- runs[[nm]]
    }
    nm <- pair_labels("ldir", ref_role)$stem
    if (nm %in% names(baseline)) {
      st <- .quietly(run_ldir_stage(result$ldir$ldir_clean, ref, ftir_run, cfg, out_dirs))
      runs[[nm]] <- list(icp = st$results$ldir_icp, match = st$results$ldir_ref_match,
                         method = st$results$ldir_aligned$align_method[1])
    }
    for (nm in names(runs)) {
      r <- runs[[nm]]; b <- baseline[[nm]]
      if (is.null(r$match)) next
      tp <- extract_transform_params(r$icp$transform)
      # The coarse (pre-ICP) transform, where the result keeps it
      pre <- if (identical(r$method, "landmark")) r$landmark$transform else r$ransac$transform
      pre_rot <- if (is.matrix(pre)) round(extract_transform_params(pre)$rotation_deg, 2) else NA
      kept <- mean(pair_keys(b$matched, b$src_label, b$ref_label) %in%
                   pair_keys(r$match$matched, b$src_label, b$ref_label))
      rows[[length(rows) + 1]] <- data.frame(
        pair = nm, with = base_label, parameter = key, value = as.character(val),
        is_baseline = length(base) == 0 && identical(config[[key]], val),
        method = r$method %||% NA_character_, coarse_rotation_deg = pre_rot,
        scale = round(tp$scale, 4), rotation_deg = round(tp$rotation_deg, 2),
        tx = round(tp$tx, 1), ty = round(tp$ty, 1),
        n_matched = r$match$match_stats$n_matched,
        median_distance_um = round(r$match$match_stats$median_distance, 1),
        baseline_pairs_kept_pct = round(100 * kept, 1))
    }
    unlink(cfg$output_dir, recursive = TRUE)
  }
  do.call(rbind, rows)
}


# ---- LDIR join check --------------------------------------------------------
# ldir_join_check() is in R/ldir_join_check.R (the pipeline logs it too).

# ---- Everything -------------------------------------------------------------

#' Run all checks on a pipeline result and write them to out_dir.
#'
#' @param result run_pipeline() result
#' @param out_dir Folder for the CSVs (NULL = do not write)
#' @param n_iter Iterations of the null model
#' @param alignment Also run the (slower) alignment sweep
#' @return list(null_model, gate, alignment, ldir_join)
validate_run <- function(result, out_dir = NULL, n_iter = 200, alignment = TRUE) {
  config <- result$config
  pairs <- validation_pairs(result)
  null <- do.call(rbind, lapply(names(pairs), function(nm)
    cbind(pair = nm, .quietly(null_model_matches(pairs[[nm]], config, n_iter)))))
  gate <- do.call(rbind, lapply(names(pairs), function(nm)
    cbind(pair = nm, .quietly(sensitivity_gate(pairs[[nm]], config)))))
  align <- if (isTRUE(alignment))
    rbind(sensitivity_alignment(result),
          sensitivity_alignment(result, tier2_alignment_sweep,
                                base = list(landmark_skip_full_ransac = FALSE)))
  join <- result$ldir$join_check %||%
    .quietly(ldir_join_check(result$ldir$ldir_ref_match$matched,
                             ref_prefix = result$ldir$labels$ref_prefix %||% "raman"))
  v <- list(null_model = null, gate = gate, alignment = align, ldir_join = join)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(null, file.path(out_dir, "null_model.csv"), row.names = FALSE)
    utils::write.csv(gate, file.path(out_dir, "sensitivity_gate.csv"), row.names = FALSE)
    if (!is.null(align))
      utils::write.csv(align, file.path(out_dir, "sensitivity_alignment.csv"), row.names = FALSE)
    if (!is.null(join)) {
      utils::write.csv(join$summary, file.path(out_dir, "ldir_join_check.csv"), row.names = FALSE)
      utils::write.csv(join$pairs, file.path(out_dir, "ldir_join_pairs.csv"), row.names = FALSE)
    }
  }
  v
}
