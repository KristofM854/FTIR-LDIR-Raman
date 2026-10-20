# =============================================================================
# 07_match.R — Particle matching (spatial + morphological)
# =============================================================================
#
# Provides two matching strategies:
#   1. Hungarian (optimal global assignment) — default
#   2. Greedy nearest-neighbour — fallback for very large datasets
#
# Both produce identical output structure. Hungarian is preferred because it
# avoids the sub-optimal greedy order-dependence in dense regions.
#
# Also provides:
#   - Ambiguity metrics (per-match risk scoring)
#   - Composite (1:many) matching for fragmented particles
# =============================================================================


#' Match particles between two instruments based on spatial proximity
#'
#' Dispatches to Hungarian or greedy matching depending on config$match_method.
#' Default is "hungarian" (optimal, for any number of particles: it is solved
#' per group of nearby particles); "greedy" only when set explicitly.
#'
#' @param src_df Source data frame with x_aligned, y_aligned
#' @param ref_df Reference data frame with x_norm, y_norm
#' @param config Configuration list
#' @param src_label Label for the source instrument (default "ftir")
#' @param ref_label Label for the reference instrument (default "raman")
#' @return List with matched, unmatched_src, unmatched_ref, match_stats
match_particles <- function(src_df, ref_df, config,
                            src_label = "ftir", ref_label = "raman") {
  log_message("Matching particles (", src_label, " -> ", ref_label, ")")

  method <- config$match_method %||% "hungarian"

  # Optimal assignment needs clue. Switching silently to greedy would change
  # the results, so a missing package is an error, not a fallback.
  if (method == "hungarian" && !requireNamespace("clue", quietly = TRUE))
    stop("match_method = 'hungarian' needs the 'clue' package: install.packages('clue')",
         call. = FALSE)

  force_complete <- isTRUE(config$ldir_force_complete_match) && src_label == "ldir"

  # LDIR spatial precision (coarse circle-calibrated centroids on large
  # particles) is inherently lower than FTIR/Raman, so genuine LDIR matches
  # sit at larger residuals — the FTIR-tuned gate rejects them. Use a
  # dedicated, looser threshold for any pairing that involves LDIR.
  is_ldir_pair <- src_label == "ldir" || ref_label == "ldir"
  if (is_ldir_pair && !is.null(config$match_dist_threshold_ldir_um)) {
    eff <- config$match_dist_threshold_ldir_um
    log_message("  LDIR pairing: using distance gate ", eff,
                " um (default ", config$match_dist_threshold_um, " um)")
    config$match_dist_threshold_um <- eff
  }

  if (method == "hungarian") {
    result <- .match_hungarian(src_df, ref_df, config, src_label, ref_label,
                               force_complete = force_complete)
  } else {
    result <- .match_greedy(src_df, ref_df, config, src_label, ref_label)
  }

  # Add ambiguity metrics
  if (nrow(result$matched) > 0) {
    result$matched <- .compute_ambiguity(src_df, ref_df, result$matched,
                                         config, src_label, ref_label)
  }

  result$match_stats <- c(result$match_stats,
                          .chance_stats(src_df, ref_df, result$match_stats$n_matched,
                                        config$match_dist_threshold_um, config,
                                        src_label, ref_label))
  result
}


# How many of the pairs a random arrangement would also give: each source
# particle inside the reference field has a chance partner within the gate
# with p = 1 - exp(-rho * pi * gate^2) (reference density rho over its convex
# hull). Expected chance pairs = p x source particles inside the field; the
# chance share is that over the pairs found. A lower bound on chance: the
# adaptive gate of large particles is wider than the base gate used here.
.chance_stats <- function(src_df, ref_df, n_matched, gate, config, src_label, ref_label) {
  cr <- chance_partner_rate(ref_df$x_norm, ref_df$y_norm, gate)
  sx <- src_df$x_aligned; sy <- src_df$y_aligned
  n_in <- if (is.finite(cr$area_um2) && length(sx) > 0) {
    ok <- is.finite(sx) & is.finite(sy)
    rx <- ref_df$x_norm; ry <- ref_df$y_norm
    fin <- is.finite(rx) & is.finite(ry)
    h <- grDevices::chull(rx[fin], ry[fin])
    sum(.points_in_polygon(sx[ok], sy[ok], rx[fin][h], ry[fin][h]))
  } else NA_integer_
  expected <- cr$p * n_in
  share <- if (isTRUE(n_matched > 0)) min(1, expected / n_matched) else NA_real_
  lvl <- if (isTRUE(share > cfg_value(config, "chance_share_warn"))) "WARN" else "INFO"
  log_message(sprintf(paste0("  Chance pairs (%s -> %s, gate %g um, %d reference particles): ",
                             "%.1f of %d expected at random = %s of the pairs"),
                      src_label, ref_label, gate, cr$n, expected, n_matched,
                      if (is.na(share)) "NA" else paste0(round(100 * share), "%")),
              level = lvl)
  if (identical(lvl, "WARN"))
    log_message("  The reference is too dense for this gate to separate true pairs ",
                "from chance: treat these pairs and their agreement with caution",
                level = "WARN")
  list(chance_gate_um = gate, chance_p_partner = cr$p, chance_expected = expected,
       chance_share = share)
}

# Ray casting: which points lie inside the polygon (px, py)
.points_in_polygon <- function(x, y, px, py) {
  n <- length(px); inside <- logical(length(x))
  j <- n
  for (i in seq_len(n)) {
    cross <- ((py[i] > y) != (py[j] > y)) &
      (x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i] + 1e-300) + px[i])
    inside <- xor(inside, cross)
    j <- i
  }
  inside
}


# ============================================================================
# Hungarian (optimal) matching
# ============================================================================

.match_hungarian <- function(src_df, ref_df, config, src_label, ref_label,
                             force_complete = FALSE) {
  log_message("  Using Hungarian (optimal) matching",
              if (force_complete) " [force_complete]" else "")

  dist_thresh     <- config$match_dist_threshold_um
  adaptive_factor <- cfg_value(config, "match_adaptive_dist_factor")
  lam_area        <- config$match_lambda_area   %||% 0.3
  lam_feret       <- config$match_lambda_feret  %||% 0.3
  lam_aspect      <- config$match_lambda_aspect %||% 0.1
  max_ratio       <- cfg_value(config, "match_max_size_ratio")

  src_x  <- src_df$x_aligned
  src_y  <- src_df$y_aligned
  ref_x  <- ref_df$x_norm
  ref_y  <- ref_df$y_norm
  n_src  <- nrow(src_df)
  n_ref  <- nrow(ref_df)

  # Per-particle adaptive threshold (ignored when force_complete = TRUE)
  src_thresh <- rep(dist_thresh, n_src)
  if (!force_complete && adaptive_factor > 0 && "major_um" %in% names(src_df)) {
    major <- src_df$major_um
    src_thresh <- pmax(dist_thresh, adaptive_factor * ifelse(is.na(major), 0, major))
  }

  # When force_complete, search all reference particles, not just k nearest
  k <- if (force_complete) n_ref else min(20, n_ref)
  nn <- RANN::nn2(cbind(ref_x, ref_y), cbind(src_x, src_y), k = k)

  # Candidate pairs (src i, ref j, cost) within the gate; force_complete:
  # every candidate regardless of distance
  BIG <- 1e12
  cand <- vector("list", n_src)

  for (i in seq_len(n_src)) {
    thr_i <- src_thresh[i]
    js <- integer(0); cs <- numeric(0)
    for (ki in seq_len(k)) {
      j <- nn$nn.idx[i, ki]
      d <- nn$nn.dists[i, ki]
      if (!force_complete && d > thr_i) next
      if (!.size_ratio_ok(src_df, ref_df, i, j, max_ratio)) next
      # Multi-feature cost
      penalty <- 0

      if (lam_area > 0 && "area_um2" %in% names(src_df) && "area_um2" %in% names(ref_df)) {
        a_src <- src_df$area_um2[i]
        a_ref <- ref_df$area_um2[j]
        if (!is.na(a_src) && !is.na(a_ref) && a_src > 0 && a_ref > 0) {
          penalty <- penalty + lam_area * abs(log(a_src / a_ref))
        }
      }

      if (lam_feret > 0 && "feret_max_um" %in% names(src_df) && "feret_max_um" %in% names(ref_df)) {
        f_src <- src_df$feret_max_um[i]
        f_ref <- ref_df$feret_max_um[j]
        if (!is.na(f_src) && !is.na(f_ref) && f_src > 0 && f_ref > 0) {
          penalty <- penalty + lam_feret * abs(log(f_src / f_ref))
        }
      }

      if (lam_aspect > 0 && "major_um" %in% names(src_df) && "minor_um" %in% names(src_df) &&
          "major_um" %in% names(ref_df) && "minor_um" %in% names(ref_df)) {
        asp_s <- src_df$major_um[i] / max(src_df$minor_um[i], 1, na.rm = TRUE)
        asp_r <- ref_df$major_um[j] / max(ref_df$minor_um[j], 1, na.rm = TRUE)
        if (!is.na(asp_s) && !is.na(asp_r)) {
          penalty <- penalty + lam_aspect * abs(asp_s - asp_r) / max(asp_s, asp_r, 1)
        }
      }

      js <- c(js, j); cs <- c(cs, d + penalty)
    }
    if (length(js)) cand[[i]] <- list(j = js, c = cs)
  }
  ei <- rep(seq_len(n_src), vapply(cand, function(x) length(x$j), integer(1)))
  ej <- as.integer(unlist(lapply(cand, `[[`, "j")))
  ec <- unlist(lapply(cand, `[[`, "c"))

  # Solve the assignment. Pairs only exist within the gate, so the problem
  # splits into independent groups of nearby particles (connected components
  # of the candidate graph); solving each group on its own gives the same
  # optimum as one dense n x n problem (which needs n^2 memory and ~n^3 time:
  # 1.7 GB and most of the run time for a 14,000-particle field sample).
  # Objective in both: as many pairs as possible (each unpaired particle
  # costs BIG), then the smallest summed cost.
  assignment <- rep(NA_integer_, n_src)
  cost_of <- rep(NA_real_, n_src)
  if (length(ei) > 0) {
    comp <- .bipartite_components(ei, ej, n_src, n_ref)
    log_message("  ", length(unique(comp$edge)), " groups of nearby particles, largest ",
                max(tabulate(match(comp$edge, unique(comp$edge)))), " candidate pairs")
    for (g in split(seq_along(ei), comp$edge)) {
      rows <- sort(unique(ei[g])); cols <- sort(unique(ej[g]))
      m <- matrix(BIG, length(rows), length(cols))
      ri <- match(ei[g], rows); ci <- match(ej[g], cols)
      m[cbind(ri, ci)] <- pmin(m[cbind(ri, ci)], ec[g])
      if (length(rows) <= length(cols)) {
        a <- as.integer(clue::solve_LSAP(m, maximum = FALSE))
        pick <- cbind(seq_along(rows), a)
      } else {
        a <- as.integer(clue::solve_LSAP(t(m), maximum = FALSE))
        pick <- cbind(a, seq_along(cols))
      }
      ok <- m[pick] < BIG
      assignment[rows[pick[ok, 1]]] <- cols[pick[ok, 2]]
      cost_of[rows[pick[ok, 1]]] <- m[pick][ok]
    }
  }

  # Extract valid matches (cost < BIG means within threshold, or force_complete)
  match_rows <- list()
  src_taken <- logical(n_src)
  ref_taken <- logical(n_ref)

  for (i in seq_len(n_src)) {
    j <- assignment[i]
    if (is.na(j)) next  # no partner within the gate

    # Compute actual distance for this pair
    d_actual <- sqrt((src_x[i] - ref_x[j])^2 + (src_y[i] - ref_y[j])^2)

    src_taken[i] <- TRUE
    ref_taken[j] <- TRUE

    match_rows[[length(match_rows) + 1]] <- data.frame(
      src_idx        = i,
      ref_idx        = j,
      match_distance = d_actual,
      match_score    = cost_of[i]
    )
  }

  .assemble_match_result(src_df, ref_df, match_rows, src_taken, ref_taken,
                         src_label, ref_label)
}


#' May source particle i and reference particle j pair, by size? Not when their
#' Feret max differ by more than max_ratio (each instrument sizes a particle
#' somewhat differently, but a 148 um particle is not a 6 um one). Missing
#' sizes, or max_ratio NULL / Inf, do not restrict.
.size_ratio_ok <- function(src_df, ref_df, i, j, max_ratio) {
  if (is.null(max_ratio) || !is.finite(max_ratio)) return(TRUE)
  a <- src_df$feret_max_um[i]; b <- ref_df$feret_max_um[j]
  if (is.null(a) || is.null(b) || is.na(a) || is.na(b) || a <= 0 || b <= 0) return(TRUE)
  max(a, b) / min(a, b) <= max_ratio
}


#' Connected components of a bipartite candidate graph.
#'
#' @param ei,ej Edge endpoints: source row i, reference row j
#' @param n_src,n_ref Number of source / reference particles
#' @return list(edge = component id of each edge)
.bipartite_components <- function(ei, ej, n_src, n_ref) {
  parent <- seq_len(n_src + n_ref)
  find <- function(x) {
    r <- x
    while (parent[r] != r) r <- parent[r]
    while (parent[x] != r) { nx <- parent[x]; parent[x] <<- r; x <- nx }
    r
  }
  for (e in seq_along(ei)) {
    a <- find(ei[e]); b <- find(n_src + ej[e])
    if (a != b) parent[max(a, b)] <- min(a, b)
  }
  list(edge = vapply(ei, find, integer(1)))
}


# ============================================================================
# Greedy (fallback) matching — original algorithm
# ============================================================================

.match_greedy <- function(src_df, ref_df, config, src_label, ref_label) {
  log_message("  Using greedy nearest-neighbour matching")

  dist_thresh      <- config$match_dist_threshold_um
  adaptive_factor  <- cfg_value(config, "match_adaptive_dist_factor")
  size_weight      <- config$match_size_weight
  size_metric      <- config$match_size_metric
  max_ratio        <- cfg_value(config, "match_max_size_ratio")

  src_x  <- src_df$x_aligned
  src_y  <- src_df$y_aligned
  ref_x  <- ref_df$x_norm
  ref_y  <- ref_df$y_norm
  n_src  <- nrow(src_df)
  n_ref  <- nrow(ref_df)

  # Per-particle adaptive threshold
  src_thresh <- rep(dist_thresh, n_src)
  if (adaptive_factor > 0 && "major_um" %in% names(src_df)) {
    major <- src_df$major_um
    src_thresh <- pmax(dist_thresh, adaptive_factor * ifelse(is.na(major), 0, major))
  }

  k <- min(5, n_ref)
  nn <- RANN::nn2(cbind(ref_x, ref_y), cbind(src_x, src_y), k = k)

  cand_list <- vector("list", n_src * k)
  cand_n <- 0L

  for (i in seq_len(n_src)) {
    thresh_i <- src_thresh[i]
    for (j in seq_len(k)) {
      d <- nn$nn.dists[i, j]
      if (d > thresh_i) next

      r_idx <- nn$nn.idx[i, j]
      if (!.size_ratio_ok(src_df, ref_df, i, r_idx, max_ratio)) next

      size_penalty <- 0
      if (size_weight > 0 && size_metric %in% names(src_df) &&
          size_metric %in% names(ref_df)) {
        s_src  <- src_df[[size_metric]][i]
        s_ref  <- ref_df[[size_metric]][r_idx]
        if (!is.na(s_src) && !is.na(s_ref) && (s_src + s_ref) > 0) {
          size_penalty <- size_weight * abs(s_src - s_ref) / ((s_src + s_ref) / 2)
        }
      }

      score <- d + size_penalty
      cand_n <- cand_n + 1L
      cand_list[[cand_n]] <- data.frame(
        src_idx  = i,
        ref_idx  = r_idx,
        distance = d,
        score    = score
      )
    }
  }

  if (cand_n > 0) {
    candidates <- do.call(rbind, cand_list[seq_len(cand_n)])
  } else {
    candidates <- data.frame(src_idx = integer(), ref_idx = integer(),
                             distance = numeric(), score = numeric())
  }

  candidates <- candidates[order(candidates$score), ]

  match_rows <- list()
  src_taken  <- logical(n_src)
  ref_taken  <- logical(n_ref)

  for (row_i in seq_len(nrow(candidates))) {
    fi <- candidates$src_idx[row_i]
    ri <- candidates$ref_idx[row_i]
    if (src_taken[fi] || ref_taken[ri]) next

    src_taken[fi] <- TRUE
    ref_taken[ri] <- TRUE

    match_rows[[length(match_rows) + 1]] <- data.frame(
      src_idx        = fi,
      ref_idx        = ri,
      match_distance = candidates$distance[row_i],
      match_score    = candidates$score[row_i]
    )
  }

  .assemble_match_result(src_df, ref_df, match_rows, src_taken, ref_taken,
                         src_label, ref_label)
}


# ============================================================================
# Shared: assemble match result structure
# ============================================================================

.assemble_match_result <- function(src_df, ref_df, match_rows,
                                   src_taken, ref_taken,
                                   src_label, ref_label) {
  n_src <- nrow(src_df)
  n_ref <- nrow(ref_df)

  if (length(match_rows) > 0) {
    match_info <- do.call(rbind, match_rows)
  } else {
    match_info <- data.frame(src_idx = integer(), ref_idx = integer(),
                             match_distance = numeric(), match_score = numeric())
  }

  # Rename idx columns for compatibility
  names(match_info)[names(match_info) == "src_idx"] <- paste0(src_label, "_idx")
  names(match_info)[names(match_info) == "ref_idx"] <- paste0(ref_label, "_idx")

  if (nrow(match_info) > 0) {
    src_matched  <- src_df[match_info[[paste0(src_label, "_idx")]], ]
    ref_matched  <- ref_df[match_info[[paste0(ref_label, "_idx")]], ]

    names(src_matched)  <- paste0(src_label, "_", names(src_matched))
    names(ref_matched)  <- paste0(ref_label, "_", names(ref_matched))

    matched <- cbind(
      match_id = seq_len(nrow(match_info)),
      match_info,
      src_matched,
      ref_matched
    )
    rownames(matched) <- NULL
  } else {
    matched <- data.frame()
  }

  unmatched_src <- src_df[!src_taken, ]
  unmatched_ref <- ref_df[!ref_taken, ]

  match_stats <- list(
    n_src            = n_src,
    n_ref            = n_ref,
    n_matched        = nrow(match_info),
    n_unmatched_src  = sum(!src_taken),
    n_unmatched_ref  = sum(!ref_taken),
    match_rate_src   = if (n_src > 0) nrow(match_info) / n_src else 0,
    match_rate_ref   = if (n_ref > 0) nrow(match_info) / n_ref else 0,
    mean_distance    = if (nrow(match_info) > 0) mean(match_info$match_distance) else NA_real_,
    median_distance  = if (nrow(match_info) > 0) median(match_info$match_distance) else NA_real_,
    max_distance     = if (nrow(match_info) > 0) max(match_info$match_distance) else NA_real_
  )

  # Back-compat aliases
  match_stats$n_ftir           <- match_stats$n_src
  match_stats$n_raman          <- match_stats$n_ref
  match_stats$n_unmatched_ftir <- match_stats$n_unmatched_src
  match_stats$n_unmatched_raman <- match_stats$n_unmatched_ref
  match_stats$match_rate_ftir  <- match_stats$match_rate_src
  match_stats$match_rate_raman <- match_stats$match_rate_ref

  log_message("  Matched: ", match_stats$n_matched, " pairs")
  log_message("  Unmatched ", src_label, ":  ", match_stats$n_unmatched_src)
  log_message("  Unmatched ", ref_label, ": ", match_stats$n_unmatched_ref)
  log_message("  Match rate (", src_label, "):  ",
              round(match_stats$match_rate_src * 100, 1), "%")
  log_message("  Mean match distance: ",
              round(match_stats$mean_distance, 2), " µm")

  # Unmatched particles under each side's label (unmatched_ftir and
  # unmatched_raman for the default FTIR -> Raman labels), and as
  # unmatched_src / unmatched_ref whatever the labels
  result <- list(
    matched         = matched,
    match_stats     = match_stats,
    labels          = c(src = src_label, ref = ref_label),
    unmatched_src   = unmatched_src,
    unmatched_ref   = unmatched_ref
  )
  result[[paste0("unmatched_", src_label)]] <- unmatched_src
  result[[paste0("unmatched_", ref_label)]] <- unmatched_ref

  result
}


# ============================================================================
# Ambiguity metrics
# ============================================================================

#' Compute per-match ambiguity scores
#'
#' For each matched pair, counts how many reference particles are within
#' an ambiguity radius and computes a risk score.
#'
#' @param src_df Source data frame
#' @param ref_df Reference data frame
#' @param matched Matched pairs data frame
#' @param config Configuration list
#' @param src_label Source instrument label
#' @param ref_label Reference instrument label
#' @return matched data frame with added ambiguity columns
.compute_ambiguity <- function(src_df, ref_df, matched, config,
                               src_label = "ftir", ref_label = "raman") {
  amb_radius <- config$ambiguity_radius_um %||% 50

  src_x_col <- paste0(src_label, "_x_aligned")
  src_y_col <- paste0(src_label, "_y_aligned")

  # Compute all ref-to-matched distances for each match
  ref_x <- ref_df$x_norm
  ref_y <- ref_df$y_norm
  n_matched <- nrow(matched)

  matched$n_ref_candidates    <- NA_integer_
  matched$next_best_distance  <- NA_real_
  matched$ambiguity_score     <- NA_real_
  matched$ambiguity_flag      <- FALSE

  for (mk in seq_len(n_matched)) {
    sx <- matched[[src_x_col]][mk]
    sy <- matched[[src_y_col]][mk]

    if (is.na(sx) || is.na(sy)) next

    dists_all <- sqrt((ref_x - sx)^2 + (ref_y - sy)^2)

    # Count candidates within ambiguity radius
    n_cand <- sum(dists_all <= amb_radius, na.rm = TRUE)

    # Find next-best distance (exclude the matched ref)
    ref_idx_col <- paste0(ref_label, "_idx")
    matched_ref_idx <- matched[[ref_idx_col]][mk]
    dists_others <- dists_all
    if (!is.na(matched_ref_idx) && matched_ref_idx <= length(dists_others)) {
      dists_others[matched_ref_idx] <- Inf
    }
    next_best <- min(dists_others, na.rm = TRUE)

    d_match <- matched$match_distance[mk]
    amb_score <- if (is.finite(next_best) && next_best > 0) {
      d_match / next_best
    } else {
      0
    }

    matched$n_ref_candidates[mk]   <- n_cand
    matched$next_best_distance[mk] <- if (is.finite(next_best)) next_best else NA_real_
    matched$ambiguity_score[mk]    <- amb_score
    matched$ambiguity_flag[mk]     <- amb_score > 0.7 || n_cand >= 3
  }

  n_flagged <- sum(matched$ambiguity_flag, na.rm = TRUE)
  if (n_flagged > 0) {
    log_message("  Ambiguity: ", n_flagged, " of ", n_matched,
                " matches flagged as ambiguous")
  }

  matched
}


# ============================================================================
# Composite (1:many) matching
# ============================================================================

#' Find composite matches — one source particle matched to multiple ref fragments
#'
#' For unmatched source particles, checks if multiple nearby ref particles
#' collectively explain the source (combined area comparison).
#'
#' @param src_df Source data frame (all particles, with x_aligned, y_aligned)
#' @param ref_unmatched Unmatched reference particles (with x_norm, y_norm)
#' @param config Configuration list
#' @param max_radius_factor Multiplier on source major_um/2 for search radius
#' @return Data frame of composite matches (one row per source particle)
find_composite_matches <- function(src_df, ref_unmatched, config,
                                   max_radius_factor = 1.5) {
  if (nrow(src_df) == 0 || nrow(ref_unmatched) == 0) {
    return(data.frame(
      src_particle_id = character(), n_ref_parts = integer(),
      ref_ids = character(), combined_area = numeric(),
      src_area = numeric(), area_ratio = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  log_message("  Searching for composite (1:many) matches")

  composites <- list()
  ref_x <- ref_unmatched$x_norm
  ref_y <- ref_unmatched$y_norm
  # Each reference fragment belongs to at most one composite: larger source
  # particles claim their fragments first (ties by ID, so the result does not
  # depend on row order).
  used <- logical(nrow(ref_unmatched))
  visit <- order(-ifelse(is.na(src_df$area_um2), 0, src_df$area_um2),
                 as.character(src_df$particle_id))

  for (i in visit) {
    sx <- src_df$x_aligned[i]
    sy <- src_df$y_aligned[i]
    if (is.na(sx) || is.na(sy)) next

    # Search radius proportional to source particle size
    r <- max_radius_factor * src_df$major_um[i] / 2
    if (is.na(r) || r <= 0) r <- config$match_dist_threshold_um

    dists <- sqrt((ref_x - sx)^2 + (ref_y - sy)^2)
    candidates <- which(dists <= r & !used)

    if (length(candidates) < 2) next

    # Combined area check
    combined_area <- sum(ref_unmatched$area_um2[candidates], na.rm = TRUE)
    src_area <- src_df$area_um2[i]
    if (is.na(src_area) || src_area <= 0 || combined_area <= 0) next

    area_ratio <- combined_area / src_area

    if (area_ratio >= 0.3 && area_ratio <= 5.0) {
      used[candidates] <- TRUE
      composites[[length(composites) + 1]] <- data.frame(
        .row            = i,
        src_particle_id = src_df$particle_id[i],
        n_ref_parts     = length(candidates),
        ref_ids         = paste(ref_unmatched$particle_id[candidates], collapse = ";"),
        combined_area   = combined_area,
        src_area        = src_area,
        area_ratio      = round(area_ratio, 3),
        stringsAsFactors = FALSE
      )
    }
  }

  result <- if (length(composites) > 0) {
    r <- do.call(rbind, composites)
    r <- r[order(r$.row), , drop = FALSE]      # source table order
    r$.row <- NULL
    rownames(r) <- NULL
    r
  } else {
    data.frame(
      src_particle_id = character(), n_ref_parts = integer(),
      ref_ids = character(), combined_area = numeric(),
      src_area = numeric(), area_ratio = numeric(),
      stringsAsFactors = FALSE
    )
  }

  log_message("  Found ", nrow(result), " composite matches")
  result
}


# Null-coalescing operator (if not already defined)
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
