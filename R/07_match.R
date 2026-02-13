# =============================================================================
# 07_match.R — Particle matching
# =============================================================================

#' Match FTIR particles to Raman particles based on spatial proximity
#'
#' For each transformed FTIR particle, finds the nearest Raman particle.
#' A match is accepted if:
#'   - spatial distance ≤ threshold
#'   - (optional) size similarity constraint is met
#'
#' Uses a greedy one-to-one assignment: once a Raman particle is matched,
#' it cannot be matched again. Pairs are accepted in order of increasing distance.
#'
#' @param ftir_df FTIR data frame with x_aligned, y_aligned (in Raman coordinate frame)
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param config Configuration list
#' @return List with:
#'   matched         — data frame of matched pairs with columns from both datasets
#'   unmatched_ftir  — data frame of FTIR particles without a match
#'   unmatched_raman — data frame of Raman particles without a match
#'   match_stats     — summary statistics
match_particles <- function(ftir_df, raman_df, config) {
  log_message("Matching particles")

  dist_thresh  <- config$match_dist_threshold_um
  size_weight  <- config$match_size_weight
  size_metric  <- config$match_size_metric

  ftir_x  <- ftir_df$x_aligned
  ftir_y  <- ftir_df$y_aligned
  raman_x <- raman_df$x_norm
  raman_y <- raman_df$y_norm
  n_ftir  <- nrow(ftir_df)
  n_raman <- nrow(raman_df)

  # Find k nearest Raman neighbors for each FTIR particle
  # Use k > 1 so we can fall back if first choice is already taken
  k <- min(5, n_raman)
  nn <- RANN::nn2(cbind(raman_x, raman_y), cbind(ftir_x, ftir_y), k = k)

  # Build candidate list: (ftir_idx, raman_idx, distance, combined_score)
  candidates <- data.frame(
    ftir_idx  = integer(),
    raman_idx = integer(),
    distance  = numeric(),
    score     = numeric()
  )

  for (i in seq_len(n_ftir)) {
    for (j in seq_len(k)) {
      d <- nn$nn.dists[i, j]
      if (d > dist_thresh) next

      r_idx <- nn$nn.idx[i, j]

      # Optional size similarity score
      size_penalty <- 0
      if (size_weight > 0 && size_metric %in% names(ftir_df) &&
          size_metric %in% names(raman_df)) {
        s_ftir  <- ftir_df[[size_metric]][i]
        s_raman <- raman_df[[size_metric]][r_idx]
        if (!is.na(s_ftir) && !is.na(s_raman) && (s_ftir + s_raman) > 0) {
          size_penalty <- size_weight * abs(s_ftir - s_raman) / ((s_ftir + s_raman) / 2)
        }
      }

      score <- d + size_penalty

      candidates <- rbind(candidates, data.frame(
        ftir_idx  = i,
        raman_idx = r_idx,
        distance  = d,
        score     = score
      ))
    }
  }

  # Greedy one-to-one matching: accept pairs in order of increasing score
  candidates <- candidates[order(candidates$score), ]

  matched_ftir  <- integer()
  matched_raman <- integer()
  match_rows    <- list()

  ftir_taken  <- logical(n_ftir)
  raman_taken <- logical(n_raman)

  for (row_i in seq_len(nrow(candidates))) {
    fi <- candidates$ftir_idx[row_i]
    ri <- candidates$raman_idx[row_i]

    if (ftir_taken[fi] || raman_taken[ri]) next

    ftir_taken[fi]  <- TRUE
    raman_taken[ri] <- TRUE

    match_rows[[length(match_rows) + 1]] <- data.frame(
      ftir_idx       = fi,
      raman_idx      = ri,
      match_distance = candidates$distance[row_i],
      match_score    = candidates$score[row_i]
    )
  }

  if (length(match_rows) > 0) {
    match_info <- do.call(rbind, match_rows)
  } else {
    match_info <- data.frame(ftir_idx = integer(), raman_idx = integer(),
                             match_distance = numeric(), match_score = numeric())
  }

  # Build matched data frame with columns from both datasets
  if (nrow(match_info) > 0) {
    ftir_matched  <- ftir_df[match_info$ftir_idx, ]
    raman_matched <- raman_df[match_info$raman_idx, ]

    # Prefix columns to avoid collisions
    names(ftir_matched)  <- paste0("ftir_", names(ftir_matched))
    names(raman_matched) <- paste0("raman_", names(raman_matched))

    matched <- cbind(
      match_id = seq_len(nrow(match_info)),
      match_info,
      ftir_matched,
      raman_matched
    )
    rownames(matched) <- NULL
  } else {
    matched <- data.frame()
  }

  # Unmatched particles
  unmatched_ftir  <- ftir_df[!ftir_taken, ]
  unmatched_raman <- raman_df[!raman_taken, ]

  # Statistics
  match_stats <- list(
    n_ftir           = n_ftir,
    n_raman          = n_raman,
    n_matched        = nrow(match_info),
    n_unmatched_ftir = sum(!ftir_taken),
    n_unmatched_raman = sum(!raman_taken),
    match_rate_ftir  = if (n_ftir > 0) nrow(match_info) / n_ftir else 0,
    match_rate_raman = if (n_raman > 0) nrow(match_info) / n_raman else 0,
    mean_distance    = if (nrow(match_info) > 0) mean(match_info$match_distance) else NA_real_,
    median_distance  = if (nrow(match_info) > 0) median(match_info$match_distance) else NA_real_,
    max_distance     = if (nrow(match_info) > 0) max(match_info$match_distance) else NA_real_
  )

  log_message("  Matched: ", match_stats$n_matched, " pairs")
  log_message("  Unmatched FTIR:  ", match_stats$n_unmatched_ftir)
  log_message("  Unmatched Raman: ", match_stats$n_unmatched_raman)
  log_message("  Match rate (FTIR):  ", round(match_stats$match_rate_ftir * 100, 1), "%")
  log_message("  Mean match distance: ", round(match_stats$mean_distance, 2), " µm")

  list(
    matched         = matched,
    unmatched_ftir  = unmatched_ftir,
    unmatched_raman = unmatched_raman,
    match_stats     = match_stats
  )
}
