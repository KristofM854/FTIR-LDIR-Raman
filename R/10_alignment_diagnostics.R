# =============================================================================
# 10_alignment_diagnostics.R — LDIR–Raman alignment residual diagnostics
# =============================================================================
#
# Produces per-run diagnostic artifacts in output/<run>/debug/ that reveal
# the nature and magnitude of any remaining alignment error:
#
#   ldir_raman_residuals.csv        — per-particle nearest-neighbour residuals
#   ldir_raman_residual_quiver.png  — quiver (vector field) plot of residuals
#   ldir_raman_residual_stats.json  — summary stats (mean/median/p90/max)
#   ldir_raman_transform_audit.json — normalization params + transform matrix
#
# All functions are designed to be called once per pipeline run (not in Shiny).
# =============================================================================


# ---------------------------------------------------------------------------
# 1. compute_alignment_residuals
# ---------------------------------------------------------------------------

#' Compute per-LDIR-particle nearest-Raman-particle residual vectors
#'
#' Uses a kd-tree (RANN::nn2) so this stays O(n log n) even for large clouds.
#'
#' @param ldir_df   Data frame with x_aligned, y_aligned columns
#' @param raman_df  Data frame with x_norm, y_norm columns
#' @return Data frame with columns:
#'   ldir_row, ldir_x_aligned, ldir_y_aligned,
#'   raman_row, raman_x_norm, raman_y_norm,
#'   dx, dy, dist
#'   or NULL if either cloud is empty.
compute_alignment_residuals <- function(ldir_df, raman_df) {
  ldir_valid  <- ldir_df[is.finite(ldir_df$x_aligned)  &
                          is.finite(ldir_df$y_aligned),  ]
  raman_valid <- raman_df[is.finite(raman_df$x_norm) &
                           is.finite(raman_df$y_norm),  ]
  if (nrow(ldir_valid) == 0 || nrow(raman_valid) == 0) return(NULL)

  nn <- RANN::nn2(
    data  = cbind(raman_valid$x_norm, raman_valid$y_norm),
    query = cbind(ldir_valid$x_aligned, ldir_valid$y_aligned),
    k     = 1
  )
  ridx <- nn$nn.idx[, 1]

  data.frame(
    ldir_row       = seq_len(nrow(ldir_valid)),
    ldir_x_aligned = ldir_valid$x_aligned,
    ldir_y_aligned = ldir_valid$y_aligned,
    raman_row      = ridx,
    raman_x_norm   = raman_valid$x_norm[ridx],
    raman_y_norm   = raman_valid$y_norm[ridx],
    dx             = raman_valid$x_norm[ridx] - ldir_valid$x_aligned,
    dy             = raman_valid$y_norm[ridx] - ldir_valid$y_aligned,
    dist           = nn$nn.dists[, 1],
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# 2. write_residual_quiver_png
# ---------------------------------------------------------------------------

#' Save a quiver (vector field) diagnostic plot of alignment residuals
#'
#' Shows LDIR (green) and Raman (blue) point clouds together with orange
#' arrows pointing from each LDIR aligned position toward its nearest Raman
#' neighbour.  Downsamples to max_arrows to avoid clutter.
#'
#' @param residuals   Output of compute_alignment_residuals()
#' @param raman_df    Full Raman data frame (x_norm, y_norm)
#' @param output_path PNG path to write
#' @param max_arrows  Maximum number of residual arrows (default 200)
write_residual_quiver_png <- function(residuals, raman_df, output_path,
                                       max_arrows = 200) {
  tryCatch({
    if (is.null(residuals) || nrow(residuals) == 0) return(invisible(NULL))

    arrow_df <- if (nrow(residuals) > max_arrows) {
      residuals[sample.int(nrow(residuals), max_arrows), ]
    } else {
      residuals
    }

    raman_valid <- raman_df[is.finite(raman_df$x_norm) &
                             is.finite(raman_df$y_norm), ]

    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = raman_valid,
        ggplot2::aes(x = x_norm, y = y_norm),
        shape = 1, colour = "steelblue", alpha = 0.45, size = 1.5
      ) +
      ggplot2::geom_point(
        data = residuals,
        ggplot2::aes(x = ldir_x_aligned, y = ldir_y_aligned),
        shape = 2, colour = "forestgreen", alpha = 0.45, size = 1.5
      ) +
      ggplot2::geom_segment(
        data = arrow_df,
        ggplot2::aes(
          x    = ldir_x_aligned, y    = ldir_y_aligned,
          xend = raman_x_norm,   yend = raman_y_norm
        ),
        colour    = "orange",
        alpha     = 0.55,
        linewidth = 0.35,
        arrow     = ggplot2::arrow(
          length = ggplot2::unit(0.06, "inches"),
          type   = "open"
        )
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title    = "LDIR\u2013Raman residual vector field",
        subtitle = paste0(
          "Green \u25b2 = LDIR aligned  |  Blue \u25cb = Raman  |  ",
          "Orange arrows = residuals  |  ",
          nrow(residuals), " LDIR pts, ",
          min(nrow(arrow_df), max_arrows), " arrows shown"
        ),
        x = "\u00b5m (normalised)", y = "\u00b5m (normalised)"
      )

    ggplot2::ggsave(output_path, p, width = 9, height = 9, dpi = 150)
    log_message("  Residual quiver saved: ", output_path)
  }, error = function(e) {
    log_message("  Could not save residual quiver: ", e$message, level = "WARN")
  })
  invisible(output_path)
}


# ---------------------------------------------------------------------------
# 3. write_residual_stats_json
# ---------------------------------------------------------------------------

#' Compute and save alignment residual summary statistics
#'
#' @param residuals     Output of compute_alignment_residuals() (or NULL)
#' @param n_ldir        Total LDIR particle count (before NN subset)
#' @param n_raman       Total Raman particle count
#' @param output_path   JSON path to write
#' @return Named list of stats (invisibly) so the caller can pass to audit JSON
write_residual_stats_json <- function(residuals, n_ldir, n_raman, output_path) {
  if (!is.null(residuals) && nrow(residuals) > 0) {
    d <- residuals$dist
    stats <- list(
      n_ldir_points    = n_ldir,
      n_raman_points   = n_raman,
      n_residual_pairs = nrow(residuals),
      mean_residual    = round(mean(d),            2),
      median_residual  = round(stats::median(d),   2),
      p90_residual     = round(stats::quantile(d, 0.9, names = FALSE), 2),
      max_residual     = round(max(d),             2)
    )
  } else {
    stats <- list(
      n_ldir_points    = n_ldir,
      n_raman_points   = n_raman,
      n_residual_pairs = 0L,
      mean_residual    = NA_real_,
      median_residual  = NA_real_,
      p90_residual     = NA_real_,
      max_residual     = NA_real_
    )
  }
  tryCatch(
    writeLines(
      jsonlite::toJSON(stats, pretty = TRUE, auto_unbox = TRUE, null = "null",
                       na = "null"),
      output_path
    ),
    error = function(e)
      log_message("  Could not write residual stats: ", e$message, level = "WARN")
  )
  invisible(stats)
}


# ---------------------------------------------------------------------------
# 4. write_transform_audit_json
# ---------------------------------------------------------------------------

#' Decompose a 3×3 similarity matrix into human-readable parameters
.decompose_transform <- function(M) {
  a  <- M[1, 1]; b <- M[2, 1]
  sc <- sqrt(a^2 + b^2)
  list(
    scale        = round(sc,             6),
    rotation_deg = round(atan2(b, a) * 180 / pi, 4),
    translation_x = round(M[1, 3],       2),
    translation_y = round(M[2, 3],       2)
  )
}

#' Save a transform audit JSON with normalization params, matrix and diagnostics
#'
#' @param ldir_norm_params  List: centroid_x, centroid_y, scale_factor, y_flip_applied
#' @param raman_norm_params List: centroid_x, centroid_y, scale_factor
#' @param transform_matrix  3×3 numeric matrix (LDIR_norm → Raman_norm)
#' @param residual_stats    List returned by write_residual_stats_json() or NULL
#' @param output_path       JSON path to write
write_transform_audit_json <- function(ldir_norm_params, raman_norm_params,
                                        transform_matrix, residual_stats,
                                        output_path) {
  tryCatch({
    decomposed <- .decompose_transform(transform_matrix)

    # inlier_ratio: fraction of pairs with residual < 500 µm
    inlier_ratio <- if (!is.null(residual_stats) &&
                         !is.null(residual_stats$n_residual_pairs) &&
                         residual_stats$n_residual_pairs > 0 &&
                         !is.null(residual_stats$p90_residual)) {
      # Approximate from p90 — we don't have the full distribution here;
      # p90 < 500 means at least 90 % are inliers.  Report as a string flag.
      if (!is.na(residual_stats$p90_residual) &&
          residual_stats$p90_residual < 500) ">0.90" else "<0.90"
    } else {
      NA_character_
    }

    audit <- list(
      ldir_normalization = list(
        centroid_x    = ldir_norm_params$centroid_x,
        centroid_y    = ldir_norm_params$centroid_y,
        scale_factor  = ldir_norm_params$scale_factor,
        y_flip_applied = isTRUE(ldir_norm_params$y_flip_applied)
      ),
      raman_normalization = list(
        centroid_x   = raman_norm_params$centroid_x,
        centroid_y   = raman_norm_params$centroid_y,
        scale_factor = raman_norm_params$scale_factor
      ),
      transform_matrix = lapply(seq_len(nrow(transform_matrix)),
                                function(i) as.numeric(transform_matrix[i, ])),
      transform_decomposed = decomposed,
      diagnostics = list(
        residual_mean   = residual_stats$mean_residual,
        residual_median = residual_stats$median_residual,
        residual_p90    = residual_stats$p90_residual,
        inlier_ratio    = inlier_ratio
      )
    )

    writeLines(
      jsonlite::toJSON(audit, pretty = TRUE, auto_unbox = TRUE, null = "null",
                       na = "null"),
      output_path
    )
    log_message("  Transform audit saved: ", output_path)
  }, error = function(e) {
    log_message("  Could not write transform audit: ", e$message, level = "WARN")
  })
  invisible(output_path)
}


# ---------------------------------------------------------------------------
# 5. write_alignment_diagnostics  (orchestrator)
# ---------------------------------------------------------------------------

#' Generate all LDIR–Raman alignment diagnostics for a pipeline run
#'
#' Produces residual CSV, quiver PNG, stats JSON, and transform audit JSON.
#' All artifacts are written to <output_dir>/debug/ (created if needed).
#' Failures are caught and logged as WARNs — never interrupts the pipeline.
#'
#' @param ldir_aligned     Data frame with x_aligned, y_aligned
#' @param raman_clean      Data frame with x_norm, y_norm
#' @param ldir_norm_params List with centroid_x/y, scale_factor, y_flip_applied
#' @param norm_result      Return value of normalize_coordinates() —
#'                         used for raman_centroid and raman_scale
#' @param ldir_icp         ICP result list with $transform (3×3 matrix)
#' @param config           Pipeline config list with $output_dir
#' @return Path to the debug directory (invisibly)
write_alignment_diagnostics <- function(ldir_aligned, raman_clean,
                                         ldir_norm_params, norm_result,
                                         ldir_icp, config) {
  diag_dir <- file.path(config$output_dir, "debug")
  if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

  residuals <- tryCatch(
    compute_alignment_residuals(ldir_aligned, raman_clean),
    error = function(e) {
      log_message("  compute_alignment_residuals failed: ", e$message, level = "WARN")
      NULL
    }
  )

  # --- Residual CSV ---
  if (!is.null(residuals)) {
    tryCatch(
      utils::write.csv(residuals,
                       file.path(diag_dir, "ldir_raman_residuals.csv"),
                       row.names = FALSE),
      error = function(e)
        log_message("  Could not write residuals CSV: ", e$message, level = "WARN")
    )
  }

  # --- Quiver PNG ---
  write_residual_quiver_png(
    residuals   = residuals,
    raman_df    = raman_clean,
    output_path = file.path(diag_dir, "ldir_raman_residual_quiver.png")
  )

  # --- Stats JSON (also returned for use in audit) ---
  n_ldir  <- sum(is.finite(ldir_aligned$x_aligned))
  n_raman <- sum(is.finite(raman_clean$x_norm))
  residual_stats <- write_residual_stats_json(
    residuals   = residuals,
    n_ldir      = n_ldir,
    n_raman     = n_raman,
    output_path = file.path(diag_dir, "ldir_raman_residual_stats.json")
  )

  # --- Transform audit JSON ---
  raman_np <- list(
    centroid_x   = norm_result$raman_centroid[1],
    centroid_y   = norm_result$raman_centroid[2],
    scale_factor = norm_result$raman_scale %||% 1
  )
  write_transform_audit_json(
    ldir_norm_params = ldir_norm_params,
    raman_norm_params = raman_np,
    transform_matrix  = ldir_icp$transform,
    residual_stats    = residual_stats,
    output_path       = file.path(diag_dir, "ldir_raman_transform_audit.json")
  )

  invisible(diag_dir)
}
