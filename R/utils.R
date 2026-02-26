# =============================================================================
# utils.R — Shared utility functions for FTIR–Raman particle matching
# =============================================================================

# ---------------------------------------------------------------------------
# Output directory with timestamped subfolders
# ---------------------------------------------------------------------------

#' Create a timestamped run subfolder inside the base output directory
#'
#' Format: output/YYYY-MM-DD_1, output/YYYY-MM-DD_2, etc.
#' Automatically increments the run number for multiple runs on the same day.
#'
#' @param base_dir Base output directory (e.g., "output")
#' @return Path to the new run-specific subfolder (already created)
make_run_dir <- function(base_dir = "output") {
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  today <- format(Sys.Date(), "%Y-%m-%d")
  existing <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

  # Find existing run numbers for today
  pattern <- paste0("^", gsub("-", "\\\\-", today), "_(\\d+)$")
  matches <- regmatches(existing, regexec(pattern, existing))
  run_numbers <- as.integer(vapply(matches, function(m) {
    if (length(m) == 2) m[2] else NA_character_
  }, character(1)))
  run_numbers <- run_numbers[!is.na(run_numbers)]

  next_run <- if (length(run_numbers) == 0) 1L else max(run_numbers) + 1L
  run_dir <- file.path(base_dir, paste0(today, "_", next_run))
  dir.create(run_dir, recursive = TRUE)

  log_message("Run output directory: ", run_dir)
  run_dir
}

# ---------------------------------------------------------------------------
# Coordinate parsing
# ---------------------------------------------------------------------------

#' Parse FTIR coordinate string "[x;y]" into numeric x, y
#' @param coord_str Character vector of strings like "[1234;5678]" or "[1234.5;5678.9]"
#' @return Data frame with columns x_um and y_um
parse_ftir_coordinates <- function(coord_str) {
  # Remove brackets and split on semicolon

  cleaned <- gsub("\\[|\\]", "", coord_str)
  parts   <- strsplit(cleaned, ";")

  x_um <- vapply(parts, function(p) as.numeric(p[1]), numeric(1))
  y_um <- vapply(parts, function(p) as.numeric(p[2]), numeric(1))

  data.frame(x_um = x_um, y_um = y_um)
}

# ---------------------------------------------------------------------------
# 2D similarity transform utilities (3x3 homogeneous matrices)
# ---------------------------------------------------------------------------

#' Build a 3x3 homogeneous similarity transform matrix
#'
#' Without reflection:
#'   | a  -b  tx |      a = s*cos(theta), b = s*sin(theta)
#'   | b   a  ty |
#'   | 0   0   1 |
#'
#' With reflection (y-flip before rotation+scale):
#'   | a   b  tx |
#'   | b  -a  ty |
#'   | 0   0   1 |
#'
#' @param a Numeric, s*cos(theta)
#' @param b Numeric, s*sin(theta)
#' @param tx Numeric, x-translation
#' @param ty Numeric, y-translation
#' @param reflect Logical, whether the transform includes reflection
#' @return 3x3 matrix
build_transform_matrix <- function(a, b, tx, ty, reflect = FALSE) {
  if (!reflect) {
    matrix(c(a, b, 0,
             -b, a, 0,
             tx, ty, 1), nrow = 3, byrow = FALSE)
  } else {
    matrix(c(a, b, 0,
             b, -a, 0,
             tx, ty, 1), nrow = 3, byrow = FALSE)
  }
}

#' Estimate a 2D similarity transform from point correspondences
#'
#' Given source points P and destination points Q, find the similarity
#' transform T such that T(P) ≈ Q, minimizing sum of squared residuals.
#' Optional per-point weights allow down-weighting unreliable correspondences
#' (e.g. elongated particles whose centroid is uncertain).
#'
#' @param src_x Numeric vector, source x-coordinates
#' @param src_y Numeric vector, source y-coordinates
#' @param dst_x Numeric vector, destination x-coordinates
#' @param dst_y Numeric vector, destination y-coordinates
#' @param allow_reflection Logical, whether to also try reflection and pick best
#' @param weights Optional numeric vector of per-point weights (positive).
#'        Higher weight = more influence. NULL means equal weights.
#' @return List with: matrix (3x3), a, b, tx, ty, scale, rotation_deg,
#'         reflected, residual_rms
estimate_similarity_transform <- function(src_x, src_y, dst_x, dst_y,
                                          allow_reflection = TRUE,
                                          weights = NULL) {
  n <- length(src_x)
  stopifnot(n >= 2, n == length(src_y), n == length(dst_x), n == length(dst_y))

  # Response vector
  q_vec <- as.numeric(rbind(dst_x, dst_y))  # interleaved: qx1, qy1, qx2, qy2, ...

  # Build per-row weight vector (each point contributes 2 rows: x and y)
  if (!is.null(weights)) {
    stopifnot(length(weights) == n)
    w_sqrt <- sqrt(pmax(weights, 0))
    w_row <- rep(w_sqrt, each = 2)  # interleaved: w1, w1, w2, w2, ...
  } else {
    w_row <- NULL
  }

  # --- No-reflection model ---
  # qx_i = a*px_i - b*py_i + tx
  # qy_i = b*px_i + a*py_i + ty
  A_no <- matrix(0, nrow = 2 * n, ncol = 4)
  for (i in seq_len(n)) {
    row_x <- 2 * i - 1
    row_y <- 2 * i
    A_no[row_x, ] <- c(src_x[i], -src_y[i], 1, 0)
    A_no[row_y, ] <- c(src_y[i],  src_x[i], 0, 1)
  }

  # Apply weights via row scaling: W*A*x = W*q  (weighted least squares)
  if (!is.null(w_row)) {
    A_no_w <- A_no * w_row
    q_no_w <- q_vec * w_row
  } else {
    A_no_w <- A_no
    q_no_w <- q_vec
  }

  fit_no   <- qr.solve(A_no_w, q_no_w)
  res_no   <- q_vec - A_no %*% fit_no  # residuals on unweighted data
  rms_no   <- sqrt(mean(res_no^2))

  results <- list()
  results$no_reflect <- list(
    a = fit_no[1], b = fit_no[2], tx = fit_no[3], ty = fit_no[4],
    rms = rms_no, reflect = FALSE
  )

  if (allow_reflection) {
    # --- Reflection model (y-flip before rotation+scale) ---
    # qx_i =  a*px_i + b*py_i + tx
    # qy_i =  b*px_i - a*py_i + ty
    A_ref <- matrix(0, nrow = 2 * n, ncol = 4)
    for (i in seq_len(n)) {
      row_x <- 2 * i - 1
      row_y <- 2 * i
      A_ref[row_x, ] <- c( src_x[i],  src_y[i], 1, 0)
      A_ref[row_y, ] <- c(-src_y[i],  src_x[i], 0, 1)
    }

    if (!is.null(w_row)) {
      A_ref_w <- A_ref * w_row
      q_ref_w <- q_vec * w_row
    } else {
      A_ref_w <- A_ref
      q_ref_w <- q_vec
    }

    fit_ref <- qr.solve(A_ref_w, q_ref_w)
    res_ref <- q_vec - A_ref %*% fit_ref
    rms_ref <- sqrt(mean(res_ref^2))

    results$reflect <- list(
      a = fit_ref[1], b = fit_ref[2], tx = fit_ref[3], ty = fit_ref[4],
      rms = rms_ref, reflect = TRUE
    )
  }

  # Pick the best
  if (allow_reflection && results$reflect$rms < results$no_reflect$rms) {
    best <- results$reflect
  } else {
    best <- results$no_reflect
  }

  scale     <- sqrt(best$a^2 + best$b^2)
  rot_rad   <- atan2(best$b, best$a)
  rot_deg   <- rot_rad * 180 / pi

  M <- build_transform_matrix(best$a, best$b, best$tx, best$ty, best$reflect)

  list(
    matrix       = M,
    a            = best$a,
    b            = best$b,
    tx           = best$tx,
    ty           = best$ty,
    scale        = scale,
    rotation_deg = rot_deg,
    reflected    = best$reflect,
    residual_rms = best$rms
  )
}

#' Apply a 3x3 homogeneous transform to 2D points
#'
#' @param x Numeric vector of x-coordinates
#' @param y Numeric vector of y-coordinates
#' @param M 3x3 homogeneous transform matrix
#' @return Data frame with columns x_transformed, y_transformed
apply_transform_points <- function(x, y, M) {
  pts <- rbind(x, y, rep(1, length(x)))  # 3 x n
  result <- M %*% pts                     # 3 x n
  data.frame(
    x_transformed = result[1, ],
    y_transformed = result[2, ]
  )
}

#' Compose two 3x3 transforms: apply T1 first, then T2
#' @param T1 3x3 matrix (applied first)
#' @param T2 3x3 matrix (applied second)
#' @return 3x3 matrix T2 %*% T1
compose_transforms <- function(T1, T2) {
  T2 %*% T1
}

#' Extract human-readable parameters from a 3x3 similarity transform matrix
#' @param M 3x3 homogeneous similarity transform matrix
#' @return List with scale, rotation_deg, tx, ty, reflected
extract_transform_params <- function(M) {
  a  <- M[1, 1]
  b  <- M[2, 1]
  tx <- M[1, 3]
  ty <- M[2, 3]

  # Check reflection: det of upper-left 2x2
  det_ul <- M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
  reflected <- det_ul < 0

  scale   <- sqrt(a^2 + b^2)
  rot_deg <- atan2(b, a) * 180 / pi

  list(
    scale        = scale,
    rotation_deg = rot_deg,
    tx           = tx,
    ty           = ty,
    reflected    = reflected
  )
}

#' Create a translation-only 3x3 matrix
make_translation_matrix <- function(tx, ty) {
  M <- diag(3)
  M[1, 3] <- tx
  M[2, 3] <- ty
  M
}

#' Create a rotation-only 3x3 matrix (rotation about origin)
#' @param angle_deg Rotation angle in degrees
make_rotation_matrix <- function(angle_deg) {
  theta <- angle_deg * pi / 180
  ct <- cos(theta)
  st <- sin(theta)
  matrix(c(ct, st, 0,
           -st, ct, 0,
           0, 0, 1), nrow = 3, byrow = FALSE)
}

#' Create a y-axis reflection matrix
make_mirror_y_matrix <- function() {
  matrix(c(1, 0, 0,
           0, -1, 0,
           0, 0, 1), nrow = 3, byrow = FALSE)
}

# ---------------------------------------------------------------------------
# Nearest-neighbor helpers
# ---------------------------------------------------------------------------

#' Find nearest neighbors using RANN
#' @param query_x, query_y Coordinates of query points
#' @param ref_x, ref_y Coordinates of reference points
#' @param k Number of neighbors to return
#' @return List with nn_idx (indices into ref) and nn_dist (distances)
find_nearest_neighbors <- function(query_x, query_y, ref_x, ref_y, k = 1) {
  query_mat <- cbind(query_x, query_y)
  ref_mat   <- cbind(ref_x, ref_y)

  nn <- RANN::nn2(data = ref_mat, query = query_mat, k = k)

  list(
    nn_idx  = nn$nn.idx,
    nn_dist = nn$nn.dists
  )
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log_message <- function(..., level = "INFO") {
  msg <- paste0("[", Sys.time(), "] [", level, "] ", paste0(..., collapse = ""))
  message(msg)
}

# ---------------------------------------------------------------------------
# Debug: A3 particle trace (Step 5)
# ---------------------------------------------------------------------------

#' Dump a single particle's coordinates at one pipeline stage to a CSV
#'
#' Writes one row per particle_id in \code{ids} to
#'   \code{<debug_dir>/<id>_stage_<stage>.csv}
#' If the particle is not found a "not_found" note is written instead.
#' Every write is verified with stopifnot(file.exists()).
#'
#' @param df Data frame at the current pipeline stage
#' @param ids Character vector of particle_id values to dump (e.g. c("A3","MP_11"))
#' @param stage Character label (used in filename and "stage" column)
#' @param debug_dir Path to debug directory
dump_particle <- function(df, ids, stage, debug_dir) {
  if (is.null(debug_dir) || !dir.exists(debug_dir)) return(invisible(NULL))

  all_cols <- c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
                "x_aligned", "y_aligned", "coord_source", "material",
                "feret_max_um", "quality")

  for (pid in ids) {
    row_idx <- which(df$particle_id == pid)
    out_path <- file.path(debug_dir,
                          paste0(gsub("[^A-Za-z0-9_-]", "_", pid),
                                 "_stage_", stage, ".csv"))

    if (length(row_idx) == 0) {
      row_df <- data.frame(
        particle_id = pid, stage = stage,
        note = "not_found_at_this_stage",
        stringsAsFactors = FALSE
      )
    } else {
      row_df <- df[row_idx[1], intersect(all_cols, names(df)), drop = FALSE]
      # Fill missing coord columns with NA
      for (col in all_cols) {
        if (!col %in% names(row_df)) row_df[[col]] <- NA
      }
      row_df$stage <- stage
    }

    tryCatch({
      write.csv(row_df, out_path, row.names = FALSE)
      stopifnot(file.exists(out_path))
    }, error = function(e) {
      log_message("  WARN: dump_particle write failed for '", pid,
                  "' stage '", stage, "': ", e$message, level = "WARN")
    })
  }

  log_message("  Dumped stage '", stage, "' for: ",
              paste(ids, collapse = ", "))
}


#' Trace a particle through pipeline stages
#'
#' Appends a snapshot row to debug_A3_trace.csv for each pipeline stage.
#' Traces ALL particles so the CSV can be filtered to any ID later.
#'
#' @param df Data frame with particle data at a given stage
#' @param stage Character label for this pipeline stage
#' @param config Config list (needs config$debug_dir)
trace_particle_snapshot <- function(df, stage, config) {
  if (!isTRUE(config$debug) || is.null(config$debug_dir)) return(invisible(NULL))

  trace_file <- file.path(config$debug_dir, "debug_A3_trace.csv")

  # Extract available coordinate columns
  cols_available <- intersect(
    c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
      "x_aligned", "y_aligned", "coord_source", "material",
      "feret_max_um"),
    names(df)
  )

  snapshot <- df[, cols_available, drop = FALSE]
  snapshot$stage <- stage

  # Fill missing columns with NA
  for (col in c("particle_id", "x_um", "y_um", "x_norm", "y_norm",
                "x_aligned", "y_aligned", "coord_source", "material",
                "feret_max_um")) {
    if (!col %in% names(snapshot)) snapshot[[col]] <- NA
  }

  tryCatch({
    if (file.exists(trace_file)) {
      write.table(snapshot, trace_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    } else {
      write.csv(snapshot, trace_file, row.names = FALSE)
    }
    stopifnot(file.exists(trace_file))
  }, error = function(e) {
    log_message("  WARN: trace_particle_snapshot write failed: ", e$message,
                level = "WARN")
  })

  log_message("  Traced ", nrow(snapshot), " particles at stage: ", stage)
}


# ---------------------------------------------------------------------------
# Debug: Branch A/B comparison + residual vectors (Steps 2 & 6)
# ---------------------------------------------------------------------------

#' Run Branch A/B Y-flip comparison and save debug artifacts
#'
#' Generates overlay plots for both Y-flip settings and computes quality
#' metrics (ICP RMS, median match distance) for each.
#'
#' @param ldir_aligned LDIR data frame with x_aligned, y_aligned
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param ldir_raman_match Match result for LDIR-Raman
#' @param ldir_icp ICP result for LDIR alignment
#' @param config Configuration list
debug_ldir_branches <- function(ldir_aligned, raman_df, ldir_raman_match,
                                 ldir_icp, config) {
  if (!isTRUE(config$debug) || is.null(config$debug_dir)) return(invisible(NULL))

  debug_dir <- config$debug_dir

  tryCatch({
    # Current branch overlay (the active flip setting)
    branch_label <- if (isTRUE(config$ldir_flip_y_for_alignment)) "A" else "B"

    # Step 5: Write overlay_plot_data_ldir.csv so we can verify plotted columns
    overlay_csv <- file.path(debug_dir, "overlay_plot_data_ldir.csv")
    plot_cols <- intersect(c("particle_id", "x_aligned", "y_aligned",
                              "x_norm", "y_norm", "x_um", "y_um",
                              "material", "feret_max_um", "coord_source"),
                           names(ldir_aligned))
    write.csv(ldir_aligned[seq_len(min(200, nrow(ldir_aligned))), plot_cols,
                            drop = FALSE],
              overlay_csv, row.names = FALSE)
    stopifnot(file.exists(overlay_csv))
    log_message("  Debug: wrote overlay_plot_data_ldir.csv (",
                min(200, nrow(ldir_aligned)), " rows, cols: ",
                paste(plot_cols, collapse = ","), ")")

    # A3 in the final overlay (what is actually plotted)
    trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
    for (pid in trace_ids) {
      row_idx <- which(ldir_aligned$particle_id == pid)
      if (length(row_idx) > 0) {
        r <- ldir_aligned[row_idx[1], plot_cols, drop = FALSE]
        log_message("  Debug overlay coords for '", pid, "': ",
                    "x_aligned=", round(r$x_aligned, 2),
                    " y_aligned=", round(r$y_aligned, 2))
      } else {
        log_message("  Debug: '", pid, "' NOT FOUND in ldir_aligned", level = "WARN")
      }
    }

    # Save current branch overlay
    icp_rms_val <- if (length(ldir_icp$rms_history) > 0)
      round(tail(ldir_icp$rms_history, 1), 1) else NA_real_

    p_current <- plot_overlay(
      ldir_aligned, raman_df, ldir_raman_match,
      ftir_color = "darkgreen", raman_color = "steelblue",
      src_label = "ldir"
    ) + ggplot2::labs(
      title = paste0("LDIR-Raman Overlay — Branch ", branch_label,
                      " (flip_y=", config$ldir_flip_y_for_alignment, ")"),
      subtitle = paste0("Green=LDIR(x_aligned/y_aligned), Blue=Raman(x_norm/y_norm) | ",
                        "ICP RMS=", icp_rms_val, " \u00b5m")
    )
    overlay_png <- file.path(debug_dir,
                              paste0("overlay_branch", branch_label, ".png"))
    ggplot2::ggsave(overlay_png, p_current, width = 10, height = 8, dpi = 150)
    stopifnot(file.exists(overlay_png))

    # Compute metrics for current branch
    matched <- ldir_raman_match$matched
    metrics <- data.frame(
      branch = branch_label,
      flip_y = config$ldir_flip_y_for_alignment,
      icp_rms = icp_rms_val,
      n_matched = nrow(matched),
      median_match_dist = if (nrow(matched) > 0)
        round(median(matched$match_distance), 2) else NA_real_,
      stringsAsFactors = FALSE
    )

    # Write branch summary
    summary_file <- file.path(debug_dir, "branch_summary.txt")
    summary_lines <- c(
      paste0("Branch ", branch_label, " (ldir_flip_y_for_alignment = ",
             config$ldir_flip_y_for_alignment, "):"),
      paste0("  ICP RMS:             ", metrics$icp_rms, " \u00b5m"),
      paste0("  Matched pairs:       ", metrics$n_matched),
      paste0("  Median match dist:   ", metrics$median_match_dist, " \u00b5m"),
      ""
    )
    writeLines(summary_lines, summary_file)
    stopifnot(file.exists(summary_file))
    log_message("  Debug: saved Branch ", branch_label, " overlay + metrics")

    # --- Residual vectors (Step 6) ---
    if (nrow(matched) > 0) {
      src_x_col <- "ldir_x_aligned"
      src_y_col <- "ldir_y_aligned"
      ref_x_col <- "raman_x_norm"
      ref_y_col <- "raman_y_norm"

      if (all(c(src_x_col, src_y_col, ref_x_col, ref_y_col) %in% names(matched))) {
        resid_df <- data.frame(
          x = matched[[src_x_col]],
          y = matched[[src_y_col]],
          dx = matched[[ref_x_col]] - matched[[src_x_col]],
          dy = matched[[ref_y_col]] - matched[[src_y_col]],
          dist = matched$match_distance
        )
        resid_df <- resid_df[complete.cases(resid_df), ]

        if (nrow(resid_df) > 0) {
          p_resid <- ggplot2::ggplot(resid_df) +
            ggplot2::geom_segment(
              ggplot2::aes(x = x, y = y,
                           xend = x + dx * 5, yend = y + dy * 5,
                           colour = dist),
              arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")),
              linewidth = 0.5, alpha = 0.7
            ) +
            ggplot2::scale_colour_viridis_c(name = "Distance (µm)") +
            ggplot2::coord_equal() +
            ggplot2::labs(
              title = "LDIR Residual Vector Field",
              subtitle = "Arrows magnified 5x | Systematic patterns = local distortion",
              x = "X (µm)", y = "Y (µm)"
            ) +
            ggplot2::theme_minimal()

          ggplot2::ggsave(file.path(debug_dir, "ldir_residual_vectors.png"),
                          p_resid, width = 10, height = 8, dpi = 150)
          log_message("  Debug: saved ldir_residual_vectors.png")
        }
      }
    }

  }, error = function(e) {
    log_message("  Debug branch comparison failed: ", e$message, level = "WARN")
  })
}
