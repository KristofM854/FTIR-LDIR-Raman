# =============================================================================
# 10_export.R — Export results for downstream use
# =============================================================================

#' Export all pipeline outputs
#'
#' @param match_result Result from match_particles()
#' @param agreement Result from analyze_agreement()
#' @param diagnostics List of ggplot objects from generate_diagnostics()
#' @param icp_result Result from icp_refine()
#' @param norm_result Result from normalize_coordinates()
#' @param config Configuration list
#' @param ftir_scan_bounds Optional scan bounds
#' @param ldir_results Optional list of LDIR pipeline results
#' @return Invisible NULL
export_results <- function(match_result, agreement, diagnostics,
                           icp_result, norm_result, config,
                           ftir_scan_bounds = NULL,
                           ldir_results = NULL) {
  out_dir <- config$output_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  log_message("Exporting results to: ", out_dir)

  # --- 1. Matched particle table ---
  if (nrow(match_result$matched) > 0) {
    write.csv(match_result$matched,
              file.path(out_dir, "matched_particles.csv"),
              row.names = FALSE)
    log_message("  Wrote matched_particles.csv (", nrow(match_result$matched), " pairs)")
  }

  # --- 2. Unmatched particles ---
  if (nrow(match_result$unmatched_ftir) > 0) {
    write.csv(match_result$unmatched_ftir,
              file.path(out_dir, "unmatched_ftir.csv"),
              row.names = FALSE)
  }
  if (nrow(match_result$unmatched_raman) > 0) {
    write.csv(match_result$unmatched_raman,
              file.path(out_dir, "unmatched_raman.csv"),
              row.names = FALSE)
  }

  # --- 3. Transform parameters ---
  params <- icp_result$params
  param_lines <- c(
    "# FTIR-to-Raman Transform Parameters",
    paste0("# Generated: ", Sys.time()),
    "",
    paste0("scale:        ", round(params$scale, 6)),
    paste0("rotation_deg: ", round(params$rotation_deg, 4)),
    paste0("tx:           ", round(params$tx, 4)),
    paste0("ty:           ", round(params$ty, 4)),
    paste0("reflected:    ", params$reflected),
    "",
    "# Normalization parameters (applied before transform)",
    paste0("ftir_centroid_x:  ", round(norm_result$ftir_centroid[1], 4)),
    paste0("ftir_centroid_y:  ", round(norm_result$ftir_centroid[2], 4)),
    paste0("raman_centroid_x: ", round(norm_result$raman_centroid[1], 4)),
    paste0("raman_centroid_y: ", round(norm_result$raman_centroid[2], 4)),
    paste0("ftir_scale:       ", round(norm_result$ftir_scale, 6)),
    paste0("raman_scale:      ", round(norm_result$raman_scale, 6)),
    "",
    "# ICP refinement info",
    paste0("icp_converged:    ", icp_result$converged),
    paste0("icp_iterations:   ", icp_result$n_iterations),
    paste0("icp_final_rms:    ",
           round(tail(icp_result$rms_history, 1), 4), " um"),
    ""
  )

  if (!is.null(ftir_scan_bounds)) {
    param_lines <- c(param_lines,
      "# FTIR scan area (physical extent of the background image)",
      paste0("ftir_scan_xmin:   ", ftir_scan_bounds$x_min),
      paste0("ftir_scan_xmax:   ", ftir_scan_bounds$x_max),
      paste0("ftir_scan_ymin:   ", ftir_scan_bounds$y_min),
      paste0("ftir_scan_ymax:   ", ftir_scan_bounds$y_max),
      ""
    )
  }

  param_lines <- c(param_lines,
    "# 3x3 Transform matrix (homogeneous, FTIR_norm -> Raman_norm)",
    paste0("matrix_row1: ", paste(round(icp_result$transform[1, ], 8), collapse = ", ")),
    paste0("matrix_row2: ", paste(round(icp_result$transform[2, ], 8), collapse = ", ")),
    paste0("matrix_row3: ", paste(round(icp_result$transform[3, ], 8), collapse = ", "))
  )
  writeLines(param_lines,
             file.path(out_dir, "transform_params.txt"))
  log_message("  Wrote transform_params.txt")

  # --- 4. Agreement summary (tiered) ---
  if (!is.null(agreement) && !is.null(agreement$agreement_detail) &&
      nrow(agreement$agreement_detail) > 0) {
    write.csv(agreement$agreement_detail,
              file.path(out_dir, "agreement_summary.csv"),
              row.names = FALSE)
    write.csv(agreement$summary_df,
              file.path(out_dir, "agreement_pairwise.csv"),
              row.names = FALSE)
    log_message("  Wrote agreement CSVs")
  }

  # Category summary
  if (!is.null(agreement$category_summary) && nrow(agreement$category_summary) > 0) {
    write.csv(agreement$category_summary,
              file.path(out_dir, "agreement_by_category.csv"),
              row.names = FALSE)
    log_message("  Wrote agreement_by_category.csv")
  }

  # --- 5. Match statistics ---
  stats <- match_result$match_stats
  stats_lines <- c(
    "# Match Statistics",
    paste0("# Generated: ", Sys.time()),
    "",
    paste0("n_ftir:            ", stats$n_ftir),
    paste0("n_raman:           ", stats$n_raman),
    paste0("n_matched:         ", stats$n_matched),
    paste0("n_unmatched_ftir:  ", stats$n_unmatched_ftir),
    paste0("n_unmatched_raman: ", stats$n_unmatched_raman),
    paste0("match_rate_ftir:   ", round(stats$match_rate_ftir * 100, 1), "%"),
    paste0("match_rate_raman:  ", round(stats$match_rate_raman * 100, 1), "%"),
    paste0("mean_distance_um:  ", round(stats$mean_distance, 3)),
    paste0("median_distance_um:", round(stats$median_distance, 3)),
    paste0("max_distance_um:   ", round(stats$max_distance, 3))
  )

  # Tiered agreement rates
  if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
    tr <- agreement$tiered_rates
    stats_lines <- c(stats_lines, "",
      "# Tiered agreement",
      paste0("agreement_n_pairs:       ", tr$n_total),
      paste0("agreement_exact_pct:     ", tr$exact_pct, "%"),
      paste0("agreement_family_pct:    ", tr$family_pct, "%"),
      paste0("agreement_family_or_better_pct: ", tr$family_or_better_pct, "%")
    )
  }

  writeLines(stats_lines,
             file.path(out_dir, "match_statistics.txt"))
  log_message("  Wrote match_statistics.txt")

  # --- 6. Triage export: top mismatched/high-risk pairs ---
  export_triage(match_result, agreement, out_dir)

  # --- 7. LDIR results ---
  if (!is.null(ldir_results)) {
    export_ldir_results(ldir_results, out_dir)
  }

  # --- 8. Diagnostic plots ---
  if (length(diagnostics) > 0) {
    plot_dir <- file.path(out_dir, "plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    for (plot_name in names(diagnostics)) {
      tryCatch({
        ggplot2::ggsave(
          filename = file.path(plot_dir, paste0(plot_name, ".png")),
          plot     = diagnostics[[plot_name]],
          width    = config$plot_width,
          height   = config$plot_height,
          dpi      = 150
        )
      }, error = function(e) {
        log_message("  Warning: Could not save plot '", plot_name, "': ", e$message,
                    level = "WARN")
      })
    }

    # Combined PDF
    tryCatch({
      pdf_path <- file.path(plot_dir, "all_diagnostics.pdf")
      grDevices::pdf(pdf_path,
                     width = config$plot_width, height = config$plot_height)
      for (plot_name in names(diagnostics)) {
        print(diagnostics[[plot_name]])
      }
      grDevices::dev.off()
      log_message("  Wrote all_diagnostics.pdf")
    }, error = function(e) {
      log_message("  Warning: Could not save combined PDF: ", e$message,
                  level = "WARN")
    })
  }

  log_message("Export complete.")
  invisible(NULL)
}


#' Export triage list: top mismatched and high-risk pairs
export_triage <- function(match_result, agreement, out_dir, n_top = 50) {
  if (is.null(agreement) || nrow(agreement$summary_df) == 0) return(invisible(NULL))
  if (nrow(match_result$matched) == 0) return(invisible(NULL))

  matched <- match_result$matched
  summary <- agreement$summary_df

  triage_df <- data.frame(match_id = summary$match_id, stringsAsFactors = FALSE)

  amb_score <- if ("ambiguity_score" %in% names(matched)) {
    matched$ambiguity_score[match(triage_df$match_id, matched$match_id)]
  } else {
    rep(0, nrow(triage_df))
  }
  amb_score[is.na(amb_score)] <- 0

  disagree_score <- ifelse(summary$tier == "Disagree", 1, 0)
  dist_score <- summary$match_distance / 100

  triage_df$triage_score <- amb_score * 2 + disagree_score + dist_score
  triage_df$tier <- summary$tier
  triage_df$material_a <- summary$material_a_raw
  triage_df$material_b <- summary$material_b_raw
  triage_df$family_a <- summary$family_a
  triage_df$family_b <- summary$family_b
  triage_df$match_distance <- summary$match_distance
  triage_df$ambiguity_score <- amb_score

  idx <- match(triage_df$match_id, matched$match_id)
  for (coord_col in c("ftir_x_um", "ftir_y_um", "raman_x_um", "raman_y_um",
                       "ftir_x_aligned", "ftir_y_aligned")) {
    if (coord_col %in% names(matched)) {
      triage_df[[coord_col]] <- matched[[coord_col]][idx]
    }
  }

  triage_df <- triage_df[order(-triage_df$triage_score), ]
  triage_df <- head(triage_df, n_top)

  if (nrow(triage_df) > 0) {
    write.csv(triage_df, file.path(out_dir, "triage_top_pairs.csv"), row.names = FALSE)
    log_message("  Wrote triage_top_pairs.csv (", nrow(triage_df), " pairs)")
  }
  invisible(NULL)
}


#' Export LDIR-specific results
export_ldir_results <- function(ldir_results, out_dir) {
  log_message("  Exporting LDIR results")

  if (!is.null(ldir_results$ldir_raman_match) &&
      nrow(ldir_results$ldir_raman_match$matched) > 0) {
    write.csv(ldir_results$ldir_raman_match$matched,
              file.path(out_dir, "ldir_raman_matched.csv"), row.names = FALSE)
    log_message("    LDIR-Raman matched: ", nrow(ldir_results$ldir_raman_match$matched))
  }

  if (!is.null(ldir_results$ldir_ftir_match) &&
      nrow(ldir_results$ldir_ftir_match$matched) > 0) {
    write.csv(ldir_results$ldir_ftir_match$matched,
              file.path(out_dir, "ldir_ftir_matched.csv"), row.names = FALSE)
    log_message("    LDIR-FTIR matched: ", nrow(ldir_results$ldir_ftir_match$matched))
  }

  if (!is.null(ldir_results$triplets) && nrow(ldir_results$triplets) > 0) {
    write.csv(ldir_results$triplets,
              file.path(out_dir, "triplets_3way.csv"), row.names = FALSE)
    log_message("    Three-way triplets: ", nrow(ldir_results$triplets))
  }

  if (!is.null(ldir_results$ldir_icp)) {
    lp <- ldir_results$ldir_icp$params
    writeLines(c(
      "# LDIR-to-Raman Transform Parameters",
      paste0("scale:        ", round(lp$scale, 6)),
      paste0("rotation_deg: ", round(lp$rotation_deg, 4)),
      paste0("tx:           ", round(lp$tx, 4)),
      paste0("ty:           ", round(lp$ty, 4)),
      paste0("reflected:    ", lp$reflected),
      paste0("icp_converged:  ", ldir_results$ldir_icp$converged),
      paste0("icp_final_rms:  ",
             round(tail(ldir_results$ldir_icp$rms_history, 1), 4), " um")
    ), file.path(out_dir, "ldir_transform_params.txt"))
  }

  if (!is.null(ldir_results$ldir_raman_agreement) &&
      nrow(ldir_results$ldir_raman_agreement$agreement_detail) > 0) {
    write.csv(ldir_results$ldir_raman_agreement$agreement_detail,
              file.path(out_dir, "ldir_raman_agreement.csv"), row.names = FALSE)
  }

  if (!is.null(ldir_results$ldir_clean)) {
    ldir_df <- data.frame(
      material = names(table(ldir_results$ldir_clean$material)),
      count    = as.integer(table(ldir_results$ldir_clean$material)),
      stringsAsFactors = FALSE
    )
    write.csv(ldir_df, file.path(out_dir, "ldir_material_distribution.csv"),
              row.names = FALSE)
  }

  if (!is.null(ldir_results$ldir_with_coords) &&
      "coord_match_cost" %in% names(ldir_results$ldir_with_coords)) {
    cols <- intersect(c("particle_id", "x_um", "y_um", "coord_match_cost",
                        "coord_source", "material", "quality", "feret_max_um"),
                      names(ldir_results$ldir_with_coords))
    write.csv(ldir_results$ldir_with_coords[, cols],
              file.path(out_dir, "ldir_coord_quality.csv"), row.names = FALSE)
  }

  invisible(NULL)
}
