# =============================================================================
# 10_export.R — Export results for downstream use
# =============================================================================

#' Export all pipeline outputs
#'
#' Writes:
#'   - matched particle table (CSV)
#'   - unmatched particle tables (CSV)
#'   - transform parameters (JSON-like text)
#'   - agreement summary (CSV)
#'   - diagnostic plots (PDF/PNG)
#'
#' @param match_result Result from match_particles()
#' @param agreement Result from analyze_agreement()
#' @param diagnostics List of ggplot objects from generate_diagnostics()
#' @param icp_result Result from icp_refine()
#' @param norm_result Result from normalize_coordinates()
#' @param config Configuration list
#' @return Invisible NULL
export_results <- function(match_result, agreement, diagnostics,
                           icp_result, norm_result, config) {
  out_dir <- config$output_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_message("Exporting results to: ", out_dir)

  # --- 1. Matched particle table ---
  if (nrow(match_result$matched) > 0) {
    write.csv(match_result$matched,
              file.path(out_dir, paste0("matched_particles_", timestamp, ".csv")),
              row.names = FALSE)
    log_message("  Wrote matched_particles.csv (", nrow(match_result$matched), " pairs)")
  }

  # --- 2. Unmatched particles ---
  if (nrow(match_result$unmatched_ftir) > 0) {
    write.csv(match_result$unmatched_ftir,
              file.path(out_dir, paste0("unmatched_ftir_", timestamp, ".csv")),
              row.names = FALSE)
  }
  if (nrow(match_result$unmatched_raman) > 0) {
    write.csv(match_result$unmatched_raman,
              file.path(out_dir, paste0("unmatched_raman_", timestamp, ".csv")),
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
           round(tail(icp_result$rms_history, 1), 4), " µm"),
    "",
    "# 3x3 Transform matrix (homogeneous, FTIR_norm -> Raman_norm)",
    paste0("matrix_row1: ", paste(round(icp_result$transform[1, ], 8), collapse = ", ")),
    paste0("matrix_row2: ", paste(round(icp_result$transform[2, ], 8), collapse = ", ")),
    paste0("matrix_row3: ", paste(round(icp_result$transform[3, ], 8), collapse = ", "))
  )
  writeLines(param_lines,
             file.path(out_dir, paste0("transform_params_", timestamp, ".txt")))
  log_message("  Wrote transform_params.txt")

  # --- 4. Agreement summary ---
  if (!is.null(agreement) && !is.null(agreement$agreement_detail) &&
      nrow(agreement$agreement_detail) > 0) {
    write.csv(agreement$agreement_detail,
              file.path(out_dir, paste0("agreement_summary_", timestamp, ".csv")),
              row.names = FALSE)

    write.csv(agreement$summary_df,
              file.path(out_dir, paste0("agreement_pairwise_", timestamp, ".csv")),
              row.names = FALSE)
    log_message("  Wrote agreement CSVs")
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
  writeLines(stats_lines,
             file.path(out_dir, paste0("match_statistics_", timestamp, ".txt")))
  log_message("  Wrote match_statistics.txt")

  # --- 6. Diagnostic plots ---
  if (length(diagnostics) > 0) {
    plot_dir <- file.path(out_dir, "plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    for (plot_name in names(diagnostics)) {
      tryCatch({
        ggplot2::ggsave(
          filename = file.path(plot_dir,
                               paste0(plot_name, "_", timestamp, ".png")),
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

    # Also save all plots as a single PDF
    tryCatch({
      pdf_path <- file.path(plot_dir, paste0("all_diagnostics_", timestamp, ".pdf"))
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
