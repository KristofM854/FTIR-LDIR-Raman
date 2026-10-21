# =============================================================================
# 10_export.R — Export results for downstream use
# =============================================================================

#' Export all pipeline outputs
#'
#' @param pairs Named list of align_and_match() results (FTIR-type source ->
#'   reference), each with $src_name (transform header, "FTIR_perkin") and
#'   optionally $scan_bounds; skipped pairs write their raw table as
#'   unmatched_<src>.csv
#' @param primary The run's primary pair, list(match, agreement, labels): the
#'   first of FTIR (PerkinElmer), FTIR (Bruker), LDIR -> reference. Its
#'   reference-side unmatched particles, match statistics, agreement summary
#'   and triage list are written; NULL when nothing was matched.
#' @param diagnostics List of ggplot objects from generate_diagnostics()
#' @param config Configuration list
#' @param ldir_results Optional list of LDIR pipeline results
#' @param aligned Named list (by file name: "raman", "ftir_perkin", ...) of
#'   every instrument's particles in the shared frame (x_aligned, y_aligned),
#'   written as 04_alignment/aligned_<name>.csv
#' @param input_paths Named list of input file paths used in this run
#'   (e.g. list(ftir = "...", raman = "...", ldir = "...", ldir_image = "..."))
#' @param ldir_image_info Result from canonicalize_ldir_image() or NULL
#' @param image_infos Named list of instrument image info objects
#' @return Invisible NULL
export_results <- function(pairs, primary, diagnostics, config,
                           ldir_results = NULL,
                           aligned = list(),
                           input_paths = list(),
                           ldir_image_info = NULL,
                           image_infos = list()) {
  out_dir <- config$output_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Create stage-based subdirectories
  dirs <- get_output_dirs(out_dir)

  log_message("Exporting results to: ", out_dir)

  # --- 0. Finalize manifest (mark stage as "export") ---
  run_id <- basename(out_dir)
  tryCatch({
    manifest_path <- file.path(dirs$manifest, "manifest.json")
    # Also check legacy location for update
    legacy_manifest <- file.path(out_dir, "manifest.json")
    if (file.exists(manifest_path)) {
      update_manifest_stage(out_dir, stage = "export_complete")
    } else if (file.exists(legacy_manifest)) {
      update_manifest_stage(out_dir, stage = "export_complete")
    } else {
      write_manifest(
        run_dir         = out_dir,
        run_id          = run_id,
        config          = config,
        input_paths     = input_paths,
        ldir_image_info = ldir_image_info,
        image_infos     = image_infos,
        stage           = "export_complete"
      )
    }
  }, error = function(e) {
    log_message("  Could not write/update manifest: ", e$message, level = "WARN")
  })

  # --- 1-3. FTIR-type pairs: matches, transform, agreement tables ---
  for (pair in pairs) export_pair(pair, dirs)

  # --- 4-6. The primary pair: statistics, agreement summary, triage ---
  if (!is.null(primary)) export_primary_pair(primary, dirs)

  # --- 7. LDIR results ---
  if (!is.null(ldir_results)) {
    export_ldir_results(ldir_results, dirs)
  }

  # Every instrument in the shared frame (the viewer draws these)
  for (nm in names(aligned)) {
    d <- aligned[[nm]]
    if (is.data.frame(d) && nrow(d) > 0 && all(c("x_aligned", "y_aligned") %in% names(d)))
      write.csv(d, file.path(dirs$alignment, paste0("aligned_", nm, ".csv")), row.names = FALSE)
  }

  # --- 8. Diagnostic plots ---
  if (length(diagnostics) > 0) {
    plot_dir <- dirs$plots

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


#' Write one FTIR-type -> reference pair: matched_<stem>.csv,
#' unmatched_<src>_vs_<ref>.csv, transform_params_<stem>.txt and the
#' agreement tables. A skipped pair (too few particles to align) writes its
#' table as unmatched_<src>.csv.
#'
#' @param pair align_and_match() result with $labels, $src_name and
#'   optionally $scan_bounds and $raw (the table, for a skipped pair)
#' @param dirs get_output_dirs() of the run
export_pair <- function(pair, dirs) {
  lb <- pair$labels
  src_file <- lb$src_file; ref_file <- lb$ref_file
  if (isTRUE(pair$skipped)) {
    if (is.data.frame(pair$raw) && nrow(pair$raw) > 0) {
      write.csv(pair$raw, file.path(dirs$matches, paste0("unmatched_", src_file, ".csv")),
                row.names = FALSE)
      log_message("  Wrote unmatched_", src_file, ".csv (", nrow(pair$raw), " particles)")
    }
    return(invisible(NULL))
  }
  m <- pair$match
  if (nrow(m$matched) > 0) {
    write.csv(m$matched, file.path(dirs$matches, paste0("matched_", lb$stem, ".csv")),
              row.names = FALSE)
    log_message("  Wrote matched_", lb$stem, ".csv (", nrow(m$matched), " pairs)")
  }
  if (nrow(m$unmatched_src) > 0) {
    f <- paste0("unmatched_", src_file, "_vs_", ref_file, ".csv")
    write.csv(m$unmatched_src, file.path(dirs$matches, f), row.names = FALSE)
    log_message("  Wrote ", f, " (", nrow(m$unmatched_src), " particles)")
  }
  write_transform_params(pair, file.path(dirs$alignment, paste0("transform_params_", lb$stem, ".txt")),
                         pair$src_name, scan_bounds = pair$scan_bounds,
                         ref_name = pair$ref_name %||% "Raman")
  if (!is.null(pair$agreement)) {
    write_agreement_tables(pair$agreement, dirs$agreement, lb$stem)
    log_message("  Wrote agreement tables (", lb$stem, ")")
  }
  invisible(NULL)
}


#' The primary pair's run-level outputs: the reference particles it left
#' unmatched (unmatched_<ref>_vs_<src>.csv), agreement_summary.csv,
#' agreement_by_category.csv, match_statistics.txt and triage_top_pairs.csv.
#'
#' @param primary list(match, agreement, labels)
#' @param dirs get_output_dirs() of the run
export_primary_pair <- function(primary, dirs) {
  match_result <- primary$match
  agreement <- primary$agreement
  lb <- primary$labels
  if (nrow(match_result$unmatched_ref) > 0)
    write.csv(match_result$unmatched_ref,
              file.path(dirs$matches, paste0("unmatched_", lb$ref_file, "_vs_", lb$src_file, ".csv")),
              row.names = FALSE)

  # Agreement summary (tiered)
  if (!is.null(agreement) && !is.null(agreement$agreement_detail) &&
      nrow(agreement$agreement_detail) > 0) {
    write.csv(agreement$agreement_detail,
              file.path(dirs$agreement, "agreement_summary.csv"),
              row.names = FALSE)
    log_message("  Wrote agreement CSVs")
  }

  # Category summary
  if (!is.null(agreement$category_summary) && nrow(agreement$category_summary) > 0) {
    write.csv(agreement$category_summary,
              file.path(dirs$agreement, "agreement_by_category.csv"),
              row.names = FALSE)
    log_message("  Wrote agreement_by_category.csv")
  }

  # Match statistics
  stats <- match_result$match_stats
  a <- lb$src_prefix; b <- lb$ref_prefix
  line <- function(key, value) paste0(sprintf("%-19s", paste0(key, ":")), value)
  stats_lines <- c(
    "# Match Statistics",
    paste0("# Generated: ", Sys.time()),
    paste0("# Pair: ", lb$stem),
    "",
    line(paste0("n_", a), stats$n_src),
    line(paste0("n_", b), stats$n_ref),
    line("n_matched", stats$n_matched),
    line(paste0("n_unmatched_", a), stats$n_unmatched_src),
    line(paste0("n_unmatched_", b), stats$n_unmatched_ref),
    line(paste0("match_rate_", a), paste0(round(stats$match_rate_src * 100, 1), "%")),
    line(paste0("match_rate_", b), paste0(round(stats$match_rate_ref * 100, 1), "%")),
    line("mean_distance_um", round(stats$mean_distance, 3)),
    line("median_distance_um", round(stats$median_distance, 3)),
    line("max_distance_um", round(stats$max_distance, 3)),
    line("chance_expected", round(stats$chance_expected %||% NA, 1)),
    line("chance_share", round(stats$chance_share %||% NA, 3))
  )

  # Tiered agreement rates
  if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
    tr <- agreement$tiered_rates
    stats_lines <- c(stats_lines, "",
      "# Agreement (rates over resolved pairs)",
      paste0("agreement_pairs_matched:        ", tr$n_pairs_matched),
      paste0("agreement_excluded_quality:     ", tr$n_excluded_quality_a + tr$n_excluded_quality_b),
      paste0("agreement_excluded_unidentified:", tr$n_excluded_unidentified),
      paste0("agreement_n_scored:             ", tr$n_total),
      paste0("agreement_n_unresolved:         ", tr$n_unresolved),
      paste0("agreement_n_resolved:           ", tr$n_resolved),
      paste0("agreement_family_pct_PRIMARY:   ", tr$family_or_better_pct, "%"),
      paste0("agreement_family_95ci:          ", tr$family_ci_low_pct, "-",
             tr$family_ci_high_pct, "%"),
      paste0("agreement_exact_pct_SECONDARY:  ", tr$exact_pct, "%"),
      paste0("agreement_exact_95ci:           ", tr$exact_ci_low_pct, "-",
             tr$exact_ci_high_pct, "%"),
      paste0("agreement_high_conf_n_resolved: ", tr$high_conf_n_resolved %||% NA,
             " (p_true >= ", tr$high_conf_min %||% NA, ")"),
      paste0("agreement_high_conf_family_pct: ", tr$high_conf_family_pct %||% NA, "%"),
      paste0("agreement_family_only_pct:      ", tr$family_pct, "%"),
      paste0("agreement_filler_pigment_pct:   ", tr$filler_pct, "%"),
      paste0("agreement_disagree_pct:         ", tr$disagree_pct, "%")
    )
  }

  writeLines(stats_lines,
             file.path(dirs$diagnostics, "match_statistics.txt"))
  log_message("  Wrote match_statistics.txt")

  # Triage export: top mismatched/high-risk pairs
  export_triage(.as_ftir_raman_columns(match_result, c(a, b)), agreement, dirs$diagnostics)
  invisible(NULL)
}


#' Write one instrument pair's agreement tables, the same for every pair:
#'   agreement_pairwise_<pair>.csv   one row per scored matched pair (tier,
#'                                   polymers, families, join warning, distance)
#'   agreement_by_material_<pair>.csv  counts and rates per material of A
#' (The LDIR and Bruker "pairwise" files used to hold the per-material table.)
write_agreement_tables <- function(agreement, out_dir, pair) {
  if (is.null(agreement)) return(invisible(NULL))
  if (!is.null(agreement$summary_df) && nrow(agreement$summary_df) > 0)
    write.csv(agreement$summary_df,
              file.path(out_dir, paste0("agreement_pairwise_", pair, ".csv")), row.names = FALSE)
  if (!is.null(agreement$agreement_detail) && nrow(agreement$agreement_detail) > 0)
    write.csv(agreement$agreement_detail,
              file.path(out_dir, paste0("agreement_by_material_", pair, ".csv")), row.names = FALSE)
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

  # A plastic-vs-filler/pigment pair is worth a look, but less than a plain
  # disagreement: both identifications can be right for the same particle.
  # Unresolved pairs need a human to look at the names, so they rank with
  # disagreements.
  disagree_score <- ifelse(summary$tier %in% c("Disagree", "Unresolved"), 1,
                           ifelse(summary$tier == "Filler/pigment", 0.5, 0))
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
#'
#' @param ldir_results List of LDIR pipeline results
#' @param dirs Named list of stage directories from get_output_dirs()
export_ldir_results <- function(ldir_results, dirs) {
  log_message("  Exporting LDIR results")
  # LDIR -> reference: ldir_raman, else ldir_ftir_perkin / ldir_ftir_bruker
  stem <- ldir_results$labels$stem %||% "ldir_raman"
  ref_file <- ldir_results$labels$ref_file %||% "raman"
  ref_name <- c(raman = "Raman", ftir_perkin = "FTIR (PerkinElmer)",
                ftir_bruker = "FTIR (Bruker)")[[ref_file]]
  lm <- ldir_results$ldir_ref_match

  if (!is.null(lm) && nrow(lm$matched) > 0) {
    write.csv(lm$matched,
              file.path(dirs$matches, paste0("matched_", stem, ".csv")), row.names = FALSE)
    log_message("    LDIR-", ref_name, " matched: ", nrow(lm$matched))
  }

  if (!is.null(ldir_results$ldir_ftir_match) &&
      nrow(ldir_results$ldir_ftir_match$matched) > 0) {
    write.csv(ldir_results$ldir_ftir_match$matched,
              file.path(dirs$matches, "matched_ldir_ftir_perkin.csv"), row.names = FALSE)
    log_message("    LDIR-FTIR matched: ", nrow(ldir_results$ldir_ftir_match$matched))
  }

  # Unmatched LDIR particles (from the LDIR-reference match)
  if (!is.null(lm$unmatched_src) && nrow(lm$unmatched_src) > 0) {
    write.csv(lm$unmatched_src,
              file.path(dirs$matches, paste0("unmatched_ldir_vs_", ref_file, ".csv")),
              row.names = FALSE)
    log_message("    Unmatched LDIR: ", nrow(lm$unmatched_src))
  }

  if (!is.null(ldir_results$join_check)) {
    write.csv(ldir_results$join_check$summary,
              file.path(dirs$diagnostics, "ldir_join_check.csv"), row.names = FALSE)
    write.csv(ldir_results$join_check$pairs,
              file.path(dirs$diagnostics, "ldir_join_pairs.csv"), row.names = FALSE)
    log_message("    LDIR join check: ldir_join_check.csv, ldir_join_pairs.csv")
  }
  if (!is.null(ldir_results$triplets) && nrow(ldir_results$triplets) > 0) {
    write.csv(ldir_results$triplets,
              file.path(dirs$matches, "triplets_3way.csv"), row.names = FALSE)
    log_message("    Three-way triplets: ", nrow(ldir_results$triplets))
  }

  if (!is.null(ldir_results$ldir_icp)) {
    lp <- ldir_results$ldir_icp$params
    writeLines(c(
      paste0("# LDIR-to-", ref_name, " Transform Parameters"),
      paste0("scale:        ", round(lp$scale, 6)),
      paste0("rotation_deg: ", round(lp$rotation_deg, 4)),
      paste0("tx:           ", round(lp$tx, 4)),
      paste0("ty:           ", round(lp$ty, 4)),
      paste0("reflected:    ", lp$reflected),
      paste0("icp_converged:  ", ldir_results$ldir_icp$converged),
      paste0("icp_final_rms:  ",
             round(tail(ldir_results$ldir_icp$rms_history, 1), 4), " um")
    ), file.path(dirs$alignment, paste0("transform_params_", stem, ".txt")))
  }

  write_agreement_tables(ldir_results$ldir_ref_agreement, dirs$agreement, stem)

  if (!is.null(ldir_results$ldir_clean)) {
    ldir_df <- data.frame(
      material = names(table(ldir_results$ldir_clean$material)),
      count    = as.integer(table(ldir_results$ldir_clean$material)),
      stringsAsFactors = FALSE
    )
    write.csv(ldir_df, file.path(dirs$diagnostics, "ldir_material_distribution.csv"),
              row.names = FALSE)
  }

  if (!is.null(ldir_results$ldir_with_coords) &&
      "coord_match_cost" %in% names(ldir_results$ldir_with_coords)) {
    cols <- intersect(c("particle_id", "x_um", "y_um", "coord_match_cost",
                        "coord_source", "material", "quality", "feret_max_um"),
                      names(ldir_results$ldir_with_coords))
    write.csv(ldir_results$ldir_with_coords[, cols],
              file.path(dirs$joined, "ldir_coord_quality.csv"), row.names = FALSE)
  }

  # Image-extracted coordinates (before Hungarian join) for Shiny viewer
  if (!is.null(ldir_results$ldir_image_extracted) &&
      nrow(ldir_results$ldir_image_extracted) > 0) {
    write.csv(ldir_results$ldir_image_extracted,
              file.path(dirs$joined, "ldir_image_extracted.csv"), row.names = FALSE)
    log_message("    LDIR image-extracted particles: ",
                nrow(ldir_results$ldir_image_extracted))
  }

  invisible(NULL)
}
