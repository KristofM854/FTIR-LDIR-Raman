# =============================================================================
# run_pipeline.R — the whole pipeline as one function of its config
# =============================================================================
#
#   config <- make_config(ftir_path = ..., raman_path = ..., output_dir = "output")
#   config$ldir_path <- ...            # optional instruments and images
#   result <- run_pipeline(config)
#
# Stages: run folder + provenance; ingest (read_instrument, R/01_readers.R);
# image placement; Raman reference frame and FTIR (PerkinElmer / Bruker)
# alignment + matching + agreement (R/align_and_match.R); LDIR
# (R/ldir_stage.R); cross-instrument matches; diagnostics; export; report.
# Everything is written to a new run folder under config$output_dir.
# =============================================================================


#' Run the pipeline.
#'
#' @param config Pipeline config (make_config() plus input paths:
#'   ftir_path, raman_path, and optionally ftir_bruker_path, ldir_path,
#'   ftir_image, ftir_bruker_image, raman_image, ldir_image)
#' @param report Also write the PDF / HTML report
#' @return Invisibly, a list: config (as run, incl. output_dir), raman (the
#'   reference frame), ftir, ftir_bruker (align_and_match() results), ldir
#'   (LDIR stage results), cross (Bruker<->PerkinElmer, Bruker<->LDIR matches),
#'   composites, tps_assessment
run_pipeline <- function(config, report = TRUE) {

  # Spectral-library index for material classification (R/08b_material_map.R).
  # Load it now so the log states up front whether trade names / synonyms from
  # the licensed library listings are being resolved in this run.
  if (!is.null(config$spectral_library_dir))
    options(ftir.spectral_library_dir = config$spectral_library_dir)
  if (is.null(get_spectral_library_index()))
    log_message("Spectral-library index: not available -- material families ",
                "come from the built-in name rules only")

  # Create a timestamped run subfolder (output/YYYY-MM-DD_1, _2, ...)
  config$output_dir <- make_run_dir(config$output_dir)

  # --- Collect input file paths for provenance ---
  .run_input_paths <- list(
    ftir       = if (!is.null(config$ftir_path)  && nzchar(config$ftir_path))  config$ftir_path  else NULL,
    raman      = if (!is.null(config$raman_path) && nzchar(config$raman_path)) config$raman_path else NULL,
    ldir       = if (!is.null(config$ldir_path)  && nzchar(config$ldir_path))  config$ldir_path  else NULL,
    ldir_image  = if (!is.null(config$ldir_image)  && nzchar(config$ldir_image))  config$ldir_image  else NULL,
    ftir_image  = if (!is.null(config$ftir_image)  && nzchar(config$ftir_image))  config$ftir_image  else NULL,
    raman_image = if (!is.null(config$raman_image) && nzchar(config$raman_image)) config$raman_image else NULL,
    ftir_bruker       = if (!is.null(config$ftir_bruker_path)  && nzchar(config$ftir_bruker_path))  config$ftir_bruker_path  else NULL,
    ftir_bruker_image = if (!is.null(config$ftir_bruker_image) && nzchar(config$ftir_bruker_image)) config$ftir_bruker_image else NULL
  )
  .run_input_paths <- .run_input_paths[!vapply(.run_input_paths, is.null, logical(1))]

  # --- Canonicalize instrument images and create previews in inputs/ ---
  .ftir_image_info <- NULL
  .raman_image_info <- NULL
  .ldir_image_info <- NULL
  inputs_dir <- file.path(config$output_dir, "inputs")

  if (!is.null(config$ftir_image) && nzchar(config$ftir_image) && file.exists(config$ftir_image)) {
    .ftir_image_info <- canonicalize_instrument_image(
      src_path = config$ftir_image,
      inputs_dir = inputs_dir,
      instrument = "ftir",
      max_preview_px = 2000L
    )
  }
  if (!is.null(config$raman_image) && nzchar(config$raman_image) && file.exists(config$raman_image)) {
    .raman_image_info <- canonicalize_instrument_image(
      src_path = config$raman_image,
      inputs_dir = inputs_dir,
      instrument = "raman",
      max_preview_px = 2000L
    )
  }
  if (!is.null(config$ldir_image) && nzchar(config$ldir_image) && file.exists(config$ldir_image)) {
    .ldir_image_info <- canonicalize_instrument_image(
      src_path = config$ldir_image,
      inputs_dir = inputs_dir,
      instrument = "ldir",
      max_preview_px = 2000L
    )
    if (!is.null(.ldir_image_info$canonical_path) && file.exists(.ldir_image_info$canonical_path)) {
      log_message("Using canonical PNG for LDIR processing: ", .ldir_image_info$canonical_path)
      config$ldir_image_canonical <- .ldir_image_info$canonical_path
      config$ldir_image_preview   <- .ldir_image_info$preview_path
    }
  }

  # --- Write run manifest (provenance) as early as possible ---
  tryCatch({
    write_manifest(
      run_dir         = config$output_dir,
      run_id          = basename(config$output_dir),
      config          = config,
      input_paths     = .run_input_paths,
      ldir_image_info = .ldir_image_info,
      image_infos     = list(
        ftir_image = .ftir_image_info,
        raman_image = .raman_image_info,
        ldir_image = .ldir_image_info
      ),
      stage           = "started"
    )
  }, error = function(e) {
    log_message("  Could not write manifest: ", e$message, level = "WARN")
  })
  # Full parameter set + R session, so a run can be reproduced from its folder.
  write_run_environment(config$output_dir, config)

  # --- Debug mode setup (Step 0: bulletproof) ---
  # config$debug = TRUE (or config_overrides = list(debug = TRUE) for main.R)
  # enables debug artifacts.
  # Writes to an absolute path so nothing can silently swallow errors.
  if (isTRUE(config$debug)) {
    run_id <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    config$run_id <- run_id

    # Build absolute path -- avoids any working-directory ambiguity
    debug_dir_abs <- normalizePath(
      file.path(config$output_dir, "debug"),
      winslash = "/", mustWork = FALSE
    )
    dir.create(debug_dir_abs, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(debug_dir_abs)) {
      stop("DEBUG: failed to create debug directory: ", debug_dir_abs)
    }

    # Heartbeat file -- if this is absent the whole debug run failed
    heartbeat <- file.path(debug_dir_abs, "DEBUG_ALIVE.txt")
    writeLines(c(paste0("run_id: ", run_id),
                 paste0("created: ", Sys.time()),
                 paste0("debug_dir: ", debug_dir_abs)),
               heartbeat)
    if (!file.exists(heartbeat)) {
      stop("DEBUG: heartbeat write failed \u2014 check permissions for ", debug_dir_abs)
    }

    config$debug_dir <- debug_dir_abs
    message("DEBUG DIR: ", debug_dir_abs)
    log_message("Debug mode ON \u2014 artifacts in: ", debug_dir_abs)
  }

  # ---------------------------------------------------------------------------
  # 2. Data ingestion
  # ---------------------------------------------------------------------------

  log_message(strrep("=", 60))
  log_message("Multi-Instrument Particle Matching Pipeline")
  log_message(strrep("=", 60))

  # Create stage-based output directories for audit exports
  .out_dirs <- get_output_dirs(config$output_dir)
  # Canonical image objects of this run (R/particle_contract.R), written to
  # 00_manifest/images.json at the end.
  .image_objects <- list()

  # --- FTIR ---
  ftir_raw  <- read_instrument("ftir", config)
  write.csv(ftir_raw, file.path(.out_dirs$ingested, "ftir_perkin_ingested.csv"), row.names = FALSE)

  # --- Raman ---
  raman_raw <- read_instrument("raman", config)
  write.csv(raman_raw, file.path(.out_dirs$ingested, "raman_ingested.csv"), row.names = FALSE)

  # --- Auto-calibrate Raman image placement ----------------------------------
  # Measure where the exported Raman micrograph actually sits under the
  # particles (the export may be a cropped/zoomed view whose footprint differs
  # from the WITec panel Width/Height).  On a confident, non-mirrored fit, patch
  # the run manifest's config_snapshot so the viewer places the image exactly --
  # no per-run tools/diagnose_raman_placement.R --apply needed.  The config
  # values stay as the scale-search seed; a weak fit leaves them untouched.
  # NOTE: this deliberately does NOT require config$raman_image_width_um /
  # _height_um. Those are per-dataset values the operator copies out of WITec,
  # and they are only the scale-search SEED -- the search spans 0.3-2.4x it and
  # derives its own seed from the particle bbox when they are absent. Gating on
  # them meant the auto-calibration skipped exactly the runs that needed it:
  # a new dataset with no WITec values entered fell straight through to
  # particle-bbox placement, which is visibly offset on sparse scans.
  if (!is.null(.raman_image_info)) {
    .raman_cal_img <- if (!is.null(.raman_image_info$canonical_path) &&
                          file.exists(.raman_image_info$canonical_path))
      .raman_image_info$canonical_path else config$raman_image
    .raman_fit <- tryCatch(
      measure_raman_image_placement(
        .raman_cal_img, raman_raw$x_um, raman_raw$y_um,
        config$raman_image_width_um, config$raman_image_height_um),
      error = function(e) { log_message("  Raman placement auto-measure error: ",
                                        conditionMessage(e), level = "WARN"); NULL })
    if (!is.null(.raman_fit) && !isTRUE(.raman_fit$mirrored)) {
      log_message(sprintf(paste0(
        "Raman image auto-calibration: extent %.0f x %.0f um at center (%.0f, %.0f), ",
        "scale %.2f x WITec, %.0f%% of particles on image ",
        "(random baseline %.0f%%, lift %.2f)."),
        .raman_fit$width_um, .raman_fit$height_um,
        .raman_fit$center_x_um, .raman_fit$center_y_um,
        .raman_fit$scale, 100 * .raman_fit$frac_bright,
        100 * (.raman_fit$baseline %||% NA_real_), .raman_fit$lift %||% NA_real_))
      config$raman_image_width_um    <- .raman_fit$width_um
      config$raman_image_height_um   <- .raman_fit$height_um
      config$raman_image_center_x_um <- .raman_fit$center_x_um
      config$raman_image_center_y_um <- .raman_fit$center_y_um
      update_manifest_config_snapshot(config$output_dir, list(
        raman_image_width_um    = .raman_fit$width_um,
        raman_image_height_um   = .raman_fit$height_um,
        raman_image_center_x_um = .raman_fit$center_x_um,
        raman_image_center_y_um = .raman_fit$center_y_um))
    } else {
      # Worth a WARN, not an INFO: the configured raman_image_* values are
      # PER-DATASET. If they were left over from another scan the viewer either
      # rejects them (containment check) and falls back to particle-bbox
      # placement, or -- worse -- accepts them and draws the image in the wrong
      # place. Either way the user needs to know the measurement did not land.
      log_message("Raman image auto-calibration: no confident non-mirrored fit; ",
                  "keeping configured raman_image_* values. If the Raman tab ",
                  "shows the image offset from the particles, re-measure with ",
                  "tools/diagnose_raman_placement.R --apply.", level = "WARN")
    }
  }

  # Raman image object: the extent the viewer uses (auto-calibrated above, else
  # the configured WITec values). Without an extent the viewer falls back to the
  # particle bounding box, which is not recorded as an image object.
  if (!is.null(.raman_image_info) &&
      all(!vapply(config[c("raman_image_width_um", "raman_image_height_um",
                           "raman_image_center_x_um", "raman_image_center_y_um")],
                  is.null, logical(1)))) {
    .w <- config$raman_image_width_um;    .h <- config$raman_image_height_um
    .cx <- config$raman_image_center_x_um; .cy <- config$raman_image_center_y_um
    .image_objects$raman <- tryCatch(image_object_from_extent(
      "raman", config$raman_image,
      .raman_image_info$orig_width, .raman_image_info$orig_height,
      c(xmin = .cx - .w / 2, xmax = .cx + .w / 2, ymin = .cy - .h / 2, ymax = .cy + .h / 2),
      provenance = if (!is.null(.raman_fit) && !isTRUE(.raman_fit$mirrored))
        "registration" else "config"),
      error = function(e) NULL)
  }

  # --- LDIR (optional) ---
  ldir_raw <- NULL
  has_ldir <- !is.null(config$ldir_path) && nzchar(config$ldir_path)
  if (has_ldir) {
    # Low-HQI identifications are relabelled "unknown" (kept, not dropped) to
    # mirror the LDIR software; see read_instrument().
    ldir_raw <- read_instrument("ldir", config)
    # Auto-load coordinate-swap corrections from the CSV sidecar next to the LDIR
    # Excel (<stem>_coord_swaps.csv), unless swaps were set explicitly in config.
    # [[ ]] and list(): `config$ldir_coord_swaps <- NULL` would DELETE the
    # element, after which config$ldir_coord_swaps partial-matches
    # ldir_coord_swaps_suffix and returns "_coord_swaps" as a bogus swap.
    if (is.null(config[["ldir_coord_swaps"]]))
      config["ldir_coord_swaps"] <- list(load_ldir_coord_swaps(config$ldir_path, config))
    write.csv(ldir_raw, file.path(.out_dirs$ingested, "ldir_ingested.csv"), row.names = FALSE)
    log_message("LDIR data loaded: ", nrow(ldir_raw), " particles",
                if (all(is.na(ldir_raw$x_um))) " (no coordinates)" else "")
  } else {
    log_message("LDIR: not provided \u2014 skipping LDIR analysis")
  }

  # --- FTIR Bruker (optional) ---
  ftir_bruker_raw <- NULL
  has_ftir_bruker <- !is.null(config$ftir_bruker_path) && nzchar(config$ftir_bruker_path)
  if (has_ftir_bruker) {
    ftir_bruker_raw <- read_instrument("ftir_bruker", config)
    write.csv(ftir_bruker_raw, file.path(.out_dirs$ingested, "ftir_bruker_ingested.csv"), row.names = FALSE)
    log_message("FTIR (Bruker) data loaded: ", nrow(ftir_bruker_raw), " particles")
  } else {
    log_message("FTIR (Bruker): not provided \u2014 skipping")
  }

  # --- FTIR image placement (for background display only) ---
  # FTIR particle um coordinates come directly from the data file; the FTIR
  # image is only drawn behind them. The PNG carries no geometry metadata, so
  # place_ftir_image() (R/ftir_image_placement.R) takes the extent from the
  # config (ftir_image_width_um / _height_um / _center_*), else registers the
  # particles visible in the image against ALL ingested particles, else falls
  # back to an approximate particle-extent fit. The result is written to the run
  # directory once and reused by the viewer and the report, so no plot derives
  # the image extent from whatever (filtered) particles it happens to draw.
  ftir_scan_bounds <- NULL
  .placement_provenance <- c(config_extent = "config", registration = "registration",
                             particle_extent_approx = "particle_extent_approx")
  .placement_image_object <- function(pl, instrument, path) {
    if (is.null(pl)) return(NULL)
    # In the file's own pixels: a downsampled raster covers the first
    # fx * width (fy * height) pixels of the file (see .raster_coverage).
    cov <- .raster_coverage(pl$raster)
    tryCatch(image_object_from_extent(
      instrument, path, cov$fx * cov$orig_w, cov$fy * cov$orig_h,
      c(xmin = pl$xmin, xmax = pl$xmax, ymin = pl$ymin, ymax = pl$ymax),
      provenance = .placement_provenance[[pl$method]],
      # A row-flipped placement: the file's first row is the lowest y.
      y_axis = if (isTRUE(pl$flip_rows)) "down" else "up"),
      error = function(e) {
        log_message("  ", instrument, " image object not recorded: ", conditionMessage(e),
                    level = "WARN"); NULL })
  }
  .place_ftir_run_image <- function(img_path, df, prefix, key, label) {
    if (is.null(img_path) || !nzchar(img_path) || is.null(df) || nrow(df) == 0)
      return(NULL)
    img_raw <- read_image_any(img_path, verbose = FALSE)
    if (is.null(img_raw)) {
      log_message("  WARNING: could not read ", label, " image \u2014 placement unavailable",
                  level = "WARN")
      return(NULL)
    }
    ok <- is.finite(df$x_um) & is.finite(df$y_um)
    pl <- place_ftir_image(img_raw, df$x_um[ok], df$y_um[ok], df$feret_max_um[ok],
                           cfg = config, prefix = prefix)
    if (is.null(pl)) return(NULL)
    write_ftir_image_placement(pl, file.path(config$output_dir, FTIR_PLACEMENT_FILES[[key]]))
    reg <- pl$registration
    log_message(sprintf("  %s image %d x %d px placed by %s: X [%.1f, %.1f] Y [%.1f, %.1f] um",
                        label, ncol(img_raw), nrow(img_raw), pl$method,
                        pl$xmin, pl$xmax, pl$ymin, pl$ymax))
    if (!is.null(reg))
      log_message(sprintf("    registration: %d of %d particles matched to %d image blobs, RMS %.1f um (%.2f px), pixel %.3f x %.3f um",
                          reg$n_matched, reg$n_particles, reg$n_blobs, reg$rms_um, reg$rms_px,
                          reg$bx, abs(reg$by)))
    if (identical(pl$method, "particle_extent_approx"))
      log_message("    APPROXIMATE: no config extent and the image could not be registered ",
                  "to the particles; set ", prefix, "_width_um / _height_um in the config ",
                  "for an exact overlay", level = "WARN")
    .image_objects[[key]] <<- .placement_image_object(pl, key, img_path)
    pl
  }
  if (!is.null(config$ftir_image) && nzchar(config$ftir_image)) {
    log_message(strrep("-", 50))
    log_message("FTIR image: placing the viewer background")
    .ftir_pl <- .place_ftir_run_image(config$ftir_image, ftir_raw, "ftir_image",
                                      "ftir", "FTIR (PerkinElmer)")
    if (!is.null(.ftir_pl))
      ftir_scan_bounds <- list(x_min = .ftir_pl$xmin, x_max = .ftir_pl$xmax,
                               y_min = .ftir_pl$ymin, y_max = .ftir_pl$ymax)
  }
  if (!is.null(config$ftir_bruker_image) && nzchar(config$ftir_bruker_image))
    invisible(.place_ftir_run_image(config$ftir_bruker_image, ftir_bruker_raw,
                                    "ftir_bruker_image", "ftir_bruker", "FTIR (Bruker)"))

  # ---------------------------------------------------------------------------
  # 3-9. Raman reference frame; FTIR (PerkinElmer) alignment, matching, agreement
  #
  # Raman is the reference frame. Every spatial step runs on the FULL particle
  # sets; quality filters apply only in the agreement analysis.
  # ---------------------------------------------------------------------------

  ref <- prepare_raman_reference(raman_raw, config)
  write.csv(ref$prefiltered, file.path(.out_dirs$prefiltered, "raman_prefiltered.csv"), row.names = FALSE)

  ftir_pair <- align_and_match(ftir_raw, ref, config, "FTIR (PerkinElmer)")
  write.csv(ftir_pair$prefiltered, file.path(.out_dirs$prefiltered, "ftir_perkin_prefiltered.csv"), row.names = FALSE)
  if (isTRUE(ftir_pair$skipped))
    stop("FTIR (PerkinElmer): too few particles to align to Raman.", call. = FALSE)

  raman_clean      <- ref$clean
  raman_for_match  <- ref$for_match
  match_result     <- ftir_pair$match
  agreement        <- ftir_pair$agreement
  icp_result       <- ftir_pair$icp
  norm_result      <- ftir_pair$norm
  ftir_aligned     <- ftir_pair$aligned
  ftir_aligned_all <- ftir_pair$aligned_all


  # ---------------------------------------------------------------------------
  # 10. TPS assessment (systematic distortion check)
  # ---------------------------------------------------------------------------

  tps_assessment <- assess_tps_need(match_result)
  log_message("TPS assessment: ", tps_assessment$message)
  if (isTRUE(tps_assessment$recommend_tps)) {
    log_message("  RECOMMENDATION: Consider thin-plate spline warp for improved alignment",
                level = "WARN")
  }

  # ---------------------------------------------------------------------------
  # 11. Composite matching (1:many) for unmatched FTIR
  # ---------------------------------------------------------------------------

  composites <- data.frame()
  if (nrow(match_result$unmatched_ftir) > 0 && nrow(match_result$unmatched_raman) > 0) {
    # Ensure unmatched FTIR particles have aligned coordinates
    unmatched_ftir_src <- match_result$unmatched_ftir
    if (!"x_aligned" %in% names(unmatched_ftir_src)) {
      unmatched_ftir_src <- apply_normalization(unmatched_ftir_src,
                                                norm_result$ftir_centroid,
                                                norm_result$ftir_scale)
      tf <- apply_transform_points(unmatched_ftir_src$x_norm,
                                   unmatched_ftir_src$y_norm,
                                   icp_result$transform)
      unmatched_ftir_src$x_aligned <- tf$x_transformed
      unmatched_ftir_src$y_aligned <- tf$y_transformed
    }

    composites <- find_composite_matches(
      unmatched_ftir_src, match_result$unmatched_raman, config
    )

    if (nrow(composites) > 0) {
      log_message("Composite matches found: ", nrow(composites),
                  " FTIR particles matched to multiple Raman fragments")
    }
  }

  # ---------------------------------------------------------------------------
  # 5b. FTIR (Bruker): the same alignment, matching and agreement
  # ---------------------------------------------------------------------------

  bruker_pair <- NULL
  if (has_ftir_bruker && !is.null(ftir_bruker_raw) && nrow(ftir_bruker_raw) > 0) {
    bruker_pair <- align_and_match(ftir_bruker_raw, ref, config, "FTIR (Bruker)")
    write.csv(bruker_pair$prefiltered,
              file.path(.out_dirs$prefiltered, "ftir_bruker_prefiltered.csv"),
              row.names = FALSE)
  }
  bruker_aligned <- if (!isTRUE(bruker_pair$skipped)) bruker_pair$aligned

  # ---------------------------------------------------------------------------
  # 12. LDIR (image -> coordinates -> alignment -> matching), R/ldir_stage.R
  # ---------------------------------------------------------------------------

  ldir_results <- NULL
  if (has_ldir && !is.null(ldir_raw)) {
    .ldir_stage <- run_ldir_stage(ldir_raw, ref, ftir_pair, config, .out_dirs)
    ldir_results <- .ldir_stage$results
    .image_objects$ldir <- .ldir_stage$image
  }

  # ---------------------------------------------------------------------------
  # 12b. Cross-instrument pairwise matches (non-Raman pairs)
  # ---------------------------------------------------------------------------
  # Every instrument is aligned into the Raman frame, so any pair can be matched
  # directly. Perkin<->LDIR is matched in the LDIR stage.
  bruker_perkin_match <- match_aligned_pair(bruker_aligned, ftir_aligned, config,
                                            "ftir_bruker", "ftir_perkin")
  if (!is.null(bruker_perkin_match))
    log_message("FTIR (Bruker) <-> FTIR (PerkinElmer): ",
                bruker_perkin_match$match_stats$n_matched, " matched")
  bruker_ldir_match <- match_aligned_pair(bruker_aligned, ldir_results$ldir_aligned, config,
                                          "ftir_bruker", "ldir")
  if (!is.null(bruker_ldir_match))
    log_message("FTIR (Bruker) <-> LDIR: ",
                bruker_ldir_match$match_stats$n_matched, " matched")


  # ---------------------------------------------------------------------------
  # 13. Diagnostics (use full datasets for overlay, matched pairs for detail)
  # ---------------------------------------------------------------------------

  diagnostics <- generate_diagnostics(
    ftir_aligned_all, raman_clean,
    match_result, icp_result, agreement,
    config
  )

  # Add LDIR diagnostics if spatial matching was performed
  if (!is.null(ldir_results) && !is.null(ldir_results$ldir_aligned) &&
      !is.null(ldir_results$ldir_raman_match)) {

    # Step 4: Assert LDIR overlay uses aligned coordinates
    ldir_plot_df <- ldir_results$ldir_aligned
    stopifnot(all(c("x_aligned", "y_aligned") %in% colnames(ldir_plot_df)))

    # Step 5: dump traced particles + full snapshot for final overlay
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_plot_df, trace_ids, "final_overlay", config$debug_dir)
      trace_particle_snapshot(ldir_plot_df, "final_overlay", config)

      # Log the exact plotted coordinates for traced particles
      for (pid in trace_ids) {
        ri <- which(ldir_plot_df$particle_id == pid)
        if (length(ri) > 0) {
          message(sprintf("OVERLAY CHECK '%s': x_aligned=%.2f  y_aligned=%.2f",
                          pid, ldir_plot_df$x_aligned[ri[1]],
                          ldir_plot_df$y_aligned[ri[1]]))
        } else {
          message(sprintf("OVERLAY CHECK '%s': NOT IN ldir_plot_df", pid))
        }
      }
    }

    ldir_diag <- generate_ldir_diagnostics(
      ldir_plot_df, raman_clean,
      ldir_results$ldir_raman_match, ldir_results$ldir_raman_agreement,
      debug_subtitle = if (isTRUE(config$debug)) {
        # Build per-particle subtitle for A3 etc.
        parts <- character(0)
        for (pid in (config$debug_trace_ids %||% c("A3", "MP_11"))) {
          ri <- which(ldir_plot_df$particle_id == pid)
          if (length(ri) > 0) {
            parts <- c(parts, sprintf("%s: (%.1f, %.1f)",
                                      pid,
                                      ldir_plot_df$x_aligned[ri[1]],
                                      ldir_plot_df$y_aligned[ri[1]]))
          }
        }
        if (length(parts) > 0) paste(parts, collapse = "  ") else NULL
      } else NULL
    )
    diagnostics <- c(diagnostics, ldir_diag)
    log_message("  Added ", length(ldir_diag), " LDIR diagnostic plots")
  }

  # ---------------------------------------------------------------------------
  # 14. Export
  # ---------------------------------------------------------------------------

  export_results(
    match_result, agreement, diagnostics,
    icp_result, norm_result, config,
    ftir_scan_bounds = ftir_scan_bounds,
    ldir_results     = ldir_results,
    input_paths      = .run_input_paths,
    ldir_image_info  = .ldir_image_info,
    image_infos      = list(
      ftir_image = .ftir_image_info,
      raman_image = .raman_image_info,
      ldir_image = .ldir_image_info
    )
  )

  # Use stage-based output directories for remaining exports
  .export_dirs <- get_output_dirs(config$output_dir)

  # Export FTIR Bruker results
  if (!is.null(bruker_pair) && !isTRUE(bruker_pair$skipped)) {
    bm <- bruker_pair$match
    if (nrow(bm$matched) > 0) {
      write.csv(bm$matched, file.path(.export_dirs$matches, "matched_ftir_bruker_raman.csv"),
                row.names = FALSE)
      log_message("  Wrote matched_ftir_bruker_raman.csv (", nrow(bm$matched), " pairs)")
    }
    if (nrow(bm$unmatched_ftir) > 0) {
      write.csv(bm$unmatched_ftir,
                file.path(.export_dirs$matches, "unmatched_ftir_bruker_vs_raman.csv"),
                row.names = FALSE)
      log_message("  Wrote unmatched_ftir_bruker_vs_raman.csv (",
                  nrow(bm$unmatched_ftir), " particles)")
    }
    write_transform_params(bruker_pair,
                           file.path(.export_dirs$alignment, "transform_params_ftir_bruker_raman.txt"),
                           "FTIR_bruker")
    if (!is.null(bruker_pair$agreement)) {
      write_agreement_tables(bruker_pair$agreement, .export_dirs$agreement, "ftir_bruker_raman")
      log_message("  Wrote agreement tables (ftir_bruker_raman)")
    }
  } else if (has_ftir_bruker && !is.null(ftir_bruker_raw) && nrow(ftir_bruker_raw) > 0) {
    # Alignment skipped (too few particles)
    write.csv(ftir_bruker_raw,
              file.path(.export_dirs$matches, "unmatched_ftir_bruker.csv"),
              row.names = FALSE)
    log_message("  Wrote unmatched_ftir_bruker.csv (", nrow(ftir_bruker_raw), " particles)")
  }

  # Export cross-instrument matches (Bruker<->Perkin, Bruker<->LDIR). Perkin<->LDIR
  # is written by export_ldir_results() as matched_ldir_ftir_perkin.csv.
  if (!is.null(bruker_perkin_match) && nrow(bruker_perkin_match$matched) > 0) {
    write.csv(bruker_perkin_match$matched,
              file.path(.export_dirs$matches, "matched_ftir_bruker_ftir_perkin.csv"),
              row.names = FALSE)
    log_message("  Wrote matched_ftir_bruker_ftir_perkin.csv (",
                nrow(bruker_perkin_match$matched), " pairs)")
  }
  if (!is.null(bruker_ldir_match) && nrow(bruker_ldir_match$matched) > 0) {
    write.csv(bruker_ldir_match$matched,
              file.path(.export_dirs$matches, "matched_ftir_bruker_ldir.csv"),
              row.names = FALSE)
    log_message("  Wrote matched_ftir_bruker_ldir.csv (",
                nrow(bruker_ldir_match$matched), " pairs)")
  }

  # Export composite matches if found
  if (nrow(composites) > 0) {
    write.csv(composites,
              file.path(.export_dirs$matches, "composite_matches.csv"),
              row.names = FALSE)
    log_message("  Wrote composite_matches.csv (", nrow(composites), " composites)")
  }

  # Export TPS assessment
  if (!is.null(tps_assessment$quadrant_residuals)) {
    write.csv(tps_assessment$quadrant_residuals,
              file.path(.export_dirs$diagnostics, "tps_quadrant_residuals.csv"),
              row.names = FALSE)
    tps_lines <- c(
      "# TPS (Thin-Plate Spline) Assessment",
      paste0("# Generated: ", Sys.time()),
      "",
      paste0("systematic_distortion: ", tps_assessment$systematic_distortion),
      paste0("recommend_tps:         ", tps_assessment$recommend_tps),
      paste0("dx_range_um:           ", tps_assessment$dx_range_um),
      paste0("dy_range_um:           ", tps_assessment$dy_range_um),
      paste0("assessment:            ", tps_assessment$message)
    )
    writeLines(tps_lines, file.path(.export_dirs$diagnostics, "tps_assessment.txt"))
  }

  # ---------------------------------------------------------------------------
  # 15. Summary
  # ---------------------------------------------------------------------------

  log_message(strrep("=", 60))
  log_message("Pipeline complete!")
  log_message(strrep("=", 60))
  log_message("  Alignment method:  ", ftir_pair$method,
              if (ftir_pair$method == "landmark")
                paste0(" (", ftir_pair$landmark$n_inliers, " landmark inliers, ",
                       round(ftir_pair$landmark$mean_residual, 1), " \u00b5m mean residual)")
              else "")
  log_message("  FTIR landmarks:    ", ftir_pair$landmark$n_ftir_landmarks,
              " (confident: ", ftir_pair$landmark$confident, ")")
  log_message("  FTIR plastic anchors (alignment): ", nrow(ftir_pair$for_align))
  log_message("  Raman targets (alignment):       ", nrow(ref$for_align))
  log_message("  FTIR particles (matching):       ", nrow(ftir_pair$for_match))
  log_message("  Raman particles (matching):      ", nrow(raman_for_match))
  log_message("  Matched pairs:      ", match_result$match_stats$n_matched)
  log_message("  Unmatched FTIR:     ", match_result$match_stats$n_unmatched_ftir)
  log_message("  Unmatched Raman:    ", match_result$match_stats$n_unmatched_raman)
  log_message("  FTIR match rate:    ",
              round(match_result$match_stats$match_rate_ftir * 100, 1), "%")
  if (!is.na(agreement$agreement_rate)) {
    log_message("  Material agreement: ", round(agreement$agreement_rate * 100, 1), "%")
  }
  if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
    tr <- agreement$tiered_rates
    log_message("  Agreement:         Family ", tr$family_or_better_pct, "% (primary), Exact ",
                tr$exact_pct, "%, Filler/pigment ", tr$filler_pct, "% of ", tr$n_resolved,
                " resolved pairs; ", tr$n_unresolved, " unresolved")
  }
  if (nrow(composites) > 0) {
    log_message("  Composite matches: ", nrow(composites))
  }
  log_message("  TPS assessment:    ", tps_assessment$message)
  if (has_ldir && !is.null(ldir_results)) {
    log_message("  --- LDIR ---")
    log_message("  LDIR particles:    ", nrow(ldir_results$ldir_clean))
    log_message("  LDIR coordinates:  ",
                if (ldir_results$has_coords) "spatial (from image)" else "none")
    if (!is.null(ldir_results$ldir_landmark)) {
      log_message("  LDIR alignment:    ", ldir_results$ldir_alignment_method,
                  " (landmarks: ", ldir_results$ldir_landmark$n_ftir_landmarks,
                  ", confident: ", ldir_results$ldir_landmark$confident, ")")
    }
    if (!is.null(ldir_results$ldir_raman_match)) {
      log_message("  LDIR-Raman matched: ",
                  ldir_results$ldir_raman_match$match_stats$n_matched)
    }
    if (!is.null(ldir_results$ldir_ftir_match)) {
      log_message("  LDIR-FTIR matched:  ",
                  ldir_results$ldir_ftir_match$match_stats$n_matched)
    }
    if (nrow(ldir_results$triplets) > 0) {
      log_message("  Three-way triplets: ", nrow(ldir_results$triplets))
    }
    if (!is.null(ldir_results$ldir_raman_agreement) &&
        !is.null(ldir_results$ldir_raman_agreement$tiered_rates) &&
        ldir_results$ldir_raman_agreement$tiered_rates$n_total > 0) {
      ltr <- ldir_results$ldir_raman_agreement$tiered_rates
      log_message("  LDIR-Raman agreement: Family ", ltr$family_or_better_pct,
                  "% (primary), Exact ", ltr$exact_pct, "%, Filler/pigment ",
                  ltr$filler_pct, "% of ", ltr$n_resolved, " resolved pairs; ",
                  ltr$n_unresolved, " unresolved, ", ltr$n_join_warning,
                  " with LDIR join warning")
    }
  }
  log_message("  Results in: ", config$output_dir)
  # Rewrite the parameter dump with the values the run ended with (the Raman
  # image auto-calibration and LDIR coordinate-swap loading update config).
  write_run_environment(config$output_dir, config)
  write_image_objects(config$output_dir, .image_objects)

  if (isTRUE(report)) write_pipeline_report(config)

  invisible(list(config = config, raman = ref, ftir = ftir_pair, ftir_bruker = bruker_pair,
                 ldir = ldir_results,
                 cross = list(ftir_bruker_ftir_perkin = bruker_perkin_match,
                              ftir_bruker_ldir = bruker_ldir_match),
                 composites = composites, tps_assessment = tps_assessment))
}


#' Write the PDF and interactive HTML report of a finished run.
#'
#' Harmonised with the Shiny viewer's "Download report": same page set, same
#' renderers, and the run is loaded back through the viewer's own loaders
#' (load_run_data + build_instrument_dfs), so the report sees exactly the
#' columns and match_status values the viewer sees. The viewer's code is
#' sourced into a private environment, not the global one.
#'
#' @param config Config of the run (config$output_dir = the run folder)
#' @param app_dir Folder of the Shiny viewer
write_pipeline_report <- function(config, app_dir = "shiny_app") {
  report_env <- new.env(parent = globalenv())
  report_env$config <- config
  tryCatch({
    log_message("Generating report...")
    source(file.path(app_dir, "global.R"), local = report_env)
    local({


      # --- Report-only quality gate --------------------------------------------
      # Everything above this point (ingest, image recognition, alignment, ICP,
      # matching, agreement) deliberately runs on the FULL particle set -- narrowing
      # it there would change which pairs the registration can find. Only the
      # report is filtered, using the same per-instrument defaults as the viewer's
      # quality sliders.
      #
      # The scales differ per instrument and are NOT interchangeable: FTIR
      # (PerkinElmer and Bruker) and LDIR carry quality on 0-1, Raman carries HQI
      # on 0-100. Confirmed against the ingested CSVs.
      REPORT_QUALITY_RANGE <- list(
        ftir        = c(0.70, 1),
        ftir_bruker = c(0.70, 1),
        raman       = c(70,   100),   # HQI scale
        ldir        = c(0.85, 1)
      )
      .rq_label <- c(ftir = "FTIR (PerkinElmer)", ftir_bruker = "FTIR (Bruker)",
                     raman = "Raman", ldir = "LDIR")

      # Load the run back through the viewer's loaders so the report and the
      # viewer agree on every column, not just the numbers.
      .run_info <- list(dir = config$output_dir, format = "subdir")
      .rdata <- load_run_data(.run_info)
      .dfs   <- build_instrument_dfs(.rdata)
      # Unfiltered frames + the stored FTIR placements: the report's FTIR pages
      # place their image from these, not from the report-filtered particles.
      .dfs_all <- .dfs
      .ftir_placements <- lapply(FTIR_PLACEMENT_FILES, function(f)
        read_ftir_image_placement(file.path(config$output_dir, f)))

      # Minimum Feret Max (um) for the report, same for every instrument and the
      # same default the viewer's size sliders start at (DEFAULT_MIN_SIZE_UM in
      # global.R, so the two cannot drift apart).
      .report_min_size <- DEFAULT_MIN_SIZE_UM

      # Mirrors filter_instrument() in the viewer: quality first, then size,
      # dropping NA in both (the viewer's sliders bound finite values only).
      for (k in names(REPORT_QUALITY_RANGE)) {
        d <- .dfs[[k]]
        if (is.null(d) || nrow(d) == 0) next
        rng <- REPORT_QUALITY_RANGE[[k]]
        n0  <- nrow(d)
        if ("quality" %in% names(d)) {
          q <- suppressWarnings(as.numeric(d$quality))
          d <- d[!is.na(q) & q >= rng[1] & q <= rng[2], , drop = FALSE]
        }
        n_after_q <- nrow(d)
        if ("feret_max" %in% names(d)) {
          s <- suppressWarnings(as.numeric(d$feret_max))
          d <- d[!is.na(s) & s >= .report_min_size, , drop = FALSE]
        }
        .dfs[[k]] <- d
        log_message(sprintf(
          "  Report filter %-19s quality %s-%s: %d of %d | size >= %s um: %d of %d",
          .rq_label[[k]], format(rng[1]), format(rng[2]), n_after_q, n0,
          format(.report_min_size), nrow(d), n_after_q))
      }

      .quality_note <- c(
        "Filters applied to this report:",
        vapply(names(REPORT_QUALITY_RANGE), function(k) {
          r <- REPORT_QUALITY_RANGE[[k]]
          sprintf("  %-19s quality %s - %s%s", .rq_label[[k]],
                  format(r[1]), format(r[2]),
                  if (identical(k, "raman")) "  (HQI)" else "")
        }, character(1)),
        sprintf("  %-19s Feret Max >= %s um (all instruments)",
                "", format(.report_min_size)),
        "",
        "Alignment and matching upstream used ALL particles; the filter",
        "affects only the figures and tables below.")

      # Background images, resolved the same way the viewer resolves them.
      .man       <- tryCatch(load_run_manifest(config$output_dir),
                             error = function(e) list())
      .meta      <- .man$config_snapshot %||% list()
      .img_paths <- tryCatch(get_run_image_paths(.man, config$output_dir),
                             error = function(e) list())

      report_pages <- build_report_pages(
        .dfs, meta = .meta, img_paths = .img_paths,
        run_label = basename(config$output_dir),
        run_id    = .man$run_id %||% basename(config$output_dir),
        quality_note = .quality_note,
        manifest  = .man,
        placements = .ftir_placements,
        place_dfs  = .dfs_all)

      report_file <- file.path(config$output_dir, "particle_report.pdf")
      n_pages <- write_report_pdf(report_pages, report_file)
      log_message("  Report written: ", report_file, " (", n_pages, " pages)")

      # Interactive HTML: same pages, plus the plotly barplot after the static one.
      html_file <- file.path(config$output_dir, "particle_report.html")
      tryCatch({
        .DEV <- c("FTIR (PerkinElmer)" = "ftir", "FTIR (Bruker)" = "ftir_bruker",
                  "Raman" = "raman", "LDIR" = "ldir")
        device_counts <- lapply(setNames(nm = names(.DEV)), function(lbl) {
          x <- .dfs[[.DEV[[lbl]]]]
          if (is.null(x) || nrow(x) == 0 || !"material" %in% names(x)) return(NULL)
          table(classify_family_vec(x$material))
        })
        .mk_fig <- function(rel) tryCatch(
          build_plotly_barplot(device_counts, rel_mode = rel),
          error = function(e) {
            log_message("  WARNING: plotly build failed: ", e$message); NULL })
        # Interactive twins of the two static barplot pages (2 = absolute,
        # 3 = relative), each inserted directly after its static counterpart.
        # Interactive twins of the two static barplot pages (2 = absolute,
        # 3 = relative), each inserted directly after its static counterpart.
        .bar_specs <- list(
          list(fig = .mk_fig(FALSE), after = 2L,
               title = "Material Comparison \u2014 Interactive (absolute counts)",
               caption = paste0("Hover over bars for exact counts. ",
                                "Use the legend to show/hide families.")),
          list(fig = .mk_fig(TRUE), after = 3L,
               title = "Material Comparison \u2014 Interactive (relative share)",
               caption = paste0("Share of each instrument's own total (%). ",
                                "Hover for the underlying count.")))
        write_report_html(report_pages, html_file, plotly_figs = .bar_specs)
        log_message("  HTML report written: ", html_file)
      }, error = function(e) {
        log_message("  WARNING: HTML report generation failed: ", e$message)
      })

    }, envir = report_env)
  }, error = function(e) {
    log_message("  WARNING: Report generation failed: ", e$message)
  })
  invisible(NULL)
}
