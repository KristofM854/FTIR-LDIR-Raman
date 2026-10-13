# =============================================================================
# ldir_stage.R — LDIR: coordinates from the image, alignment to Raman, matching
# =============================================================================
#
# The LDIR export has no particle coordinates (unless a column-mapping reader
# supplied them). This stage
#   12a  extracts particle positions from the LDIR image (analyzed overlay if
#        present, else the optical image with its scan circle) and joins them
#        to the table by size (join_ldir_coords) plus manual swaps
#   12b  prefilters
#   12c  normalises (normalize_coords_ldir)
#   12d  aligns to Raman: Tier 0 Procrustes from a landmark map, Tier 1
#        landmarks, Tier 2 material (or descriptor) RANSAC vs global
#        registration with reflection forbidden; ICP accepted only if it does
#        not lose one-to-one pairs
#   12g  matches to Raman (optional TPS refinement), scores agreement, matches
#        to FTIR (PerkinElmer), and builds three-way triplets
# =============================================================================


#' Run the LDIR stage.
#'
#' @param ldir_raw Canonical LDIR particle table (read_instrument("ldir", ...))
#' @param ref Raman frame (prepare_raman_reference())
#' @param ftir_pair FTIR (PerkinElmer) result of align_and_match(): its aligned
#'   particles are the LDIR<->FTIR reference, its matches build the triplets
#' @param config Pipeline config
#' @param out_dirs get_output_dirs() of the run
#' @return list(results = the ldir_results list consumed by export and the
#'   viewer, image = the LDIR image object or NULL)
run_ldir_stage <- function(ldir_raw, ref, ftir_pair, config, out_dirs) {
  raman_clean         <- ref$clean
  raman_for_transform <- ref$for_transform
  raman_norm_align    <- ref$norm_align
  raman_for_match     <- ref$for_match
  ftir_aligned        <- ftir_pair$aligned
  match_result        <- ftir_pair$match
  icp_result          <- ftir_pair$icp
  norm_result         <- ftir_pair$norm
  .out_dirs           <- out_dirs
  ldir_image_object    <- NULL
  ldir_image_extracted <- NULL
  .ldir_circle_info    <- NULL
  ldir_landmark_result <- NULL
  use_ldir_landmark    <- FALSE

  log_message(strrep("-", 50))
  log_message("LDIR Spatial Pipeline")

  # 12a. Extract coordinates from LDIR image BEFORE filtering (join on raw data)
  ldir_with_coords <- ldir_raw
  has_ldir_coords  <- FALSE
  ldir_table_has_coords <- sum(!is.na(ldir_raw$x_um) & !is.na(ldir_raw$y_um)) >= 10

  if (ldir_table_has_coords) {
    # Coordinates came with the table (generic column-mapping reader): no
    # image extraction or image join.
    log_message("LDIR coordinates taken from the particle table (no image join)")
    has_ldir_coords <- TRUE
  } else if (!is.null(config$ldir_image) && nzchar(config$ldir_image)) {
    log_message("Extracting LDIR coordinates from companion image")
    # Python is required for LDIR image analysis: fail here, with what is
    # missing, rather than partway through the detection.
    setup_python_detector()

    # Prefer the canonical PNG (format-guaranteed) produced by canonicalize_ldir_image().
    # Fall back to the original path if canonicalization was skipped.
    ldir_img_for_extraction <- if (!is.null(config$ldir_image_canonical) &&
                                    file.exists(config$ldir_image_canonical)) {
      log_message("  Using canonical PNG: ", config$ldir_image_canonical)
      config$ldir_image_canonical
    } else {
      log_message("  Using original image (no canonical available): ", config$ldir_image)
      config$ldir_image
    }

    # Compute scan bounds from configured scan diameter (circular filter)
    ldir_scan_diam <- config$ldir_scan_diameter_um
    if (is.null(ldir_scan_diam)) ldir_scan_diam <- 13000
    ldir_scan_bounds <- list(
      x_min = 0, x_max = ldir_scan_diam,
      y_min = 0, y_max = ldir_scan_diam
    )

    # Prefer the LDIR software's analyzed particle-overlay image when present.
    # Its coloured blobs are the machine's own segmentation and carry
    # correctly-scaled sizes (the optical image over-sizes large particles ~3x),
    # so feeding them to join_ldir_coords fixes size-based matching for large
    # particles.  Falls back to optical-image circle-calibrated extraction.
    #
    # Search next to the ORIGINAL selected image (config$ldir_image), where the
    # user's "<name>_analyzed.<ext>" export lives -- not ldir_img_for_extraction,
    # which may point at the canonicalized copy in the run's inputs/ folder.
    processed_img <- find_ldir_processed_image(config$ldir_image, config)
    if (!is.null(processed_img)) {
      log_message("Using LDIR processed image for coordinate extraction: ",
                  basename(processed_img))
      ldir_proc_result     <- extract_ldir_processed_image_coords(
        processed_img, scan_bounds = ldir_scan_bounds,
        expected_count = nrow(ldir_raw),
        expected_total_area_um2 = sum(ldir_raw$area_um2, na.rm = TRUE),
        config = config
      )
      ldir_image_particles <- ldir_proc_result$particles
      .ldir_circle_info    <- ldir_proc_result$circle_info
      .ldir_image_object   <- ldir_proc_result$image
    } else {
      # Circle-calibrated extraction -- expected_count uses RAW count for best matching
      ldir_extract_result  <- extract_ldir_image_coords(
        ldir_img_for_extraction,
        scan_bounds    = ldir_scan_bounds,
        expected_count = nrow(ldir_raw),
        config         = config
      )
      ldir_image_particles <- ldir_extract_result$particles
      .ldir_circle_info    <- ldir_extract_result$circle_info
      .ldir_image_object   <- ldir_extract_result$image
    }

    # Persist circle calibration to manifest
    tryCatch(
      update_manifest_ldir_circle(config$output_dir, .ldir_circle_info),
      error = function(e)
        log_message("  Could not update manifest ldir_circle: ", e$message, level = "WARN")
    )

    ldir_image_object <- .ldir_image_object

    # Save raw image-extracted coordinates (before join) for Shiny viewer
    ldir_image_extracted <- ldir_image_particles

    # Join image coordinates with RAW Excel data via size-based Hungarian matching
    ldir_with_coords <- join_ldir_coords(ldir_raw, ldir_image_particles, config = config)
    # Apply any manual coordinate-join corrections (config$ldir_coord_swaps)
    ldir_with_coords <- apply_ldir_coord_swaps(ldir_with_coords, config)
    write.csv(ldir_with_coords, file.path(.out_dirs$joined, "ldir_joined_raw.csv"), row.names = FALSE)

    # Validate join quality via scan-order correlation
    scan_order <- validate_ldir_scan_order(ldir_with_coords)
    log_message("  Scan order validation: ", scan_order$message)

    # Debug traces
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_raw, trace_ids, "after_ingest", config$debug_dir)
      dump_particle(ldir_with_coords, trace_ids, "after_coords_join", config$debug_dir)
      trace_particle_snapshot(ldir_raw, "after_ingest", config)
      trace_particle_snapshot(ldir_with_coords, "after_coords_join", config)
    }

    n_with_coords <- sum(!is.na(ldir_with_coords$x_um))
    has_ldir_coords <- n_with_coords >= 10

    # Guard: skip alignment if the circle was not reliably detected
    if (has_ldir_coords && !isTRUE(.ldir_circle_info$detected)) {
      log_message("  LDIR circle not reliably detected \u2014 skipping LDIR alignment",
                  level = "WARN")
      has_ldir_coords <- FALSE
    }

    log_message("  LDIR particles with coordinates: ", n_with_coords,
                " of ", nrow(ldir_with_coords))
  }

  # 12b. Pre-filter LDIR particles AFTER coordinate join
  ldir_clean <- prefilter_ldir(
    ldir_with_coords,
    min_quality = config$ldir_quality_threshold,
    min_size_um = config$min_particle_size_um
  )
  write.csv(ldir_clean, file.path(.out_dirs$prefiltered, "ldir_prefiltered.csv"), row.names = FALSE)

  # Update ldir_with_coords to the filtered version for downstream use
  ldir_with_coords <- ldir_clean

  # 12c-j. Spatial alignment & matching (only if coordinates available)
  ldir_raman_match     <- NULL
  ldir_ftir_match      <- NULL
  ldir_raman_agreement <- NULL
  ldir_icp             <- NULL
  ldir_aligned         <- NULL
  triplets             <- data.frame()

  if (has_ldir_coords) {
    log_message("LDIR spatial matching enabled")

    # 12c. Normalize LDIR coordinates using explicit normalize_coords_ldir()
    #
    # This returns an auditable norm_params object (centroid, scale, y_flip)
    # and writes ldir_norm_params.json to debug_dir if debug=TRUE.
    # The Y-flip toggle resolves the reflection ambiguity between LDIR and Raman.
    ldir_norm_result <- normalize_coords_ldir(
      ldir_with_coords,
      flip_y       = isTRUE(config$ldir_flip_y_for_alignment),
      scale_coords = isTRUE(config$normalize_scale),
      rotate_deg   = cfg_value(config, "ldir_rotate_deg_for_alignment"),
      debug_dir    = if (isTRUE(config$debug)) config$debug_dir else NULL
    )
    ldir_with_coords <- ldir_norm_result$df
    ldir_norm_params <- ldir_norm_result$norm_params

    # Step 5: dump traced particles after normalization
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_with_coords, trace_ids, "after_normalization", config$debug_dir)
      trace_particle_snapshot(ldir_with_coords, "after_normalization", config)
    }

    # Step 3: ldir_valid is now extracted AFTER x_norm/y_norm exist
    ldir_valid <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                    !is.na(ldir_with_coords$y_um), ]

    # Step 1: dump traced particles after coords join (before normalization)
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_with_coords, trace_ids, "after_coords_join", config$debug_dir)
    }

    # 12d. Tiered LDIR->Raman alignment
    #
    # Tier 0 (Step 3): Explicit Procrustes -- if config$ldir_landmark_map is set,
    #   use named LDIR<->Raman correspondences to fit via SVD. This guarantees
    #   A3 (and other named landmarks) have minimal residuals by construction.
    #   When ldir_procrustes_lock=TRUE (default), this is the FINAL transform.
    #   ICP still runs but only for diagnostics (its output is not used).
    #
    # Tier 1: Size-based landmark RANSAC (large particles & fibers).
    # Tier 2: Material-based RANSAC: PET/PP/PC anchors (fallback).

    ldir_procrustes  <- NULL
    ldir_anchor_pairs <- NULL   # passed to ICP for pinned-weight refinement
    use_procrustes_final <- FALSE

    if (!is.null(config$ldir_landmark_map) && length(config$ldir_landmark_map) > 0) {
      log_message(strrep("-", 50))
      log_message("  LDIR Tier 0: Explicit Procrustes alignment from landmark map")

      ldir_procrustes <- fit_similarity_from_landmarks(
        src_df       = ldir_valid,
        tgt_df       = raman_for_transform,
        landmark_map = config$ldir_landmark_map,
        debug_dir    = if (isTRUE(config$debug)) config$debug_dir else NULL
      )

      if (ldir_procrustes$success) {
        log_message("  Procrustes SUCCESS: ", ldir_procrustes$message)

        # Build anchor_pairs (src_idx in ldir_valid, tgt_idx in raman_for_transform)
        # so ICP can pin them with high weight even when refining
        lp <- ldir_procrustes$landmark_pairs
        if (nrow(lp) > 0) {
          a_src <- match(lp$ldir_id,  ldir_valid$particle_id)
          a_tgt <- match(lp$raman_id, raman_for_transform$particle_id)
          ok_a  <- !is.na(a_src) & !is.na(a_tgt)
          if (sum(ok_a) > 0) {
            ldir_anchor_pairs <- data.frame(
              src_idx = a_src[ok_a],
              tgt_idx = a_tgt[ok_a]
            )
          }
        }

        use_procrustes_final <- isTRUE(config$ldir_procrustes_lock)
      } else {
        log_message("  Procrustes failed (", ldir_procrustes$message,
                    ") \u2014 falling through to RANSAC", level = "WARN")
      }
    }

    ldir_landmark_result <- tryCatch({
      landmark_align(ldir_valid, raman_for_transform, config, src_label = "LDIR")
    }, error = function(e) {
      log_message("  LDIR landmark alignment failed: ", e$message, level = "WARN")
      list(success = FALSE, confident = FALSE,
           n_ftir_landmarks = 0, n_raman_landmarks = 0)
    })

    use_ldir_landmark <- !use_procrustes_final &&
                         ldir_landmark_result$confident &&
                         config$landmark_skip_full_ransac

    if (use_ldir_landmark) {
      log_message("  Using LDIR landmark transform (Tier 1) \u2014 skipping material RANSAC")
      ldir_ransac <- list(
        transform = ldir_landmark_result$transform,
        params    = ldir_landmark_result$params,
        n_inliers = ldir_landmark_result$n_inliers
      )
    } else if (!use_procrustes_final) {
      # Tier 2: Material-based RANSAC (default) or descriptor RANSAC (opt-in)
      log_message(strrep("-", 50))

      if (isTRUE(config$ldir_use_descriptor_ransac)) {
        log_message("  LDIR Tier 2: Descriptor RANSAC alignment")
        ldir_for_align <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
        log_message("  LDIR alignment anchors: ", nrow(ldir_for_align), " particles")

        ldir_ransac <- tryCatch({
          descriptor_ransac_align(
            ldir_for_align, raman_norm_align, config,
            scale_min     = config$icp_min_scale        %||% 0.5,
            scale_max     = config$icp_max_scale        %||% 2.0,
            rot_limit_deg = 180   # any orientation (filters can be rotated)
          )
        }, error = function(e) {
          log_message("  Descriptor RANSAC failed: ", e$message, level = "WARN")
          NULL
        })

        if (is.null(ldir_ransac)) {
          log_message("  Descriptor RANSAC returned NULL \u2014 falling back to material RANSAC",
                      level = "WARN")
          config$ldir_use_descriptor_ransac <- FALSE  # force fallback in logging
        }
      }

      if (!isTRUE(config$ldir_use_descriptor_ransac)) {
        log_message("  LDIR Tier 2: Material-based RANSAC alignment")
        ldir_for_align <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
        if (!is.null(config$align_ldir_materials)) {
          ldir_mat_mask <- grepl(
            paste(config$align_ldir_materials, collapse = "|"),
            ldir_for_align$material, ignore.case = TRUE
          )
          if (sum(ldir_mat_mask) >= 4) {
            ldir_for_align <- ldir_for_align[ldir_mat_mask, ]
          } else {
            log_message("  Insufficient LDIR anchor materials (",
                        sum(ldir_mat_mask),
                        ") \u2014 using all LDIR particles for alignment")
          }
        }
        log_message("  LDIR alignment anchors: ", nrow(ldir_for_align), " particles")

        # LDIR and Raman image the same filter from the same side, so there is
        # no physical reflection between them. Allowing mirror lets the aligner
        # lock onto a spurious REFLECTED optimum that matches sparse particles
        # to the wrong neighbours (observed: reflected transform, wrong
        # scale/rotation, ~19 spurious matches). Forbid reflection unless the
        # user explicitly enables it.
        .ldir_reflect <- isTRUE(config$ldir_allow_reflection)
        .ldir_align_cfg <- config
        .ldir_align_cfg$ransac_allow_mirror <- .ldir_reflect
        if (!.ldir_reflect)
          log_message("  LDIR alignment: reflection disabled (no physical ",
                      "mirror between LDIR and Raman)")

        ldir_ransac <- tryCatch({
          ransac_align(ldir_for_align, raman_norm_align, .ldir_align_cfg,
                       src_label = "LDIR")
        }, error = function(e) {
          log_message("  LDIR RANSAC failed: ", e$message, level = "WARN")
          NULL
        })

        # Global registration: robust for sparse LDIR<->Raman where the
        # coarse RANSAC locks onto a poor local optimum (few inliers despite
        # many achievable). Run it and keep whichever transform pairs more
        # particles. Disable with config$ldir_use_global_register = FALSE.
        if (!isFALSE(config$ldir_use_global_register)) {
          ldir_global <- tryCatch(
            global_register_align(ldir_for_align, raman_norm_align,
                                  .ldir_align_cfg, allow_mirror = .ldir_reflect),
            error = function(e) {
              log_message("  LDIR global registration failed: ",
                          e$message, level = "WARN"); NULL })
          if (!is.null(ldir_global)) {
            r_in <- if (!is.null(ldir_ransac)) ldir_ransac$n_inliers else -1
            if (ldir_global$n_inliers > r_in) {
              log_message("  LDIR alignment: using global registration (",
                          ldir_global$n_inliers, " inliers vs RANSAC ",
                          max(r_in, 0), ")")
              ldir_ransac <- ldir_global
            } else {
              log_message("  LDIR alignment: keeping RANSAC (", r_in,
                          " inliers vs global ", ldir_global$n_inliers, ")")
            }
          }
        }
      }
    } else {
      # Procrustes locked -- still set ldir_ransac from Procrustes for logging
      ldir_ransac <- list(
        transform = ldir_procrustes$matrix,
        params    = ldir_procrustes$params,
        n_inliers = ldir_procrustes$n_pairs
      )
    }

    if (!is.null(ldir_ransac)) {
      if (!use_procrustes_final) {
        log_message("  LDIR-Raman RANSAC: scale=",
                    round(ldir_ransac$params$scale, 4),
                    ", rot=", round(ldir_ransac$params$rotation_deg, 2), " deg",
                    ", reflected=", ldir_ransac$params$reflected,
                    ", inliers=", ldir_ransac$n_inliers)
      }

      # Initial transform for ICP (Procrustes if locked, RANSAC otherwise)
      icp_initial_transform <- if (use_procrustes_final) {
        ldir_procrustes$matrix
      } else {
        ldir_ransac$transform
      }

      # ICP refinement on all LDIR particles with coordinates
      # Anchor pairs are pinned with 100x weight (Step 4)
      ldir_for_icp <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                        !is.na(ldir_with_coords$y_um), ]

      # Resolve anchor_pairs indices into ldir_for_icp row space
      icp_anchor <- NULL
      if (!is.null(ldir_anchor_pairs) && nrow(ldir_anchor_pairs) > 0) {
        # ldir_anchor_pairs$src_idx points into ldir_valid; remap to ldir_for_icp
        valid_ids <- ldir_valid$particle_id[ldir_anchor_pairs$src_idx]
        icp_src   <- match(valid_ids, ldir_for_icp$particle_id)
        ok_remap  <- !is.na(icp_src)
        if (sum(ok_remap) > 0) {
          icp_anchor <- data.frame(
            src_idx = icp_src[ok_remap],
            tgt_idx = ldir_anchor_pairs$tgt_idx[ok_remap]
          )
          log_message("  ICP anchor pairs (pinned 100x weight): ",
                      sum(ok_remap), " landmarks")
        }
      }

      ldir_icp <- tryCatch({
        icp_refine(ldir_for_icp, raman_for_transform,
                   icp_initial_transform, config,
                   anchor_pairs = icp_anchor)
      }, error = function(e) {
        log_message("  LDIR ICP failed: ", e$message,
                    " \u2014 using initial transform", level = "WARN")
        list(
          transform    = icp_initial_transform,
          params       = extract_transform_params(icp_initial_transform),
          converged    = FALSE,
          n_iterations = 0,
          rms_history  = numeric(0)
        )
      })

      log_message("  LDIR-Raman ICP: scale=",
                  round(ldir_icp$params$scale, 4),
                  ", rot=", round(ldir_icp$params$rotation_deg, 2), " deg",
                  ", converged=", ldir_icp$converged)

      # --- Accept ICP only if it did not make the pairing worse --------------
      # ICP minimises nearest-neighbour RMS, and that objective rewards
      # collapse (a shrunken cloud sits inside dense regions of the target and
      # scores a LOWER RMS than the correct pose). So RMS cannot be used to
      # decide whether the refinement helped. One-to-one inlier count can: it
      # cannot be inflated by shrinking, because each reference particle is
      # consumed at most once. If ICP lost pairs, keep the pose it started from.
      .icp_tol <- config$ransac_inlier_dist_um %||% 200
      .n_init <- align_transform_inliers(icp_initial_transform,
                                         ldir_for_icp$x_norm, ldir_for_icp$y_norm,
                                         raman_for_transform$x_norm,
                                         raman_for_transform$y_norm, .icp_tol)
      .n_icp  <- align_transform_inliers(ldir_icp$transform,
                                         ldir_for_icp$x_norm, ldir_for_icp$y_norm,
                                         raman_for_transform$x_norm,
                                         raman_for_transform$y_norm, .icp_tol)
      if (is.finite(.n_init) && is.finite(.n_icp)) {
        log_message("  LDIR-Raman ICP pairing: ", .n_init, " -> ", .n_icp,
                    " one-to-one inliers at ", .icp_tol, " um")
        if (.n_icp < .n_init) {
          log_message("  LDIR-Raman ICP REJECTED: refinement lost ",
                      .n_init - .n_icp, " pairs. Reverting to the pre-ICP ",
                      "transform (scale=",
                      round(extract_transform_params(icp_initial_transform)$scale, 4),
                      ", rot=",
                      round(extract_transform_params(icp_initial_transform)$rotation_deg, 2),
                      " deg).", level = "WARN")
          ldir_icp$transform <- icp_initial_transform
          ldir_icp$params    <- extract_transform_params(icp_initial_transform)
          ldir_icp$icp_rejected <- TRUE
        }
      }

      # When Procrustes is locked, discard ICP transform and keep Procrustes
      if (use_procrustes_final) {
        log_message("  Procrustes lock active \u2014 retaining Procrustes transform ",
                    "(ICP ran for diagnostics only)")
        ldir_icp$transform <- ldir_procrustes$matrix
        ldir_icp$params    <- ldir_procrustes$params
      }

      # Quality check
      ldir_final_rms <- if (length(ldir_icp$rms_history) > 0)
        tail(ldir_icp$rms_history, 1) else NA_real_
      ftir_final_rms  <- if (length(icp_result$rms_history) > 0)
        tail(icp_result$rms_history, 1) else NA_real_

      if (!is.na(ldir_final_rms) && ldir_final_rms > 100) {
        log_message("  LDIR alignment quality: POOR (ICP RMS = ",
                    round(ldir_final_rms, 1), " \u00b5m). Possible causes:",
                    level = "WARN")
        log_message("    1. LDIR scan area (ldir_scan_diameter_um=",
                    config$ldir_scan_diameter_um, " \u00b5m) may not match actual scan",
                    level = "WARN")
        log_message("    2. Too few anchor material particles for robust RANSAC",
                    level = "WARN")
        log_message("    3. Systematic coordinate flip or offset in image extraction",
                    level = "WARN")
        log_message("    Recommendation: set config$ldir_landmark_map with explicit",
                    " A3\u2194Raman correspondences, or verify ldir_scan_diameter_um",
                    level = "WARN")
      } else if (!is.na(ldir_final_rms) && !is.na(ftir_final_rms) &&
                 ldir_final_rms > 3 * ftir_final_rms) {
        log_message("  LDIR alignment quality: MARGINAL (RMS ",
                    round(ldir_final_rms, 1), " \u00b5m vs FTIR ",
                    round(ftir_final_rms, 1), " \u00b5m)", level = "WARN")
      }

      # 12f. Apply transform to all LDIR particles with coordinates
      ldir_aligned <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                        !is.na(ldir_with_coords$y_um), ]
      ldir_tf <- apply_transform_points(
        ldir_aligned$x_norm, ldir_aligned$y_norm, ldir_icp$transform
      )
      ldir_aligned$x_aligned <- ldir_tf$x_transformed
      ldir_aligned$y_aligned <- ldir_tf$y_transformed
      # Tag which transform path was used (for provenance in debug CSV)
      ldir_aligned$align_method <- if (use_procrustes_final) "procrustes" else
                                    if (use_ldir_landmark) "landmark_ransac" else
                                    if (isTRUE(config$ldir_use_descriptor_ransac))
                                      "descriptor_ransac" else "ransac_icp"

      # --- Transform guardrails ---
      tryCatch({
        .M      <- ldir_icp$transform
        .a      <- .M[1, 1]; .b <- .M[2, 1]
        .tf_sc  <- sqrt(.a^2 + .b^2)
        .tf_rot <- atan2(.b, .a) * 180 / pi
        .tf_tx  <- .M[1, 3]; .tf_ty <- .M[2, 3]
        .fov_um <- config$ldir_scan_diameter_um %||% 13000
        if (.tf_sc < (config$icp_min_scale %||% 0.5) ||
            .tf_sc > (config$icp_max_scale %||% 2.0))
          log_message("  ALIGNMENT GUARDRAIL: scale=", round(.tf_sc, 3),
                      " outside [", config$icp_min_scale %||% 0.5, ", ",
                      config$icp_max_scale %||% 2.0,
                      "] \u2014 transform may be unreliable", level = "WARN")
        if (sqrt(.tf_tx^2 + .tf_ty^2) > .fov_um)
          log_message("  ALIGNMENT GUARDRAIL: translation=",
                      round(sqrt(.tf_tx^2 + .tf_ty^2)), " \u00b5m > FOV (",
                      .fov_um, " \u00b5m) \u2014 possible offset error", level = "WARN")
      }, error = function(e)
        log_message("  Guardrail check failed: ", e$message, level = "WARN"))

      # --- Alignment diagnostics (residuals, quiver, stats, audit) ---
      .residuals_for_overlay <- NULL
      tryCatch({
        .residuals_for_overlay <- compute_alignment_residuals(ldir_aligned, raman_clean)
        write_alignment_diagnostics(
          ldir_aligned     = ldir_aligned,
          raman_clean      = raman_clean,
          ldir_norm_params = ldir_norm_params,
          norm_result      = norm_result,
          ldir_icp         = ldir_icp,
          config           = config
        )
      }, error = function(e) {
        log_message("  Alignment diagnostics failed: ", e$message, level = "WARN")
        .residuals_for_overlay <- NULL
      })

      # --- Overlay diagnostic PNG (with optional residual arrows) ---
      tryCatch({
        diag_dir   <- file.path(config$output_dir, "debug")
        if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
        ldir_diag  <- ldir_aligned[is.finite(ldir_aligned$x_aligned) &
                                    is.finite(ldir_aligned$y_aligned), ]
        raman_diag <- raman_clean[is.finite(raman_clean$x_norm) &
                                   is.finite(raman_clean$y_norm), ]

        p_overlay <- ggplot2::ggplot() +
          ggplot2::geom_point(data = raman_diag,
                              ggplot2::aes(x = x_norm, y = y_norm),
                              shape = 1, colour = "steelblue", alpha = 0.5, size = 1.5) +
          ggplot2::geom_point(data = ldir_diag,
                              ggplot2::aes(x = x_aligned, y = y_aligned),
                              shape = 2, colour = "forestgreen", alpha = 0.5, size = 1.5) +
          ggplot2::coord_fixed() + ggplot2::theme_minimal() +
          ggplot2::labs(
            title    = "LDIR-Raman overlay diagnostic",
            subtitle = paste0("LDIR: green \u25b2 (x_aligned/y_aligned)  |  ",
                              "Raman: blue \u25cb (x_norm/y_norm)  |  ",
                              "Orange arrows: residuals (up to 150)"),
            x = "\u00b5m", y = "\u00b5m"
          )

        # Add residual arrows if available (downsampled to 150)
        if (!is.null(.residuals_for_overlay) &&
            nrow(.residuals_for_overlay) > 0) {
          arrow_df <- if (nrow(.residuals_for_overlay) > 150)
            # evenly spaced, not random: the same run must draw the same plot
            .residuals_for_overlay[unique(round(seq(1, nrow(.residuals_for_overlay),
                                                    length.out = 150))), ]
          else .residuals_for_overlay
          p_overlay <- p_overlay +
            ggplot2::geom_segment(
              data = arrow_df,
              ggplot2::aes(x = ldir_x_aligned, y = ldir_y_aligned,
                           xend = raman_x_norm,  yend = raman_y_norm),
              colour = "orange", alpha = 0.45, linewidth = 0.35,
              arrow  = ggplot2::arrow(length = ggplot2::unit(0.06, "inches"),
                                      type   = "open")
            )
        }

        ggplot2::ggsave(file.path(diag_dir, "ldir_raman_overlay_diag.png"),
                        p_overlay, width = 8, height = 8, dpi = 150)
        log_message("  Overlay diagnostic saved: ", diag_dir, "/ldir_raman_overlay_diag.png")
      }, error = function(e)
        log_message("  Could not save overlay diagnostic: ", e$message, level = "WARN"))

      # --- Persist alignment method + matrix to manifest ---
      tryCatch({
        .align_method_str <- ldir_aligned$align_method[1] %||% "ransac_icp"
        if (!is.null(.ldir_circle_info))   # absent when the table had coordinates
          update_manifest_ldir_circle(config$output_dir, .ldir_circle_info)  # refresh
        m_path <- resolve_manifest_path(config$output_dir)
        if (file.exists(m_path) && requireNamespace("jsonlite", quietly = TRUE)) {
          m_upd <- jsonlite::fromJSON(m_path, simplifyVector = FALSE)
          m_upd$ldir_raman_alignment_method <- .align_method_str
          m_upd$ldir_raman_similarity_matrix <- lapply(
            seq_len(nrow(ldir_icp$transform)),
            function(i) as.numeric(ldir_icp$transform[i, ])
          )
          writeLines(jsonlite::toJSON(m_upd, pretty = TRUE, auto_unbox = TRUE,
                                       null = "null", na = "null"), m_path)
        }
      }, error = function(e)
        log_message("  Could not persist alignment to manifest: ", e$message,
                    level = "WARN"))

      # Step 4: Assert aligned coordinates exist before any plotting
      stopifnot(all(c("x_aligned", "y_aligned") %in% colnames(ldir_aligned)))

      # Step 5: dump traced particles + full snapshot after alignment
      if (isTRUE(config$debug)) {
        trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
        dump_particle(ldir_aligned, trace_ids, "after_alignment", config$debug_dir)
        trace_particle_snapshot(ldir_aligned, "after_alignment", config)
      }

      # 12g. LDIR<->Raman matching
      ldir_raman_match <- match_particles(
        ldir_aligned, raman_for_match, config,
        src_label = "ldir", ref_label = "raman"
      )

      # 12g-bis. TPS local refinement: a single global similarity can't overlay
      # every particle when the two coordinate systems differ by a small
      # non-rigid distortion, leaving peripheral particles unmatched despite
      # clearly corresponding. Fit a regularized thin-plate spline to the
      # residual displacement at the confident matches, warp all LDIR
      # coordinates locally, and re-match -- adopting the result only if it
      # increases matches. Displacement is capped so non-overlapping debris
      # cannot be flung onto spurious partners. (No-op under force-complete
      # matching, which already matches everything.)
      if (isTRUE(config$ldir_tps_refine) &&
          !is.null(ldir_raman_match$matched) &&
          nrow(ldir_raman_match$matched) >= (config$ldir_tps_min_controls %||% 6)) {
        mm <- ldir_raman_match$matched
        warp <- tryCatch(
          tps_fit_warp(mm$ldir_x_aligned, mm$ldir_y_aligned,
                       mm$raman_x_norm,   mm$raman_y_norm,
                       lambda = config$ldir_tps_lambda %||% 0.5),
          error = function(e) { log_message("  TPS fit failed: ",
                                            e$message, level = "WARN"); NULL })
        if (!is.null(warp)) {
          w <- tps_apply(warp, ldir_aligned$x_aligned, ldir_aligned$y_aligned)
          ldir_tps <- ldir_aligned
          ldir_tps$x_aligned <- w$x; ldir_tps$y_aligned <- w$y
          match2 <- try_or(
            match_particles(ldir_tps, raman_for_match, config,
                            src_label = "ldir", ref_label = "raman"),
            default = NULL, what = "LDIR TPS re-match")
          n0 <- ldir_raman_match$match_stats$n_matched
          n1 <- if (!is.null(match2)) match2$match_stats$n_matched else -1L
          if (!is.null(match2) && n1 > n0) {
            log_message("  LDIR TPS refinement: matches ", n0, " -> ", n1,
                        " (control residual ", round(warp$ctrl_res), " um, ",
                        warp$n, " controls)")
            ldir_aligned     <- ldir_tps
            ldir_raman_match <- match2
          } else {
            log_message("  LDIR TPS refinement: no gain (", n0, " vs ",
                        max(n1, 0L), ") \u2014 keeping global alignment")
          }
        }
      }

      # 12h. LDIR<->Raman agreement analysis
      ldir_raman_agreement <- analyze_agreement(
        ldir_raman_match, config,
        instrument_a = "LDIR", instrument_b = "Raman"
      )

      # 12i. LDIR<->FTIR matching (both already in Raman coordinate frame)
      # Set up aligned FTIR as reference: matcher expects ref with x_norm/y_norm
      ftir_as_ref <- ftir_aligned
      ftir_as_ref$x_norm <- ftir_as_ref$x_aligned
      ftir_as_ref$y_norm <- ftir_as_ref$y_aligned
      ldir_ftir_match <- match_particles(
        ldir_aligned, ftir_as_ref, config,
        src_label = "ldir", ref_label = "ftir"
      )
      log_message("  LDIR-FTIR direct matching: ",
                  ldir_ftir_match$match_stats$n_matched, " pairs")

      # 12j. Three-way triplets (FTIR<->Raman and LDIR<->Raman sharing a Raman
      # particle), scored with the agreement quality gate
      triplets <- build_triplets(match_result, ldir_raman_match, raman_for_match, config)

    } else {
      log_message("  LDIR spatial alignment failed \u2014 ",
                  "material comparison only", level = "WARN")
    }

  } else {
    log_message("  LDIR: no spatial coordinates \u2014 material comparison only")
  }

  # Non-spatial material distribution (always available)
  ldir_mats <- table(ldir_clean$material)
  log_message("  LDIR material distribution (top 5):")
  for (m in names(head(sort(ldir_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", ldir_mats[m])
  }

  ldir_results <- list(
    ldir_clean           = ldir_clean,
    ldir_with_coords     = ldir_with_coords,
    ldir_aligned         = ldir_aligned,
    ldir_raman_match     = ldir_raman_match,
    ldir_ftir_match      = ldir_ftir_match,
    ldir_raman_agreement = ldir_raman_agreement,
    ldir_icp             = ldir_icp,
    ldir_landmark        = ldir_landmark_result,
    ldir_alignment_method = if (isTRUE(use_ldir_landmark)) "landmark" else "ransac",
    triplets             = triplets,
    has_coords           = has_ldir_coords,
    ldir_material_dist   = ldir_mats,
    ldir_image_extracted = ldir_image_extracted
  )

  # --- Debug: Branch A/B comparison + residual vectors ---
  if (isTRUE(config$debug) && has_ldir_coords && !is.null(ldir_aligned) &&
      !is.null(ldir_raman_match)) {
    debug_ldir_branches(
      ldir_aligned, raman_clean, ldir_raman_match,
      ldir_icp, config
    )
  }

  list(results = ldir_results, image = ldir_image_object)
}
