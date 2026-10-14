# =============================================================================
# align_and_match.R — the shared Raman frame and the per-pair spatial steps
# =============================================================================
#
# Raman is the reference: every other instrument is aligned into the
# Raman-normalised frame (Raman coordinates minus the Raman full-cloud
# centroid). prepare_raman_reference() builds that frame once;
# align_and_match() aligns one FTIR-type particle table (PerkinElmer, Bruker)
# to it, matches and scores agreement; match_aligned_pair() matches two
# instruments that are both already in the frame. The LDIR path, which first
# needs coordinates from its image, is in R/ldir_stage.R.
# =============================================================================


#' Centroid and scale of a particle cloud (the normalisation of one instrument).
#'
#' @param df Particle table with x_um, y_um
#' @param normalize_scale Also divide by the RMS distance from the centroid
#' @param what Label for messages
#' @return list(centroid = c(x, y), scale)
cloud_normalization <- function(df, normalize_scale = FALSE, what = "particle cloud") {
  if (nrow(df) == 0)
    stop("Cannot normalize coordinates: ", what, " is empty. ",
         "A centroid needs at least one particle.", call. = FALSE)
  if (all(is.na(df$x_um)) || all(is.na(df$y_um)))
    stop("Cannot normalize coordinates: a coordinate column of ", what,
         " is entirely NA.", call. = FALSE)
  cx <- mean(df$x_um, na.rm = TRUE)
  cy <- mean(df$y_um, na.rm = TRUE)
  scale <- 1
  if (isTRUE(normalize_scale)) {
    rms <- sqrt(mean((df$x_um - cx)^2 + (df$y_um - cy)^2))
    if (is.finite(rms) && rms > 0) scale <- rms
  }
  log_message("  ", what, " centroid: (", round(cx, 1), ", ", round(cy, 1), ")",
              if (scale != 1) paste0(", RMS scale ", round(scale, 1)) else "")
  list(centroid = c(cx, cy), scale = scale)
}


#' Raman particle sets in the shared (Raman-normalised) frame.
#'
#' @param raman_raw Canonical Raman particle table
#' @param config Pipeline config
#' @return list with
#'   prefiltered    Raman after the coordinate prefilter, before normalisation
#'   clean          the same, normalised (x_norm, y_norm): the full cloud
#'   norm_align     material/size/HQI anchors for the coarse aligners
#'   for_transform  particles >= align_raman_min_size_um, for landmarks and ICP
#'   for_match      all particles >= min_particle_size_um, for matching
#'   centroid, scale  the frame's normalisation
prepare_raman_reference <- function(raman_raw, config) {
  anchor_min <- config$align_min_anchor_count %||% 4
  clean <- prefilter_raman(raman_raw, min_hqi = 0, min_size_um = 0)

  # Alignment anchors: material, then size, then HQI. Each filter is honoured
  # only while it leaves enough anchors: three filters that are individually
  # reasonable can intersect to nothing.
  for_align <- select_material_anchors(clean, config$align_raman_materials,
                                       anchor_min, "Raman")
  min_size <- config$align_raman_min_size_um
  if (!is.null(min_size) && min_size > 0 && any(!is.na(for_align$feret_max_um))) {
    size_mask <- is.na(for_align$feret_max_um) | for_align$feret_max_um >= min_size
    if (sum(size_mask) >= anchor_min) {
      n_before <- nrow(for_align)
      for_align <- for_align[size_mask, ]
      log_message("Raman alignment: removed ", n_before - nrow(for_align),
                  " particles < ", min_size, " um")
    } else {
      log_message("Raman alignment: size filter (>= ", min_size, " um) would ",
                  "leave only ", sum(size_mask), " anchors \u2014 skipping it",
                  level = "WARN")
    }
  }
  # Only reliably identified particles should drive material-based alignment
  hqi <- config$raman_hqi_threshold
  if (!is.null(hqi) && hqi > 0 && any(!is.na(for_align$quality))) {
    hqi_mask <- !is.na(for_align$quality) & for_align$quality >= hqi
    if (sum(hqi_mask) >= anchor_min) {
      n_before <- nrow(for_align)
      for_align <- for_align[hqi_mask, ]
      log_message("Raman alignment: HQI filter (>= ", hqi, "): kept ",
                  nrow(for_align), " of ", n_before)
    } else {
      log_message("Raman alignment: HQI filter (>= ", hqi, ") would leave only ",
                  sum(hqi_mask), " anchors \u2014 skipping it", level = "WARN")
    }
  }
  log_message("Raman alignment target particles: ", nrow(for_align),
              " (material + >= ", min_size, " um + HQI >= ", hqi, ")")

  # All particles for spatial matching (no HQI filter: HQI is applied only in
  # the agreement analysis, so matching has full coverage).
  for_match <- prefilter_raman(raman_raw, min_hqi = 0,
                               min_size_um = config$min_particle_size_um)
  log_message("Raman for matching: ", nrow(for_match),
              " particles (all HQI, spatial only; HQI filter applied in agreement)")

  # The frame: centroid of the FULL cleaned cloud, not of an anchor subset (a
  # subset centroid depends on which particles each instrument detected, and
  # would start the aligners from a biased translation).
  log_message("Normalizing coordinates (Raman reference frame)")
  nm <- cloud_normalization(clean, config$normalize_scale, "Raman")
  norm <- function(df) apply_normalization(df, nm$centroid, nm$scale)
  clean_n <- norm(clean)

  # Landmark / ICP set: Raman particles large enough to be visible to every
  # instrument; only size matters for geometry.
  for_transform <- clean_n
  if (!is.null(min_size) && min_size > 0 && any(!is.na(clean_n$feret_max_um)))
    for_transform <- clean_n[is.na(clean_n$feret_max_um) | clean_n$feret_max_um >= min_size, ]
  log_message("Raman for spatial transform (>= ", min_size, " um): ",
              nrow(for_transform), " particles")

  list(prefiltered = clean, clean = clean_n, norm_align = norm(for_align),
       for_align = for_align, for_transform = for_transform,
       for_match = norm(for_match), centroid = nm$centroid, scale = nm$scale)
}


#' Align an FTIR-type particle table to the Raman frame, match and score.
#'
#' Tiered alignment, then matching and agreement:
#'   Tier 1  landmark alignment on large particles and fibres; used when
#'           confident (and landmark_skip_full_ransac is TRUE)
#'   Tier 2  coarse RANSAC and global registration on the material anchors;
#'           the one pairing more particles wins
#'   ICP     always, to polish the transform
#'
#' @param src_raw Canonical particle table (ftir or ftir_bruker role)
#' @param ref Result of prepare_raman_reference()
#' @param config Pipeline config
#' @param label Instrument label for messages ("FTIR (PerkinElmer)")
#' @return list: prefiltered, clean, for_align, for_match, aligned, aligned_all,
#'   landmark, ransac, method, icp, match, agreement, norm. With too few
#'   particles to align (warned) only label, prefiltered and skipped = TRUE.
align_and_match <- function(src_raw, ref, config, label) {
  log_message(strrep("-", 50))
  log_message(label, " \u2194 Raman: alignment and matching")
  anchor_min <- config$align_min_anchor_count %||% 4

  prefiltered <- prefilter_ftir(src_raw, min_quality = 0, min_size_um = 0)
  # Material anchors for the coarse aligners; select_material_anchors() falls
  # back to the full cloud when too few carry an anchor material.
  for_align <- select_material_anchors(prefiltered, config$align_ftir_materials,
                                       anchor_min, label)
  for_match <- prefilter_ftir(src_raw, min_quality = config$ftir_quality_threshold,
                              min_size_um = config$min_particle_size_um)
  if (nrow(for_align) < 3) {
    log_message(label, ": too few particles (", nrow(for_align),
                ") for alignment \u2014 skipping", level = "WARN")
    return(list(label = label, prefiltered = prefiltered, skipped = TRUE))
  }

  # Normalisation: the source's full-cloud centroid; anchor and match subsets
  # use the same one so they stay in one frame.
  nm <- cloud_normalization(prefiltered, config$normalize_scale, label)
  norm <- function(df) apply_normalization(df, nm$centroid, nm$scale)
  clean <- norm(prefiltered)
  norm_align <- norm(for_align)
  for_match <- norm(for_match)

  # Tier 1: landmarks
  landmark <- landmark_align(clean, ref$for_transform, config, src_label = label)
  ransac <- NULL
  if (landmark$confident && config$landmark_skip_full_ransac) {
    log_message("Using landmark transform (Tier 1) \u2014 skipping full RANSAC")
    transform <- landmark$transform
    method <- "landmark"
  } else {
    # Tier 2: coarse RANSAC vs global registration (rotation x scale sweep,
    # translation by voting; robust when the clouds overlap only partly).
    # Disable the latter with config$ftir_use_global_register = FALSE.
    log_message(strrep("-", 50))
    log_message("Tier 2: Full alignment (anchor sets)")
    ransac <- ransac_align(norm_align, ref$norm_align, config, src_label = label)
    log_message("RANSAC transform: scale = ", round(ransac$params$scale, 4),
                ", rotation = ", round(ransac$params$rotation_deg, 2), " deg",
                ", reflected = ", ransac$params$reflected,
                ", inliers = ", ransac$n_inliers)
    transform <- ransac$transform
    method <- "ransac"
    if (!isFALSE(config$ftir_use_global_register)) {
      global <- tryCatch(
        global_register_align(norm_align, ref$norm_align, config,
                              allow_mirror = config$ransac_allow_mirror),
        error = function(e) {
          log_message("Global registration failed: ", e$message, level = "WARN")
          NULL })
      if (!is.null(global)) {
        n_r <- tier2_pairs(ransac$transform, norm_align, ref$norm_align, config)
        n_g <- tier2_pairs(global$transform, norm_align, ref$norm_align, config)
        if (n_g > n_r) {
          log_message("Tier 2: using global registration (", n_g,
                      " one-to-one pairs vs RANSAC ", n_r, ")")
          ransac <- global
          transform <- global$transform
          method <- "global_register"
        } else {
          log_message("Tier 2: keeping RANSAC (", n_r,
                      " one-to-one pairs vs global ", n_g, ")")
        }
      }
    }
  }

  icp <- icp_refine(clean, ref$for_transform, transform, config)
  log_message("ICP refined transform: scale = ", round(icp$params$scale, 4),
              ", rotation = ", round(icp$params$rotation_deg, 2), " deg",
              ", converged = ", icp$converged)

  aligned     <- apply_ftir_transform(for_match, icp$transform)
  aligned_all <- apply_ftir_transform(clean, icp$transform)

  match <- match_particles(aligned, ref$for_match, config)
  ms <- match$match_stats
  log_message(label, " \u2194 Raman: ", ms$n_matched, " matched, ", ms$n_unmatched_ftir,
              " unmatched ", label, ", ", ms$n_unmatched_raman, " unmatched Raman")
  agreement <- analyze_agreement(match, config)

  list(label = label, prefiltered = prefiltered, skipped = FALSE, clean = clean,
       for_align = for_align, for_match = for_match,
       aligned = aligned, aligned_all = aligned_all,
       landmark = landmark, ransac = ransac, method = method, icp = icp,
       match = match, agreement = agreement,
       norm = list(ftir_centroid = nm$centroid, ftir_scale = nm$scale,
                   raman_centroid = ref$centroid, raman_scale = ref$scale))
}


#' Score a Tier 2 candidate transform: one-to-one pairs between the
#' transformed source anchors and the reference anchors within
#' ransac_inlier_dist_um. The aligners' own inlier counts are not comparable
#' (RANSAC counts many-to-one nearest neighbours), so the race between RANSAC
#' and global registration is decided on this common score; ties keep RANSAC.
#' @return Integer (-1 when the transform is missing)
tier2_pairs <- function(M, src, ref, config) {
  if (is.null(M)) return(-1L)
  n <- align_transform_inliers(M, src$x_norm, src$y_norm, ref$x_norm, ref$y_norm,
                               config$ransac_inlier_dist_um %||% 200)
  if (is.na(n)) -1L else n
}


#' Match two instruments that are both aligned into the Raman frame.
#'
#' The matcher reads the reference's x_norm / y_norm, so the reference's
#' aligned coordinates are put there. The acceptance gate follows from the
#' labels (a pair involving LDIR uses the LDIR gate).
#'
#' @param src_aligned,ref_aligned Particle tables with x_aligned / y_aligned
#' @param src_label,ref_label Labels ("ftir_bruker", "ftir_perkin", "ldir", ...)
#' @return match_particles() result, or NULL when either side is empty
match_aligned_pair <- function(src_aligned, ref_aligned, config, src_label, ref_label) {
  ok <- function(d) is.data.frame(d) && nrow(d) > 0 &&
    all(c("x_aligned", "y_aligned") %in% names(d))
  if (!ok(src_aligned) || !ok(ref_aligned)) return(NULL)
  ref_aligned$x_norm <- ref_aligned$x_aligned
  ref_aligned$y_norm <- ref_aligned$y_aligned
  match_particles(src_aligned, ref_aligned, config,
                  src_label = src_label, ref_label = ref_label)
}


#' Three-way triplets: FTIR-Raman and LDIR-Raman pairs sharing a Raman particle.
#'
#' Each triplet's three pairwise family comparisons use the same rules as the
#' pairwise agreement: a side below its instrument's agreement quality cutoff
#' (ftir_agreement_min_quality, raman_hqi_threshold, ldir_hqi_unknown_threshold)
#' or without an identification is not scored, and "Unknown" / "Other polymer"
#' never agree. A triplet with fewer than two scored sides has no consensus.
#'
#' @param ftir_match match_particles() result FTIR -> Raman
#' @param ldir_match match_particles() result LDIR -> Raman
#' @param raman_for_match Raman particles the matches refer to
#' @param config Pipeline config
#' @param ldir_ftir_match Optional match_particles() result LDIR -> FTIR. When
#'   given, `direct_ldir_ftir_match` says whether the triplet's LDIR and FTIR
#'   particles are also each other's direct match (FALSE: the triplet links
#'   them only through the Raman particle, so one of the two Raman pairs may
#'   be wrong; NA: the LDIR particle has no direct FTIR match).
#' @return data.frame, one row per triplet (zero rows when there are none)
build_triplets <- function(ftir_match, ldir_match, raman_for_match, config,
                           ldir_ftir_match = NULL) {
  fm <- ftir_match$matched; lm <- ldir_match$matched
  if (is.null(fm) || is.null(lm) || nrow(fm) == 0 || nrow(lm) == 0) return(data.frame())
  col <- function(d, k) if (k %in% names(d)) d[[k]] else rep(NA, nrow(d))
  triplets <- merge(
    data.frame(raman_particle_id = fm$raman_particle_id,
               ftir_particle_id = fm$ftir_particle_id,
               ftir_material = fm$ftir_material,
               ftir_quality = col(fm, "ftir_quality"),
               ftir_raman_distance = fm$match_distance, stringsAsFactors = FALSE),
    data.frame(raman_particle_id = lm$raman_particle_id,
               ldir_particle_id = lm$ldir_particle_id,
               ldir_material = lm$ldir_material,
               ldir_quality = col(lm, "ldir_quality"),
               ldir_raman_distance = lm$match_distance, stringsAsFactors = FALSE),
    by = "raman_particle_id")
  if (nrow(triplets) == 0) {
    log_message("  Three-way triplets: 0 particles matched")
    return(triplets)
  }
  dm <- ldir_ftir_match$matched
  if (!is.null(dm) && all(c("ldir_particle_id", "ftir_particle_id") %in% names(dm))) {
    partner <- as.character(dm$ftir_particle_id)[match(triplets$ldir_particle_id,
                                                       dm$ldir_particle_id)]
    triplets$direct_ldir_ftir_match <- ifelse(is.na(partner), NA,
                                              partner == as.character(triplets$ftir_particle_id))
  }
  ri <- match(triplets$raman_particle_id, raman_for_match$particle_id)
  triplets$raman_material <- raman_for_match$material[ri]
  triplets$raman_quality  <- raman_for_match$quality[ri]

  scored <- function(prefix) {
    q <- suppressWarnings(as.numeric(triplets[[paste0(prefix, "_quality")]]))
    cut <- .agreement_quality_cutoff(prefix, config)
    q_ok <- if (cut <= 0 || all(is.na(q))) rep(TRUE, length(q)) else !is.na(q) & q >= cut
    q_ok & !is_unidentified(triplets[[paste0(prefix, "_material")]])
  }
  s_f <- scored("ftir"); s_r <- scored("raman"); s_l <- scored("ldir")

  ftir_fam  <- classify_family_vec(triplets$ftir_material)
  raman_fam <- classify_family_vec(triplets$raman_material)
  ldir_fam  <- classify_family_vec(triplets$ldir_material)
  fr_agree <- s_f & s_r & family_agrees(ftir_fam, raman_fam)
  lr_agree <- s_l & s_r & family_agrees(ldir_fam, raman_fam)
  fl_agree <- s_f & s_l & family_agrees(ftir_fam, ldir_fam)
  n_agree <- as.integer(fr_agree) + as.integer(lr_agree) + as.integer(fl_agree)
  n_scored <- as.integer(s_f) + as.integer(s_r) + as.integer(s_l)

  triplets$n_instruments_scored   <- n_scored
  triplets$n_instrument_agreement <- n_agree   # 0-3 pairwise family agreements
  triplets$ftir_family  <- ftir_fam
  triplets$raman_family <- raman_fam
  triplets$ldir_family  <- ldir_fam
  triplets$material_consensus <- ifelse(
    n_scored < 2, "Not scored (quality / no identification)",
    ifelse(n_agree == 3, "Full agreement",
    ifelse(fr_agree, "FTIR+Raman agree",
    ifelse(lr_agree, "LDIR+Raman agree",
    ifelse(fl_agree, "FTIR+LDIR agree", "No agreement")))))

  n_full <- sum(n_agree == 3)
  n_partial <- sum(n_agree >= 1 & n_agree < 3)
  n_unscored <- sum(n_scored < 2)
  log_message("  Three-way triplets: ", nrow(triplets), " total (", n_full,
              " full agreement, ", n_partial, " partial, ", n_unscored,
              " not scored, ", nrow(triplets) - n_full - n_partial - n_unscored,
              " no agreement)")
  if ("direct_ldir_ftir_match" %in% names(triplets)) {
    n_off <- sum(triplets$direct_ldir_ftir_match %in% FALSE)
    log_message("  Triplets whose LDIR and FTIR particles are not each other's direct ",
                "match: ", n_off, " of ", nrow(triplets),
                level = if (n_off > 0) "WARN" else "INFO")
  }
  triplets
}


#' Write an FTIR-type -> Raman transform to a text file.
#'
#' @param pair Result of align_and_match()
#' @param path Output file
#' @param src_name Name in the header ("FTIR_perkin", "FTIR_bruker")
#' @param scan_bounds Optional list(x_min, x_max, y_min, y_max) of the
#'   instrument's background image
write_transform_params <- function(pair, path, src_name, scan_bounds = NULL) {
  p <- pair$icp$params; nr <- pair$norm; M <- pair$icp$transform
  lines <- c(
    paste0("# ", src_name, "-to-Raman Transform Parameters"),
    paste0("# Generated: ", Sys.time()),
    "",
    paste0("scale:        ", round(p$scale, 6)),
    paste0("rotation_deg: ", round(p$rotation_deg, 4)),
    paste0("tx:           ", round(p$tx, 4)),
    paste0("ty:           ", round(p$ty, 4)),
    paste0("reflected:    ", p$reflected),
    "",
    "# Normalization parameters (applied before transform)",
    paste0("ftir_centroid_x:  ", round(nr$ftir_centroid[1], 4)),
    paste0("ftir_centroid_y:  ", round(nr$ftir_centroid[2], 4)),
    paste0("raman_centroid_x: ", round(nr$raman_centroid[1], 4)),
    paste0("raman_centroid_y: ", round(nr$raman_centroid[2], 4)),
    paste0("ftir_scale:       ", round(nr$ftir_scale, 6)),
    paste0("raman_scale:      ", round(nr$raman_scale, 6)),
    "",
    "# ICP refinement info",
    paste0("icp_converged:    ", pair$icp$converged),
    paste0("icp_iterations:   ", pair$icp$n_iterations),
    paste0("icp_final_rms:    ", round(tail(pair$icp$rms_history, 1), 4), " um"),
    "")
  if (!is.null(scan_bounds))
    lines <- c(lines,
      paste0("# ", src_name, " scan area (physical extent of the background image)"),
      paste0("ftir_scan_xmin:   ", scan_bounds$x_min),
      paste0("ftir_scan_xmax:   ", scan_bounds$x_max),
      paste0("ftir_scan_ymin:   ", scan_bounds$y_min),
      paste0("ftir_scan_ymax:   ", scan_bounds$y_max),
      "")
  lines <- c(lines,
    paste0("# 3x3 Transform matrix (homogeneous, ", src_name, "_norm -> Raman_norm)"),
    paste0("matrix_row1: ", paste(round(M[1, ], 8), collapse = ", ")),
    paste0("matrix_row2: ", paste(round(M[2, ], 8), collapse = ", ")),
    paste0("matrix_row3: ", paste(round(M[3, ], 8), collapse = ", ")))
  writeLines(lines, path)
  log_message("  Wrote ", basename(path))
  invisible(path)
}
