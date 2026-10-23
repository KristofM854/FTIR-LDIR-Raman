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


#' Instrument roles, in the order the reference frame is chosen from.
#' LDIR is never the reference: its coordinates usually come from its image.
INSTRUMENT_ROLES <- c("raman", "ftir", "ftir_bruker", "ldir")
REFERENCE_ORDER  <- c("raman", "ftir", "ftir_bruker")

#' Display name, file name and match-column prefix of each role. The file name
#' is used in output file names (matched_<a>_<b>.csv); the prefix is the
#' column prefix of the role in a match table with the reference.
ROLE_NAME   <- c(raman = "Raman", ftir = "FTIR (PerkinElmer)",
                 ftir_bruker = "FTIR (Bruker)", ldir = "LDIR")
ROLE_FILE   <- c(raman = "raman", ftir = "ftir_perkin",
                 ftir_bruker = "ftir_bruker", ldir = "ldir")

#' The instruments a config provides (input path set), in INSTRUMENT_ROLES order.
present_roles <- function(config) {
  key <- c(raman = "raman_path", ftir = "ftir_path",
           ftir_bruker = "ftir_bruker_path", ldir = "ldir_path")
  has <- vapply(INSTRUMENT_ROLES, function(r) {
    v <- config[[key[[r]]]]
    !is.null(v) && length(v) == 1 && !is.na(v) && nzchar(v)
  }, logical(1))
  INSTRUMENT_ROLES[has]
}

#' The reference frame of a run: Raman if present, else FTIR (PerkinElmer),
#' else FTIR (Bruker). Needs at least two instruments.
choose_reference <- function(roles) {
  if (length(roles) < 2)
    stop("At least two instruments are needed; got ",
         if (length(roles)) paste(ROLE_NAME[roles], collapse = ", ") else "none",
         ". Set two or more of raman_path, ftir_path, ftir_bruker_path, ldir_path.",
         call. = FALSE)
  ref <- REFERENCE_ORDER[REFERENCE_ORDER %in% roles][1]
  if (is.na(ref))
    stop("LDIR alone cannot be paired: add Raman or FTIR data.", call. = FALSE)
  ref
}

#' Column prefixes and output file stem of a pair source -> reference.
#'
#' Each output file always has the same column prefixes, whichever instrument
#' is the run's reference, so the viewer reads every file one way. The
#' prefixes are those of the original Raman-reference runs (FTIR types are
#' "ftir" against Raman, LDIR -> FTIR (PerkinElmer) uses "ldir" / "ftir", the
#' FTIR (Bruker) <-> FTIR (PerkinElmer) file "ftir_bruker" / "ftir_perkin").
#' @return list(src_prefix, ref_prefix, stem, src_file, ref_file)
pair_labels <- function(src_role, ref_role) {
  stem <- paste0(ROLE_FILE[[src_role]], "_", ROLE_FILE[[ref_role]])
  pfx <- switch(paste(src_role, ref_role),
    "ftir raman"        = c("ftir", "raman"),
    "ftir_bruker raman" = c("ftir", "raman"),
    "ldir raman"        = c("ldir", "raman"),
    "ldir ftir"         = c("ldir", "ftir"),
    "ftir_bruker ftir"  = c("ftir_bruker", "ftir_perkin"),
    c(ROLE_FILE[[src_role]], ROLE_FILE[[ref_role]]))
  list(src_prefix = pfx[1], ref_prefix = pfx[2], stem = stem,
       src_file = ROLE_FILE[[src_role]], ref_file = ROLE_FILE[[ref_role]])
}


#' Reference particle sets in the shared (reference-normalised) frame.
#'
#' The reference instrument's coordinates minus its full-cloud centroid are
#' the frame every other instrument is aligned into.
#'
#' @param raw Canonical particle table of the reference instrument
#' @param config Pipeline config
#' @param role Reference role: "raman", "ftir" or "ftir_bruker"
#' @return list with
#'   prefiltered    after the coordinate prefilter, before normalisation
#'   clean          the same, normalised (x_norm, y_norm): the full cloud
#'   norm_align     material (and for Raman size / HQI) anchors for the
#'                  coarse aligners
#'   for_transform  the landmark and ICP target (Raman: particles >=
#'                  align_raman_min_size_um; FTIR: all)
#'   for_match      the particles matched against (>= min_particle_size_um)
#'   centroid, scale  the frame's normalisation
#'   role, name     the reference instrument
#'   An FTIR reference's sets also carry x_aligned / y_aligned (= x_norm /
#'   y_norm), the columns an FTIR table has when it is aligned to Raman, so
#'   its match tables have the same columns whichever instrument is the
#'   reference.
prepare_reference <- function(raw, config, role = "raman") {
  name <- ROLE_NAME[[role]]
  is_raman <- identical(role, "raman")
  anchor_min <- config$align_min_anchor_count %||% 4
  prefilter <- if (is_raman) prefilter_raman else prefilter_ftir
  clean <- prefilter(raw, 0, 0)

  # Alignment anchors: material, then (Raman) size, then HQI. Each filter is
  # honoured only while it leaves enough anchors: three filters that are
  # individually reasonable can intersect to nothing.
  materials <- if (is_raman) config$align_raman_materials else config$align_ftir_materials
  for_align <- select_alignment_anchors(clean, materials, role, config, name)
  min_size <- if (is_raman) config$align_raman_min_size_um else NULL
  if (!is.null(min_size) && min_size > 0 && any(!is.na(for_align$feret_max_um))) {
    size_mask <- is.na(for_align$feret_max_um) | for_align$feret_max_um >= min_size
    if (sum(size_mask) >= anchor_min) {
      n_before <- nrow(for_align)
      for_align <- for_align[size_mask, ]
      log_message(name, " alignment: removed ", n_before - nrow(for_align),
                  " particles < ", min_size, " um")
    } else {
      log_message(name, " alignment: size filter (>= ", min_size, " um) would ",
                  "leave only ", sum(size_mask), " anchors \u2014 skipping it",
                  level = "WARN")
    }
  }
  # Only reliably identified particles should drive material-based alignment
  hqi <- if (is_raman) config$raman_hqi_threshold else NULL
  if (!is.null(hqi) && hqi > 0 && any(!is.na(for_align$quality))) {
    hqi_mask <- !is.na(for_align$quality) & for_align$quality >= hqi
    if (sum(hqi_mask) >= anchor_min) {
      n_before <- nrow(for_align)
      for_align <- for_align[hqi_mask, ]
      log_message(name, " alignment: HQI filter (>= ", hqi, "): kept ",
                  nrow(for_align), " of ", n_before)
    } else {
      log_message(name, " alignment: HQI filter (>= ", hqi, ") would leave only ",
                  sum(hqi_mask), " anchors \u2014 skipping it", level = "WARN")
    }
  }
  if (is_raman)
    log_message("Raman alignment target particles: ", nrow(for_align),
                " (material + >= ", min_size, " um + HQI >= ", hqi, ")")
  else
    log_message(name, " alignment target particles: ", nrow(for_align))

  # The particles matched against: those passing the instrument's quality
  # cutoff (match_quality_filter), the same set the agreement scores.
  for_match <- prefilter(raw, if (is_raman) 0 else config$ftir_quality_threshold,
                         config$min_particle_size_um)
  for_match <- .match_quality_subset(for_match, role, config, name)
  log_message(name, " for matching: ", nrow(for_match), " particles",
              if (isTRUE(cfg_value(config, "match_quality_filter")))
                " (passing the quality cutoff)" else "")

  # The frame: centroid of the FULL cleaned cloud, not of an anchor subset (a
  # subset centroid depends on which particles each instrument detected, and
  # would start the aligners from a biased translation).
  log_message("Normalizing coordinates (", name, " reference frame)")
  nm <- cloud_normalization(clean, config$normalize_scale, name)
  norm <- function(df) {
    df <- apply_normalization(df, nm$centroid, nm$scale)
    if (!is_raman) { df$x_aligned <- df$x_norm; df$y_aligned <- df$y_norm }
    df
  }
  clean_n <- norm(clean)

  # Landmark / ICP set: reference particles large enough to be visible to
  # every instrument; only size matters for geometry.
  for_transform <- clean_n
  if (!is.null(min_size) && min_size > 0 && any(!is.na(clean_n$feret_max_um)))
    for_transform <- clean_n[is.na(clean_n$feret_max_um) | clean_n$feret_max_um >= min_size, ]
  if (is_raman)
    log_message("Raman for spatial transform (>= ", min_size, " um): ",
                nrow(for_transform), " particles")

  list(prefiltered = clean, clean = clean_n, norm_align = norm(for_align),
       for_align = for_align, for_transform = for_transform,
       for_match = norm(for_match), centroid = nm$centroid, scale = nm$scale,
       role = role, name = name)
}

#' The Raman reference frame (prepare_reference() with role "raman").
prepare_raman_reference <- function(raman_raw, config)
  prepare_reference(raman_raw, config, "raman")


#' Align an FTIR-type particle table to the reference frame, match and score.
#'
#' Tiered alignment, then matching and agreement:
#'   Tier 1  landmark alignment on large particles and fibres; used when
#'           confident (and landmark_skip_full_ransac is TRUE)
#'   Tier 2  coarse RANSAC and global registration on the material anchors;
#'           the one pairing more particles wins
#'   ICP     always, to polish the transform
#'
#' @param src_raw Canonical particle table (ftir or ftir_bruker role)
#' @param ref Result of prepare_reference()
#' @param config Pipeline config
#' @param label Instrument label for messages ("FTIR (PerkinElmer)")
#' @param src_role Role of the source ("ftir", "ftir_bruker"); sets the match
#'   column prefixes and file stem with pair_labels()
#' @return list: prefiltered, clean, for_align, for_match, aligned, aligned_all,
#'   landmark, ransac, method, icp, match, agreement, norm, labels
#'   (pair_labels()). With too few particles to align (warned) only label,
#'   prefiltered and skipped = TRUE.
align_and_match <- function(src_raw, ref, config, label, src_role = "ftir") {
  ref_role <- ref$role %||% "raman"
  ref_name <- ref$name %||% "Raman"
  labels <- pair_labels(src_role, ref_role)
  log_message(strrep("-", 50))
  log_message(label, " \u2194 ", ref_name, ": alignment and matching")
  anchor_min <- config$align_min_anchor_count %||% 4

  prefiltered <- prefilter_ftir(src_raw, min_quality = 0, min_size_um = 0)
  # Material anchors for the coarse aligners; select_material_anchors() falls
  # back to the full cloud when too few carry an anchor material.
  for_align <- select_alignment_anchors(prefiltered, config$align_ftir_materials,
                                        src_role, config, label)
  for_match <- prefilter_ftir(src_raw, min_quality = config$ftir_quality_threshold,
                              min_size_um = config$min_particle_size_um)
  for_match <- .match_quality_subset(for_match, src_role, config, label)
  if (nrow(for_align) < 3) {
    log_message(label, ": too few particles (", nrow(for_align),
                ") for alignment \u2014 skipping", level = "WARN")
    return(list(label = label, prefiltered = prefiltered, skipped = TRUE,
                labels = labels))
  }

  # Normalisation: the source's full-cloud centroid; anchor and match subsets
  # use the same one so they stay in one frame.
  nm <- cloud_normalization(prefiltered, config$normalize_scale, label)
  norm <- function(df) apply_normalization(df, nm$centroid, nm$scale)
  clean <- norm(prefiltered)
  norm_align <- norm(for_align)
  for_match <- norm(for_match)

  # Tier 1: landmarks
  landmark <- landmark_align(clean, ref$for_transform, config, src_label = label,
                             ref_label = ref_name)
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

  # Is the pose fixed by the data? Landmarks accepted above chance are; a
  # Tier 2 pose (after ICP) must pair more anchors within the match gate than
  # the search would find by chance (pose_significance()). A short particle
  # list (e.g. one already filtered to plastics) against a dense reference
  # cannot fix a pose.
  if (identical(method, "landmark")) {
    determined <- TRUE; pose_p <- landmark$p_chance %||% NA_real_
  } else {
    r_pose <- config$match_dist_threshold_um
    n_pose <- align_transform_inliers(icp$transform, norm_align$x_norm, norm_align$y_norm,
                                      ref$norm_align$x_norm, ref$norm_align$y_norm, r_pose)
    if (is.na(n_pose)) n_pose <- 0L
    ps <- pose_significance(n_pose, nrow(norm_align), ref$norm_align$x_norm,
                            ref$norm_align$y_norm, r_pose)
    pose_p <- ps$p
    determined <- isTRUE(pose_p < cfg_value(config, "align_determined_p"))
    log_message(sprintf(paste0("  Pose check: %d anchor pairs of %d anchors at %g um; ",
                               "%.2f per pose by chance; p over the search = %.2g"),
                        n_pose, nrow(norm_align), r_pose, ps$expected, pose_p),
                level = if (determined) "INFO" else "WARN")
    if (!determined)
      log_message("  ", label, " \u2192 ", ref_name, ": ALIGNMENT NOT DETERMINED (too few ",
                  "anchors for this reference density). Its pairs are written for ",
                  "inspection but not scored, and pairs through it are skipped.",
                  level = "WARN")
  }

  aligned     <- apply_ftir_transform(for_match, icp$transform)
  aligned_all <- apply_ftir_transform(clean, icp$transform)

  match <- match_particles(aligned, ref$for_match, config,
                           src_label = labels$src_prefix, ref_label = labels$ref_prefix)
  ms <- match$match_stats
  log_message(label, " \u2194 ", ref_name, ": ", ms$n_matched, " matched, ",
              ms$n_unmatched_src, " unmatched ", label, ", ", ms$n_unmatched_ref,
              " unmatched ", ref_name)
  scored <- match
  if (!determined) scored$matched <- match$matched[0, , drop = FALSE]
  agreement <- analyze_agreement(scored, config,
                                 instrument_a = .agreement_name(labels$src_prefix),
                                 instrument_b = .agreement_name(labels$ref_prefix),
                                 prefix_a = labels$src_prefix, prefix_b = labels$ref_prefix)

  list(label = label, prefiltered = prefiltered, skipped = FALSE, clean = clean,
       for_align = for_align, for_match = for_match,
       aligned = aligned, aligned_all = aligned_all,
       landmark = landmark, ransac = ransac, method = method, icp = icp,
       match = match, agreement = agreement, labels = labels,
       alignment_determined = determined, pose_p = pose_p,
       norm = list(ftir_centroid = nm$centroid, ftir_scale = nm$scale,
                   raman_centroid = ref$centroid, raman_scale = ref$scale))
}


#' Instrument name in agreement tables and plots for a match-column prefix
#' ("ftir" -> "FTIR", "ftir_bruker" -> "FTIR (Bruker)").
.agreement_name <- function(prefix)
  switch(prefix, ftir = "FTIR", raman = "Raman", ldir = "LDIR",
         ftir_perkin = "FTIR (PerkinElmer)", ftir_bruker = "FTIR (Bruker)", prefix)


#' The particles that take part in matching: with match_quality_filter, those
#' passing the instrument's agreement quality cutoff (quality_passes()).
.match_quality_subset <- function(df, prefix, config, label = prefix, quiet = FALSE) {
  if (!isTRUE(cfg_value(config, "match_quality_filter")) || is.null(df) || nrow(df) == 0)
    return(df)
  keep <- quality_passes(df, prefix, config)
  if (!quiet && !all(keep))
    log_message("  ", label, " for matching: ", sum(keep), " of ", nrow(df),
                " particles pass the quality cutoff (", .agreement_quality_cutoff(prefix, config),
                ")")
  df[keep, , drop = FALSE]
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
#' labels; a pair involving LDIR uses the wider of the LDIR and base gates.
#'
#' @param src_aligned,ref_aligned Particle tables with x_aligned / y_aligned
#' @param src_label,ref_label Labels ("ftir_bruker", "ftir_perkin", "ldir", ...)
#' @return match_particles() result, or NULL when either side is empty
match_aligned_pair <- function(src_aligned, ref_aligned, config, src_label, ref_label) {
  ok <- function(d) is.data.frame(d) && nrow(d) > 0 &&
    all(c("x_aligned", "y_aligned") %in% names(d))
  if (!ok(src_aligned) || !ok(ref_aligned)) return(NULL)
  src_aligned <- .match_quality_subset(src_aligned, src_label, config, src_label, quiet = TRUE)
  ref_aligned <- .match_quality_subset(ref_aligned, ref_label, config, ref_label, quiet = TRUE)
  if (!ok(src_aligned) || !ok(ref_aligned)) return(NULL)
  ref_aligned$x_norm <- ref_aligned$x_aligned
  ref_aligned$y_norm <- ref_aligned$y_aligned
  # Both sides carry their own alignment error, so the pair takes the wider of
  # the two gates (the LDIR gate is set for LDIR against the reference only;
  # dataset A: FTIR (Bruker) sigma 12 um, LDIR 7 um).
  if (!is.null(config$match_dist_threshold_ldir_um))
    config$match_dist_threshold_ldir_um <- max(config$match_dist_threshold_ldir_um,
                                               config$match_dist_threshold_um)
  config$match_gate_from_sigma <- FALSE    # the per-run LDIR gate is for LDIR vs the reference
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


#' Write a source -> reference transform to a text file.
#'
#' The ftir_* / raman_* keys are the source's and the reference's
#' normalisation, whatever the instruments (the viewer reads these keys).
#'
#' @param pair Result of align_and_match() (or list(icp, norm))
#' @param path Output file
#' @param src_name Name in the header ("FTIR_perkin", "FTIR_bruker")
#' @param ref_name Name of the reference in the header
#' @param scan_bounds Optional list(x_min, x_max, y_min, y_max) of the
#'   instrument's background image
write_transform_params <- function(pair, path, src_name, scan_bounds = NULL,
                                   ref_name = "Raman") {
  p <- pair$icp$params; nr <- pair$norm; M <- pair$icp$transform
  lines <- c(
    paste0("# ", src_name, "-to-", ref_name, " Transform Parameters"),
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
    paste0("# 3x3 Transform matrix (homogeneous, ", src_name, "_norm -> ", ref_name, "_norm)"),
    paste0("matrix_row1: ", paste(round(M[1, ], 8), collapse = ", ")),
    paste0("matrix_row2: ", paste(round(M[2, ], 8), collapse = ", ")),
    paste0("matrix_row3: ", paste(round(M[3, ], 8), collapse = ", ")))
  writeLines(lines, path)
  log_message("  Wrote ", basename(path))
  invisible(path)
}
