# =============================================================================
# detection_rates.R — how often one instrument finds what another found
# =============================================================================
# For every matched instrument pair, in both directions: of the particles the
# base instrument found (passing its quality cutoff, >= min_particle_size_um,
# and inside the area the other instrument covered), the share the other
# instrument also has at the same place, per size class of the base particle,
# in two senses:
#   identified  paired with a particle that passes the other instrument's
#               quality cutoff (only those are matched): detected AND
#               confidently identified by the other instrument;
#   detected    any particle of the other instrument, of any quality and size,
#               within the pair's distance gate. On a dense field part of that
#               is chance; the chance-corrected rate is
#               (detected - expected by chance) / (n - expected by chance).
#
# Each rate has a 95% Wilson interval; the identified rate also has a
# chance-corrected value: the pairs weighted by p_true (pair_confidence()).
# =============================================================================

DETECTION_SIZE_BREAKS <- c(20, 50, 100, 200, Inf)
DETECTION_SIZE_LABELS <- c("20-50", "50-100", "100-200", ">200")

#' Instrument roles of a pair stem ("ftir_bruker_ldir" -> c("ftir_bruker", "ldir"))
.pair_roles <- function(stem) {
  roles <- c("ftir_perkin", "ftir_bruker", "raman", "ldir")
  for (a in roles) if (startsWith(stem, paste0(a, "_"))) {
    b <- sub(paste0("^", a, "_"), "", stem)
    if (b %in% roles) return(c(a, b))
  }
  NULL
}

#' Detection rates for every pair of the run.
#'
#' @param aligned Named list, instrument role ("ftir_perkin", "ftir_bruker",
#'   "raman", "ldir") -> all its particles in the shared frame (x_aligned,
#'   y_aligned, particle_id, feret_max_um, quality, material)
#' @param matches Named list, pair stem ("ldir_raman", ...) -> matched table;
#'   its first two *_particle_id columns are the stem's first and second
#'   instrument, an optional p_true column the pair confidence
#' @param config Pipeline config (quality cutoffs, min_particle_size_um)
#' @param gates Named list, pair stem -> the pair's distance gate (um); the
#'   detected rate is computed for the pairs that have one
#' @return data.frame: found_by, of, group ("all" / "plastics"), size_class
#'   (incl. "all sizes"), n_base; identified: n_found, rate_pct, ci_low_pct,
#'   ci_high_pct, n_found_true (sum of p_true), rate_true_pct; detected:
#'   gate_um, n_detected, detected_pct, det_ci_low_pct, det_ci_high_pct,
#'   det_chance_expected, detected_corrected_pct
detection_rates <- function(aligned, matches, config, gates = list()) {
  adapt <- cfg_value(config, "match_adaptive_dist_factor")
  min_size <- cfg_value(config, "min_particle_size_um")
  plastic_fams <- c(synthetic_families, semi_synthetic_families)
  rows <- list()
  for (stem in names(matches)) {
    roles <- .pair_roles(stem)
    m <- matches[[stem]]
    if (is.null(roles) || is.null(m) || !all(roles %in% names(aligned))) next
    id_cols <- grep("_particle_id$", names(m), value = TRUE)
    if (length(id_cols) < 2) next
    for (dir in 1:2) {
      base_role <- roles[dir]; other_role <- roles[3 - dir]
      base <- aligned[[base_role]]; other <- aligned[[other_role]]
      if (is.null(base) || is.null(other) || nrow(base) == 0 || nrow(other) < 3) next
      prefix <- if (base_role == "ftir_perkin") "ftir" else base_role
      # (a particle without a size is matched, so it is in the base as well;
      # it has no size class and counts in "all sizes" only)
      keep <- quality_passes(base, prefix, config) &
        (is.na(base$feret_max_um) | base$feret_max_um >= min_size) &
        is.finite(base$x_aligned) & is.finite(base$y_aligned)
      # Inside the area the other instrument covered (its convex hull)
      ox <- other$x_aligned; oy <- other$y_aligned
      fin <- is.finite(ox) & is.finite(oy)
      h <- grDevices::chull(ox[fin], oy[fin])
      inside <- .points_in_polygon(base$x_aligned, base$y_aligned, ox[fin][h], oy[fin][h])
      # A particle the other instrument matched was within its reach, even when
      # it lies just outside the hull of the other's particles
      ids <- as.character(m[[id_cols[dir]]])
      inside <- inside | as.character(base$particle_id) %in% ids
      b <- base[keep & inside, , drop = FALSE]
      if (nrow(b) == 0) next

      # Detected: any particle of the other instrument within the pair's gate
      # (widened for large particles as in matching), with its chance level
      gate <- gates[[stem]]
      det <- p_chance <- rep(NA, nrow(b))
      if (!is.null(gate) && is.finite(gate) && sum(fin) >= 3) {
        g_i <- pmax(gate, adapt * ifelse(is.finite(b$major_um %||% NA), b$major_um, 0))
        d_nn <- RANN::nn2(cbind(ox[fin], oy[fin]), cbind(b$x_aligned, b$y_aligned), k = 1)$nn.dists[, 1]
        det <- d_nn <= g_i
        rho <- chance_partner_rate(ox[fin], oy[fin], 1)$density
        p_chance <- 1 - exp(-rho * pi * g_i^2)
      }

      hit <- match(as.character(b$particle_id), ids)
      p_true <- if (!is.null(m$p_true)) m$p_true[hit] else rep(NA_real_, nrow(b))
      found <- !is.na(hit)
      # A paired particle is detected (its partner may lie in a gate widened for
      # the partner's size rather than this particle's)
      if (!all(is.na(det))) det <- det | found
      fam <- b$material_family %||% classify_family_vec(b$material)
      size_class <- as.character(cut(b$feret_max_um, DETECTION_SIZE_BREAKS,
                                     DETECTION_SIZE_LABELS, right = FALSE))
      for (grp in c("all", "plastics")) {
        g <- if (grp == "all") rep(TRUE, nrow(b)) else fam %in% plastic_fams
        for (sc in c(DETECTION_SIZE_LABELS, "all sizes")) {
          # A particle below the first size class (min_particle_size_um < 20)
          # has no class: it counts in "all sizes" only
          s <- g & (sc == "all sizes" | (!is.na(size_class) & size_class == sc))
          n <- sum(s); k <- sum(found & s)
          ci <- wilson_interval(k, n)
          k_true <- if (all(is.na(p_true[found & s])) && k > 0) NA_real_
                    else sum(p_true[found & s], na.rm = TRUE)
          rows[[length(rows) + 1]] <- data.frame(
            found_by = other_role, of = base_role, group = grp, size_class = sc,
            n_base = n, n_found = k,
            rate_pct = if (n > 0) round(100 * k / n, 1) else NA_real_,
            ci_low_pct = round(100 * ci$low, 1), ci_high_pct = round(100 * ci$high, 1),
            n_found_true = round(k_true, 1),
            rate_true_pct = if (n > 0) round(100 * k_true / n, 1) else NA_real_,
            .detected_row(det[s], p_chance[s], gate),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Partner status of every particle, every pair both ways.
#'
#' The matching pairs only particles that pass both instruments' quality
#' cutoffs, so "unmatched" alone does not say whether the other instrument
#' found nothing at that place or found something it could not identify
#' confidently. For every particle of the base instrument (any quality,
#' >= min_particle_size_um) this gives the other instrument's nearest particle
#' within the pair's distance gate (any quality) and a status:
#'   matched               paired (both pass their cutoffs)
#'   partner_below_cutoff  the other found a particle there, below its cutoff
#'   this_below_cutoff     this particle is below its own cutoff; the other
#'                         identified a particle there
#'   both_below_cutoff     both found a particle there, neither identified it
#'   paired_elsewhere      both pass, but the other's particle is paired with
#'                         another particle
#'   partner_not_matchable both pass their quality cutoff, but the other's
#'                         particle is in no pair: it is below the size limit
#'                         (min_particle_size_um) or the two sizes differ by
#'                         more than match_max_size_ratio
#'   not_found             nothing of the other within the gate
#'   outside_area          outside the area the other instrument covered
#' with p_chance, the chance of an unrelated particle of the other within the
#' gate (its density): on a dense field part of the "found" statuses is chance.
#'
#' @inheritParams detection_rates
#' @return data.frame(of, found_by, particle_id, passes, status, partner_id,
#'   partner_material, partner_quality, partner_passes, partner_distance_um,
#'   gate_um, p_chance)
partner_status <- function(aligned, matches, config, gates = list()) {
  adapt <- cfg_value(config, "match_adaptive_dist_factor")
  min_size <- cfg_value(config, "min_particle_size_um")
  pfx <- function(role) if (role == "ftir_perkin") "ftir" else role
  rows <- list()
  for (stem in names(matches)) {
    roles <- .pair_roles(stem); m <- matches[[stem]]; gate <- gates[[stem]]
    if (is.null(roles) || is.null(m) || !all(roles %in% names(aligned)) ||
        is.null(gate) || !is.finite(gate)) next
    id_cols <- grep("_particle_id$", names(m), value = TRUE)
    if (length(id_cols) < 2) next
    for (dir in 1:2) {
      base <- aligned[[roles[dir]]]; other <- aligned[[roles[3 - dir]]]
      if (is.null(base) || is.null(other) || nrow(base) == 0 || nrow(other) < 3) next
      b <- base[(is.na(base$feret_max_um) | base$feret_max_um >= min_size) &
                is.finite(base$x_aligned) & is.finite(base$y_aligned), , drop = FALSE]
      o <- other[is.finite(other$x_aligned) & is.finite(other$y_aligned), , drop = FALSE]
      if (nrow(b) == 0 || nrow(o) < 3) next
      b_pass <- quality_passes(b, pfx(roles[dir]), config)
      o_pass <- quality_passes(o, pfx(roles[3 - dir]), config)
      h <- grDevices::chull(o$x_aligned, o$y_aligned)
      ids_b <- as.character(m[[id_cols[dir]]]); ids_o <- as.character(m[[id_cols[3 - dir]]])
      inside <- .points_in_polygon(b$x_aligned, b$y_aligned, o$x_aligned[h], o$y_aligned[h]) |
        as.character(b$particle_id) %in% ids_b
      g_i <- pmax(gate, adapt * ifelse(is.finite(b$major_um %||% NA), b$major_um, 0))
      nn <- RANN::nn2(cbind(o$x_aligned, o$y_aligned), cbind(b$x_aligned, b$y_aligned), k = 1)
      j <- nn$nn.idx[, 1]; d <- nn$nn.dists[, 1]
      # A matched particle's partner is its pair partner, not the nearest one
      hit <- match(as.character(b$particle_id), ids_b)
      jm <- match(ids_o[hit], as.character(o$particle_id))
      j[!is.na(jm)] <- jm[!is.na(jm)]
      d[!is.na(jm)] <- sqrt((o$x_aligned[jm[!is.na(jm)]] - b$x_aligned[!is.na(jm)])^2 +
                            (o$y_aligned[jm[!is.na(jm)]] - b$y_aligned[!is.na(jm)])^2)
      near <- d <= g_i | !is.na(hit)
      st <- ifelse(!is.na(hit), "matched",
            ifelse(!inside, "outside_area",
            ifelse(!near, "not_found",
            ifelse(b_pass & !o_pass[j], "partner_below_cutoff",
            ifelse(!b_pass & o_pass[j], "this_below_cutoff",
            ifelse(!b_pass & !o_pass[j], "both_below_cutoff",
            # both pass: the other's particle is paired with another particle,
            # or it cannot be paired at all (below the size limit, or the two
            # sizes differ by more than match_max_size_ratio)
            ifelse(as.character(o$particle_id[j]) %in% ids_o, "paired_elsewhere",
                   "partner_not_matchable")))))))
      rho <- chance_partner_rate(o$x_aligned, o$y_aligned, 1)$density
      keep_p <- near & st != "outside_area"
      rows[[length(rows) + 1]] <- data.frame(
        of = roles[dir], found_by = roles[3 - dir], particle_id = as.character(b$particle_id),
        passes = b_pass, status = st,
        partner_id = ifelse(keep_p, as.character(o$particle_id[j]), NA_character_),
        partner_material = ifelse(keep_p, as.character(o$material[j]), NA_character_),
        partner_quality = ifelse(keep_p, suppressWarnings(as.numeric(o$quality[j])), NA_real_),
        partner_passes = ifelse(keep_p, o_pass[j], NA),
        partner_distance_um = ifelse(keep_p, round(d, 1), NA_real_),
        gate_um = round(g_i, 1), p_chance = round(1 - exp(-rho * pi * g_i^2), 3),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}


# Detected-rate columns of one class: raw rate with its Wilson interval, the
# expected chance detections and the chance-corrected rate (clamped to 0-100).
.detected_row <- function(det, p_chance, gate) {
  n <- length(det)
  if (n == 0 || all(is.na(det)))
    return(data.frame(gate_um = gate %||% NA_real_, n_detected = NA_integer_,
                      detected_pct = NA_real_, det_ci_low_pct = NA_real_,
                      det_ci_high_pct = NA_real_, det_chance_expected = NA_real_,
                      detected_corrected_pct = NA_real_))
  k <- sum(det); e <- sum(p_chance)
  ci <- wilson_interval(k, n)
  corr <- if (n - e > 0) min(100, max(0, 100 * (k - e) / (n - e))) else NA_real_
  data.frame(gate_um = round(gate, 1), n_detected = k, detected_pct = round(100 * k / n, 1),
             det_ci_low_pct = round(100 * ci$low, 1), det_ci_high_pct = round(100 * ci$high, 1),
             det_chance_expected = round(e, 1), detected_corrected_pct = round(corr, 1))
}

#' Display name of an instrument by its file name ("ftir_perkin" -> "FTIR (PerkinElmer)")
.file_role_name <- function(x) {
  r <- names(ROLE_FILE)[ROLE_FILE == x]
  if (length(r)) ROLE_NAME[[r[1]]] else x
}

#' Log one line per direction: the rate over all sizes, all particles.
log_detection_rates <- function(det) {
  if (is.null(det) || nrow(det) == 0) return(invisible(NULL))
  a <- det[det$group == "all" & det$size_class == "all sizes" & det$n_base > 0, ]
  for (i in seq_len(nrow(a)))
    log_message(sprintf("  Of the %d %s particles in its area, %s identified %s%s%s",
                        a$n_base[i], .file_role_name(a$of[i]), .file_role_name(a$found_by[i]),
                        format_rate_ci(a$n_found[i], a$n_base[i]),
                        if (is.na(a$rate_true_pct[i])) ""
                        else sprintf(" (chance-corrected %.0f%%)", a$rate_true_pct[i]),
                        if (is.null(a$n_detected) || is.na(a$n_detected[i])) ""
                        else sprintf(" and detected %s (chance-corrected %.0f%%; any particle within %g um)",
                                     format_rate_ci(a$n_detected[i], a$n_base[i]),
                                     a$detected_corrected_pct[i], a$gate_um[i])))
  invisible(NULL)
}
