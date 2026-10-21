# =============================================================================
# detection_rates.R — how often one instrument finds what another found
# =============================================================================
# For every matched instrument pair, in both directions: of the particles the
# base instrument found (passing its quality cutoff, >= min_particle_size_um,
# and inside the area the other instrument covered), the share the other
# instrument also found at the same place, per size class of the base
# particle. "Found" means paired with a particle that passes the other
# instrument's quality cutoff (only those are matched), so the rate is the
# detection AND confident identification by the other instrument.
#
# Each rate has a 95% Wilson interval, and a chance-corrected rate: the found
# pairs weighted by p_true (pair_confidence()), the expected number of true
# partners among them.
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
#' @return data.frame: found_by, of, group ("all" / "plastics"), size_class
#'   (incl. "all sizes"), n_base, n_found, rate_pct, ci_low_pct, ci_high_pct,
#'   n_found_true (sum of p_true), rate_true_pct
detection_rates <- function(aligned, matches, config) {
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
      keep <- quality_passes(base, prefix, config) &
        !is.na(base$feret_max_um) & base$feret_max_um >= min_size &
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

      hit <- match(as.character(b$particle_id), ids)
      p_true <- if (!is.null(m$p_true)) m$p_true[hit] else rep(NA_real_, nrow(b))
      found <- !is.na(hit)
      fam <- b$material_family %||% classify_family_vec(b$material)
      size_class <- as.character(cut(b$feret_max_um, DETECTION_SIZE_BREAKS,
                                     DETECTION_SIZE_LABELS, right = FALSE))
      for (grp in c("all", "plastics")) {
        g <- if (grp == "all") rep(TRUE, nrow(b)) else fam %in% plastic_fams
        for (sc in c(DETECTION_SIZE_LABELS, "all sizes")) {
          s <- g & (sc == "all sizes" | size_class == sc)
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
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
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
    log_message(sprintf("  %s finds %s of the %d %s particles in its area%s",
                        .file_role_name(a$found_by[i]),
                        format_rate_ci(a$n_found[i], a$n_base[i]), a$n_base[i],
                        .file_role_name(a$of[i]),
                        if (is.na(a$rate_true_pct[i])) ""
                        else sprintf(" (chance-corrected %.0f%%)", a$rate_true_pct[i])))
  invisible(NULL)
}
