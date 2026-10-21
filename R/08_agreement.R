# =============================================================================
# 08_agreement.R — Agreement analysis between instrument identifications
# =============================================================================
#
# Supports tiered agreement scoring (R/08b_material_map.R):
#   Exact          — same polymer (secondary agreement)
#   Family         — different polymers of the same family; Exact + Family is
#                    the PRIMARY agreement ("family hit")
#   Filler/pigment — one instrument reports a plastic, the other a filler,
#                    pigment or additive compounded into plastics (TiO2,
#                    CaCO3, talc, ...): compatible, but the polymer is only
#                    confirmed by one instrument
#   Unresolved     — names that do not resolve to a comparable family; needs
#                    inspection, excluded from the rates
#   Disagree       — anything else
#
# Also provides category-level summaries (Synthetic / Semi-synthetic / Natural /
# Additive/Pigment / Inorganic / Other — see material_category_levels).
# =============================================================================

#' Quality cutoff for agreement scoring, per instrument prefix.
#'
#' Alignment and matching use every particle; only the agreement statistics
#' are restricted to identifications the instrument itself trusts:
#'   ftir  (PerkinElmer and Bruker): AAU score >= ftir_agreement_min_quality
#'   raman: HQI >= raman_hqi_threshold
#'   ldir : quality >= ldir_hqi_unknown_threshold (the Agilent "unknown" rule)
.agreement_quality_cutoff <- function(prefix, config) {
  if (startsWith(prefix, "ftir")) prefix <- "ftir"   # ftir_perkin, ftir_bruker
  key <- switch(prefix,
                ftir  = "ftir_agreement_min_quality",
                raman = "raman_hqi_threshold",
                ldir  = "ldir_hqi_unknown_threshold",
                NULL)
  if (is.null(key)) return(0)
  v <- if (is.null(config)) make_config()[[key]] else cfg_value(config, key)
  if (is.null(v) || !is.numeric(v)) 0 else v
}

# Pair-level tier counts and rates. Rates are over RESOLVED pairs (Unresolved
# pairs are neither agreements nor disagreements); unresolved_pct is over all
# scored pairs.
.tier_rates <- function(tier) {
  n_scored   <- length(tier)
  n_unres    <- sum(tier == "Unresolved")
  n_resolved <- n_scored - n_unres
  pct <- function(n, d) if (d > 0) round(n / d * 100, 1) else NA_real_
  n_exact <- sum(tier == "Exact"); n_family <- sum(tier == "Family")
  n_filler <- sum(tier == "Filler/pigment"); n_dis <- sum(tier == "Disagree")
  ci_fam <- wilson_interval(n_exact + n_family, n_resolved)
  ci_exa <- wilson_interval(n_exact, n_resolved)
  list(n_total = n_scored, n_resolved = n_resolved, n_unresolved = n_unres,
       # 95% Wilson intervals of the primary and secondary rates (%)
       family_ci_low_pct  = round(100 * ci_fam$low, 1),
       family_ci_high_pct = round(100 * ci_fam$high, 1),
       exact_ci_low_pct   = round(100 * ci_exa$low, 1),
       exact_ci_high_pct  = round(100 * ci_exa$high, 1),
       n_exact = n_exact, n_family = n_family, n_filler = n_filler,
       n_disagree = n_dis,
       family_or_better_pct = pct(n_exact + n_family, n_resolved),  # PRIMARY
       exact_pct    = pct(n_exact, n_resolved),                     # secondary
       family_pct   = pct(n_family, n_resolved),
       filler_pct   = pct(n_filler, n_resolved),
       disagree_pct = pct(n_dis, n_resolved),
       unresolved_pct = pct(n_unres, n_scored))
}


#' Analyze agreement between matched instrument pairs
#'
#' Scores each matched pair's materials in tiers (see
#' score_tiered_agreement_vec()): the primary agreement is the family hit
#' (Exact + Family), the secondary the exact polymer hit. Before scoring, a
#' pair is removed when either side's identification is below its
#' instrument's quality cutoff (.agreement_quality_cutoff()) or carries no
#' identification ("unknown", NA): low-quality identifications never enter an
#' agreement denominator. Unresolved pairs are counted but excluded from the
#' rates. LDIR pairs whose coordinates came from a low-confidence image join
#' are kept and counted (n_join_warning).
#'
#' @param match_result List returned by match_particles()
#' @param config Configuration list (quality cutoffs; NULL = defaults)
#' @param instrument_a Label for the first instrument (default "FTIR")
#' @param instrument_b Label for the second instrument (default "Raman")
#' @param prefix_a,prefix_b Match-column prefixes of the two instruments
#'   (default: the lowercased labels)
#' @return List with agreement metrics, tiered scoring, and summaries
analyze_agreement <- function(match_result, config = NULL,
                              instrument_a = "FTIR",
                              instrument_b = "Raman",
                              prefix_a = tolower(instrument_a),
                              prefix_b = tolower(instrument_b)) {
  log_message("Analyzing material agreement (", instrument_a, " vs ", instrument_b, ")")

  matched <- match_result$matched

  if (is.null(matched) || nrow(matched) == 0) {
    log_message("  No matched pairs to analyze.")
    return(.empty_agreement())
  }
  n_pairs <- nrow(matched)

  mat_col_a <- paste0(prefix_a, "_material")
  mat_col_b <- paste0(prefix_b, "_material")

  # ---- Quality: identifications each instrument trusts ----
  keep <- rep(TRUE, n_pairs)
  n_low_q <- c(a = 0L, b = 0L)
  for (side in c("a", "b")) {
    prefix <- if (side == "a") prefix_a else prefix_b
    qcol   <- paste0(prefix, "_quality")
    cutoff <- .agreement_quality_cutoff(prefix, config)
    q <- matched[[qcol]]
    if (cutoff <= 0 || is.null(q)) next
    q <- suppressWarnings(as.numeric(q))
    if (all(is.na(q))) {
      log_message("  ", toupper(prefix), " has no quality scores -- its quality ",
                  "cutoff (", cutoff, ") cannot be applied", level = "WARN")
      next
    }
    ok <- !is.na(q) & q >= cutoff
    n_low_q[side] <- sum(!ok & keep)
    keep <- keep & ok
  }
  # ---- Identifications: "unknown"/NA labels are not scored ----
  mat_a_all <- as.character(matched[[mat_col_a]])
  mat_b_all <- as.character(matched[[mat_col_b]])
  unident <- keep & (is_unidentified(mat_a_all) | is_unidentified(mat_b_all))
  keep <- keep & !unident

  log_message("  Pairs: ", n_pairs, " matched; excluded ", n_low_q[["a"]], " below ",
              instrument_a, " quality cutoff, ", n_low_q[["b"]], " below ",
              instrument_b, " quality cutoff, ", sum(unident),
              " without an identification; ", sum(keep), " scored")

  matched_for_agree <- matched[keep, , drop = FALSE]
  if (nrow(matched_for_agree) == 0) {
    log_message("  No pairs with sufficient quality for agreement analysis.")
    out <- .empty_agreement()
    out$tiered_rates$n_pairs_matched <- n_pairs
    return(out)
  }

  mat_a_raw <- mat_a_all[keep]
  mat_b_raw <- mat_b_all[keep]

  # Exact polymer names
  mat_a <- canonical_polymer_vec(mat_a_raw)
  mat_b <- canonical_polymer_vec(mat_b_raw)

  # ---- Tiered agreement scoring ----
  tier <- score_tiered_agreement_vec(mat_a_raw, mat_b_raw)
  fam_a <- classify_family_vec(mat_a_raw)
  fam_b <- classify_family_vec(mat_b_raw)
  cat_a <- classify_category_vec(fam_a)
  cat_b <- classify_category_vec(fam_b)

  rates <- .tier_rates(tier)
  rates$n_pairs_matched          <- n_pairs
  rates$n_excluded_quality_a     <- n_low_q[["a"]]
  rates$n_excluded_quality_b     <- n_low_q[["b"]]
  rates$n_excluded_unidentified  <- sum(unident)

  # LDIR coordinates from a low-confidence image join: kept, but counted
  jw_cols <- intersect(c(paste0(prefix_a, "_coord_join_warning"),
                         paste0(prefix_b, "_coord_join_warning")),
                       names(matched_for_agree))
  join_warning <- if (length(jw_cols) > 0)
    Reduce(`|`, lapply(jw_cols, function(cn) {
      v <- as.logical(matched_for_agree[[cn]]); !is.na(v) & v }))
  else rep(FALSE, length(tier))
  rates$n_join_warning <- sum(join_warning)

  log_message("  Agreement (n = ", rates$n_total, " scored, ", rates$n_unresolved,
              " unresolved, ", rates$n_resolved, " resolved):")
  log_message("    Family agreement (primary): ",
              format_rate_ci(rates$n_exact + rates$n_family, rates$n_resolved), "  ",
              "[Exact ", rates$n_exact, " + Family ", rates$n_family, "; 95% interval]")
  log_message("    Exact polymer (secondary):  ", format_rate_ci(rates$n_exact, rates$n_resolved))
  log_message("    Filler/pigment: ", rates$n_filler, " (", rates$filler_pct, "%), ",
              "Disagree: ", rates$n_disagree, " (", rates$disagree_pct, "%)")
  if (rates$n_unresolved > 0)
    log_message("    Unresolved (needs inspection): ", rates$n_unresolved, " (",
                rates$unresolved_pct, "% of scored)", level = "WARN")
  if (rates$n_join_warning > 0)
    log_message("    ", rates$n_join_warning, " scored pairs use LDIR coordinates ",
                "from a low-confidence image join (coord_join_warning)", level = "WARN")

  # The same over the pairs that are very likely true partners (p_true from
  # the positional-error fit, pair_confidence()): on a dense sample, chance
  # pairs dilute the rate; this one is nearly free of them.
  p_true <- matched_for_agree$p_true %||% rep(NA_real_, length(tier))
  hc_min <- cfg_value(config %||% list(), "pair_confidence_high")
  hc <- !is.na(p_true) & p_true >= hc_min
  r_hc <- .tier_rates(tier[hc])
  rates$high_conf_min <- hc_min
  rates$high_conf_n_resolved <- if (any(!is.na(p_true))) r_hc$n_resolved else NA_integer_
  rates$high_conf_family_pct <- r_hc$family_or_better_pct
  rates$high_conf_family_ci_low_pct  <- r_hc$family_ci_low_pct
  rates$high_conf_family_ci_high_pct <- r_hc$family_ci_high_pct
  if (any(!is.na(p_true)))
    log_message("    Pairs with p_true >= ", hc_min, ": ", r_hc$n_resolved, " resolved; ",
                "family agreement ", format_rate_ci(r_hc$n_exact + r_hc$n_family, r_hc$n_resolved))

  # ---- Confusion matrices (exact polymer and family) ----
  conf_mat <- table(setNames(mat_a, NULL), setNames(mat_b, NULL))
  names(dimnames(conf_mat)) <- c(instrument_a, instrument_b)
  conf_fam <- table(setNames(fam_a, NULL), setNames(fam_b, NULL))
  names(dimnames(conf_fam)) <- c(instrument_a, instrument_b)

  # ---- Breakdowns: by polymer of A, by category of A, by size of A ----
  breakdown <- function(group, levels = sort(unique(group))) {
    rows <- lapply(levels, function(g) {
      m <- group == g
      if (!any(m, na.rm = TRUE)) return(NULL)
      r <- .tier_rates(tier[m & !is.na(m)])
      data.frame(group = g, n_pairs = r$n_total, n_resolved = r$n_resolved,
                 n_exact = r$n_exact, n_family = r$n_family, n_filler = r$n_filler,
                 n_disagree = r$n_disagree, n_unresolved = r$n_unresolved,
                 family_or_better_pct = r$family_or_better_pct,
                 exact_pct = r$exact_pct, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }

  agreement_detail <- breakdown(mat_a)
  names(agreement_detail)[1] <- "material_a"
  agreement_detail$family_a   <- vapply(agreement_detail$material_a, function(m)
    fam_a[match(m, mat_a)], character(1))
  agreement_detail$category_a <- classify_category_vec(agreement_detail$family_a)

  category_summary <- breakdown(cat_a, intersect(material_category_levels, unique(cat_a)))
  names(category_summary)[1] <- "category"

  size_col <- paste0(prefix_a, "_feret_max_um")
  if (size_col %in% names(matched_for_agree) && any(!is.na(matched_for_agree[[size_col]]))) {
    size_class <- cut(matched_for_agree[[size_col]],
                      breaks = c(0, 50, 100, 200, 300, Inf),
                      labels = c("<50", "50-100", "100-200", "200-300", ">300"),
                      right = FALSE)
    size_analysis <- breakdown(as.character(size_class), levels(size_class))
    names(size_analysis)[1] <- "size_class"
  } else {
    size_analysis <- data.frame()
  }

  # ---- Pair-level summary ----
  summary_df <- data.frame(
    match_id           = matched_for_agree$match_id,
    material_a_raw     = mat_a_raw,
    material_b_raw     = mat_b_raw,
    polymer_a          = mat_a,
    polymer_b          = mat_b,
    family_a           = fam_a,
    family_b           = fam_b,
    category_a         = cat_a,
    category_b         = cat_b,
    tier               = tier,
    agree_family       = tier %in% c("Exact", "Family"),
    agree_exact        = tier == "Exact",
    coord_join_warning = join_warning,
    match_distance     = matched_for_agree$match_distance,
    p_true             = p_true,
    stringsAsFactors = FALSE
  )

  list(
    confusion_matrix  = conf_mat,
    confusion_family  = conf_fam,
    # Primary agreement (family hit) as a fraction of resolved pairs
    agreement_rate    = if (rates$n_resolved > 0)
                          (rates$n_exact + rates$n_family) / rates$n_resolved else NA_real_,
    tiered_rates      = rates,
    agreement_detail  = agreement_detail,
    category_summary  = category_summary,
    size_analysis     = size_analysis,
    summary_df        = summary_df,
    instrument_a      = instrument_a,
    instrument_b      = instrument_b
  )
}


#' Return empty agreement structure (for no-match cases)
.empty_agreement <- function() {
  list(
    confusion_matrix = table(NULL),
    agreement_rate   = NA_real_,
    tiered_rates     = .tier_rates(character(0)),
    agreement_detail = data.frame(),
    category_summary = data.frame(),
    size_analysis    = data.frame(),
    summary_df       = data.frame()
  )
}
