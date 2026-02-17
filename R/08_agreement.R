# =============================================================================
# 08_agreement.R — Agreement analysis between instrument identifications
# =============================================================================
#
# Supports tiered agreement scoring:
#   Exact    — both instruments agree on the same canonical name
#   Family   — same polymer family (e.g., Cellulose + CAB both → Cellulose fam)
#   Disagree — different families
#
# Also provides category-level summaries (Synthetic / Semi-synthetic / Natural).
# =============================================================================

#' Normalize a material name to a canonical short form
#'
#' Maps verbose instrument-specific names to a common short abbreviation.
#' E.g., "Polyethylene terephtalate (PET)" -> "PET",
#'       "Polypropylene (PP) Lab bench" -> "PP"
#'
#' @param x Character vector of material names
#' @return Character vector of normalized names
normalize_material <- function(x) {
  x <- trimws(x)
  result <- character(length(x))

  for (i in seq_along(x)) {
    name <- x[i]

    # Try to extract abbreviation from parentheses, e.g., "Polycarbonate (PC)" -> "PC"
    m <- regmatches(name, regexpr("\\(([A-Za-z0-9 ]+)\\)", name))
    if (length(m) == 1 && nchar(m) > 0) {
      result[i] <- gsub("[()]", "", m)
    } else {
      result[i] <- name
    }
  }

  # Normalize common FTIR short names
  result <- gsub("^Polypro$", "PP", result)

  # Strip suffixes like "Lab bench", "Nylon 12" for broader matching
  result <- gsub("\\s+(Lab bench|Nylon \\d+)$", "", result, ignore.case = TRUE)

  toupper(trimws(result))
}


#' Analyze agreement between matched instrument pairs
#'
#' For matched particle pairs, compares material/polymer labels using both
#' exact-match and tiered (Exact/Family/Disagree) scoring.
#'
#' @param match_result List returned by match_particles()
#' @param config Configuration list (optional; uses raman_hqi_threshold)
#' @param instrument_a Label for the first instrument (default "FTIR")
#' @param instrument_b Label for the second instrument (default "Raman")
#' @return List with agreement metrics, tiered scoring, and summaries
analyze_agreement <- function(match_result, config = NULL,
                              instrument_a = "FTIR",
                              instrument_b = "Raman") {
  log_message("Analyzing material agreement (", instrument_a, " vs ", instrument_b, ")")

  matched <- match_result$matched
  prefix_a <- tolower(instrument_a)
  prefix_b <- tolower(instrument_b)

  if (nrow(matched) == 0) {
    log_message("  No matched pairs to analyze.")
    return(.empty_agreement())
  }

  # Column names depend on instrument prefixes
  mat_col_a <- paste0(prefix_a, "_material")
  mat_col_b <- paste0(prefix_b, "_material")
  qual_col_b <- paste0(prefix_b, "_quality")

  # Only evaluate agreement for pairs where instrument B quality is trusted
  hqi_threshold <- if (!is.null(config)) config$raman_hqi_threshold else 0
  # Only apply HQI filter if instrument B is Raman
  if (instrument_b == "Raman" && !is.null(matched[[qual_col_b]]) && hqi_threshold > 0) {
    raman_hqi <- matched[[qual_col_b]]
    hqi_ok <- !is.na(raman_hqi) & raman_hqi >= hqi_threshold
    n_excluded <- sum(!hqi_ok)
    if (n_excluded > 0) {
      log_message("  Excluding ", n_excluded, " pairs with Raman HQI < ",
                  hqi_threshold, " from agreement scoring")
    }
    matched_for_agree <- matched[hqi_ok, ]
  } else {
    matched_for_agree <- matched
  }

  if (nrow(matched_for_agree) == 0) {
    log_message("  No pairs with sufficient quality for agreement analysis.")
    return(.empty_agreement())
  }

  # Extract material labels
  mat_a_raw <- matched_for_agree[[mat_col_a]]
  mat_b_raw <- matched_for_agree[[mat_col_b]]

  # Handle NAs
  mat_a_raw[is.na(mat_a_raw)] <- "Unknown"
  mat_b_raw[is.na(mat_b_raw)] <- "Unknown"

  # Use material equivalence mapping if configured, else abbreviation normalization
  mapped <- map_paired_materials(mat_a_raw, mat_b_raw, config)
  mat_a <- mapped$ftir_canonical
  mat_b <- mapped$raman_canonical

  # ---- Tiered agreement scoring ----
  tier <- score_tiered_agreement_vec(mat_a_raw, mat_b_raw)

  # Family + category enrichment
  fam_a <- classify_family_vec(mat_a_raw)
  fam_b <- classify_family_vec(mat_b_raw)
  cat_a <- classify_category_vec(fam_a)
  cat_b <- classify_category_vec(fam_b)

  # ---- Classic exact agreement (on normalized names) ----
  agree_mask     <- mat_a == mat_b
  agreement_rate <- mean(agree_mask)

  # Tiered rates
  n_total   <- length(tier)
  n_exact   <- sum(tier == "Exact")
  n_family  <- sum(tier == "Family")
  n_disagree <- sum(tier == "Disagree")

  log_message("  Tiered agreement (n=", n_total, "):")
  log_message("    Exact:    ", n_exact, " (", round(n_exact / n_total * 100, 1), "%)")
  log_message("    Family:   ", n_family, " (", round(n_family / n_total * 100, 1), "%)")
  log_message("    Disagree: ", n_disagree, " (", round(n_disagree / n_total * 100, 1), "%)")

  # ---- Confusion matrix (use normalized names) ----
  conf_mat <- table(setNames(mat_a, NULL), setNames(mat_b, NULL))
  names(dimnames(conf_mat)) <- c(instrument_a, instrument_b)

  # ---- Per-material agreement detail ----
  materials_a <- unique(mat_a)
  agreement_detail <- data.frame(
    material_a     = character(),
    family_a       = character(),
    category_a     = character(),
    n_pairs        = integer(),
    n_exact        = integer(),
    n_family       = integer(),
    n_disagree     = integer(),
    exact_pct      = numeric(),
    family_or_better_pct = numeric(),
    stringsAsFactors = FALSE
  )

  for (mat in sort(materials_a)) {
    mask <- mat_a == mat
    n_tot <- sum(mask)
    tier_sub <- tier[mask]
    n_ex <- sum(tier_sub == "Exact")
    n_fm <- sum(tier_sub == "Family")
    n_ds <- sum(tier_sub == "Disagree")
    fam_val <- classify_family(mat)
    cat_val <- classify_category(fam_val)

    agreement_detail <- rbind(agreement_detail, data.frame(
      material_a     = mat,
      family_a       = fam_val,
      category_a     = cat_val,
      n_pairs        = n_tot,
      n_exact        = n_ex,
      n_family       = n_fm,
      n_disagree     = n_ds,
      exact_pct      = round(n_ex / n_tot * 100, 1),
      family_or_better_pct = round((n_ex + n_fm) / n_tot * 100, 1),
      stringsAsFactors = FALSE
    ))
  }

  # ---- Category-level summary ----
  category_summary <- data.frame(
    category       = character(),
    n_pairs        = integer(),
    n_exact        = integer(),
    n_family       = integer(),
    n_disagree     = integer(),
    exact_pct      = numeric(),
    family_or_better_pct = numeric(),
    stringsAsFactors = FALSE
  )
  for (cat_val in c("Synthetic", "Semi-synthetic", "Natural/Organic", "Unknown")) {
    mask <- cat_a == cat_val
    if (!any(mask)) next
    n_tot <- sum(mask)
    tier_sub <- tier[mask]
    n_ex <- sum(tier_sub == "Exact")
    n_fm <- sum(tier_sub == "Family")
    n_ds <- sum(tier_sub == "Disagree")
    category_summary <- rbind(category_summary, data.frame(
      category       = cat_val,
      n_pairs        = n_tot,
      n_exact        = n_ex,
      n_family       = n_fm,
      n_disagree     = n_ds,
      exact_pct      = round(n_ex / n_tot * 100, 1),
      family_or_better_pct = round((n_ex + n_fm) / n_tot * 100, 1),
      stringsAsFactors = FALSE
    ))
  }

  # ---- Size-dependent analysis ----
  size_col <- paste0(prefix_a, "_feret_max_um")
  if (size_col %in% names(matched_for_agree) && any(!is.na(matched_for_agree[[size_col]]))) {
    size_vals <- matched_for_agree[[size_col]]
    size_class <- cut(
      size_vals,
      breaks = c(0, 50, 100, 200, 300, Inf),
      labels = c("<50", "50-100", "100-200", "200-300", ">300"),
      right = FALSE
    )

    size_analysis <- data.frame(
      size_class     = levels(size_class),
      n_pairs        = as.integer(table(size_class)),
      n_agree        = as.integer(tapply(agree_mask, size_class, sum)),
      n_exact        = as.integer(tapply(tier == "Exact", size_class, sum)),
      n_family       = as.integer(tapply(tier == "Family", size_class, sum)),
      stringsAsFactors = FALSE
    )
    size_analysis$n_disagree    <- size_analysis$n_pairs - size_analysis$n_exact - size_analysis$n_family
    size_analysis$agreement_pct <- ifelse(
      size_analysis$n_pairs > 0,
      round(size_analysis$n_agree / size_analysis$n_pairs * 100, 1),
      NA_real_
    )
    size_analysis$family_or_better_pct <- ifelse(
      size_analysis$n_pairs > 0,
      round((size_analysis$n_exact + size_analysis$n_family) / size_analysis$n_pairs * 100, 1),
      NA_real_
    )
  } else {
    size_analysis <- data.frame(
      size_class = character(), n_pairs = integer(),
      n_agree = integer(), n_exact = integer(), n_family = integer(),
      n_disagree = integer(), agreement_pct = numeric(),
      family_or_better_pct = numeric()
    )
  }

  # ---- Pair-level summary ----
  summary_df <- data.frame(
    match_id           = matched_for_agree$match_id,
    material_a_raw     = mat_a_raw,
    material_b_raw     = mat_b_raw,
    material_a_norm    = mat_a,
    material_b_norm    = mat_b,
    family_a           = fam_a,
    family_b           = fam_b,
    category_a         = cat_a,
    category_b         = cat_b,
    tier               = tier,
    agree_exact        = agree_mask,
    match_distance     = matched_for_agree$match_distance,
    stringsAsFactors = FALSE
  )

  # Back-compat column names
  if ("ftir_material" %in% names(matched_for_agree)) {
    summary_df$ftir_material_raw  <- mat_a_raw
    summary_df$raman_material_raw <- mat_b_raw
    summary_df$ftir_material_norm <- mat_a
    summary_df$raman_material_norm <- mat_b
  }

  log_message("  Agreement detail by ", instrument_a, " material:")
  for (i in seq_len(nrow(agreement_detail))) {
    log_message("    ", agreement_detail$material_a[i], " [", agreement_detail$family_a[i], "]: ",
                "Exact ", agreement_detail$exact_pct[i], "%, ",
                "Family+ ", agreement_detail$family_or_better_pct[i], "% (",
                agreement_detail$n_pairs[i], " pairs)")
  }

  list(
    confusion_matrix  = conf_mat,
    agreement_rate    = agreement_rate,
    tiered_rates      = list(
      n_total  = n_total,
      n_exact  = n_exact,
      n_family = n_family,
      n_disagree = n_disagree,
      exact_pct   = round(n_exact / n_total * 100, 1),
      family_pct  = round(n_family / n_total * 100, 1),
      family_or_better_pct = round((n_exact + n_family) / n_total * 100, 1)
    ),
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
    tiered_rates     = list(n_total = 0, n_exact = 0, n_family = 0, n_disagree = 0,
                            exact_pct = NA_real_, family_pct = NA_real_,
                            family_or_better_pct = NA_real_),
    agreement_detail = data.frame(),
    category_summary = data.frame(),
    size_analysis    = data.frame(),
    summary_df       = data.frame()
  )
}
