# =============================================================================
# 08_agreement.R — Agreement analysis between FTIR and Raman identifications
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


#' Analyze agreement between FTIR and Raman material identifications
#'
#' For matched particle pairs, compares the material/polymer labels from
#' both instruments. Material names are normalized to canonical abbreviations
#' (e.g., "Polyethylene terephtalate (PET)" and "PET" both become "PET").
#'
#' @param match_result List returned by match_particles()
#' @return List with:
#'   confusion_matrix — table of FTIR material vs Raman material
#'   agreement_rate   — fraction of matched pairs with matching material
#'   agreement_detail — per-material agreement rates
#'   size_analysis    — disagreement broken down by particle size class
#'   summary_df       — tidy data frame of pair-level agreement
analyze_agreement <- function(match_result) {
  log_message("Analyzing material agreement")

  matched <- match_result$matched

  if (nrow(matched) == 0) {
    log_message("  No matched pairs to analyze.")
    return(list(
      confusion_matrix = table(NULL),
      agreement_rate   = NA_real_,
      agreement_detail = data.frame(),
      size_analysis    = data.frame(),
      summary_df       = data.frame()
    ))
  }

  # Extract material labels
  ftir_mat_raw  <- matched$ftir_material
  raman_mat_raw <- matched$raman_material

  # Handle NAs
  ftir_mat_raw[is.na(ftir_mat_raw)]   <- "Unknown"
  raman_mat_raw[is.na(raman_mat_raw)] <- "Unknown"

  # Normalize to canonical abbreviations for comparison
  ftir_mat  <- normalize_material(ftir_mat_raw)
  raman_mat <- normalize_material(raman_mat_raw)

  # Confusion matrix (use normalized names)
  conf_mat <- table(FTIR = ftir_mat, Raman = raman_mat)

  # Overall agreement (normalized comparison)
  agree_mask     <- ftir_mat == raman_mat
  agreement_rate <- mean(agree_mask)
  log_message("  Overall agreement rate: ", round(agreement_rate * 100, 1), "%")

  # Per-FTIR-material agreement
  ftir_materials <- unique(ftir_mat)
  agreement_detail <- data.frame(
    ftir_material   = character(),
    n_pairs         = integer(),
    n_agree         = integer(),
    agreement_pct   = numeric(),
    stringsAsFactors = FALSE
  )

  for (mat in sort(ftir_materials)) {
    mask   <- ftir_mat == mat
    n_tot  <- sum(mask)
    n_ok   <- sum(agree_mask[mask])
    agreement_detail <- rbind(agreement_detail, data.frame(
      ftir_material = mat,
      n_pairs       = n_tot,
      n_agree       = n_ok,
      agreement_pct = round(n_ok / n_tot * 100, 1),
      stringsAsFactors = FALSE
    ))
  }

  # Size-dependent disagreement analysis
  # Use FTIR feret_max_um for size classes
  size_col <- "ftir_feret_max_um"
  if (size_col %in% names(matched) && any(!is.na(matched[[size_col]]))) {
    size_vals <- matched[[size_col]]
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
      stringsAsFactors = FALSE
    )
    size_analysis$n_disagree    <- size_analysis$n_pairs - size_analysis$n_agree
    size_analysis$agreement_pct <- ifelse(
      size_analysis$n_pairs > 0,
      round(size_analysis$n_agree / size_analysis$n_pairs * 100, 1),
      NA_real_
    )
  } else {
    size_analysis <- data.frame(
      size_class = character(), n_pairs = integer(),
      n_agree = integer(), n_disagree = integer(),
      agreement_pct = numeric()
    )
  }

  # Pair-level summary
  summary_df <- data.frame(
    match_id             = matched$match_id,
    ftir_material_raw    = ftir_mat_raw,
    raman_material_raw   = raman_mat_raw,
    ftir_material_norm   = ftir_mat,
    raman_material_norm  = raman_mat,
    agree                = agree_mask,
    match_distance       = matched$match_distance,
    stringsAsFactors = FALSE
  )

  log_message("  Agreement detail by FTIR material:")
  for (i in seq_len(nrow(agreement_detail))) {
    log_message("    ", agreement_detail$ftir_material[i], ": ",
                agreement_detail$agreement_pct[i], "% (",
                agreement_detail$n_agree[i], "/",
                agreement_detail$n_pairs[i], ")")
  }

  list(
    confusion_matrix = conf_mat,
    agreement_rate   = agreement_rate,
    agreement_detail = agreement_detail,
    size_analysis    = size_analysis,
    summary_df       = summary_df
  )
}
