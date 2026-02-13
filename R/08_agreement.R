# =============================================================================
# 08_agreement.R — Agreement analysis between FTIR and Raman identifications
# =============================================================================

#' Analyze agreement between FTIR and Raman material identifications
#'
#' For matched particle pairs, compares the material/polymer labels from
#' both instruments. Generates confusion matrix, match statistics, and
#' size-dependent disagreement analysis.
#'
#' @param match_result List returned by match_particles()
#' @return List with:
#'   confusion_matrix — table of FTIR material vs Raman material
#'   agreement_rate   — fraction of matched pairs with identical material
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
  ftir_mat  <- matched$ftir_material
  raman_mat <- matched$raman_material

  # Handle NAs
  ftir_mat[is.na(ftir_mat)]   <- "Unknown"
  raman_mat[is.na(raman_mat)] <- "Unknown"

  # Confusion matrix
  conf_mat <- table(FTIR = ftir_mat, Raman = raman_mat)

  # Overall agreement
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
    match_id       = matched$match_id,
    ftir_material  = ftir_mat,
    raman_material = raman_mat,
    agree          = agree_mask,
    match_distance = matched$match_distance,
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
