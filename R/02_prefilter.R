# =============================================================================
# 02_prefilter.R — Pre-filtering of particle datasets
# =============================================================================

#' Pre-filter FTIR particle data
#'
#' Removes rows with missing coordinates and applies optional quality/size filters.
#'
#' @param df Data frame from ingest_ftir()
#' @param min_quality Minimum AAU quality score (0 = keep all)
#' @param min_size_um Minimum particle size in µm (applied to feret_max_um)
#' @return Filtered data frame
prefilter_ftir <- function(df, min_quality = 0, min_size_um = 0) {
  n_start <- nrow(df)
  log_message("Pre-filtering FTIR: ", n_start, " particles")

  # Remove rows with missing coordinates
  df <- df[!is.na(df$x_um) & !is.na(df$y_um), ]
  n_coord <- nrow(df)
  if (n_coord < n_start) {
    log_message("  Removed ", n_start - n_coord, " particles with missing coordinates")
  }

  # Apply quality threshold
  if (min_quality > 0 && any(!is.na(df$quality))) {
    keep <- is.na(df$quality) | df$quality >= min_quality
    df <- df[keep, ]
    log_message("  Quality filter (>= ", min_quality, "): kept ", nrow(df),
                " of ", n_coord, " particles")
  }

  # Apply minimum size filter
  if (min_size_um > 0 && any(!is.na(df$feret_max_um))) {
    keep <- is.na(df$feret_max_um) | df$feret_max_um >= min_size_um
    df <- df[keep, ]
    log_message("  Size filter (>= ", min_size_um, " µm): kept ", nrow(df), " particles")
  }

  log_message("  FTIR after filtering: ", nrow(df), " particles")
  df
}


#' Pre-filter Raman particle data
#'
#' Removes rows with missing coordinates and applies HQI / size filters.
#'
#' @param df Data frame from ingest_raman()
#' @param min_hqi Minimum HQI threshold
#' @param min_size_um Minimum particle size in µm (applied to feret_max_um)
#' @return Filtered data frame
prefilter_raman <- function(df, min_hqi = 0, min_size_um = 0) {
  n_start <- nrow(df)
  log_message("Pre-filtering Raman: ", n_start, " particles")

  # Remove rows with missing coordinates
  df <- df[!is.na(df$x_um) & !is.na(df$y_um), ]
  n_coord <- nrow(df)
  if (n_coord < n_start) {
    log_message("  Removed ", n_start - n_coord, " particles with missing coordinates")
  }

  # Apply HQI threshold
  if (min_hqi > 0 && any(!is.na(df$quality))) {
    keep <- !is.na(df$quality) & df$quality >= min_hqi
    df <- df[keep, ]
    log_message("  HQI filter (>= ", min_hqi, "): kept ", nrow(df),
                " of ", n_coord, " particles")
  }

  # Apply minimum size filter
  if (min_size_um > 0 && any(!is.na(df$feret_max_um))) {
    keep <- is.na(df$feret_max_um) | df$feret_max_um >= min_size_um
    df <- df[keep, ]
    log_message("  Size filter (>= ", min_size_um, " µm): kept ", nrow(df), " particles")
  }

  log_message("  Raman after filtering: ", nrow(df), " particles")
  df
}
