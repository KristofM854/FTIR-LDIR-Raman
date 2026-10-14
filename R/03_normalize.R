# =============================================================================
# 03_normalize.R — Initial coordinate normalization
# =============================================================================

#' Apply an existing centroid (and scale) to another data frame
#'
#' Every particle set that takes part in a spatial step must be centred on the
#' SAME centroid as the set the transform was estimated from, otherwise the
#' transform is applied in the wrong frame. This re-applies the normalization
#' of a cloud (cloud_normalization(), R/align_and_match.R) to any other data
#' frame.
#'
#' @param df Data frame with x_um, y_um
#' @param centroid Numeric c(cx, cy) from cloud_normalization()
#' @param scale Scale divisor (1 when normalize_scale was FALSE)
#' @return df with x_norm, y_norm added / overwritten
apply_normalization <- function(df, centroid, scale = 1) {
  if (is.null(scale) || is.na(scale) || scale == 0) scale <- 1
  df$x_norm <- (df$x_um - centroid[1]) / scale
  df$y_norm <- (df$y_um - centroid[2]) / scale
  df
}
