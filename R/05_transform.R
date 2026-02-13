# =============================================================================
# 05_transform.R — Apply spatial transform to FTIR coordinates
# =============================================================================

#' Apply estimated transform to FTIR particle coordinates
#'
#' Transforms the normalized FTIR coordinates into the Raman coordinate frame.
#'
#' @param ftir_df FTIR data frame with x_norm, y_norm
#' @param transform_matrix 3x3 homogeneous transform matrix
#' @return FTIR data frame with added columns x_aligned, y_aligned
apply_ftir_transform <- function(ftir_df, transform_matrix) {
  log_message("Applying transform to FTIR coordinates")

  result <- apply_transform_points(ftir_df$x_norm, ftir_df$y_norm, transform_matrix)

  ftir_df$x_aligned <- result$x_transformed
  ftir_df$y_aligned <- result$y_transformed

  log_message("  Transformed ", nrow(ftir_df), " FTIR particles")
  ftir_df
}
