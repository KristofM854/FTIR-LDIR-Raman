# =============================================================================
# 03_normalize.R — Initial coordinate normalization
# =============================================================================

#' Normalize coordinates of both point clouds
#'
#' Centers both point clouds at their respective centroids.
#' Optionally normalizes scale (divides by RMS distance from centroid).
#'
#' @param ftir_df FTIR data frame with x_um, y_um
#' @param raman_df Raman data frame with x_um, y_um
#' @param normalize_scale Logical, whether to normalize scale
#' @return List with:
#'   ftir  — FTIR df with added x_norm, y_norm
#'   raman — Raman df with added x_norm, y_norm
#'   ftir_centroid  — c(cx, cy) of FTIR
#'   raman_centroid — c(cx, cy) of Raman
#'   ftir_scale     — scale factor applied to FTIR (1 if not normalized)
#'   raman_scale    — scale factor applied to Raman (1 if not normalized)
normalize_coordinates <- function(ftir_df, raman_df, normalize_scale = FALSE) {
  log_message("Normalizing coordinates")

  # Centroids
  ftir_cx  <- mean(ftir_df$x_um, na.rm = TRUE)
  ftir_cy  <- mean(ftir_df$y_um, na.rm = TRUE)
  raman_cx <- mean(raman_df$x_um, na.rm = TRUE)
  raman_cy <- mean(raman_df$y_um, na.rm = TRUE)

  log_message("  FTIR centroid:  (", round(ftir_cx, 1), ", ", round(ftir_cy, 1), ")")
  log_message("  Raman centroid: (", round(raman_cx, 1), ", ", round(raman_cy, 1), ")")

  # Center
  ftir_df$x_norm  <- ftir_df$x_um - ftir_cx
  ftir_df$y_norm  <- ftir_df$y_um - ftir_cy
  raman_df$x_norm <- raman_df$x_um - raman_cx
  raman_df$y_norm <- raman_df$y_um - raman_cy

  ftir_scale  <- 1
  raman_scale <- 1

  if (normalize_scale) {
    # RMS distance from centroid
    ftir_rms  <- sqrt(mean(ftir_df$x_norm^2 + ftir_df$y_norm^2))
    raman_rms <- sqrt(mean(raman_df$x_norm^2 + raman_df$y_norm^2))

    if (ftir_rms > 0 && raman_rms > 0) {
      ftir_scale  <- ftir_rms
      raman_scale <- raman_rms

      ftir_df$x_norm  <- ftir_df$x_norm / ftir_scale
      ftir_df$y_norm  <- ftir_df$y_norm / ftir_scale
      raman_df$x_norm <- raman_df$x_norm / raman_scale
      raman_df$y_norm <- raman_df$y_norm / raman_scale

      log_message("  Scale normalization: FTIR RMS = ", round(ftir_rms, 1),
                  ", Raman RMS = ", round(raman_rms, 1))
    }
  }

  list(
    ftir           = ftir_df,
    raman          = raman_df,
    ftir_centroid  = c(ftir_cx, ftir_cy),
    raman_centroid = c(raman_cx, raman_cy),
    ftir_scale     = ftir_scale,
    raman_scale    = raman_scale
  )
}
