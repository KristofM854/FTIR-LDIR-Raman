# =============================================================================
# 09_diagnostics.R — QC plots and diagnostic visualizations
# =============================================================================

#' Generate all diagnostic plots
#'
#' @param ftir_df FTIR data frame with x_aligned, y_aligned
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param match_result Result from match_particles()
#' @param icp_result Result from icp_refine()
#' @param agreement Result from analyze_agreement()
#' @param config Configuration list
#' @return List of ggplot objects
generate_diagnostics <- function(ftir_df, raman_df, match_result,
                                 icp_result, agreement, config) {
  log_message("Generating diagnostic plots")

  plots <- list()

  # 1. Overlay plot: FTIR vs Raman in aligned coordinate frame
  plots$overlay <- plot_overlay(ftir_df, raman_df, match_result)

  # 2. Match distance distribution
  plots$distance_hist <- plot_distance_distribution(match_result)

  # 3. ICP convergence
  plots$icp_convergence <- plot_icp_convergence(icp_result)

  # 4. Residual scatter
  plots$residual_scatter <- plot_residual_scatter(ftir_df, raman_df, match_result)

  # 5. Confusion matrix heatmap
  if (!is.null(agreement) && nrow(match_result$matched) > 0) {
    plots$confusion <- plot_confusion_matrix(agreement)
  }

  # 6. Size comparison for matched pairs
  if (nrow(match_result$matched) > 0) {
    plots$size_comparison <- plot_size_comparison(match_result)
  }

  log_message("  Generated ", length(plots), " diagnostic plots")
  plots
}


#' Overlay plot of aligned FTIR and Raman particles
plot_overlay <- function(ftir_df, raman_df, match_result) {
  matched <- match_result$matched

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = raman_df,
      ggplot2::aes(x = x_norm, y = y_norm),
      color = "steelblue", alpha = 0.4, size = 1.5, shape = 16
    ) +
    ggplot2::geom_point(
      data = ftir_df,
      ggplot2::aes(x = x_aligned, y = y_aligned),
      color = "tomato", alpha = 0.6, size = 2, shape = 17
    )

  # Draw lines connecting matched pairs
  if (nrow(matched) > 0) {
    segments_df <- data.frame(
      x    = matched$ftir_x_aligned,
      y    = matched$ftir_y_aligned,
      xend = matched$raman_x_norm,
      yend = matched$raman_y_norm
    )
    p <- p + ggplot2::geom_segment(
      data = segments_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey40", alpha = 0.3, linewidth = 0.3
    )
  }

  p <- p +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "FTIR–Raman Particle Overlay (aligned coordinates)",
      subtitle = paste0("Red triangles = FTIR, Blue circles = Raman, ",
                        "Lines = matched pairs"),
      x = "X (µm, normalized)", y = "Y (µm, normalized)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Distance distribution of matched particle pairs
plot_distance_distribution <- function(match_result) {
  matched <- match_result$matched

  if (nrow(matched) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No matches") +
             ggplot2::theme_void())
  }

  p <- ggplot2::ggplot(matched, ggplot2::aes(x = match_distance)) +
    ggplot2::geom_histogram(
      bins = 30, fill = "steelblue", color = "white", alpha = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = median(matched$match_distance),
      linetype = "dashed", color = "tomato", linewidth = 0.8
    ) +
    ggplot2::labs(
      title = "Match Distance Distribution",
      subtitle = paste0("Median = ",
                        round(median(matched$match_distance), 2), " µm, ",
                        "n = ", nrow(matched), " pairs"),
      x = "Match distance (µm)", y = "Count"
    ) +
    ggplot2::theme_minimal()

  p
}


#' ICP convergence plot (RMS error over iterations)
plot_icp_convergence <- function(icp_result) {
  if (length(icp_result$rms_history) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No ICP iterations") +
             ggplot2::theme_void())
  }

  conv_df <- data.frame(
    iteration = seq_along(icp_result$rms_history),
    rms       = icp_result$rms_history
  )

  p <- ggplot2::ggplot(conv_df, ggplot2::aes(x = iteration, y = rms)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      title = "ICP Convergence",
      subtitle = paste0("Converged: ", icp_result$converged,
                        ", Iterations: ", icp_result$n_iterations),
      x = "Iteration", y = "RMS error (µm)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Residual scatter: distance and direction of match offsets
plot_residual_scatter <- function(ftir_df, raman_df, match_result) {
  matched <- match_result$matched

  if (nrow(matched) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No matches") +
             ggplot2::theme_void())
  }

  residuals_df <- data.frame(
    dx = matched$raman_x_norm - matched$ftir_x_aligned,
    dy = matched$raman_y_norm - matched$ftir_y_aligned
  )

  p <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = dx, y = dy)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Alignment Residuals",
      subtitle = "Offset from FTIR (transformed) to Raman (nearest match)",
      x = "dX (µm)", y = "dY (µm)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Confusion matrix heatmap
plot_confusion_matrix <- function(agreement) {
  conf <- agreement$confusion_matrix

  if (length(conf) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No confusion data") +
             ggplot2::theme_void())
  }

  # Convert to data frame for ggplot
  conf_df <- as.data.frame(conf)
  names(conf_df) <- c("FTIR", "Raman", "Count")

  p <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Raman, y = FTIR, fill = Count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(
      title = "Material Confusion Matrix",
      subtitle = paste0("Overall agreement: ",
                        round(agreement$agreement_rate * 100, 1), "%"),
      x = "Raman material", y = "FTIR material"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' Scatter plot comparing particle sizes between FTIR and Raman
plot_size_comparison <- function(match_result) {
  matched <- match_result$matched

  # Try to find comparable size columns
  ftir_size  <- matched$ftir_feret_max_um
  raman_size <- matched$raman_feret_max_um

  if (all(is.na(ftir_size)) || all(is.na(raman_size))) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "Size data not available") +
             ggplot2::theme_void())
  }

  size_df <- data.frame(ftir_size = ftir_size, raman_size = raman_size)
  size_df <- size_df[!is.na(size_df$ftir_size) & !is.na(size_df$raman_size), ]

  p <- ggplot2::ggplot(size_df, ggplot2::aes(x = ftir_size, y = raman_size)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "tomato") +
    ggplot2::labs(
      title = "Particle Size Comparison (matched pairs)",
      x = "FTIR Feret max (µm)", y = "Raman Feret max (µm)"
    ) +
    ggplot2::theme_minimal()

  p
}
