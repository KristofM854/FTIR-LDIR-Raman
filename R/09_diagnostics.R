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

  # 5. Residual vector field (spatial pattern of offsets)
  plots$residual_vectors <- plot_residual_vectors(match_result)

  # 6. Confusion matrix heatmap
  if (!is.null(agreement) && nrow(match_result$matched) > 0) {
    plots$confusion <- plot_confusion_matrix(agreement)
  }

  # 7. Size comparison for matched pairs
  if (nrow(match_result$matched) > 0) {
    plots$size_comparison <- plot_size_comparison(match_result)
  }

  # 8. Ambiguity distribution
  if (nrow(match_result$matched) > 0 && "ambiguity_score" %in% names(match_result$matched)) {
    plots$ambiguity <- plot_ambiguity(match_result)
  }

  # 9. Tiered agreement bar chart
  if (!is.null(agreement) && !is.null(agreement$tiered_rates) &&
      agreement$tiered_rates$n_total > 0) {
    plots$tiered_agreement <- plot_tiered_agreement(agreement)
  }

  log_message("  Generated ", length(plots), " diagnostic plots")
  plots
}


#' Overlay plot of aligned particles from two instruments
plot_overlay <- function(ftir_df, raman_df, match_result,
                         ftir_color = "tomato", raman_color = "steelblue",
                         src_label = "ftir") {
  matched <- match_result$matched

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = raman_df,
      ggplot2::aes(x = x_norm, y = y_norm),
      color = raman_color, alpha = 0.4, size = 1.5, shape = 16
    ) +
    ggplot2::geom_point(
      data = ftir_df,
      ggplot2::aes(x = x_aligned, y = y_aligned),
      color = ftir_color, alpha = 0.6, size = 2, shape = 17
    )

  # Draw lines connecting matched pairs (instrument-generic column names)
  if (nrow(matched) > 0) {
    src_x_col <- paste0(src_label, "_x_aligned")
    src_y_col <- paste0(src_label, "_y_aligned")
    ref_x_col <- "raman_x_norm"
    ref_y_col <- "raman_y_norm"
    if (all(c(src_x_col, src_y_col, ref_x_col, ref_y_col) %in% names(matched))) {
      segments_df <- data.frame(
        x    = matched[[src_x_col]],
        y    = matched[[src_y_col]],
        xend = matched[[ref_x_col]],
        yend = matched[[ref_y_col]]
      )
      p <- p + ggplot2::geom_segment(
        data = segments_df,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "grey40", alpha = 0.3, linewidth = 0.3
      )
    }
  }

  p <- p +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "FTIR\u2013Raman Particle Overlay (aligned coordinates)",
      subtitle = paste0("Red triangles = FTIR, Blue circles = Raman, ",
                        "Lines = matched pairs"),
      x = "X (\u00b5m, normalized)", y = "Y (\u00b5m, normalized)"
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
                        round(median(matched$match_distance), 2), " \u00b5m, ",
                        "n = ", nrow(matched), " pairs"),
      x = "Match distance (\u00b5m)", y = "Count"
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
      x = "Iteration", y = "RMS error (\u00b5m)"
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
      x = "dX (\u00b5m)", y = "dY (\u00b5m)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Residual vector field — spatial plot showing direction and magnitude of offsets
#'
#' Shows arrows at each matched FTIR particle pointing toward its Raman match.
#' Systematic patterns (all arrows pointing one direction in a region) indicate
#' local distortion that a rigid transform cannot capture.
plot_residual_vectors <- function(match_result, magnification = 5) {
  matched <- match_result$matched

  if (nrow(matched) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No matches") +
             ggplot2::theme_void())
  }

  df <- data.frame(
    x   = matched$ftir_x_aligned,
    y   = matched$ftir_y_aligned,
    dx  = matched$raman_x_norm - matched$ftir_x_aligned,
    dy  = matched$raman_y_norm - matched$ftir_y_aligned,
    dist = matched$match_distance
  )
  df <- df[!is.na(df$x) & !is.na(df$y), ]

  if (nrow(df) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid residuals") +
             ggplot2::theme_void())
  }

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y,
                   xend = x + dx * magnification,
                   yend = y + dy * magnification,
                   colour = dist),
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")),
      linewidth = 0.5, alpha = 0.7
    ) +
    ggplot2::scale_colour_viridis_c(name = "Distance (\u00b5m)") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Residual Vector Field",
      subtitle = paste0("Arrows magnified ", magnification, "x | ",
                        "Systematic patterns indicate local distortion"),
      x = "X (\u00b5m)", y = "Y (\u00b5m)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Assess whether thin-plate spline (TPS) warp is justified
#'
#' Splits the matched pairs into spatial quadrants and checks if residual
#' vectors differ systematically between quadrants.
#'
#' @param match_result Result from match_particles()
#' @return List with systematic_distortion flag, quadrant residuals, recommend_tps
assess_tps_need <- function(match_result) {
  matched <- match_result$matched
  if (nrow(matched) < 20) {
    return(list(systematic_distortion = FALSE, recommend_tps = FALSE,
                quadrant_residuals = NULL, message = "Too few matches"))
  }

  df <- data.frame(
    x  = matched$ftir_x_aligned,
    y  = matched$ftir_y_aligned,
    dx = matched$raman_x_norm - matched$ftir_x_aligned,
    dy = matched$raman_y_norm - matched$ftir_y_aligned
  )
  df <- df[complete.cases(df), ]

  df$qx <- ifelse(df$x > median(df$x), "R", "L")
  df$qy <- ifelse(df$y > median(df$y), "T", "B")
  df$quad <- paste0(df$qy, df$qx)

  quad_means <- aggregate(cbind(dx, dy) ~ quad, data = df, FUN = mean)
  dx_range <- diff(range(quad_means$dx))
  dy_range <- diff(range(quad_means$dy))

  list(
    systematic_distortion = dx_range > 10 | dy_range > 10,
    quadrant_residuals    = quad_means,
    dx_range_um           = round(dx_range, 2),
    dy_range_um           = round(dy_range, 2),
    recommend_tps         = dx_range > 20 | dy_range > 20,
    message               = paste0("Quadrant dx range: ", round(dx_range, 1),
                                   " \u00b5m, dy range: ", round(dy_range, 1), " \u00b5m")
  )
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

  conf_df <- as.data.frame(conf)
  names(conf_df) <- c("Instrument_A", "Instrument_B", "Count")

  # Labels from agreement object
  lab_a <- agreement$instrument_a %||% "FTIR"
  lab_b <- agreement$instrument_b %||% "Raman"

  # Add tiered info to subtitle if available
  subtitle_text <- paste0("Overall agreement: ",
                          round(agreement$agreement_rate * 100, 1), "%")
  if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
    tr <- agreement$tiered_rates
    subtitle_text <- paste0(
      "Exact: ", tr$exact_pct, "% | Family: ", tr$family_pct,
      "% | Family+: ", tr$family_or_better_pct, "%"
    )
  }

  p <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Instrument_B, y = Instrument_A,
                                              fill = Count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(
      title = paste0("Material Confusion Matrix (", lab_a, " vs ", lab_b, ")"),
      subtitle = subtitle_text,
      x = paste0(lab_b, " material"), y = paste0(lab_a, " material")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' Scatter plot comparing particle sizes between instruments
plot_size_comparison <- function(match_result) {
  matched <- match_result$matched

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
      x = "FTIR Feret max (\u00b5m)", y = "Raman Feret max (\u00b5m)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Ambiguity score distribution
plot_ambiguity <- function(match_result) {
  matched <- match_result$matched

  p <- ggplot2::ggplot(matched, ggplot2::aes(x = ambiguity_score)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white",
                            alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0.7, linetype = "dashed", color = "tomato") +
    ggplot2::labs(
      title = "Match Ambiguity Score Distribution",
      subtitle = paste0("Flagged (>0.7): ",
                        sum(matched$ambiguity_flag, na.rm = TRUE),
                        " of ", nrow(matched), " pairs"),
      x = "Ambiguity score (match dist / next best dist)", y = "Count"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Tiered agreement bar chart
plot_tiered_agreement <- function(agreement) {
  if (is.null(agreement$agreement_detail) || nrow(agreement$agreement_detail) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") +
             ggplot2::theme_void())
  }

  detail <- agreement$agreement_detail

  # Reshape for stacked bar
  bar_data <- data.frame(
    material = rep(detail$material_a, 3),
    tier     = rep(c("Exact", "Family", "Disagree"), each = nrow(detail)),
    count    = c(detail$n_exact, detail$n_family, detail$n_disagree),
    stringsAsFactors = FALSE
  )
  bar_data$tier <- factor(bar_data$tier, levels = c("Disagree", "Family", "Exact"))

  p <- ggplot2::ggplot(bar_data, ggplot2::aes(x = material, y = count, fill = tier)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = c(
      Exact = "#2166ac", Family = "#92c5de", Disagree = "#f4a582"
    )) +
    ggplot2::labs(
      title = "Tiered Agreement by Material",
      subtitle = paste0("Overall: Exact ",
                        agreement$tiered_rates$exact_pct, "%, Family+ ",
                        agreement$tiered_rates$family_or_better_pct, "%"),
      x = "Material", y = "Count", fill = "Agreement"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' Generate LDIR-specific diagnostic plots
#'
#' @param ldir_aligned LDIR data frame with x_aligned, y_aligned
#' @param raman_df Raman data frame with x_norm, y_norm
#' @param ldir_raman_match Match result for LDIR-Raman
#' @param ldir_agreement Agreement result for LDIR-Raman
#' @return List of ggplot objects
generate_ldir_diagnostics <- function(ldir_aligned, raman_df,
                                      ldir_raman_match, ldir_agreement) {
  plots <- list()

  # LDIR-Raman overlay
  if (nrow(ldir_aligned) > 0 && nrow(raman_df) > 0) {
    plots$ldir_overlay <- plot_overlay(
      ldir_aligned, raman_df, ldir_raman_match,
      ftir_color = "darkgreen", raman_color = "steelblue",
      src_label = "ldir"
    )
    plots$ldir_overlay <- plots$ldir_overlay +
      ggplot2::labs(
        title = "LDIR\u2013Raman Particle Overlay (aligned coordinates)",
        subtitle = "Green triangles = LDIR, Blue circles = Raman"
      )
  }

  # LDIR-Raman confusion
  if (!is.null(ldir_agreement) && !is.null(ldir_agreement$confusion_matrix)) {
    plots$ldir_confusion <- plot_confusion_matrix(ldir_agreement)
  }

  plots
}


# Null-coalescing operator (if not already defined)
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
