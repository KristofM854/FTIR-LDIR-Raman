# =============================================================================
# rate_intervals.R — confidence intervals for proportions
# =============================================================================
# Agreement and detection rates come from tens of pairs, sometimes fewer: a
# rate without its interval cannot be compared between instruments or samples.
# Shared by the pipeline and the viewer (shiny_app/global.R sources it).
# =============================================================================

#' Wilson score interval for a proportion k / n.
#'
#' Unlike the normal approximation it stays inside 0-1 and is usable for small
#' n and for rates of 0% or 100%.
#'
#' @param k Successes; @param n Trials (vectors allowed)
#' @param level Confidence level
#' @return data.frame(low, high) as fractions; NA where n is 0
wilson_interval <- function(k, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- ifelse(n > 0, k / n, NA_real_)
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  data.frame(low = pmax(0, centre - half), high = pmin(1, centre + half))
}

#' "58% (45-70%)": a rate with its 95% Wilson interval, for logs and tables.
#' "-" when n is 0. ASCII only: it goes to logs, CSVs and the PDF in any locale.
format_rate_ci <- function(k, n) {
  ci <- wilson_interval(k, n)
  ifelse(n > 0,
         sprintf("%.0f%% (%.0f-%.0f%%)", 100 * k / n, 100 * ci$low, 100 * ci$high),
         "-")
}
