# =============================================================================
# precision_study.R — repeatability and reproducibility of particle counts
# =============================================================================
# A precision study after ISO 5725-2: the same filter measured n times (e.g.
# 3) on each of p instruments (e.g. 4). The measurand is a particle count per
# filter -- in total, per polymer family and per size class -- so no particle
# has to be matched between measurements: a remounted, rotated filter makes no
# difference, and field samples work as well as reference standards.
#
# Per stratum (group x size class), with the instruments in the role of the
# ISO 5725 "laboratories" and the replicates as the repeat measurements:
#   s_r   repeatability standard deviation (within instrument)
#   s_L   between-instrument standard deviation
#   s_R   reproducibility standard deviation, s_R^2 = s_L^2 + s_r^2
#   r, R  repeatability / reproducibility limits, 2.8 s_r and 2.8 s_R
#   95% intervals for s_r (chi-square) and s_R (Satterthwaite)
#   Mandel's h (between) and k (within) consistency statistics with their
#   5% and 1% indicator levels (ISO 5725-2, 7.3.1)
# and, as counts are discrete, each instrument's replicate variance against
# Poisson counting noise (dispersion index s^2 / mean, chi-square test).
#
# Counting follows the main pipeline: particles >= min_particle_size_um; the
# polymer groups count only particles passing the instrument's quality cutoff
# (quality_passes()); size classes as in the detection rates.
#
# Depends on R/01_readers.R (read_instrument), R/align_helpers.R
# (quality_passes, with R/08_agreement.R), R/08b_material_map.R (families)
# and R/ftir_segments.R.
# =============================================================================

PRECISION_SIZE_BREAKS <- c(20, 50, 100, 200, Inf)
PRECISION_SIZE_LABELS <- c("20-50", "50-100", "100-200", ">200")
PRECISION_ROLES <- c(ftir = "FTIR (PerkinElmer)", ftir_bruker = "FTIR (Bruker)",
                     raman = "Raman", ldir = "LDIR")

#' The study design from a CSV (columns instrument, replicate, file) or a
#' data frame. Instrument = role key (ftir, ftir_bruker, raman, ldir; also
#' accepted: ftir_perkin, ftir_pe, perkinelmer, bruker) or display name.
read_precision_design <- function(design) {
  d <- if (is.character(design) && length(design) == 1)
    utils::read.csv(design, stringsAsFactors = FALSE, strip.white = TRUE) else design
  names(d) <- tolower(names(d))
  miss <- setdiff(c("instrument", "replicate", "file"), names(d))
  if (length(miss)) stop("Precision design needs the column(s) ", paste(miss, collapse = ", "),
                         call. = FALSE)
  key <- tolower(gsub("[^a-z]", "", tolower(d$instrument)))
  d$instrument <- ifelse(grepl("bruker", key), "ftir_bruker",
                  ifelse(grepl("raman", key), "raman",
                  ifelse(grepl("ldir", key), "ldir",
                  ifelse(grepl("ftir|perkin", key), "ftir", NA_character_))))
  if (anyNA(d$instrument)) stop("Precision design: unknown instrument(s) in rows ",
                                paste(which(is.na(d$instrument)), collapse = ", "), call. = FALSE)
  if (anyDuplicated(paste(d$instrument, d$replicate)))
    stop("Precision design: an instrument lists the same replicate twice", call. = FALSE)
  d$replicate <- as.character(d$replicate)
  d
}

#' Particles of one measurement. A PerkinElmer segment table is read with all
#' segments of its scan (coordinates do not matter for counts).
.precision_read <- function(role, file, config) {
  files <- if (role == "ftir") ftir_segment_files(file, FTIR_SEGMENT_TABLE) else NULL
  if (length(files) < 2) files <- file
  do.call(rbind, lapply(files, function(f) {
    d <- read_instrument(role, config, path = f)
    d[, intersect(c("particle_id", "feret_max_um", "material", "material_family", "quality"),
                  names(d)), drop = FALSE]
  }))
}

#' Counts per measurement and stratum.
#'
#' Groups: "All particles" (>= the minimum size, any quality), "Plastics"
#' (synthetic + semi-synthetic families) and every polymer family found, the
#' last two among particles passing the quality cutoff. Size classes: "all
#' sizes" and PRECISION_SIZE_LABELS. Strata missing in a measurement count 0.
#'
#' @param design read_precision_design() result
#' @param particles Optional named list (one entry per design row, named
#'   "<instrument>|<replicate>") of particle tables, in place of reading files
#' @return data.frame: instrument, replicate, group, size_class, count
precision_counts <- function(design, config, particles = NULL) {
  min_size <- cfg_value(config, "min_particle_size_um")
  plastic_fams <- c(synthetic_families, semi_synthetic_families)
  rows <- list()
  for (i in seq_len(nrow(design))) {
    role <- design$instrument[i]; rep_id <- design$replicate[i]
    key <- paste(role, rep_id, sep = "|")
    df <- if (!is.null(particles)) particles[[key]] else {
      log_message("Precision study: ", PRECISION_ROLES[[role]], " replicate ", rep_id,
                  ": ", basename(design$file[i]))
      .precision_read(role, design$file[i], config)
    }
    if (is.null(df)) stop("No particles for ", key, call. = FALSE)
    df <- df[!is.na(df$feret_max_um) & df$feret_max_um >= min_size, , drop = FALSE]
    fam <- df$material_family %||% classify_family_vec(df$material)
    qok <- quality_passes(df, role, config)
    size <- as.character(cut(df$feret_max_um, PRECISION_SIZE_BREAKS, PRECISION_SIZE_LABELS,
                             right = FALSE))
    groups <- c(list("All particles" = rep(TRUE, nrow(df)),
                     "Plastics" = qok & fam %in% plastic_fams),
                stats::setNames(lapply(sort(unique(fam[qok & !fam %in% "Unknown"])),
                                       function(f) qok & fam == f),
                                sort(unique(fam[qok & !fam %in% "Unknown"]))))
    for (g in names(groups)) for (sc in c("all sizes", PRECISION_SIZE_LABELS))
      rows[[length(rows) + 1]] <- data.frame(
        instrument = role, replicate = rep_id, group = g, size_class = sc,
        count = sum(groups[[g]] & (sc == "all sizes" | size == sc), na.rm = TRUE),
        stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  # A family found in one measurement only still counts 0 in the others
  grid <- merge(unique(out[, c("group", "size_class")]),
                unique(design[, c("instrument", "replicate")]), by = NULL)
  out <- merge(grid, out, all.x = TRUE)
  out$count[is.na(out$count)] <- 0L
  out[order(match(out$group, unique(c("All particles", "Plastics", sort(unique(out$group))))),
            match(out$size_class, c("all sizes", PRECISION_SIZE_LABELS)),
            match(out$instrument, names(PRECISION_ROLES)), out$replicate), ]
}

#' Mandel's h and k indicator levels (ISO 5725-2, 7.3.1).
#' @param p Number of instruments; @param n Replicates per instrument
mandel_critical <- function(p, n, alpha = c(0.05, 0.01)) {
  h <- vapply(alpha, function(a) {
    if (p < 3) return(NA_real_)
    t <- stats::qt(1 - a / 2, p - 2)
    (p - 1) * t / sqrt(p * (t^2 + p - 2))
  }, 0)
  # k^2 = p / (1 + (p - 1) / F), F = one instrument's variance over the pooled
  # others ~ F(n - 1, (p - 1)(n - 1)); p = 4, n = 3: 1.59 (5%), 1.77 (1%),
  # confirmed by simulation (tests/testthat/test-precision-study.R)
  k <- vapply(alpha, function(a) {
    if (p < 2 || n < 2) return(NA_real_)
    f <- stats::qf(1 - a, n - 1, (p - 1) * (n - 1))
    sqrt(p / (1 + (p - 1) / f))
  }, 0)
  list(h_5 = h[1], h_1 = h[2], k_5 = k[1], k_1 = k[2])
}

#' ISO 5725-2 statistics for one stratum.
#' @param y Counts; @param lab Instrument of each count
#' @return list(stats = one-row data.frame, labs = per-instrument data.frame)
iso5725_stratum <- function(y, lab) {
  lab <- factor(lab, levels = unique(lab))
  ni <- tapply(y, lab, length); mi <- tapply(y, lab, mean)
  si <- tapply(y, lab, function(v) if (length(v) > 1) stats::sd(v) else NA_real_)
  ok <- ni >= 2
  p <- sum(ok); N <- sum(ni[ok])
  m <- sum(ni[ok] * mi[ok]) / N
  na_stats <- data.frame(n_instruments = p, n_measurements = N, mean = m,
                         s_r = NA_real_, s_r_low = NA_real_, s_r_high = NA_real_,
                         s_L = NA_real_, s_R = NA_real_, s_R_low = NA_real_, s_R_high = NA_real_,
                         r_limit = NA_real_, R_limit = NA_real_,
                         cv_r_pct = NA_real_, cv_L_pct = NA_real_, cv_R_pct = NA_real_)
  labs <- data.frame(instrument = names(ni), n = as.integer(ni), mean = as.numeric(mi),
                     sd = as.numeric(si), stringsAsFactors = FALSE)
  if (p < 2 || N - p < 1) return(list(stats = na_stats, labs = labs))

  df_r <- N - p
  s_r2 <- sum((ni[ok] - 1) * si[ok]^2) / df_r
  ms_b <- sum(ni[ok] * (mi[ok] - m)^2) / (p - 1)          # = n * s_d^2 when balanced
  n_bar <- (N - sum(ni[ok]^2) / N) / (p - 1)
  s_L2 <- max(0, (ms_b - s_r2) / n_bar)
  s_R2 <- s_L2 + s_r2
  ci_sd <- function(s2, df) sqrt(s2 * df / stats::qchisq(c(0.975, 0.025), df))
  ci_r <- if (s_r2 > 0) ci_sd(s_r2, df_r) else c(0, 0)
  # Satterthwaite: s_R^2 = ms_b / n_bar + (1 - 1/n_bar) s_r^2
  a <- ms_b / n_bar; b <- (1 - 1 / n_bar) * s_r2
  df_R <- if (s_L2 > 0) (a + b)^2 / (a^2 / (p - 1) + b^2 / df_r) else df_r
  ci_R <- if (s_R2 > 0) ci_sd(s_R2, df_R) else c(0, 0)
  pct <- function(s) if (m > 0) 100 * s / m else NA_real_
  stats <- data.frame(n_instruments = p, n_measurements = N, mean = m,
                      s_r = sqrt(s_r2), s_r_low = ci_r[1], s_r_high = ci_r[2],
                      s_L = sqrt(s_L2), s_R = sqrt(s_R2), s_R_low = ci_R[1], s_R_high = ci_R[2],
                      r_limit = 2.8 * sqrt(s_r2), R_limit = 2.8 * sqrt(s_R2),
                      cv_r_pct = pct(sqrt(s_r2)), cv_L_pct = pct(sqrt(s_L2)),
                      cv_R_pct = pct(sqrt(s_R2)))

  # Mandel's h and k (balanced forms; the mean n for unequal replicates)
  sd_means <- stats::sd(mi[ok])
  labs$h <- ifelse(ok & sd_means > 0, (mi - mean(mi[ok])) / sd_means, NA_real_)
  labs$k <- ifelse(ok & s_r2 > 0, si / sqrt(mean(si[ok]^2)), NA_real_)
  crit <- mandel_critical(p, round(mean(ni[ok])))
  flag <- function(v, c5, c1) ifelse(is.na(v) | is.na(c5), "",
                               ifelse(abs(v) > c1, "1%", ifelse(abs(v) > c5, "5%", "")))
  labs$h_flag <- flag(labs$h, crit$h_5, crit$h_1)
  labs$k_flag <- flag(labs$k, crit$k_5, crit$k_1)
  # Counting noise: a Poisson count has variance = mean
  labs$dispersion <- ifelse(labs$mean > 0, labs$sd^2 / labs$mean, NA_real_)
  labs$dispersion_p <- ifelse(labs$mean > 0 & labs$n > 1,
                              stats::pchisq((labs$n - 1) * labs$dispersion, labs$n - 1,
                                            lower.tail = FALSE), NA_real_)
  attr(labs, "critical") <- crit
  list(stats = stats, labs = labs)
}

#' The ISO 5725 statistics for every stratum of precision_counts().
#' @return list(iso = per-stratum table, instruments = per-stratum x instrument
#'   table, critical = Mandel indicator levels)
precision_iso5725 <- function(counts) {
  keys <- unique(counts[, c("group", "size_class")])
  iso <- list(); labs <- list(); crit <- NULL
  for (i in seq_len(nrow(keys))) {
    s <- counts[counts$group == keys$group[i] & counts$size_class == keys$size_class[i], ]
    r <- iso5725_stratum(s$count, s$instrument)
    iso[[i]] <- cbind(keys[i, , drop = FALSE], r$stats)
    labs[[i]] <- cbind(keys[rep(i, nrow(r$labs)), , drop = FALSE], r$labs)
    crit <- crit %||% attr(r$labs, "critical")
  }
  iso <- do.call(rbind, iso); rownames(iso) <- NULL
  labs <- do.call(rbind, labs); rownames(labs) <- NULL
  list(iso = iso, instruments = labs, critical = crit)
}

#' One-call driver: counts and statistics.
run_precision_study <- function(design, config, particles = NULL) {
  design <- read_precision_design(design)
  counts <- precision_counts(design, config, particles)
  st <- precision_iso5725(counts)
  list(design = design, counts = counts, iso = st$iso, instruments = st$instruments,
       critical = st$critical)
}
