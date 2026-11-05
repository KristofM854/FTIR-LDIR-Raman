# =============================================================================
# ldir_join_check.R — were the LDIR table rows attached to the right positions?
# =============================================================================
# The LDIR export has no coordinates; the pipeline attaches each table row
# (material, size) to a particle found in the image by size (join_ldir_coords).
# A wrong join moves materials and sizes to the wrong place without any error.
# Raman measured the same particles independently, so after alignment and
# matching, a correct join shows up as LDIR table values that agree with the
# Raman partner far more often than a random join would.
# Used by the LDIR stage (logged, written to 07_diagnostics/) and by
# tools/validation.R.
# =============================================================================


#' Were the LDIR table rows attached to the right image positions?
#'
#' For each LDIR particle matched to Raman, compare what the LDIR table says
#' (family, size) with what Raman found at that position. A correct join puts
#' agreeing particles together; re-assigning the LDIR rows at random among the
#' joined positions is the null. Two statistics:
#'   family   share of scorable pairs whose families agree (a sample of a
#'            single polymer makes this uninformative: the null is ~100% too)
#'   size     Spearman correlation of LDIR (table) and Raman Feret max. The
#'            join itself matches sizes, so this shows the rows sit on
#'            particles of the right size, not that they sit on the right
#'            particle (a swap between similar sizes passes)
#'   circularity  Spearman correlation of LDIR (table) and Raman circularity:
#'            not used by the join, so it tests the identity of the row;
#'            noisy on small particles (few pixels)
#' The family statistic is informative only when a perfect join could reach
#' clearly more than a random one (max_possible at least 0.2 above the
#' random mean).
#'
#' @param m Matched LDIR->Raman pairs (match_particles()$matched)
#' @param n_iter Permutations
#' @param seed RNG seed (the caller's RNG state is restored)
#' @param ref_prefix Column prefix of the partner in m ("raman"; "ftir" or
#'   "ftir_bruker" in a run without Raman)
#' @return list(summary = data.frame, pairs = per-pair table), or NULL with
#'   fewer than 5 pairs
ldir_join_check <- function(m, n_iter = 1000, seed = 2026, ref_prefix = "raman") {
  if (is.null(m) || nrow(m) < 5) return(NULL)
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  r <- function(k) m[[paste0(ref_prefix, "_", k)]]
  fam_l <- classify_family_vec(m$ldir_material); fam_r <- classify_family_vec(r("material"))
  size_l <- m$ldir_feret_max_um; size_r <- r("feret_max_um")
  ok_f <- !is_unidentified(m$ldir_material) & !is_unidentified(r("material")) &
          !fam_l %in% unresolved_families & !fam_r %in% unresolved_families
  ok_s <- is.finite(size_l) & is.finite(size_r)
  fam_stat  <- function(fl) mean(fl[ok_f] == fam_r[ok_f])
  size_stat <- function(sl) suppressWarnings(stats::cor(sl[ok_s], size_r[ok_s], method = "spearman"))
  # Circularity is not used by the join (which uses size rank and, gated,
  # eccentricity), so unlike size it cannot agree merely because the join
  # matched sizes; it needs a circularity on both sides (Raman exports one)
  col_or_na <- function(v) if (is.null(v)) rep(NA_real_, nrow(m)) else suppressWarnings(as.numeric(v))
  circ_l <- col_or_na(m$ldir_circularity); circ_r <- col_or_na(r("circularity"))
  ok_c <- is.finite(circ_l) & is.finite(circ_r)
  circ_stat <- function(cl) suppressWarnings(stats::cor(cl[ok_c], circ_r[ok_c], method = "spearman"))
  # Null: the LDIR rows (family, size, shape together) re-assigned to the positions
  perm <- replicate(n_iter, sample.int(nrow(m)))
  null_f <- apply(perm, 2, function(i) fam_stat(fam_l[i]))
  null_s <- apply(perm, 2, function(i) size_stat(size_l[i]))
  null_c <- if (sum(ok_c) > 9) apply(perm, 2, function(i) circ_stat(circ_l[i])) else NA_real_
  p_val <- function(obs, null) (1 + sum(null >= obs, na.rm = TRUE)) / (1 + sum(!is.na(null)))
  obs_f <- if (sum(ok_f) > 0) fam_stat(fam_l) else NA_real_
  obs_s <- if (sum(ok_s) > 2) size_stat(size_l) else NA_real_
  obs_c <- if (sum(ok_c) > 9) circ_stat(circ_l) else NA_real_
  # The highest family agreement any join could reach with these families:
  # each family agrees at most as often as the rarer side names it. When the
  # two instruments name different families (Raman reference chemicals
  # against LDIR natural materials), even a perfect join stays near the
  # random baseline and the statistic can neither confirm nor reject the join.
  max_f <- if (sum(ok_f) > 0) {
    tl <- table(fam_l[ok_f]); tr <- table(fam_r[ok_f]); f <- intersect(names(tl), names(tr))
    sum(pmin(tl[f], tr[f])) / sum(ok_f)
  } else NA_real_
  mean_f <- mean(null_f, na.rm = TRUE)
  summary <- data.frame(
    statistic = c("family_agreement", "size_spearman", "circularity_spearman"),
    n = c(sum(ok_f), sum(ok_s), sum(ok_c)),
    observed = round(c(obs_f, obs_s, obs_c), 3),
    null_mean = round(c(mean_f, mean(null_s, na.rm = TRUE), mean(null_c, na.rm = TRUE)), 3),
    null_q95 = round(c(stats::quantile(null_f, 0.95, na.rm = TRUE),
                       stats::quantile(null_s, 0.95, na.rm = TRUE),
                       if (all(is.na(null_c))) NA else stats::quantile(null_c, 0.95, na.rm = TRUE)), 3),
    max_possible = round(c(max_f, 1, 1), 3),
    p_value = signif(c(p_val(obs_f, null_f), p_val(obs_s, null_s),
                       if (is.na(obs_c)) NA else p_val(obs_c, null_c)), 3),
    informative = c(isTRUE(max_f - mean_f >= 0.2), TRUE, !is.na(obs_c)))
  pairs <- data.frame(
    ldir_particle_id = m$ldir_particle_id, ref_particle_id = r("particle_id"),
    match_distance_um = m$match_distance,
    ldir_family = fam_l, ref_family = fam_r,
    family_agrees = ifelse(ok_f, fam_l == fam_r, NA),
    ldir_feret_max_um = size_l, ref_feret_max_um = size_r,
    size_ratio = round(size_l / size_r, 2),
    ldir_circularity = circ_l, ref_circularity = circ_r,
    coord_join_warning = if ("ldir_coord_join_warning" %in% names(m)) m$ldir_coord_join_warning else NA,
    stringsAsFactors = FALSE)
  names(pairs) <- sub("^ref_", paste0(ref_prefix, "_"), names(pairs))
  list(summary = summary, pairs = pairs)
}


#' Log the join check in one or two lines; WARN when the join looks random.
log_ldir_join_check <- function(chk, ref_name = "Raman") {
  if (is.null(chk)) {
    log_message("  LDIR join check: fewer than 5 LDIR-", ref_name, " pairs \u2014 not run")
    return(invisible(NULL))
  }
  s <- chk$summary
  sz <- s[s$statistic == "size_spearman", ]
  fm <- s[s$statistic == "family_agreement", ]
  ci <- s[s$statistic == "circularity_spearman", ]
  p_txt <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else format(signif(p, 2))
  circ_txt <- if (nrow(ci) && isTRUE(ci$informative))
    sprintf("; circularity Spearman %.2f (random join 95%% %.2f, p %s)", ci$observed, ci$null_q95,
            p_txt(ci$p_value)) else ""
  fam_txt <- if (isTRUE(fm$informative)) sprintf(", p %s", p_txt(fm$p_value)) else
    if (isTRUE(fm$null_mean >= 0.8)) ": not informative for this sample (one family dominates: a random join agrees as well)" else
    sprintf(paste0(": not informative for this sample (a perfect join could reach only %.0f%%: ",
                   "the two instruments name different families)"), 100 * fm$max_possible)
  log_message(sprintf(paste0(
    "  LDIR join check (LDIR table vs %s partner, %d pairs): size Spearman %.2f ",
    "(random join %.2f, 95%% %.2f, p %s)%s; family agreement %.0f%% (random join %.0f%%%s)"),
    ref_name, sz$n, sz$observed, sz$null_mean, sz$null_q95, p_txt(sz$p_value), circ_txt,
    100 * fm$observed, 100 * fm$null_mean, fam_txt))
  if (nrow(ci) && isTRUE(ci$informative) && isTRUE(sz$p_value < 0.05) && !isTRUE(ci$p_value < 0.05))
    log_message("  LDIR join check: sizes agree but circularity (not used by the join) does not: ",
                "the rows may sit on particles of the right size but not always the right ",
                "particle -- interpret LDIR materials of small particles with care", level = "WARN")
  sig <- function(r) isTRUE(r$p_value < 0.05)
  if (!sig(sz) && !(isTRUE(fm$informative) && sig(fm)) && !(nrow(ci) && isTRUE(ci$informative) && sig(ci)))
    log_message("  LDIR join check: the LDIR table values agree with ", ref_name, " no better ",
                "than a random join \u2014 check the image join (ldir_joined_raw.csv, ",
                "07_diagnostics/ldir_join_pairs.csv) before interpreting LDIR results",
                level = "WARN")
  invisible(chk)
}
