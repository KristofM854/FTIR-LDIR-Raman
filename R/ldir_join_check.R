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
#'            matcher's size term favours size-similar partners, so the
#'            observed value is biased upward; the null keeps the matched
#'            positions and only re-assigns the table rows.
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
  # Null: the LDIR rows (family and size together) re-assigned to the positions
  perm <- replicate(n_iter, sample.int(nrow(m)))
  null_f <- apply(perm, 2, function(i) fam_stat(fam_l[i]))
  null_s <- apply(perm, 2, function(i) size_stat(size_l[i]))
  p_val <- function(obs, null) (1 + sum(null >= obs, na.rm = TRUE)) / (1 + sum(!is.na(null)))
  obs_f <- if (sum(ok_f) > 0) fam_stat(fam_l) else NA_real_
  obs_s <- if (sum(ok_s) > 2) size_stat(size_l) else NA_real_
  summary <- data.frame(
    statistic = c("family_agreement", "size_spearman"),
    n = c(sum(ok_f), sum(ok_s)),
    observed = round(c(obs_f, obs_s), 3),
    null_mean = round(c(mean(null_f, na.rm = TRUE), mean(null_s, na.rm = TRUE)), 3),
    null_q95 = round(c(stats::quantile(null_f, 0.95, na.rm = TRUE),
                       stats::quantile(null_s, 0.95, na.rm = TRUE)), 3),
    p_value = signif(c(p_val(obs_f, null_f), p_val(obs_s, null_s)), 3),
    informative = c(mean(null_f, na.rm = TRUE) < 0.9, TRUE))
  pairs <- data.frame(
    ldir_particle_id = m$ldir_particle_id, ref_particle_id = r("particle_id"),
    match_distance_um = m$match_distance,
    ldir_family = fam_l, ref_family = fam_r,
    family_agrees = ifelse(ok_f, fam_l == fam_r, NA),
    ldir_feret_max_um = size_l, ref_feret_max_um = size_r,
    size_ratio = round(size_l / size_r, 2),
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
  p_txt <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else format(signif(p, 2))
  log_message(sprintf(paste0(
    "  LDIR join check (LDIR table vs %s partner, %d pairs): size Spearman %.2f ",
    "(random join %.2f, 95%% %.2f, p %s); family agreement %.0f%% (random join %.0f%%%s)"),
    ref_name, sz$n, sz$observed, sz$null_mean, sz$null_q95, p_txt(sz$p_value),
    100 * fm$observed, 100 * fm$null_mean,
    if (isTRUE(fm$informative)) paste0(", p ", p_txt(fm$p_value))
    else ": not informative for this sample"))
  sig <- function(r) isTRUE(r$p_value < 0.05)
  if (!sig(sz) && !(isTRUE(fm$informative) && sig(fm)))
    log_message("  LDIR join check: the LDIR table values agree with ", ref_name, " no better ",
                "than a random join \u2014 check the image join (ldir_joined_raw.csv, ",
                "07_diagnostics/ldir_join_pairs.csv) before interpreting LDIR results",
                level = "WARN")
  invisible(chk)
}
