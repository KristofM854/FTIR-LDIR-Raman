# =============================================================================
# regression_fingerprint.R -- compact, comparable summary of a pipeline run
# =============================================================================
# Used by tests/testthat/test-example-dataset-A.R to compare a fresh run of the
# example dataset with the stored reference, and to (re)generate that reference:
#
#   Rscript tools/regression_fingerprint.R <run_dir> tests/testthat/fixtures/dataset_A_fingerprint.csv
#
# Rows: kind, key, value
#   n_rows     <stage file>                      number of rows
#   pair       <match file>|<id a>|<id b>        match distance (um)
#   transform  <pair>:<field>                    scale, rotation_deg, tx, ty
#   agreement  <agreement file>:<tier>           number of pairs in the tier
#   p_true, n_additional   per pair / match file   pair confidence, additional pairs
#   chance     <pair>:<field>                    gate, expected chance pairs, sigma, determined
#   detection, detection_pct  <of>|<by>|<group>|<size>:<field>   counts, chance-corrected rate
#   partner_status  <of>|<by>|<status>           number of particles
#   alignment, alignment_log10p  <instrument>    determined, log10 of the pose p
#   ldir_position  <LDIR id>:x_um / y_um         where the join placed each table row
# =============================================================================

regression_fingerprint <- function(run_dir) {
  rows <- list()
  add <- function(kind, key, value)
    rows[[length(rows) + 1]] <<- data.frame(kind = kind, key = key, value = value,
                                            stringsAsFactors = FALSE)
  for (sub in c("01_ingested", "03_prefiltered", "05_matches")) {
    for (f in sort(list.files(file.path(run_dir, sub), pattern = "\\.csv$"))) {
      d <- utils::read.csv(file.path(run_dir, sub, f), stringsAsFactors = FALSE)
      add("n_rows", file.path(sub, f), nrow(d))
      if (startsWith(f, "matched_") && nrow(d) > 0) {
        ids <- grep("_particle_id$", names(d), value = TRUE)[1:2]
        for (i in seq_len(nrow(d)))
          add("pair", paste(f, d[[ids[1]]][i], d[[ids[2]]][i], sep = "|"), d$match_distance[i])
      }
    }
  }
  for (f in sort(list.files(file.path(run_dir, "04_alignment"), pattern = "^transform_params_.*\\.txt$"))) {
    L <- readLines(file.path(run_dir, "04_alignment", f))
    for (k in c("scale", "rotation_deg", "tx", "ty")) {
      v <- sub(paste0("^", k, ":\\s*"), "", grep(paste0("^", k, ":"), L, value = TRUE))
      if (length(v)) add("transform", paste0(sub("\\.txt$", "", f), ":", k), as.numeric(v[1]))
    }
  }
  for (f in sort(list.files(file.path(run_dir, "06_agreement"), pattern = "^agreement_pairwise_.*\\.csv$"))) {
    d <- utils::read.csv(file.path(run_dir, "06_agreement", f), stringsAsFactors = FALSE)
    for (t in c("Exact", "Family", "Filler/pigment", "Unresolved", "Disagree"))
      add("agreement", paste0(f, ":", t), sum(d$tier == t))
  }
  # --- Statistics and the LDIR join (not only the pairs) ---
  rd <- function(...) {
    f <- file.path(run_dir, ...)
    if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  }
  # Pair confidence of every pair
  for (f in sort(list.files(file.path(run_dir, "05_matches"), pattern = "^matched_.*\\.csv$"))) {
    d <- rd("05_matches", f)
    if (is.null(d) || !nrow(d) || is.null(d$p_true)) next
    ids <- grep("_particle_id$", names(d), value = TRUE)[1:2]
    ok <- !is.na(d$p_true)
    for (i in which(ok))
      add("p_true", paste(f, d[[ids[1]]][i], d[[ids[2]]][i], sep = "|"), d$p_true[i])
    if ("pair_type" %in% names(d))
      add("n_additional", f, sum(d$pair_type == "additional"))
  }
  # Chance pairs, gate and positional error per pair
  ch <- rd("07_diagnostics", "chance_pairs.csv")
  if (!is.null(ch)) for (i in seq_len(nrow(ch)))
    for (k in c("gate_um", "chance_pairs_expected", "position_sigma_um", "alignment_determined")) {
      v <- suppressWarnings(as.numeric(ch[[k]][i]))
      if (!is.null(ch[[k]]) && is.finite(v)) add("chance", paste0(ch$pair[i], ":", k), v)
    }
  # Detection rates: counts and the chance-corrected identified rate
  det <- rd("07_diagnostics", "detection_rates.csv")
  if (!is.null(det)) for (i in seq_len(nrow(det))) {
    key <- paste(det$of[i], det$found_by[i], det$group[i], det$size_class[i], sep = "|")
    for (k in c("n_base", "n_found", "n_detected"))
      if (is.finite(det[[k]][i])) add("detection", paste0(key, ":", k), det[[k]][i])
    if (is.finite(det$rate_true_pct[i])) add("detection_pct", paste0(key, ":rate_true_pct"), det$rate_true_pct[i])
  }
  # Partner status counts
  ps <- rd("05_matches", "partner_status.csv")
  if (!is.null(ps) && nrow(ps)) {
    tb <- as.data.frame(table(of = ps$of, found_by = ps$found_by, status = ps$status),
                        stringsAsFactors = FALSE)
    tb <- tb[tb$Freq > 0, ]
    for (i in seq_len(nrow(tb)))
      add("partner_status", paste(tb$of[i], tb$found_by[i], tb$status[i], sep = "|"), tb$Freq[i])
  }
  # Alignment: determined, pose p (as log10: it spans many decades)
  al <- rd("04_alignment", "alignment_summary.csv")
  if (!is.null(al)) for (i in seq_len(nrow(al))) {
    add("alignment", paste0(al$instrument[i], ":determined"), as.numeric(isTRUE(al$determined[i])))
    if (is.finite(al$pose_p[i]) && al$pose_p[i] > 0)
      add("alignment_log10p", paste0(al$instrument[i], ":pose_p"), log10(al$pose_p[i]))
  }
  # LDIR join: where each table row was placed
  jn <- rd("02_joined", "ldir_joined_raw.csv")
  if (!is.null(jn) && all(c("particle_id", "x_um", "y_um") %in% names(jn)))
    for (i in which(is.finite(jn$x_um) & is.finite(jn$y_um))) {
      add("ldir_position", paste0(jn$particle_id[i], ":x_um"), jn$x_um[i])
      add("ldir_position", paste0(jn$particle_id[i], ":y_um"), jn$y_um[i])
    }

  out <- do.call(rbind, rows)
  out[order(out$kind, out$key), , drop = FALSE]
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) != 2) stop("usage: Rscript tools/regression_fingerprint.R <run_dir> <out.csv>")
  utils::write.csv(regression_fingerprint(a[1]), a[2], row.names = FALSE)
}
