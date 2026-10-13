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
  out <- do.call(rbind, rows)
  out[order(out$kind, out$key), , drop = FALSE]
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) != 2) stop("usage: Rscript tools/regression_fingerprint.R <run_dir> <out.csv>")
  utils::write.csv(regression_fingerprint(a[1]), a[2], row.names = FALSE)
}
