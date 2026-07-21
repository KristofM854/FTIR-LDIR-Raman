# =============================================================================
# diagnose_multirun_placement.R — plan v3 diagnostic (docs/multirun_image_placement_plan.md)
# =============================================================================
# Compares the three places a Raman WITec calibration can live, for the
# dataset behind a Multi-Run (reproducibility) view:
#
#   [A] what the Multi-Run viewer actually uses:
#       output/reproducibility/<run>/reproducibility_meta.csv
#   [B] what the single Raman tab uses for the matching pipeline run:
#       output/<run_dir>/00_manifest/manifest.json -> config_snapshot.raman_image_*
#       (matched to the replicate raw file via inputs.raman.md5 / basename)
#   [C] what is currently in R/00_config.R (make_config()) — the value a NEW
#       pipeline run would bake in today.
#
# Usage:
#   Rscript tools/diagnose_multirun_placement.R [replicate_run1_file] [repro_output_dir]
#
#   replicate_run1_file  raw Raman file used as CONFIG$runs[[1]]$file in
#                        tools/reproducibility.R (enables exact md5 matching)
#   repro_output_dir     a specific output/reproducibility/<id> folder
#                        (default: the newest one with a reproducibility_meta.csv)
#
# Read-only: writes nothing, changes nothing.
# =============================================================================

.script_path <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) normalizePath(f) else NA_character_
}, error = function(e) NA_character_)
REPO_ROOT <- if (!is.na(.script_path)) dirname(dirname(.script_path)) else getwd()

if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("jsonlite is required to read manifest.json (install.packages('jsonlite'))")

args <- commandArgs(trailingOnly = TRUE)
replicate_file <- if (length(args) >= 1) args[1] else NULL
repro_dir      <- if (length(args) >= 2) args[2] else NULL

WITEC_FIELDS <- c("raman_image_width_um", "raman_image_height_um",
                  "raman_image_center_x_um", "raman_image_center_y_um")

`%||%` <- function(a, b) if (is.null(a)) b else a

fmt <- function(v) {
  if (is.null(v) || length(v) == 0 || all(is.na(v))) "  <absent/NA>"
  else format(as.numeric(v)[1], nsmall = 1, trim = TRUE)
}

hr <- function() cat(strrep("-", 72), "\n")

# --- [A] reproducibility_meta.csv (what Multi-Run uses) ----------------------
hr()
cat("[A] Multi-Run metadata (reproducibility_meta.csv)\n")
hr()
if (is.null(repro_dir)) {
  base <- file.path(REPO_ROOT, "output", "reproducibility")
  cands <- if (dir.exists(base)) list.dirs(base, recursive = FALSE) else character(0)
  cands <- cands[file.exists(file.path(cands, "reproducibility_meta.csv"))]
  if (length(cands) > 0) {
    repro_dir <- cands[order(file.mtime(file.path(cands, "reproducibility_meta.csv")),
                             decreasing = TRUE)][1]
    cat("Using newest reproducibility output:", repro_dir, "\n")
  }
}
meta <- NULL
if (!is.null(repro_dir) && file.exists(file.path(repro_dir, "reproducibility_meta.csv"))) {
  meta <- utils::read.csv(file.path(repro_dir, "reproducibility_meta.csv"),
                          stringsAsFactors = FALSE)
  cat("instrument:", meta$instrument[1], "  n_runs:", meta$n_runs[1], "\n")
  for (f in c(WITEC_FIELDS, "raman_um_per_px"))
    cat(sprintf("  %-26s %s\n", f, fmt(meta[[f]])))
  if ("raman_calibration_source" %in% names(meta))
    cat("  raman_calibration_source:  ", meta$raman_calibration_source[1], "\n", sep = "")
  all_na <- all(vapply(WITEC_FIELDS, function(f)
    is.null(meta[[f]]) || all(is.na(meta[[f]])), logical(1)))
  if (all_na)
    cat("  => all WITec fields are NA: Multi-Run P1 is DEAD for this output;\n",
        "    it falls through to um/px (P2) or the particle-extent fit.\n")
} else {
  cat("No reproducibility_meta.csv found",
      if (!is.null(repro_dir)) paste0("in ", repro_dir) else
        "under output/reproducibility/", "\n")
  cat("Pass the folder explicitly as the 2nd argument.\n")
}

# --- [B] pipeline run manifests ----------------------------------------------
hr()
cat("[B] Pipeline runs (output/*/manifest.json -> config_snapshot)\n")
hr()
rep_md5 <- NULL
if (!is.null(replicate_file)) {
  if (file.exists(replicate_file)) {
    rep_md5 <- unname(tools::md5sum(replicate_file))
    cat("Replicate run-1 file:", replicate_file, "\n  md5:", rep_md5, "\n\n")
  } else {
    cat("WARNING: replicate file not found:", replicate_file,
        "- matching by basename only.\n\n")
  }
}
# Reproducibility outputs written after the auto-resolve fix record run1's
# MD5 directly in the meta - use it when no file path was given (or the given
# path no longer exists), so this script needs no arguments on a re-run.
if (is.null(rep_md5) && !is.null(meta) && "run1_md5" %in% names(meta) &&
    !is.na(meta$run1_md5[1]) && nzchar(meta$run1_md5[1])) {
  rep_md5 <- meta$run1_md5[1]
  cat("Using run1 MD5 recorded in reproducibility_meta.csv:", rep_md5,
      " (run1_file: ", meta$run1_file[1] %||% "?", ")\n\n", sep = "")
}
run_dirs <- list.dirs(file.path(REPO_ROOT, "output"), recursive = FALSE)
run_dirs <- run_dirs[basename(run_dirs) != "reproducibility"]
matched <- NULL
for (rd in run_dirs) {
  mp <- file.path(rd, "00_manifest", "manifest.json")
  if (!file.exists(mp)) mp <- file.path(rd, "manifest.json")
  if (!file.exists(mp)) next
  m <- tryCatch(jsonlite::fromJSON(mp, simplifyVector = TRUE),
                error = function(e) NULL)
  if (is.null(m)) { cat(basename(rd), ": unreadable manifest\n"); next }
  ram <- m$inputs$raman
  cs  <- m$config_snapshot
  cat(basename(rd), " (", m$timestamp %||% "?", ")\n", sep = "")
  cat("  inputs.raman: ", ram$basename %||% "<none>",
      "  md5: ", ram$md5 %||% "<none>", "\n", sep = "")
  for (f in WITEC_FIELDS)
    cat(sprintf("  config_snapshot.%-26s %s\n", f, fmt(cs[[f]])))
  is_match <- FALSE
  if (!is.null(rep_md5) && identical(ram$md5, rep_md5)) {
    is_match <- TRUE
    cat("  *** MD5 MATCH: this run processed the replicate run-1 file ***\n")
  } else if (!is.null(replicate_file) &&
             identical(ram$basename, basename(replicate_file))) {
    is_match <- TRUE
    cat("  *** basename match (md5 not compared/differs) ***\n")
  }
  if (is_match && is.null(matched)) matched <- list(dir = rd, cs = cs)
  cat("\n")
}
if (length(run_dirs) == 0) cat("No run directories under output/.\n")

# --- [C] current R/00_config.R ------------------------------------------------
hr()
cat("[C] Current R/00_config.R (what a NEW run would record)\n")
hr()
cfg_now <- tryCatch({
  source(file.path(REPO_ROOT, "R", "utils.R"), local = TRUE)
  source(file.path(REPO_ROOT, "R", "00_config.R"), local = TRUE)
  make_config()
}, error = function(e) { cat("  (could not load config:", conditionMessage(e), ")\n"); NULL })
if (!is.null(cfg_now))
  for (f in c(WITEC_FIELDS, "raman_um_per_px"))
    cat(sprintf("  %-26s %s\n", f, fmt(cfg_now[[f]])))

# --- verdict -------------------------------------------------------------------
hr()
cat("VERDICT (plan v3 step 4: meta [A] vs matched run's snapshot [B])\n")
hr()
if (is.null(meta)) {
  cat("Cannot compare: no reproducibility_meta.csv. Re-run tools/reproducibility.R.\n")
} else if (is.null(matched)) {
  cat("No pipeline run matched the replicate file",
      if (is.null(replicate_file)) "(none given - pass it as the 1st argument)"
      else "(no manifest has matching inputs.raman md5/basename)", "\n")
  cat("If no pipeline run ever processed this exact file, the single-tab\n",
      "comparison was apples-to-oranges (plan v3, last bullet of the lead\n",
      "hypothesis) - the two tabs are showing different physical scans.\n")
} else {
  any_mismatch <- FALSE
  for (f in WITEC_FIELDS) {
    a <- suppressWarnings(as.numeric(meta[[f]][1]))
    b <- suppressWarnings(as.numeric(matched$cs[[f]]))
    b <- if (length(b) == 0) NA_real_ else b[1]
    same <- (is.na(a) && is.na(b)) ||
            (!is.na(a) && !is.na(b) && isTRUE(all.equal(a, b, tolerance = 1e-6)))
    if (!same) any_mismatch <- TRUE
    cat(sprintf("  %-26s meta=%-14s manifest=%-14s %s\n", f,
                fmt(a), fmt(b), if (same) "OK" else "<<< MISMATCH"))
  }
  if (any_mismatch) {
    src <- if ("raman_calibration_source" %in% names(meta)) meta$raman_calibration_source[1] else NA
    if (!is.na(src) && startsWith(src, "manifest:")) {
      cat("\n=> Unexpected: meta was auto-resolved (", src, ") yet still differs\n",
          "  from the matched run's CURRENT manifest - re-run tools/reproducibility.R\n",
          "  to pick up the latest calibration, then re-check.\n", sep = "")
    } else {
      cat("\n=> MISMATCH CONFIRMED: this reproducibility_meta.csv predates the\n",
          "  auto-resolve fix (or fell back to manual/none). Re-run\n",
          "  tools/reproducibility.R - it now reads the calibration directly\n",
          "  from the matching pipeline run's manifest instead of requiring\n",
          "  hand-entry (plan v3 step 5).\n")
    }
  } else {
    cat("\n=> Values match exactly. The stale-values hypothesis is RULED OUT;\n",
        "  proceed to plan v2 §4/§5 (instrument both paths, then extract the\n",
        "  shared placement function).\n")
  }
}
