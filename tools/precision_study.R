# =============================================================================
# precision_study.R — repeatability and reproducibility across instruments
# =============================================================================
# The same filter measured n times on each of up to four instruments (e.g. 3
# replicates x 4 instruments): particle counts per group and size class, and
# the ISO 5725-2 statistics (R/precision_study.R). No particle is matched
# between measurements, so remounted / rotated filters and field samples are
# fine. Only the instruments' particle tables are needed -- no images.
#
# Usage:
#   Rscript tools/precision_study.R design.csv [output_dir]
#     design.csv: columns instrument, replicate, file
#       instrument: ftir (PerkinElmer), ftir_bruker, raman, ldir (or the names)
#       file: the export as for the main pipeline; for a segmented PerkinElmer
#             scan any one segment table (its segments are read together)
#   Rscript tools/precision_study.R
#     (interactive) file pickers per instrument; Cancel moves to the next one
#
# Output (output/precision/<timestamp>/ unless output_dir is given):
#   precision_design.csv       the measurements, with file MD5s
#   precision_counts.csv       count per measurement x group x size class
#   precision_iso5725.csv      s_r, s_L, s_R (+ 95% CI), r, R, CVs per stratum
#   precision_instruments.csv  per instrument: mean, SD, Mandel h / k, dispersion
#   precision_critical.csv     Mandel h / k indicator levels (5%, 1%)
#   precision_report.pdf       tables and figures
# The viewer's Precision tab reads the same folder.
# =============================================================================

.script_path <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) normalizePath(f) else NA_character_
}, error = function(e) NA_character_)
REPO_ROOT <- if (!is.na(.script_path)) dirname(dirname(.script_path)) else getwd()

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
for (m in c("utils.R", "particle_contract.R", "00_config.R", "01_ingest.R", "01c_ingest_ldir.R",
            "01_readers.R", "08b_material_map.R", "08_agreement.R", "align_helpers.R",
            "ftir_segments.R", "precision_study.R"))
  source(file.path(REPO_ROOT, "R", m))

#' Pick the replicate files per instrument with file dialogs
collect_precision_design_interactive <- function() {
  if (!interactive())
    stop("Interactive mode needs an interactive R session; otherwise pass a design CSV:\n",
         "  Rscript tools/precision_study.R design.csv", call. = FALSE)
  rows <- list()
  for (role in names(PRECISION_ROLES)) {
    message("\n=== ", PRECISION_ROLES[[role]], ": select the replicate files one at a time ",
            "(Cancel when done, or at once to skip this instrument) ===")
    k <- 0
    repeat {
      f <- tryCatch(file.choose(), error = function(e) NULL)
      if (is.null(f)) break
      k <- k + 1
      message("  replicate ", k, ": ", basename(f))
      rows[[length(rows) + 1]] <- data.frame(instrument = role, replicate = k, file = f,
                                             stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) stop("No files selected.", call. = FALSE)
  do.call(rbind, rows)
}

args <- commandArgs(trailingOnly = TRUE)
design <- if (length(args) >= 1) read_precision_design(args[1]) else
  read_precision_design(collect_precision_design_interactive())
reps <- table(factor(design$instrument, levels = names(PRECISION_ROLES)))
reps <- reps[reps > 0]
if (any(reps < 2)) stop("Every instrument needs at least 2 replicates: ",
                        paste(names(reps)[reps < 2], collapse = ", "), call. = FALSE)
if (length(reps) < 2) log_message("One instrument only: repeatability is computed, ",
                                  "between-instrument terms need at least two", level = "WARN")

cfg <- make_config()
res <- run_precision_study(design, cfg)

out <- if (length(args) >= 2) args[2] else
  file.path(REPO_ROOT, "output", "precision", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
des <- res$design
des$md5 <- vapply(des$file, function(f) file_md5(f) %||% NA_character_, "")
w <- function(d, n) utils::write.csv(d, file.path(out, n), row.names = FALSE)
w(des, "precision_design.csv")
w(res$counts, "precision_counts.csv")
w(res$iso, "precision_iso5725.csv")
w(res$instruments, "precision_instruments.csv")
w(as.data.frame(res$critical), "precision_critical.csv")

# Report: the viewer's page renderers (shiny_app/global.R), in a private env
rep_env <- new.env(parent = globalenv())
old_wd <- setwd(REPO_ROOT)
invisible(tryCatch({
  source(file.path("shiny_app", "global.R"), local = rep_env)
  pages <- rep_env$precision_report_pages(rep_env$read_precision_results(out))
  rep_env$write_report_pdf(pages, file.path(out, "precision_report.pdf"))
}, error = function(e) log_message("Precision report not written: ", conditionMessage(e),
                                   level = "WARN"),
finally = setwd(old_wd)))

# Console overview: all particles and plastics, all sizes
log_message(strrep("=", 60))
log_message("Precision study: ", paste(sprintf("%s x%d", PRECISION_ROLES[names(reps)], as.integer(reps)),
                                       collapse = ", "))
iso <- res$iso[res$iso$size_class == "all sizes", ]
for (g in intersect(c("All particles", "Plastics"), iso$group)) {
  r <- iso[iso$group == g, ]
  log_message(sprintf("  %-14s mean %.1f | CV_r %.0f%% (s_r %.1f) | CV_L %.0f%% | CV_R %.0f%% (s_R %.1f) | r %.1f, R %.1f",
                      g, r$mean, r$cv_r_pct, r$s_r, r$cv_L_pct, r$cv_R_pct, r$s_R, r$r_limit, r$R_limit))
}
flags <- res$instruments[res$instruments$size_class == "all sizes" &
                         (nzchar(res$instruments$h_flag) | nzchar(res$instruments$k_flag)), ]
for (i in seq_len(nrow(flags)))
  log_message(sprintf("  Mandel: %s, %s: h %.2f %s, k %.2f %s", PRECISION_ROLES[[flags$instrument[i]]],
                      flags$group[i], flags$h[i], flags$h_flag[i], flags$k[i], flags$k_flag[i]),
              level = "WARN")
log_message("  Results: ", out)
