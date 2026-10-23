# =============================================================================
# stitch_ftir_segments.R -- join a segmented PerkinElmer export (temporary)
# =============================================================================
# The pipeline stitches segmented exports by itself when the selected FTIR
# table is one segment (R/ftir_segments.R). This script does the same on its
# own, e.g. to look at the stitched image first:
#
#   Rscript tools/stitch_ftir_segments.R <folder> [output folder]
#
# <folder> holds the segment tables (<name>M1.csv ...) and the six segment
# images (<name>_M1.bmp ...). Writes <name>_stitched.csv, .png and .txt (the
# image placement, read by the pipeline when it sits next to the image).
# Remove with R/ftir_segments.R after the Spotlight software update.
# =============================================================================

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 1) stop("Usage: Rscript tools/stitch_ftir_segments.R <folder> [output folder]")
  if (!file.exists(file.path("R", "ftir_segments.R")))
    stop("Run from the repository root.")
  source(file.path("R", "ftir_segments.R"))
  f <- list.files(a[1], full.names = TRUE)
  tables <- f[grepl(FTIR_SEGMENT_TABLE, basename(f), ignore.case = TRUE)]
  images <- f[grepl(FTIR_SEGMENT_IMAGE, basename(f), ignore.case = TRUE)]
  st <- stitch_ftir_segments(tables, images, if (length(a) >= 2) a[2] else a[1])
  cat("Wrote", st$csv, "\nWrote", st$png, "\nWrote", st$placement,
      "(image placement; keep it next to the image)\n")
}
