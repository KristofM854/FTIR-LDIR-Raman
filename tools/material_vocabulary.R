# =============================================================================
# material_vocabulary.R -- how every identification label of a run is classified
# =============================================================================
# Lists, per instrument, each distinct label in a run's ingested tables with
# its count, the name it was translated to by the material table (if any),
# the exact polymer and the family. Labels that land in "Unknown" or "Other
# polymer" are the candidates for a row in inst/materials/<role>.csv.
#
#   Rscript tools/material_vocabulary.R <run_dir> [out.csv]
# =============================================================================

material_vocabulary <- function(run_dir) {
  files <- list.files(file.path(run_dir, "01_ingested"), "_ingested\\.csv$", full.names = TRUE)
  out <- lapply(files, function(f) {
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    raw <- if ("identification_raw" %in% names(d)) d$identification_raw else d$material
    key <- paste(raw, d$material, d$material_polymer, d$material_family, sep = "\r")
    n <- table(key)
    k <- strsplit(names(n), "\r", fixed = TRUE)
    data.frame(instrument = sub("_ingested\\.csv$", "", basename(f)),
               label = vapply(k, `[`, "", 1), material = vapply(k, `[`, "", 2),
               polymer = vapply(k, `[`, "", 3), family = vapply(k, `[`, "", 4),
               n = as.integer(n), stringsAsFactors = FALSE)
  })
  v <- do.call(rbind, out)
  v$needs_table_entry <- v$family %in% c("Unknown", "Other polymer") &
    !toupper(trimws(v$material)) %in% c("", "UNKNOWN", "NA")
  v[order(v$instrument, -v$n, v$label), ]
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript tools/material_vocabulary.R <run_dir> [out.csv]")
  v <- material_vocabulary(args[1])
  if (length(args) >= 2) utils::write.csv(v, args[2], row.names = FALSE)
  print(v, row.names = FALSE)
}
