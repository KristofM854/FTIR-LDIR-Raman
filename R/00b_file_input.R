# =============================================================================
# 00b_file_input.R — Multi-instrument file input framework
# =============================================================================
#
# Two input modes:
#   1. Hardcoded test mode  — paths provided programmatically
#   2. Interactive mode     — file picker dialog, auto-detect instrument type
#
# In interactive mode the user selects files one at a time (via file.choose())
# until they press Cancel. Each file is auto-classified by instrument type
# based on filename patterns. The user can override if detection is wrong.
# =============================================================================

# Instrument detection patterns (case-insensitive)
.INSTRUMENT_PATTERNS <- list(
  FTIR  = "FTIR|Spotlight|infrared|IR\\b|FT-IR",
  Raman = "Raman",
  LDIR  = "LDIR|8700"
)


#' Detect instrument type from a filename
#'
#' Checks the filename against known patterns for each instrument.
#' Returns the first match, or "UNKNOWN" if nothing matches.
#'
#' @param filepath Path (or just filename) to check
#' @return Character: "FTIR", "Raman", "LDIR", or "UNKNOWN"
detect_instrument <- function(filepath) {
  fname <- basename(filepath)
  for (inst in names(.INSTRUMENT_PATTERNS)) {
    if (grepl(.INSTRUMENT_PATTERNS[[inst]], fname, ignore.case = TRUE)) {
      return(inst)
    }
  }
  "UNKNOWN"
}


#' Detect file type (tabular vs image)
#'
#' @param filepath Path to file
#' @return Character: "tabular" or "image"
detect_file_type <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "bmp")) {
    return("image")
  }
  "tabular"
}


#' Collect input files interactively
#'
#' Opens a file picker dialog repeatedly. The user selects files one at a
#' time until they press Cancel. Each file is auto-classified by instrument
#' type based on filename patterns.
#'
#' @return Named list of file records:
#'   list(list(path="...", instrument="FTIR", type="tabular"), ...)
collect_files_interactive <- function() {
  if (!interactive()) {
    stop("Interactive file selection requires an interactive R session.\n",
         "In batch mode, set file paths directly in main.R.")
  }

  files <- list()
  message("=== Multi-Instrument File Input ===")
  message("Select data files one at a time. Press Cancel when done.")
  message("Supported formats: .csv, .xlsx, .xls (tabular), .png (image)")
  message("")

  repeat {
    message("Select file #", length(files) + 1, " (or Cancel to finish)...")
    path <- tryCatch(file.choose(), error = function(e) NULL)

    if (is.null(path)) {
      message("File selection complete.")
      break
    }

    instrument <- detect_instrument(path)
    ftype      <- detect_file_type(path)

    message("  File:       ", basename(path))
    message("  Detected:   ", instrument, " (", ftype, ")")

    if (instrument == "UNKNOWN") {
      message("  Could not auto-detect instrument type.")
      message("  Please specify: 1=FTIR, 2=Raman, 3=LDIR")
      choice <- readline("  Enter choice (1/2/3): ")
      instrument <- switch(choice,
                           "1" = "FTIR",
                           "2" = "Raman",
                           "3" = "LDIR",
                           "UNKNOWN")
      message("  Assigned: ", instrument)
    }

    files[[length(files) + 1]] <- list(
      path       = path,
      instrument = instrument,
      type       = ftype,
      filename   = basename(path)
    )
  }

  if (length(files) == 0) {
    stop("No files selected. Exiting.")
  }

  # Summary
  message("\n=== Selected files ===")
  for (i in seq_along(files)) {
    f <- files[[i]]
    message("  ", i, ". [", f$instrument, "] (", f$type, ") ", f$filename)
  }

  files
}


#' Build a file manifest from hardcoded paths
#'
#' Creates the same structure as collect_files_interactive() but from
#' explicit paths. Used for test/batch mode.
#'
#' @param ... Named arguments where names are instrument types and values
#'   are file paths. E.g., FTIR = "path/to/ftir.csv", Raman = "path/to/raman.csv"
#' @return List of file records (same format as collect_files_interactive)
build_file_manifest <- function(...) {
  args <- list(...)
  files <- list()

  for (inst in names(args)) {
    path <- args[[inst]]
    if (is.null(path) || is.na(path)) next
    if (!file.exists(path)) {
      warning("File not found for ", inst, ": ", path)
      next
    }
    files[[length(files) + 1]] <- list(
      path       = path,
      instrument = toupper(inst),
      type       = detect_file_type(path),
      filename   = basename(path)
    )
  }

  files
}


#' Group a file manifest by instrument type
#'
#' @param manifest List of file records from collect_files_interactive()
#'   or build_file_manifest()
#' @return Named list: FTIR = list(tabular=..., image=...), Raman = ..., LDIR = ...
group_files_by_instrument <- function(manifest) {
  result <- list(
    FTIR  = list(tabular = NULL, image = NULL),
    Raman = list(tabular = NULL, image = NULL),
    LDIR  = list(tabular = NULL, image = NULL)
  )

  for (f in manifest) {
    inst <- f$instrument
    if (!inst %in% names(result)) {
      log_message("  Unknown instrument type: ", inst, " — skipping", level = "WARN")
      next
    }
    result[[inst]][[f$type]] <- f$path
  }

  result
}
