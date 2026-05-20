# =============================================================================
# 00b_file_input.R — Multi-instrument file input framework
# =============================================================================
#
# Two input modes:
#   1. Hardcoded test mode  — paths provided programmatically
#   2. Interactive mode     — file picker dialog, auto-detect instrument type
#
# In interactive mode on Windows: a single multi-select dialog opens so the user
# can pick all files at once (Ctrl+click / Shift+click). On other platforms the
# old single-file loop (file.choose() until Cancel) is used as a fallback.
# Each file is auto-classified by instrument type from the filename; unknown types
# prompt the user for a manual choice.
# =============================================================================

# Instrument detection patterns (case-insensitive)
# FTIR_perkin: PerkinElmer Spotlight keywords (legacy "FTIR" also maps here)
# FTIR_bruker: Bruker OPUS / ALPHA keywords
.INSTRUMENT_PATTERNS <- list(
  FTIR_perkin = "FTIR|Spotlight|infrared|FT-IR|PerkinElmer",
  FTIR_bruker = "Bruker|OPUS|bruker|Lumos",
  Raman       = "Raman",
  LDIR        = "LDIR|8700"
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
  if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "bmp", "webp")) {
    return("image")
  }
  "tabular"
}


#' Collect input files interactively
#'
#' On Windows: opens a single multi-select file picker dialog
#' (choose.files()) so the user can select all files at once.
#' On other platforms: falls back to repeated file.choose() calls
#' (one file at a time until Cancel).
#'
#' Each file is auto-classified by instrument type based on filename
#' patterns. Files with unknown instrument type prompt for manual entry.
#'
#' @return Named list of file records:
#'   list(list(path="...", instrument="FTIR_perkin", type="tabular"), ...)
collect_files_interactive <- function() {
  if (!interactive()) {
    stop("Interactive file selection requires an interactive R session.\n",
         "In batch mode, set file paths directly in main.R.")
  }

  is_windows <- tolower(.Platform$OS.type) == "windows"

  message("=== Multi-Instrument File Input ===")
  message("Supported formats: .csv, .xlsx, .xls (tabular); .png/.jpg/.tif/.bmp/.webp (image)")
  message("")

  raw_paths <- character(0)

  if (is_windows) {
    # Windows: single dialog, Ctrl+click / Shift+click to select multiple files
    message("Select all input files at once (Ctrl+click or Shift+click for multiple)...")
    filter_str <- paste0(
      "All supported files|*.csv;*.xlsx;*.xls;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp;*.webp|",
      "Data files (csv/xlsx)|*.csv;*.xlsx;*.xls|",
      "Image files|*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp;*.webp|",
      "All files|*.*"
    )
    raw_paths <- tryCatch(
      choose.files(caption = "Select all FTIR/Raman/LDIR data + image files",
                   filters = matrix(strsplit(filter_str, "\\|")[[1]],
                                    ncol = 2, byrow = TRUE),
                   multi   = TRUE),
      error = function(e) character(0)
    )
    if (length(raw_paths) == 0) stop("No files selected. Exiting.")
    message("  Selected ", length(raw_paths), " file(s).")
  } else {
    # Non-Windows: one file at a time until Cancel
    message("Select data files one at a time. Press Cancel when done.")
    repeat {
      message("Select file #", length(raw_paths) + 1, " (or Cancel to finish)...")
      p <- tryCatch(file.choose(), error = function(e) NULL)
      if (is.null(p)) { message("File selection complete."); break }
      raw_paths <- c(raw_paths, p)
    }
    if (length(raw_paths) == 0) stop("No files selected. Exiting.")
  }

  files <- list()
  for (path in raw_paths) {
    instrument <- detect_instrument(path)
    ftype      <- detect_file_type(path)

    message("  File:       ", basename(path))
    message("  Detected:   ", instrument, " (", ftype, ")")

    if (instrument == "UNKNOWN") {
      message("  Could not auto-detect instrument type.")
      message("  Please specify: 1=FTIR (PerkinElmer), 2=FTIR (Bruker), 3=Raman, 4=LDIR")
      choice <- readline("  Enter choice (1/2/3/4): ")
      instrument <- switch(choice,
                           "1" = "FTIR_perkin",
                           "2" = "FTIR_bruker",
                           "3" = "Raman",
                           "4" = "LDIR",
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
    stop("No valid files after classification. Exiting.")
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
#' @return Named list: FTIR_perkin = list(tabular=..., image=...), FTIR_bruker = ...,
#'   Raman = ..., LDIR = ...
group_files_by_instrument <- function(manifest) {
  result <- list(
    FTIR_perkin = list(tabular = NULL, image = NULL),
    FTIR_bruker = list(tabular = NULL, image = NULL),
    Raman       = list(tabular = NULL, image = NULL),
    LDIR        = list(tabular = NULL, image = NULL)
  )

  for (f in manifest) {
    inst <- f$instrument
    # Backwards compatibility: bare "FTIR" maps to FTIR_perkin
    if (identical(inst, "FTIR")) inst <- "FTIR_perkin"
    if (!inst %in% names(result)) {
      log_message("  Unknown instrument type: ", inst, " — skipping", level = "WARN")
      next
    }
    result[[inst]][[f$type]] <- f$path
  }

  result
}
