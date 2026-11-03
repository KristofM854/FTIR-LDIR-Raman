# =============================================================================
# main.R -- Multi-instrument particle matching pipeline
# =============================================================================
#
# Orchestrates the full pipeline for aligning and matching particles
# detected by FTIR, Raman, and LDIR microspectroscopy on the same filter.
#
# Input modes:
#   Mode "explicit"   (default): one labeled file-picker dialog per instrument slot
#   Mode "hardcoded": set ftir_file / raman_file / etc. before sourcing, then
#                     set input_mode <- "hardcoded"
#   Mode "interactive": legacy pattern-based auto-detection (disabled)
#
# Usage:
#   source("main.R")
#
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Setup: load packages and source modules
# ---------------------------------------------------------------------------

source("R/load.R")
load_pipeline()

# ---------------------------------------------------------------------------
# 1. File input and configuration
# ---------------------------------------------------------------------------

# --- Input mode selection ---
# Default is "explicit": one labeled file-picker dialog per instrument slot.
# Set input_mode <- "hardcoded" before sourcing to supply paths directly.
# The old pattern-based auto-detection block is preserved below (commented out)
# in case it is needed again.
if (!exists("input_mode")) input_mode <- "explicit"

if (input_mode == "explicit") {
  # ------ Explicit per-slot file picker dialogs ------
  # A list of instruments first, then one single-file picker per slot of the
  # chosen instruments. Image slots return NULL on Cancel (skipped); a data
  # slot re-prompts once, a second Cancel ends the run. A confirmation lists
  # the files before the run starts (No ends it).

  .is_windows <- tolower(.Platform$OS.type) == "windows"

  # Remembers the folder of the most recently selected file, across every
  # dialog in this run (data files and images alike) -- including slots the
  # user cancels, since the folder is captured before the next dialog opens.
  # Without this, choose.files() re-opens in the R working directory every
  # time, forcing a re-navigate for each of the up to 8 per-instrument
  # dialogs when the dataset lives elsewhere.
  .last_dir <- NULL

  # Abort the run cleanly. Raised as a condition rather than a bare stop() so
  # a caller that sources main.R can distinguish "the operator cancelled" from
  # "the pipeline failed"; uncaught it simply halts sourcing, which is what an
  # operator at the console wants.
  .cancel_pipeline <- function(reason = "cancelled at the file chooser") {
    message("\n=== Pipeline cancelled (", reason, ") ===")
    stop(structure(
      class = c("pipeline_cancelled", "error", "condition"),
      list(message = paste0("Pipeline cancelled by the user: ", reason),
           call = NULL)))
  }

  .pick_file <- function(caption, required = FALSE,
                         filter = "All files|*.*") {
    # A required slot used to re-prompt forever on Cancel, so the only way out
    # was to kill the R process. Cancel twice in a row to abort instead.
    cancels <- 0L
    repeat {
      path <- if (.is_windows) {
        tryCatch(
          choose.files(caption = caption,
                       filters = matrix(strsplit(filter, "\\|")[[1]],
                                        ncol = 2, byrow = TRUE),
                       multi = FALSE,
                       default = if (!is.null(.last_dir))
                         file.path(.last_dir, "*.*") else ""),
          error = function(e) character(0)
        )
      } else {
        message(caption, " (Cancel to skip):")
        tryCatch(file.choose(), error = function(e) NULL)
      }
      # Normalize to NULL when nothing was selected
      if (is.null(path) || length(path) == 0 || !nzchar(path)) {
        if (required) {
          cancels <- cancels + 1L
          # A non-interactive session cannot offer a second chance -- there is
          # nobody to press Cancel again, so re-prompting would spin forever.
          if (!interactive() || cancels >= 2L)
            .cancel_pipeline("required file not selected")
          message("  This file is required \u2014 please select it.")
          message("  (Press Cancel again to abort the pipeline.)")
          next
        }
        return(NULL)
      }
      cancels <- 0L
      .last_dir <<- dirname(path)
      return(path)
    }
  }

  # Last check before the run: the chosen files, Yes to run, No to end.
  # A Tk window as wide as the longest file name (the Windows message box
  # wraps long names at a fixed width); the message box, or the console,
  # when Tk is not available.
  .confirm_files <- function(lines) {
    tk <- tryCatch(.confirm_files_tk(lines), error = function(e) NA)
    if (!is.na(tk)) return(tk)
    if (.is_windows)
      return(identical(utils::winDialog("yesno", paste0(
        "Start the run with these files?\n\n", paste(lines, collapse = "\n"))), "YES"))
    ans <- readline("Start the run with these files? [Y/n] ")
    tolower(trimws(ans)) %in% c("", "y", "yes")
  }
  .confirm_files_tk <- function(lines) {
    if (!isTRUE(capabilities("tcltk")) || !requireNamespace("tcltk", quietly = TRUE)) return(NA)
    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "Start the run?")
    ok <- tcltk::tclVar(0)
    tcltk::tkgrid(tcltk::tklabel(tt, text = "Start the run with these files?"),
                  sticky = "w", padx = 14, pady = c(12, 6))
    tcltk::tkgrid(tcltk::tklabel(tt, text = paste(lines, collapse = "\n"), justify = "left",
                                 font = "TkFixedFont"), sticky = "w", padx = 14)
    btn <- tcltk::tkframe(tt)
    close <- function(v) function() { tcltk::tclvalue(ok) <- v; tcltk::tkdestroy(tt) }
    yes <- tcltk::tkbutton(btn, text = "Yes", width = 10, command = close(1))
    tcltk::tkgrid(yes, tcltk::tkbutton(btn, text = "No", width = 10, command = close(0)), padx = 6)
    tcltk::tkgrid(btn, sticky = "e", padx = 10, pady = 12)
    tcltk::tkbind(tt, "<Return>", close(1)); tcltk::tkbind(tt, "<Escape>", close(0))
    tcltk::tcl("wm", "attributes", tt, "-topmost", 1)
    tcltk::tkfocus(yes)
    tcltk::tkwait.window(tt)
    identical(tcltk::tclvalue(ok), "1")
  }

  DATA_FILTER  <- "Data files|*.csv;*.xlsx;*.xls|All files|*.*"
  IMAGE_FILTER <- "Image files|*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp|All files|*.*"

  # Which instruments take part: one list up front, so only their dialogs
  # open (Cancel here ends the run at once). A wrong file is caught by the
  # confirmation after the last dialog -- one click instead of cancelling
  # through every remaining slot.
  .INSTRUMENTS <- c(ftir = "FTIR (PerkinElmer)", ftir_bruker = "FTIR Bruker (Lumos)",
                    raman = "Raman", ldir = "LDIR")
  .chosen <- if (interactive())
    names(.INSTRUMENTS)[.INSTRUMENTS %in% utils::select.list(
      unname(.INSTRUMENTS), multiple = TRUE, graphics = .is_windows,
      title = "Instruments in this run (Cancel ends the run)")] else names(.INSTRUMENTS)
  if (length(.chosen) == 0) .cancel_pipeline("no instruments selected")
  if (length(.chosen) < 2) .cancel_pipeline("at least two instruments are needed")

  message("=== File Input: ", paste(.INSTRUMENTS[.chosen], collapse = ", "), " ===")
  message("Data file: Cancel twice ends the run. Image: Cancel skips it.\n")

  .slot <- function(key, what, filter, image = FALSE) {
    if (!key %in% .chosen) return(NULL)
    .pick_file(paste0(.INSTRUMENTS[[key]], " \u2014 ", what), required = !image, filter = filter)
  }
  ftir_file         <- .slot("ftir", "data file (.csv/.xlsx)", DATA_FILTER)
  ftir_image        <- .slot("ftir", "microscope image (optional)", IMAGE_FILTER, image = TRUE)
  ftir_bruker_file  <- .slot("ftir_bruker", "data file (.csv/.xlsx)", DATA_FILTER)
  ftir_bruker_image <- .slot("ftir_bruker", "microscope image (optional)", IMAGE_FILTER, image = TRUE)
  raman_file        <- .slot("raman", "data file (.csv/.xlsx)", DATA_FILTER)
  raman_image       <- .slot("raman", "microscope image, full size is fine (optional)",
                             IMAGE_FILTER, image = TRUE)
  ldir_file         <- .slot("ldir", "data file (.csv/.xlsx)", DATA_FILTER)
  ldir_image        <- .slot("ldir", "companion image (optional)", IMAGE_FILTER, image = TRUE)

  .picked <- c(ftir_file = ftir_file %||% "", ftir_image = ftir_image %||% "",
               ftir_bruker_file = ftir_bruker_file %||% "",
               ftir_bruker_image = ftir_bruker_image %||% "",
               raman_file = raman_file %||% "", raman_image = raman_image %||% "",
               ldir_file = ldir_file %||% "", ldir_image = ldir_image %||% "")
  .summary <- sprintf("  %-18s %s", names(.picked),
                      ifelse(nzchar(.picked), basename(.picked), "(skipped)"))
  message("\n=== Selected files ===")
  message(paste(sprintf("  %-24s %s", names(.picked),
                        ifelse(nzchar(.picked), .picked, "(skipped)")), collapse = "\n"))
  if (interactive() && !.confirm_files(.summary))
    .cancel_pipeline("files not confirmed")
  .n_data <- sum(!vapply(list(ftir_file, ftir_bruker_file, raman_file, ldir_file),
                         is.null, logical(1)))
  if (.n_data < 2)
    .cancel_pipeline(paste0(.n_data, " instrument data file(s) selected; at least two ",
                            "are needed"))

} else if (input_mode == "interactive") {
  # ------ [LEGACY] Pattern-based auto-detection ------
  # Commented out because instrument type cannot always be reliably inferred
  # from the filename alone (e.g. Bruker Lumos files have no fixed keyword).
  # Kept here for reference; use input_mode <- "explicit" instead.
  #
  # file_manifest <- collect_files_interactive()
  # grouped       <- group_files_by_instrument(file_manifest)
  #
  # ftir_file         <- grouped$FTIR_perkin$tabular
  # ftir_image        <- grouped$FTIR_perkin$image
  # ftir_bruker_file  <- grouped$FTIR_bruker$tabular
  # ftir_bruker_image <- grouped$FTIR_bruker$image
  # raman_file        <- grouped$Raman$tabular
  # raman_image       <- grouped$Raman$image
  # ldir_file         <- grouped$LDIR$tabular
  # ldir_image        <- grouped$LDIR$image

  stop("input_mode 'interactive' is disabled. Use input_mode <- 'explicit' instead.")

} else {
  # ------ Mode: Hardcoded paths ------
  # Set these variables before sourcing main.R, then set input_mode <- "hardcoded".
  if (!exists("ftir_file"))         ftir_file         <- NULL
  if (!exists("raman_file"))        raman_file        <- NULL
  if (!exists("ldir_file"))         ldir_file         <- NULL
  if (!exists("ftir_image"))        ftir_image        <- NULL
  if (!exists("raman_image"))       raman_image       <- NULL
  if (!exists("ldir_image"))        ldir_image        <- NULL
  if (!exists("ftir_bruker_file"))  ftir_bruker_file  <- NULL
  if (!exists("ftir_bruker_image")) ftir_bruker_image <- NULL
}

# output_dir may be preset before sourcing (hardcoded mode, examples, tests).
config <- make_config(
  ftir_path  = ftir_file,
  raman_path = raman_file,
  output_dir = if (exists("output_dir")) output_dir else "output"
)

# Store optional device paths in config
config$ldir_path         <- if (exists("ldir_file"))         ldir_file         else NULL
config$ftir_bruker_path  <- if (exists("ftir_bruker_file"))  ftir_bruker_file  else NULL
config$ftir_image        <- if (exists("ftir_image"))        ftir_image        else NULL
config$raman_image       <- if (exists("raman_image"))       raman_image       else NULL
config$ldir_image        <- if (exists("ldir_image"))        ldir_image        else NULL
config$ftir_bruker_image <- if (exists("ftir_bruker_image")) ftir_bruker_image else NULL

# Parameters to change for this run, e.g.
#   config_overrides <- list(debug = TRUE, match_dist_threshold_um = 80)
# set before sourcing main.R (or build the config yourself and call
# run_pipeline(config) directly; see R/run_pipeline.R).
if (exists("config_overrides")) config <- utils::modifyList(config, config_overrides)

# ---------------------------------------------------------------------------
# 2-16. Run (R/run_pipeline.R)
# ---------------------------------------------------------------------------
pipeline_result <- run_pipeline(config)
