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
  # Each call opens a single-file picker with a descriptive caption.
  # Optional slots (images, LDIR, Bruker) return NULL when the user
  # presses Cancel -- the pipeline skips those instruments/images.
  # Mandatory slots (Raman) re-prompt until a file is chosen.

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

  DATA_FILTER  <- "Data files|*.csv;*.xlsx;*.xls|All files|*.*"
  IMAGE_FILTER <- "Image files|*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp|All files|*.*"

  message("=== File Input: one dialog per instrument slot ===")
  message("Press Cancel on any optional slot to skip it.")
  message("Raman is required; to abort the whole run, press Cancel twice ",
          "in a row on that dialog.\n")

  ftir_file         <- .pick_file("FTIR (PerkinElmer) \u2014 data file (.csv/.xlsx)",
                                   required = FALSE, filter = DATA_FILTER)
  ftir_image        <- .pick_file("FTIR (PerkinElmer) \u2014 microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  ftir_bruker_file  <- .pick_file("FTIR Bruker (Lumos) \u2014 data file (.csv/.xlsx) (optional)",
                                   required = FALSE, filter = DATA_FILTER)
  ftir_bruker_image <- .pick_file("FTIR Bruker (Lumos) \u2014 microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  raman_file        <- .pick_file("Raman \u2014 data file (.csv/.xlsx) [REQUIRED]",
                                   required = TRUE, filter = DATA_FILTER)
  raman_image       <- .pick_file("Raman \u2014 microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  ldir_file         <- .pick_file("LDIR \u2014 data file (.csv/.xlsx) (optional)",
                                   required = FALSE, filter = DATA_FILTER)
  ldir_image        <- .pick_file("LDIR \u2014 companion image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  message("\n=== Selected files ===")
  for (nm in c("ftir_file", "ftir_image", "ftir_bruker_file", "ftir_bruker_image",
               "raman_file", "raman_image", "ldir_file", "ldir_image")) {
    val <- get(nm)
    message(sprintf("  %-24s %s", nm, if (is.null(val)) "(skipped)" else val))
  }

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
