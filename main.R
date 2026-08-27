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

required_packages <- c("readxl", "ggplot2", "RANN", "ggrepel")

missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace,
                                          logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

library(ggplot2)

# Source all modules (relative to project root)
source("R/read_bmp.R")
source("R/measure_raman_placement.R")
source("R/tps_refine.R")
source("R/utils.R")
source("R/00_config.R")
source("R/00b_file_input.R")
source("R/01_ingest.R")
source("R/01b_ingest_image.R")
source("R/utils_python.R")
source("R/01c_ingest_ldir.R")
source("R/02_prefilter.R")
source("R/03_normalize.R")
source("R/03b_landmark_align.R")
source("R/03c_procrustes_align.R")
source("R/align_helpers.R")
source("R/04_ransac.R")
source("R/05_transform.R")
source("R/06_icp_refine.R")
source("R/07_match.R")
source("R/08_agreement.R")
source("R/08b_material_map.R")
source("R/09_diagnostics.R")
source("R/10_alignment_diagnostics.R")
source("R/10_export.R")

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

  .pick_file <- function(caption, required = FALSE,
                         filter = "All files|*.*") {
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
          message("  This file is required — please select it.")
          next
        }
        return(NULL)
      }
      .last_dir <<- dirname(path)
      return(path)
    }
  }

  DATA_FILTER  <- "Data files|*.csv;*.xlsx;*.xls|All files|*.*"
  IMAGE_FILTER <- "Image files|*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp|All files|*.*"

  message("=== File Input: one dialog per instrument slot ===")
  message("Press Cancel on any optional slot to skip it.\n")

  ftir_file         <- .pick_file("FTIR (PerkinElmer) — data file (.csv/.xlsx)",
                                   required = FALSE, filter = DATA_FILTER)
  ftir_image        <- .pick_file("FTIR (PerkinElmer) — microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  ftir_bruker_file  <- .pick_file("FTIR Bruker (Lumos) — data file (.csv/.xlsx) (optional)",
                                   required = FALSE, filter = DATA_FILTER)
  ftir_bruker_image <- .pick_file("FTIR Bruker (Lumos) — microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  raman_file        <- .pick_file("Raman — data file (.csv/.xlsx) [REQUIRED]",
                                   required = TRUE, filter = DATA_FILTER)
  raman_image       <- .pick_file("Raman — microscope image (optional)",
                                   required = FALSE, filter = IMAGE_FILTER)

  ldir_file         <- .pick_file("LDIR — data file (.csv/.xlsx) (optional)",
                                   required = FALSE, filter = DATA_FILTER)
  ldir_image        <- .pick_file("LDIR — companion image (optional)",
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

config <- make_config(
  ftir_path  = ftir_file,
  raman_path = raman_file,
  output_dir = "output"
)

# Store optional device paths in config
config$ldir_path         <- if (exists("ldir_file"))         ldir_file         else NULL
config$ftir_bruker_path  <- if (exists("ftir_bruker_file"))  ftir_bruker_file  else NULL
config$ftir_image        <- if (exists("ftir_image"))        ftir_image        else NULL
config$raman_image       <- if (exists("raman_image"))       raman_image       else NULL
config$ldir_image        <- if (exists("ldir_image"))        ldir_image        else NULL
config$ftir_bruker_image <- if (exists("ftir_bruker_image")) ftir_bruker_image else NULL

# Create a timestamped run subfolder (output/YYYY-MM-DD_1, _2, ...)
config$output_dir <- make_run_dir(config$output_dir)

# --- Collect input file paths for provenance ---
.run_input_paths <- list(
  ftir       = if (!is.null(config$ftir_path)  && nzchar(config$ftir_path))  config$ftir_path  else NULL,
  raman      = if (!is.null(config$raman_path) && nzchar(config$raman_path)) config$raman_path else NULL,
  ldir       = if (!is.null(config$ldir_path)  && nzchar(config$ldir_path))  config$ldir_path  else NULL,
  ldir_image  = if (!is.null(config$ldir_image)  && nzchar(config$ldir_image))  config$ldir_image  else NULL,
  ftir_image  = if (!is.null(config$ftir_image)  && nzchar(config$ftir_image))  config$ftir_image  else NULL,
  raman_image = if (!is.null(config$raman_image) && nzchar(config$raman_image)) config$raman_image else NULL
)
.run_input_paths <- .run_input_paths[!vapply(.run_input_paths, is.null, logical(1))]

# --- Canonicalize instrument images and create previews in inputs/ ---
.ftir_image_info <- NULL
.raman_image_info <- NULL
.ldir_image_info <- NULL
inputs_dir <- file.path(config$output_dir, "inputs")

if (!is.null(config$ftir_image) && nzchar(config$ftir_image) && file.exists(config$ftir_image)) {
  .ftir_image_info <- canonicalize_instrument_image(
    src_path = config$ftir_image,
    inputs_dir = inputs_dir,
    instrument = "ftir",
    max_preview_px = 2000L
  )
}
if (!is.null(config$raman_image) && nzchar(config$raman_image) && file.exists(config$raman_image)) {
  .raman_image_info <- canonicalize_instrument_image(
    src_path = config$raman_image,
    inputs_dir = inputs_dir,
    instrument = "raman",
    max_preview_px = 2000L
  )
}
if (!is.null(config$ldir_image) && nzchar(config$ldir_image) && file.exists(config$ldir_image)) {
  .ldir_image_info <- canonicalize_instrument_image(
    src_path = config$ldir_image,
    inputs_dir = inputs_dir,
    instrument = "ldir",
    max_preview_px = 2000L
  )
  if (!is.null(.ldir_image_info$canonical_path) && file.exists(.ldir_image_info$canonical_path)) {
    log_message("Using canonical PNG for LDIR processing: ", .ldir_image_info$canonical_path)
    config$ldir_image_canonical <- .ldir_image_info$canonical_path
    config$ldir_image_preview   <- .ldir_image_info$preview_path
  }
}

# --- Write run manifest (provenance) as early as possible ---
tryCatch({
  write_manifest(
    run_dir         = config$output_dir,
    run_id          = basename(config$output_dir),
    config          = config,
    input_paths     = .run_input_paths,
    ldir_image_info = .ldir_image_info,
    image_infos     = list(
      ftir_image = .ftir_image_info,
      raman_image = .raman_image_info,
      ldir_image = .ldir_image_info
    ),
    stage           = "started"
  )
}, error = function(e) {
  log_message("  Could not write manifest: ", e$message, level = "WARN")
})

# --- Debug mode setup (Step 0: bulletproof) ---
# Set config$debug <- TRUE before sourcing to enable debug artifacts.
# Writes to an absolute path so nothing can silently swallow errors.
if (isTRUE(config$debug)) {
  run_id <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  config$run_id <- run_id

  # Build absolute path -- avoids any working-directory ambiguity
  debug_dir_abs <- normalizePath(
    file.path(config$output_dir, "debug"),
    winslash = "/", mustWork = FALSE
  )
  dir.create(debug_dir_abs, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(debug_dir_abs)) {
    stop("DEBUG: failed to create debug directory: ", debug_dir_abs)
  }

  # Heartbeat file -- if this is absent the whole debug run failed
  heartbeat <- file.path(debug_dir_abs, "DEBUG_ALIVE.txt")
  writeLines(c(paste0("run_id: ", run_id),
               paste0("created: ", Sys.time()),
               paste0("debug_dir: ", debug_dir_abs)),
             heartbeat)
  if (!file.exists(heartbeat)) {
    stop("DEBUG: heartbeat write failed — check permissions for ", debug_dir_abs)
  }

  config$debug_dir <- debug_dir_abs
  message("DEBUG DIR: ", debug_dir_abs)
  log_message("Debug mode ON — artifacts in: ", debug_dir_abs)
}

# Override any defaults as needed:
# config$raman_hqi_threshold     <- 75
# config$match_dist_threshold_um <- 25
# config$ransac_allow_mirror     <- TRUE

# ---------------------------------------------------------------------------
# 2. Data ingestion
# ---------------------------------------------------------------------------

log_message(strrep("=", 60))
log_message("Multi-Instrument Particle Matching Pipeline")
log_message(strrep("=", 60))

# Create stage-based output directories for audit exports
.out_dirs <- get_output_dirs(config$output_dir)

# --- FTIR ---
ftir_raw  <- ingest_ftir(config$ftir_path, sheet = config$ftir_sheet)
write.csv(ftir_raw, file.path(.out_dirs$ingested, "ftir_perkin_ingested.csv"), row.names = FALSE)

# --- Raman ---
raman_raw <- ingest_raman(config$raman_path, sheet = config$raman_sheet)
write.csv(raman_raw, file.path(.out_dirs$ingested, "raman_ingested.csv"), row.names = FALSE)

# --- Auto-calibrate Raman image placement ----------------------------------
# Measure where the exported Raman micrograph actually sits under the
# particles (the export may be a cropped/zoomed view whose footprint differs
# from the WITec panel Width/Height).  On a confident, non-mirrored fit, patch
# the run manifest's config_snapshot so the viewer places the image exactly --
# no per-run tools/diagnose_raman_placement.R --apply needed.  The config
# values stay as the scale-search seed; a weak fit leaves them untouched.
# NOTE: this deliberately does NOT require config$raman_image_width_um /
# _height_um. Those are per-dataset values the operator copies out of WITec,
# and they are only the scale-search SEED -- the search spans 0.3-2.4x it and
# derives its own seed from the particle bbox when they are absent. Gating on
# them meant the auto-calibration skipped exactly the runs that needed it:
# a new dataset with no WITec values entered fell straight through to
# particle-bbox placement, which is visibly offset on sparse scans.
if (!is.null(.raman_image_info)) {
  .raman_cal_img <- if (!is.null(.raman_image_info$canonical_path) &&
                        file.exists(.raman_image_info$canonical_path))
    .raman_image_info$canonical_path else config$raman_image
  .raman_fit <- tryCatch(
    measure_raman_image_placement(
      .raman_cal_img, raman_raw$x_um, raman_raw$y_um,
      config$raman_image_width_um, config$raman_image_height_um),
    error = function(e) { log_message("  Raman placement auto-measure error: ",
                                      conditionMessage(e), level = "WARN"); NULL })
  if (!is.null(.raman_fit) && !isTRUE(.raman_fit$mirrored)) {
    log_message(sprintf(paste0(
      "Raman image auto-calibration: extent %.0f x %.0f um at center (%.0f, %.0f), ",
      "scale %.2f x WITec, %.0f%% of particles on image ",
      "(random baseline %.0f%%, lift %.2f)."),
      .raman_fit$width_um, .raman_fit$height_um,
      .raman_fit$center_x_um, .raman_fit$center_y_um,
      .raman_fit$scale, 100 * .raman_fit$frac_bright,
      100 * (.raman_fit$baseline %||% NA_real_), .raman_fit$lift %||% NA_real_))
    config$raman_image_width_um    <- .raman_fit$width_um
    config$raman_image_height_um   <- .raman_fit$height_um
    config$raman_image_center_x_um <- .raman_fit$center_x_um
    config$raman_image_center_y_um <- .raman_fit$center_y_um
    update_manifest_config_snapshot(config$output_dir, list(
      raman_image_width_um    = .raman_fit$width_um,
      raman_image_height_um   = .raman_fit$height_um,
      raman_image_center_x_um = .raman_fit$center_x_um,
      raman_image_center_y_um = .raman_fit$center_y_um))
  } else {
    # Worth a WARN, not an INFO: the configured raman_image_* values are
    # PER-DATASET. If they were left over from another scan the viewer either
    # rejects them (containment check) and falls back to particle-bbox
    # placement, or -- worse -- accepts them and draws the image in the wrong
    # place. Either way the user needs to know the measurement did not land.
    log_message("Raman image auto-calibration: no confident non-mirrored fit; ",
                "keeping configured raman_image_* values. If the Raman tab ",
                "shows the image offset from the particles, re-measure with ",
                "tools/diagnose_raman_placement.R --apply.", level = "WARN")
  }
}

# --- LDIR (optional) ---
ldir_raw <- NULL
has_ldir <- !is.null(config$ldir_path) && nzchar(config$ldir_path)
if (has_ldir) {
  ldir_raw <- ingest_ldir(config$ldir_path)
  # Relabel low-HQI identifications as "unknown" (kept, not dropped) to mirror
  # the LDIR software; flows into the viewer and the agreement analysis.
  ldir_raw <- relabel_ldir_low_hqi(ldir_raw, config)
  # Auto-load coordinate-swap corrections from the CSV sidecar next to the LDIR
  # Excel (<stem>_coord_swaps.csv), unless swaps were set explicitly in config.
  if (is.null(config$ldir_coord_swaps))
    config$ldir_coord_swaps <- load_ldir_coord_swaps(config$ldir_path, config)
  write.csv(ldir_raw, file.path(.out_dirs$ingested, "ldir_ingested.csv"), row.names = FALSE)
  log_message("LDIR data loaded: ", nrow(ldir_raw), " particles (no coordinates)")
} else {
  log_message("LDIR: not provided — skipping LDIR analysis")
}

# --- FTIR Bruker (optional) ---
ftir_bruker_raw <- NULL
has_ftir_bruker <- !is.null(config$ftir_bruker_path) && nzchar(config$ftir_bruker_path)
if (has_ftir_bruker) {
  ftir_bruker_raw <- ingest_ftir_bruker(config$ftir_bruker_path)
  write.csv(ftir_bruker_raw, file.path(.out_dirs$ingested, "ftir_bruker_ingested.csv"), row.names = FALSE)
  log_message("FTIR (Bruker) data loaded: ", nrow(ftir_bruker_raw), " particles")
} else {
  log_message("FTIR (Bruker): not provided — skipping")
}

# --- FTIR image scan bounds (for background display only) ---
# FTIR particle um coordinates come directly from the Excel data file.
# The FTIR image is ONLY used as a background in the single FTIR viewer.
# No particle extraction is performed on the FTIR image.
ftir_scan_bounds <- NULL
if (!is.null(config$ftir_image) && nzchar(config$ftir_image)) {
  log_message(strrep("-", 50))
  log_message("FTIR image: computing scan bounds for viewer background")

  # Estimate scan bounds from image dimensions and 25 um grid step.
  # The PerkinElmer Spotlight renders ~6 image pixels per 25 um grid cell.
  ftir_img_raw <- read_image_any(config$ftir_image, verbose = TRUE)
  if (is.null(ftir_img_raw)) {
    log_message("  WARNING: could not read FTIR image — scan bounds unavailable",
                level = "WARN")
  } else {
    ftir_grid_nx <- round((ncol(ftir_img_raw) + 1) / 6)
    ftir_grid_ny <- round((nrow(ftir_img_raw) + 1) / 6)
    ftir_scan_bounds <- list(
      x_min = 0,
      x_max = ftir_grid_nx * 25,
      y_min = 0,
      y_max = ftir_grid_ny * 25
    )
    rm(ftir_img_raw)
    log_message("  Scan bounds: [0, ", ftir_scan_bounds$x_max, "] x [0, ",
                ftir_scan_bounds$y_max, "] µm")
  }
}

# ---------------------------------------------------------------------------
# 3. Pre-filtering
#
# Strategy: use only PLASTIC particles for spatial alignment (they are the
# ones that genuinely overlap between instruments). Then apply quality
# filters for the final matching & agreement steps.
# ---------------------------------------------------------------------------

# Minimal filtering (just remove invalid coords)
ftir_clean  <- prefilter_ftir(ftir_raw, min_quality = 0, min_size_um = 0)
write.csv(ftir_clean, file.path(.out_dirs$prefiltered, "ftir_perkin_prefiltered.csv"), row.names = FALSE)

raman_clean <- prefilter_raman(raman_raw, min_hqi = 0, min_size_um = 0)
write.csv(raman_clean, file.path(.out_dirs$prefiltered, "raman_prefiltered.csv"), row.names = FALSE)

# Extract PLASTIC particles from FTIR for alignment anchoring.
# select_material_anchors() falls back to the full cloud when too few particles
# carry an anchor material, so a sample dominated by an off-list polymer (or a
# field sample with no dominant polymer at all) still aligns.
.anchor_min <- config$align_min_anchor_count %||% 4
ftir_for_align <- select_material_anchors(
  ftir_clean, config$align_ftir_materials, .anchor_min, "FTIR")

# Filter Raman alignment targets by material (if configured)
raman_for_align <- select_material_anchors(
  raman_clean, config$align_raman_materials, .anchor_min, "Raman")

# Remove Raman particles below FTIR detection limit.
# Each successive filter is honoured only while it leaves enough anchors: three
# filters that are individually reasonable can intersect to nothing, and an
# empty anchor set produces a NaN centroid rather than an error.
min_size <- config$align_raman_min_size_um
if (!is.null(min_size) && min_size > 0 && any(!is.na(raman_for_align$feret_max_um))) {
  size_mask <- is.na(raman_for_align$feret_max_um) | raman_for_align$feret_max_um >= min_size
  if (sum(size_mask) >= .anchor_min) {
    n_before <- nrow(raman_for_align)
    raman_for_align <- raman_for_align[size_mask, ]
    log_message("Raman alignment: removed ", n_before - nrow(raman_for_align),
                " particles < ", min_size, " um")
  } else {
    log_message("Raman alignment: size filter (>= ", min_size, " um) would ",
                "leave only ", sum(size_mask), " anchors — skipping it",
                level = "WARN")
  }
}
# Apply HQI threshold to material alignment anchors (only reliably
# identified particles should drive material-based alignment)
if (!is.null(config$raman_hqi_threshold) && config$raman_hqi_threshold > 0 &&
    any(!is.na(raman_for_align$quality))) {
  hqi_mask <- !is.na(raman_for_align$quality) &
              raman_for_align$quality >= config$raman_hqi_threshold
  if (sum(hqi_mask) >= .anchor_min) {
    n_before_hqi <- nrow(raman_for_align)
    raman_for_align <- raman_for_align[hqi_mask, ]
    log_message("Raman alignment: HQI filter (>= ", config$raman_hqi_threshold,
                "): kept ", nrow(raman_for_align), " of ", n_before_hqi)
  } else {
    log_message("Raman alignment: HQI filter (>= ", config$raman_hqi_threshold,
                ") would leave only ", sum(hqi_mask), " anchors — skipping it",
                level = "WARN")
  }
}
log_message("Raman alignment target particles: ", nrow(raman_for_align),
            " (material + >= ", min_size, " um + HQI >= ",
            config$raman_hqi_threshold, ")")

# Quality-filtered sets for matching (applied after alignment)
# FTIR: apply quality filter as usual
ftir_for_match <- prefilter_ftir(
  ftir_raw,
  min_quality = config$ftir_quality_threshold,
  min_size_um = config$min_particle_size_um
)
# Raman: use ALL particles for spatial matching (no HQI filter here).
# HQI filtering is applied only inside analyze_agreement() so that
# spatial matching has maximum coverage while agreement scoring
# only considers particles with reliable Raman identifications.
raman_for_match <- prefilter_raman(
  raman_raw,
  min_hqi     = 0,
  min_size_um = config$min_particle_size_um
)
log_message("Raman for matching: ", nrow(raman_for_match),
            " particles (all HQI, spatial only; HQI filter applied in agreement)")

# ---------------------------------------------------------------------------
# 4. Coordinate normalization
#
# Centroids come from the FULL cleaned clouds, not the material anchor subsets.
# A material-subset centroid centres FTIR on *FTIR's* PET particles and Raman on
# *Raman's* PET particles; those coincide only if both instruments detected the
# same particles. With different detection limits -- the normal case, and the
# rule on field samples -- the two centroids point at different things and the
# aligner starts from a biased translation. The full-cloud centroid is only a
# common origin: the aligners recover the real translation themselves.
# ---------------------------------------------------------------------------

norm_result <- normalize_coordinates(
  ftir_clean, raman_clean,
  normalize_scale = config$normalize_scale
)

ftir_clean  <- norm_result$ftir
raman_clean <- norm_result$raman

# Anchor subsets carry the SAME centroid/scale as the clouds they came from
ftir_norm_align  <- apply_normalization(ftir_for_align,
                                        norm_result$ftir_centroid,
                                        norm_result$ftir_scale)
raman_norm_align <- apply_normalization(raman_for_align,
                                        norm_result$raman_centroid,
                                        norm_result$raman_scale)

# Build spatial transform set: all Raman >= 20 um (visible to FTIR).
# This set is used for ALL spatial transform steps (landmarks, ICP)
# regardless of material or HQI -- only size matters for geometry.
min_size_spatial <- config$align_raman_min_size_um
raman_for_transform <- raman_clean
if (!is.null(min_size_spatial) && min_size_spatial > 0 &&
    any(!is.na(raman_clean$feret_max_um))) {
  raman_for_transform <- raman_clean[
    is.na(raman_clean$feret_max_um) | raman_clean$feret_max_um >= min_size_spatial, ]
}
log_message("Raman for spatial transform (>= ", min_size_spatial, " um): ",
            nrow(raman_for_transform), " particles")
log_message("ICP refinement sets: FTIR ", nrow(ftir_clean),
            ", Raman ", nrow(raman_for_transform))

# ---------------------------------------------------------------------------
# 5. Tiered alignment (FTIR ↔ Raman)
#
# Tier 1 -- Landmark alignment: use large particles & fibers to quickly
#   determine the spatial transform. If confident, skip Tier 2.
# Tier 2 -- Global registration (rotation x scale sweep, translation recovered by
#   voting) plus the legacy coarse RANSAC; whichever pairs more particles wins.
#   Only runs if Tier 1 was not confident enough.
# ICP refinement always runs to polish the transform.
# ---------------------------------------------------------------------------

# --- Tier 1: Landmark alignment ---
# Use size-filtered Raman (>= 20 um) -- landmarks are selected by size inside
landmark_result <- landmark_align(ftir_clean, raman_for_transform, config)

use_landmark_transform <- landmark_result$confident && config$landmark_skip_full_ransac

if (use_landmark_transform) {
  log_message("Using landmark transform (Tier 1) — skipping full RANSAC")
  alignment_transform <- landmark_result$transform
  alignment_method    <- "landmark"
} else {
  # --- Tier 2: Full alignment on the anchor sets ---
  log_message(strrep("-", 50))
  log_message("Tier 2: Full alignment (anchor sets)")

  ransac_result <- ransac_align(ftir_norm_align, raman_norm_align, config)

  log_message("RANSAC transform: scale = ", round(ransac_result$params$scale, 4),
              ", rotation = ", round(ransac_result$params$rotation_deg, 2), " deg",
              ", reflected = ", ransac_result$params$reflected,
              ", inliers = ", ransac_result$n_inliers)

  alignment_transform <- ransac_result$transform
  alignment_method    <- "ransac"

  # Global registration: the coarse RANSAC anchors translation on single
  # nearest-neighbour guesses, which is fragile when the clouds overlap only
  # partially -- exactly the field-sample case. Global registration sweeps
  # rotation x scale and recovers translation by voting over all pairwise
  # offsets, scoring one-to-one so a degenerate collapse cannot win. Already
  # the default on the LDIR path; run both and keep whichever pairs more.
  # Disable with config$ftir_use_global_register = FALSE.
  if (!isFALSE(config$ftir_use_global_register)) {
    ftir_global <- tryCatch(
      global_register_align(ftir_norm_align, raman_norm_align, config,
                            allow_mirror = config$ransac_allow_mirror),
      error = function(e) {
        log_message("Global registration failed: ", e$message, level = "WARN")
        NULL })
    if (!is.null(ftir_global)) {
      if (ftir_global$n_inliers > ransac_result$n_inliers) {
        log_message("Tier 2: using global registration (",
                    ftir_global$n_inliers, " inliers vs RANSAC ",
                    ransac_result$n_inliers, ")")
        ransac_result       <- ftir_global
        alignment_transform <- ftir_global$transform
        alignment_method    <- "global_register"
      } else {
        log_message("Tier 2: keeping RANSAC (", ransac_result$n_inliers,
                    " inliers vs global ", ftir_global$n_inliers, ")")
      }
    }
  }
}

# --- ICP refinement (always runs to polish the transform) ---
icp_result <- icp_refine(ftir_clean, raman_for_transform, alignment_transform, config)

log_message("ICP refined transform: scale = ", round(icp_result$params$scale, 4),
            ", rotation = ", round(icp_result$params$rotation_deg, 2), " deg",
            ", converged = ", icp_result$converged)

# ---------------------------------------------------------------------------
# 7. Apply transform to particles for matching
#
# Normalize using the SAME centroids from step 4 (full-cloud based),
# then apply the ICP-refined transform.
# ---------------------------------------------------------------------------

# Normalize filtered FTIR using alignment centroids
ftir_for_match <- apply_normalization(ftir_for_match,
                                      norm_result$ftir_centroid,
                                      norm_result$ftir_scale)

# Normalize filtered Raman using alignment centroids
raman_for_match <- apply_normalization(raman_for_match,
                                       norm_result$raman_centroid,
                                       norm_result$raman_scale)

# Apply ICP-refined transform to filtered FTIR
ftir_aligned <- apply_ftir_transform(ftir_for_match, icp_result$transform)

# Also align the full FTIR set for diagnostics
ftir_aligned_all <- apply_ftir_transform(ftir_clean, icp_result$transform)

# ---------------------------------------------------------------------------
# 8. Particle matching (FTIR ↔ Raman, spatial)
# ---------------------------------------------------------------------------

match_result <- match_particles(ftir_aligned, raman_for_match, config)

# ---------------------------------------------------------------------------
# 9. Agreement analysis (FTIR ↔ Raman)
# ---------------------------------------------------------------------------

agreement <- analyze_agreement(match_result, config)

# ---------------------------------------------------------------------------
# 10. TPS assessment (systematic distortion check)
# ---------------------------------------------------------------------------

tps_assessment <- assess_tps_need(match_result)
log_message("TPS assessment: ", tps_assessment$message)
if (isTRUE(tps_assessment$recommend_tps)) {
  log_message("  RECOMMENDATION: Consider thin-plate spline warp for improved alignment",
              level = "WARN")
}

# ---------------------------------------------------------------------------
# 11. Composite matching (1:many) for unmatched FTIR
# ---------------------------------------------------------------------------

composites <- data.frame()
if (nrow(match_result$unmatched_ftir) > 0 && nrow(match_result$unmatched_raman) > 0) {
  # Ensure unmatched FTIR particles have aligned coordinates
  unmatched_ftir_src <- match_result$unmatched_ftir
  if (!"x_aligned" %in% names(unmatched_ftir_src)) {
    unmatched_ftir_src <- apply_normalization(unmatched_ftir_src,
                                              norm_result$ftir_centroid,
                                              norm_result$ftir_scale)
    tf <- apply_transform_points(unmatched_ftir_src$x_norm,
                                 unmatched_ftir_src$y_norm,
                                 icp_result$transform)
    unmatched_ftir_src$x_aligned <- tf$x_transformed
    unmatched_ftir_src$y_aligned <- tf$y_transformed
  }

  composites <- find_composite_matches(
    unmatched_ftir_src, match_result$unmatched_raman, config
  )

  if (nrow(composites) > 0) {
    log_message("Composite matches found: ", nrow(composites),
                " FTIR particles matched to multiple Raman fragments")
  }
}

# ---------------------------------------------------------------------------
# 5b. FTIR Bruker spatial pipeline (alignment → matching → agreement)
# ---------------------------------------------------------------------------

bruker_match_result     <- NULL
bruker_agreement        <- NULL
bruker_icp_result       <- NULL
bruker_norm_result      <- NULL
bruker_aligned_all      <- NULL
bruker_alignment_method <- NULL

if (has_ftir_bruker && !is.null(ftir_bruker_raw) && nrow(ftir_bruker_raw) > 0) {
  log_message(strrep("-", 50))
  log_message("FTIR (Bruker) Spatial Pipeline")

  # Pre-filtering
  ftir_bruker_clean <- prefilter_ftir(ftir_bruker_raw, min_quality = 0, min_size_um = 0)
  write.csv(ftir_bruker_clean,
            file.path(.out_dirs$prefiltered, "ftir_bruker_prefiltered.csv"),
            row.names = FALSE)

  ftir_bruker_for_align <- select_material_anchors(
    ftir_bruker_clean, config$align_ftir_materials, .anchor_min, "FTIR (Bruker)")

  ftir_bruker_for_match <- prefilter_ftir(
    ftir_bruker_raw,
    min_quality = config$ftir_quality_threshold,
    min_size_um = config$min_particle_size_um
  )

  if (nrow(ftir_bruker_for_align) < 3) {
    log_message("FTIR (Bruker): too few anchors (", nrow(ftir_bruker_for_align),
                ") for alignment — skipping", level = "WARN")
  } else {

    # Coordinate normalization. Centroids from the FULL clouds (see step 4);
    # the anchor subsets inherit the same centroid so they stay in one frame.
    bruker_norm_result     <- normalize_coordinates(
      ftir_bruker_clean, raman_clean,
      normalize_scale = config$normalize_scale
    )
    ftir_bruker_clean      <- bruker_norm_result$ftir
    ftir_bruker_norm_align <- apply_normalization(ftir_bruker_for_align,
                                                  bruker_norm_result$ftir_centroid,
                                                  bruker_norm_result$ftir_scale)

    # Apply centroid to the for-match set as well
    ftir_bruker_for_match <- apply_normalization(ftir_bruker_for_match,
                                                 bruker_norm_result$ftir_centroid,
                                                 bruker_norm_result$ftir_scale)

    # Tiered alignment: Bruker → Raman
    log_message("Tiered alignment: FTIR (Bruker) ↔ Raman")

    # Tier 1: Landmark
    bruker_landmark_result <- landmark_align(ftir_bruker_clean, raman_for_transform, config)
    use_bruker_landmark    <- bruker_landmark_result$confident && config$landmark_skip_full_ransac

    if (use_bruker_landmark) {
      log_message("  Using Bruker landmark transform (Tier 1) — skipping full RANSAC")
      bruker_alignment_transform <- bruker_landmark_result$transform
      bruker_alignment_method    <- "landmark"
    } else {
      # Tier 2: Full alignment (same two-aligner race as the PerkinElmer path)
      log_message(strrep("-", 50))
      log_message("  Tier 2: Full alignment (Bruker)")
      bruker_ransac_result <- ransac_align(ftir_bruker_norm_align, raman_norm_align, config)
      log_message("  Bruker RANSAC: scale=", round(bruker_ransac_result$params$scale, 4),
                  ", rotation=", round(bruker_ransac_result$params$rotation_deg, 2), " deg",
                  ", reflected=", bruker_ransac_result$params$reflected,
                  ", inliers=", bruker_ransac_result$n_inliers)
      bruker_alignment_transform <- bruker_ransac_result$transform
      bruker_alignment_method    <- "ransac"

      if (!isFALSE(config$ftir_use_global_register)) {
        bruker_global <- tryCatch(
          global_register_align(ftir_bruker_norm_align, raman_norm_align, config,
                                allow_mirror = config$ransac_allow_mirror),
          error = function(e) {
            log_message("  Bruker global registration failed: ", e$message,
                        level = "WARN")
            NULL })
        if (!is.null(bruker_global) &&
            bruker_global$n_inliers > bruker_ransac_result$n_inliers) {
          log_message("  Bruker Tier 2: using global registration (",
                      bruker_global$n_inliers, " inliers vs RANSAC ",
                      bruker_ransac_result$n_inliers, ")")
          bruker_ransac_result       <- bruker_global
          bruker_alignment_transform <- bruker_global$transform
          bruker_alignment_method    <- "global_register"
        }
      }
    }

    # ICP refinement
    bruker_icp_result <- icp_refine(
      ftir_bruker_clean, raman_for_transform, bruker_alignment_transform, config
    )
    log_message("  Bruker ICP: scale=", round(bruker_icp_result$params$scale, 4),
                ", rotation=", round(bruker_icp_result$params$rotation_deg, 2), " deg",
                ", converged=", bruker_icp_result$converged)

    # Apply transform
    bruker_aligned     <- apply_ftir_transform(ftir_bruker_for_match, bruker_icp_result$transform)
    bruker_aligned_all <- apply_ftir_transform(ftir_bruker_clean,     bruker_icp_result$transform)

    # Particle matching: Bruker ↔ Raman
    bruker_match_result <- match_particles(bruker_aligned, raman_for_match, config)
    bms <- bruker_match_result$match_stats
    log_message("FTIR (Bruker) ↔ Raman: ",
                bms$n_matched, " matched, ",
                bms$n_unmatched_ftir, " unmatched Bruker, ",
                bms$n_unmatched_raman, " unmatched Raman")

    # Agreement analysis
    bruker_agreement <- analyze_agreement(bruker_match_result, config)
  }
}

# ---------------------------------------------------------------------------
# 12. LDIR spatial pipeline (image → coordinates → alignment → matching)
# ---------------------------------------------------------------------------

ldir_results <- NULL
if (has_ldir && !is.null(ldir_raw)) {
  log_message(strrep("-", 50))
  log_message("LDIR Spatial Pipeline")

  # 12a. Extract coordinates from LDIR image BEFORE filtering (join on raw data)
  ldir_with_coords <- ldir_raw
  has_ldir_coords  <- FALSE

  if (!is.null(config$ldir_image) && nzchar(config$ldir_image)) {
    log_message("Extracting LDIR coordinates from companion image")

    # Prefer the canonical PNG (format-guaranteed) produced by canonicalize_ldir_image().
    # Fall back to the original path if canonicalization was skipped.
    ldir_img_for_extraction <- if (!is.null(config$ldir_image_canonical) &&
                                    file.exists(config$ldir_image_canonical)) {
      log_message("  Using canonical PNG: ", config$ldir_image_canonical)
      config$ldir_image_canonical
    } else {
      log_message("  Using original image (no canonical available): ", config$ldir_image)
      config$ldir_image
    }

    # Compute scan bounds from configured scan diameter (circular filter)
    ldir_scan_diam <- config$ldir_scan_diameter_um
    if (is.null(ldir_scan_diam)) ldir_scan_diam <- 13000
    ldir_scan_bounds <- list(
      x_min = 0, x_max = ldir_scan_diam,
      y_min = 0, y_max = ldir_scan_diam
    )

    # Prefer the LDIR software's analyzed particle-overlay image when present.
    # Its coloured blobs are the machine's own segmentation and carry
    # correctly-scaled sizes (the optical image over-sizes large particles ~3x),
    # so feeding them to join_ldir_coords fixes size-based matching for large
    # particles.  Falls back to optical-image circle-calibrated extraction.
    #
    # Search next to the ORIGINAL selected image (config$ldir_image), where the
    # user's "<name>_analyzed.<ext>" export lives -- not ldir_img_for_extraction,
    # which may point at the canonicalized copy in the run's inputs/ folder.
    processed_img <- find_ldir_processed_image(config$ldir_image, config)
    if (!is.null(processed_img)) {
      log_message("Using LDIR processed image for coordinate extraction: ",
                  basename(processed_img))
      ldir_proc_result     <- extract_ldir_processed_image_coords(
        processed_img, scan_bounds = ldir_scan_bounds,
        expected_count = nrow(ldir_raw),
        expected_total_area_um2 = sum(ldir_raw$area_um2, na.rm = TRUE),
        config = config
      )
      ldir_image_particles <- ldir_proc_result$particles
      .ldir_circle_info    <- ldir_proc_result$circle_info
    } else {
      # Circle-calibrated extraction -- expected_count uses RAW count for best matching
      ldir_extract_result  <- extract_ldir_image_coords(
        ldir_img_for_extraction,
        scan_bounds    = ldir_scan_bounds,
        expected_count = nrow(ldir_raw),
        config         = config
      )
      ldir_image_particles <- ldir_extract_result$particles
      .ldir_circle_info    <- ldir_extract_result$circle_info
    }

    # Persist circle calibration to manifest
    tryCatch(
      update_manifest_ldir_circle(config$output_dir, .ldir_circle_info),
      error = function(e)
        log_message("  Could not update manifest ldir_circle: ", e$message, level = "WARN")
    )

    # Save raw image-extracted coordinates (before join) for Shiny viewer
    ldir_image_extracted <- ldir_image_particles

    # Join image coordinates with RAW Excel data via size-based Hungarian matching
    ldir_with_coords <- join_ldir_coords(ldir_raw, ldir_image_particles, config = config)
    # Apply any manual coordinate-join corrections (config$ldir_coord_swaps)
    ldir_with_coords <- apply_ldir_coord_swaps(ldir_with_coords, config)
    write.csv(ldir_with_coords, file.path(.out_dirs$joined, "ldir_joined_raw.csv"), row.names = FALSE)

    # Validate join quality via scan-order correlation
    scan_order <- validate_ldir_scan_order(ldir_with_coords)
    log_message("  Scan order validation: ", scan_order$message)

    # Debug traces
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_raw, trace_ids, "after_ingest", config$debug_dir)
      dump_particle(ldir_with_coords, trace_ids, "after_coords_join", config$debug_dir)
      trace_particle_snapshot(ldir_raw, "after_ingest", config)
      trace_particle_snapshot(ldir_with_coords, "after_coords_join", config)
    }

    n_with_coords <- sum(!is.na(ldir_with_coords$x_um))
    has_ldir_coords <- n_with_coords >= 10

    # Guard: skip alignment if the circle was not reliably detected
    if (has_ldir_coords && !isTRUE(.ldir_circle_info$detected)) {
      log_message("  LDIR circle not reliably detected — skipping LDIR alignment",
                  level = "WARN")
      has_ldir_coords <- FALSE
    }

    log_message("  LDIR particles with coordinates: ", n_with_coords,
                " of ", nrow(ldir_with_coords))
  }

  # 12b. Pre-filter LDIR particles AFTER coordinate join
  ldir_clean <- prefilter_ldir(
    ldir_with_coords,
    min_quality = config$ldir_quality_threshold,
    min_size_um = config$min_particle_size_um
  )
  write.csv(ldir_clean, file.path(.out_dirs$prefiltered, "ldir_prefiltered.csv"), row.names = FALSE)

  # Update ldir_with_coords to the filtered version for downstream use
  ldir_with_coords <- ldir_clean

  # 12c–j. Spatial alignment & matching (only if coordinates available)
  ldir_raman_match     <- NULL
  ldir_ftir_match      <- NULL
  ldir_raman_agreement <- NULL
  ldir_icp             <- NULL
  ldir_aligned         <- NULL
  triplets             <- data.frame()

  if (has_ldir_coords) {
    log_message("LDIR spatial matching enabled")

    # 12c. Normalize LDIR coordinates using explicit normalize_coords_ldir()
    #
    # This returns an auditable norm_params object (centroid, scale, y_flip)
    # and writes ldir_norm_params.json to debug_dir if debug=TRUE.
    # The Y-flip toggle resolves the reflection ambiguity between LDIR and Raman.
    ldir_norm_result <- normalize_coords_ldir(
      ldir_with_coords,
      flip_y       = isTRUE(config$ldir_flip_y_for_alignment),
      scale_coords = isTRUE(config$normalize_scale),
      rotate_deg   = config$ldir_rotate_deg_for_alignment %||% 0,
      debug_dir    = if (isTRUE(config$debug)) config$debug_dir else NULL
    )
    ldir_with_coords <- ldir_norm_result$df
    ldir_norm_params <- ldir_norm_result$norm_params

    # Step 5: dump traced particles after normalization
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_with_coords, trace_ids, "after_normalization", config$debug_dir)
      trace_particle_snapshot(ldir_with_coords, "after_normalization", config)
    }

    # Step 3: ldir_valid is now extracted AFTER x_norm/y_norm exist
    ldir_valid <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                    !is.na(ldir_with_coords$y_um), ]

    # Step 1: dump traced particles after coords join (before normalization)
    if (isTRUE(config$debug)) {
      trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
      dump_particle(ldir_with_coords, trace_ids, "after_coords_join", config$debug_dir)
    }

    # 12d. Tiered LDIR→Raman alignment
    #
    # Tier 0 (Step 3): Explicit Procrustes -- if config$ldir_landmark_map is set,
    #   use named LDIR↔Raman correspondences to fit via SVD. This guarantees
    #   A3 (and other named landmarks) have minimal residuals by construction.
    #   When ldir_procrustes_lock=TRUE (default), this is the FINAL transform.
    #   ICP still runs but only for diagnostics (its output is not used).
    #
    # Tier 1: Size-based landmark RANSAC (large particles & fibers).
    # Tier 2: Material-based RANSAC: PET/PP/PC anchors (fallback).

    ldir_procrustes  <- NULL
    ldir_anchor_pairs <- NULL   # passed to ICP for pinned-weight refinement
    use_procrustes_final <- FALSE

    if (!is.null(config$ldir_landmark_map) && length(config$ldir_landmark_map) > 0) {
      log_message(strrep("-", 50))
      log_message("  LDIR Tier 0: Explicit Procrustes alignment from landmark map")

      ldir_procrustes <- fit_similarity_from_landmarks(
        src_df       = ldir_valid,
        tgt_df       = raman_for_transform,
        landmark_map = config$ldir_landmark_map,
        debug_dir    = if (isTRUE(config$debug)) config$debug_dir else NULL
      )

      if (ldir_procrustes$success) {
        log_message("  Procrustes SUCCESS: ", ldir_procrustes$message)

        # Build anchor_pairs (src_idx in ldir_valid, tgt_idx in raman_for_transform)
        # so ICP can pin them with high weight even when refining
        lp <- ldir_procrustes$landmark_pairs
        if (nrow(lp) > 0) {
          a_src <- match(lp$ldir_id,  ldir_valid$particle_id)
          a_tgt <- match(lp$raman_id, raman_for_transform$particle_id)
          ok_a  <- !is.na(a_src) & !is.na(a_tgt)
          if (sum(ok_a) > 0) {
            ldir_anchor_pairs <- data.frame(
              src_idx = a_src[ok_a],
              tgt_idx = a_tgt[ok_a]
            )
          }
        }

        use_procrustes_final <- isTRUE(config$ldir_procrustes_lock)
      } else {
        log_message("  Procrustes failed (", ldir_procrustes$message,
                    ") — falling through to RANSAC", level = "WARN")
      }
    }

    ldir_landmark_result <- tryCatch({
      landmark_align(ldir_valid, raman_for_transform, config, src_label = "LDIR")
    }, error = function(e) {
      log_message("  LDIR landmark alignment failed: ", e$message, level = "WARN")
      list(success = FALSE, confident = FALSE,
           n_ftir_landmarks = 0, n_raman_landmarks = 0)
    })

    use_ldir_landmark <- !use_procrustes_final &&
                         ldir_landmark_result$confident &&
                         config$landmark_skip_full_ransac

    if (use_ldir_landmark) {
      log_message("  Using LDIR landmark transform (Tier 1) — skipping material RANSAC")
      ldir_ransac <- list(
        transform = ldir_landmark_result$transform,
        params    = ldir_landmark_result$params,
        n_inliers = ldir_landmark_result$n_inliers
      )
    } else if (!use_procrustes_final) {
      # Tier 2: Material-based RANSAC (default) or descriptor RANSAC (opt-in)
      log_message(strrep("-", 50))

      if (isTRUE(config$ldir_use_descriptor_ransac)) {
        log_message("  LDIR Tier 2: Descriptor RANSAC alignment")
        ldir_for_align <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
        log_message("  LDIR alignment anchors: ", nrow(ldir_for_align), " particles")

        ldir_ransac <- tryCatch({
          descriptor_ransac_align(
            ldir_for_align, raman_norm_align, config,
            scale_min     = config$icp_min_scale        %||% 0.5,
            scale_max     = config$icp_max_scale        %||% 2.0,
            rot_limit_deg = config$icp_max_rotation_deg %||% 90
          )
        }, error = function(e) {
          log_message("  Descriptor RANSAC failed: ", e$message, level = "WARN")
          NULL
        })

        if (is.null(ldir_ransac)) {
          log_message("  Descriptor RANSAC returned NULL — falling back to material RANSAC",
                      level = "WARN")
          config$ldir_use_descriptor_ransac <- FALSE  # force fallback in logging
        }
      }

      if (!isTRUE(config$ldir_use_descriptor_ransac)) {
        log_message("  LDIR Tier 2: Material-based RANSAC alignment")
        ldir_for_align <- ldir_with_coords[!is.na(ldir_with_coords$x_um), ]
        if (!is.null(config$align_ldir_materials)) {
          ldir_mat_mask <- grepl(
            paste(config$align_ldir_materials, collapse = "|"),
            ldir_for_align$material, ignore.case = TRUE
          )
          if (sum(ldir_mat_mask) >= 4) {
            ldir_for_align <- ldir_for_align[ldir_mat_mask, ]
          } else {
            log_message("  Insufficient LDIR anchor materials (",
                        sum(ldir_mat_mask),
                        ") — using all LDIR particles for alignment")
          }
        }
        log_message("  LDIR alignment anchors: ", nrow(ldir_for_align), " particles")

        # LDIR and Raman image the same filter from the same side, so there is
        # no physical reflection between them. Allowing mirror lets the aligner
        # lock onto a spurious REFLECTED optimum that matches sparse particles
        # to the wrong neighbours (observed: reflected transform, wrong
        # scale/rotation, ~19 spurious matches). Forbid reflection unless the
        # user explicitly enables it.
        .ldir_reflect <- isTRUE(config$ldir_allow_reflection)
        .ldir_align_cfg <- config
        .ldir_align_cfg$ransac_allow_mirror <- .ldir_reflect
        if (!.ldir_reflect)
          log_message("  LDIR alignment: reflection disabled (no physical ",
                      "mirror between LDIR and Raman)")

        ldir_ransac <- tryCatch({
          ransac_align(ldir_for_align, raman_norm_align, .ldir_align_cfg)
        }, error = function(e) {
          log_message("  LDIR RANSAC failed: ", e$message, level = "WARN")
          NULL
        })

        # Global registration: robust for sparse LDIR<->Raman where the
        # coarse RANSAC locks onto a poor local optimum (few inliers despite
        # many achievable). Run it and keep whichever transform pairs more
        # particles. Disable with config$ldir_use_global_register = FALSE.
        if (!isFALSE(config$ldir_use_global_register)) {
          ldir_global <- tryCatch(
            global_register_align(ldir_for_align, raman_norm_align,
                                  .ldir_align_cfg, allow_mirror = .ldir_reflect),
            error = function(e) {
              log_message("  LDIR global registration failed: ",
                          e$message, level = "WARN"); NULL })
          if (!is.null(ldir_global)) {
            r_in <- if (!is.null(ldir_ransac)) ldir_ransac$n_inliers else -1
            if (ldir_global$n_inliers > r_in) {
              log_message("  LDIR alignment: using global registration (",
                          ldir_global$n_inliers, " inliers vs RANSAC ",
                          max(r_in, 0), ")")
              ldir_ransac <- ldir_global
            } else {
              log_message("  LDIR alignment: keeping RANSAC (", r_in,
                          " inliers vs global ", ldir_global$n_inliers, ")")
            }
          }
        }
      }
    } else {
      # Procrustes locked -- still set ldir_ransac from Procrustes for logging
      ldir_ransac <- list(
        transform = ldir_procrustes$matrix,
        params    = ldir_procrustes$params,
        n_inliers = ldir_procrustes$n_pairs
      )
    }

    if (!is.null(ldir_ransac)) {
      if (!use_procrustes_final) {
        log_message("  LDIR-Raman RANSAC: scale=",
                    round(ldir_ransac$params$scale, 4),
                    ", rot=", round(ldir_ransac$params$rotation_deg, 2), " deg",
                    ", reflected=", ldir_ransac$params$reflected,
                    ", inliers=", ldir_ransac$n_inliers)
      }

      # Initial transform for ICP (Procrustes if locked, RANSAC otherwise)
      icp_initial_transform <- if (use_procrustes_final) {
        ldir_procrustes$matrix
      } else {
        ldir_ransac$transform
      }

      # ICP refinement on all LDIR particles with coordinates
      # Anchor pairs are pinned with 100x weight (Step 4)
      ldir_for_icp <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                        !is.na(ldir_with_coords$y_um), ]

      # Resolve anchor_pairs indices into ldir_for_icp row space
      icp_anchor <- NULL
      if (!is.null(ldir_anchor_pairs) && nrow(ldir_anchor_pairs) > 0) {
        # ldir_anchor_pairs$src_idx points into ldir_valid; remap to ldir_for_icp
        valid_ids <- ldir_valid$particle_id[ldir_anchor_pairs$src_idx]
        icp_src   <- match(valid_ids, ldir_for_icp$particle_id)
        ok_remap  <- !is.na(icp_src)
        if (sum(ok_remap) > 0) {
          icp_anchor <- data.frame(
            src_idx = icp_src[ok_remap],
            tgt_idx = ldir_anchor_pairs$tgt_idx[ok_remap]
          )
          log_message("  ICP anchor pairs (pinned 100x weight): ",
                      sum(ok_remap), " landmarks")
        }
      }

      ldir_icp <- tryCatch({
        icp_refine(ldir_for_icp, raman_for_transform,
                   icp_initial_transform, config,
                   anchor_pairs = icp_anchor)
      }, error = function(e) {
        log_message("  LDIR ICP failed: ", e$message,
                    " — using initial transform", level = "WARN")
        list(
          transform    = icp_initial_transform,
          params       = extract_transform_params(icp_initial_transform),
          converged    = FALSE,
          n_iterations = 0,
          rms_history  = numeric(0)
        )
      })

      log_message("  LDIR-Raman ICP: scale=",
                  round(ldir_icp$params$scale, 4),
                  ", rot=", round(ldir_icp$params$rotation_deg, 2), " deg",
                  ", converged=", ldir_icp$converged)

      # When Procrustes is locked, discard ICP transform and keep Procrustes
      if (use_procrustes_final) {
        log_message("  Procrustes lock active — retaining Procrustes transform ",
                    "(ICP ran for diagnostics only)")
        ldir_icp$transform <- ldir_procrustes$matrix
        ldir_icp$params    <- ldir_procrustes$params
      }

      # Quality check
      ldir_final_rms <- if (length(ldir_icp$rms_history) > 0)
        tail(ldir_icp$rms_history, 1) else NA_real_
      ftir_final_rms  <- if (length(icp_result$rms_history) > 0)
        tail(icp_result$rms_history, 1) else NA_real_

      if (!is.na(ldir_final_rms) && ldir_final_rms > 100) {
        log_message("  LDIR alignment quality: POOR (ICP RMS = ",
                    round(ldir_final_rms, 1), " \u00b5m). Possible causes:",
                    level = "WARN")
        log_message("    1. LDIR scan area (ldir_scan_diameter_um=",
                    config$ldir_scan_diameter_um, " \u00b5m) may not match actual scan",
                    level = "WARN")
        log_message("    2. Too few anchor material particles for robust RANSAC",
                    level = "WARN")
        log_message("    3. Systematic coordinate flip or offset in image extraction",
                    level = "WARN")
        log_message("    Recommendation: set config$ldir_landmark_map with explicit",
                    " A3\u2194Raman correspondences, or verify ldir_scan_diameter_um",
                    level = "WARN")
      } else if (!is.na(ldir_final_rms) && !is.na(ftir_final_rms) &&
                 ldir_final_rms > 3 * ftir_final_rms) {
        log_message("  LDIR alignment quality: MARGINAL (RMS ",
                    round(ldir_final_rms, 1), " \u00b5m vs FTIR ",
                    round(ftir_final_rms, 1), " \u00b5m)", level = "WARN")
      }

      # 12f. Apply transform to all LDIR particles with coordinates
      ldir_aligned <- ldir_with_coords[!is.na(ldir_with_coords$x_um) &
                                        !is.na(ldir_with_coords$y_um), ]
      ldir_tf <- apply_transform_points(
        ldir_aligned$x_norm, ldir_aligned$y_norm, ldir_icp$transform
      )
      ldir_aligned$x_aligned <- ldir_tf$x_transformed
      ldir_aligned$y_aligned <- ldir_tf$y_transformed
      # Tag which transform path was used (for provenance in debug CSV)
      ldir_aligned$align_method <- if (use_procrustes_final) "procrustes" else
                                    if (use_ldir_landmark) "landmark_ransac" else
                                    if (isTRUE(config$ldir_use_descriptor_ransac))
                                      "descriptor_ransac" else "ransac_icp"

      # --- Transform guardrails ---
      tryCatch({
        .M      <- ldir_icp$transform
        .a      <- .M[1, 1]; .b <- .M[2, 1]
        .tf_sc  <- sqrt(.a^2 + .b^2)
        .tf_rot <- atan2(.b, .a) * 180 / pi
        .tf_tx  <- .M[1, 3]; .tf_ty <- .M[2, 3]
        .fov_um <- config$ldir_scan_diameter_um %||% 13000
        if (.tf_sc < (config$icp_min_scale %||% 0.5) ||
            .tf_sc > (config$icp_max_scale %||% 2.0))
          log_message("  ALIGNMENT GUARDRAIL: scale=", round(.tf_sc, 3),
                      " outside [", config$icp_min_scale %||% 0.5, ", ",
                      config$icp_max_scale %||% 2.0,
                      "] — transform may be unreliable", level = "WARN")
        if (abs(.tf_rot) > (config$icp_max_rotation_deg %||% 90))
          log_message("  ALIGNMENT GUARDRAIL: rotation=", round(.tf_rot, 1),
                      "\u00b0 > ", config$icp_max_rotation_deg %||% 90,
                      "\u00b0 — likely spurious rotation", level = "WARN")
        if (sqrt(.tf_tx^2 + .tf_ty^2) > .fov_um)
          log_message("  ALIGNMENT GUARDRAIL: translation=",
                      round(sqrt(.tf_tx^2 + .tf_ty^2)), " \u00b5m > FOV (",
                      .fov_um, " \u00b5m) — possible offset error", level = "WARN")
      }, error = function(e)
        log_message("  Guardrail check failed: ", e$message, level = "WARN"))

      # --- Alignment diagnostics (residuals, quiver, stats, audit) ---
      tryCatch({
        .residuals_for_overlay <- compute_alignment_residuals(ldir_aligned, raman_clean)
        write_alignment_diagnostics(
          ldir_aligned     = ldir_aligned,
          raman_clean      = raman_clean,
          ldir_norm_params = ldir_norm_params,
          norm_result      = norm_result,
          ldir_icp         = ldir_icp,
          config           = config
        )
      }, error = function(e) {
        log_message("  Alignment diagnostics failed: ", e$message, level = "WARN")
        .residuals_for_overlay <- NULL
      })

      # --- Overlay diagnostic PNG (with optional residual arrows) ---
      tryCatch({
        diag_dir   <- file.path(config$output_dir, "debug")
        if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
        ldir_diag  <- ldir_aligned[is.finite(ldir_aligned$x_aligned) &
                                    is.finite(ldir_aligned$y_aligned), ]
        raman_diag <- raman_clean[is.finite(raman_clean$x_norm) &
                                   is.finite(raman_clean$y_norm), ]

        p_overlay <- ggplot2::ggplot() +
          ggplot2::geom_point(data = raman_diag,
                              ggplot2::aes(x = x_norm, y = y_norm),
                              shape = 1, colour = "steelblue", alpha = 0.5, size = 1.5) +
          ggplot2::geom_point(data = ldir_diag,
                              ggplot2::aes(x = x_aligned, y = y_aligned),
                              shape = 2, colour = "forestgreen", alpha = 0.5, size = 1.5) +
          ggplot2::coord_fixed() + ggplot2::theme_minimal() +
          ggplot2::labs(
            title    = "LDIR-Raman overlay diagnostic",
            subtitle = paste0("LDIR: green \u25b2 (x_aligned/y_aligned)  |  ",
                              "Raman: blue \u25cb (x_norm/y_norm)  |  ",
                              "Orange arrows: residuals (up to 150)"),
            x = "\u00b5m", y = "\u00b5m"
          )

        # Add residual arrows if available (downsampled to 150)
        if (exists(".residuals_for_overlay") &&
            !is.null(.residuals_for_overlay) &&
            nrow(.residuals_for_overlay) > 0) {
          arrow_df <- if (nrow(.residuals_for_overlay) > 150)
            .residuals_for_overlay[sample.int(nrow(.residuals_for_overlay), 150), ]
          else .residuals_for_overlay
          p_overlay <- p_overlay +
            ggplot2::geom_segment(
              data = arrow_df,
              ggplot2::aes(x = ldir_x_aligned, y = ldir_y_aligned,
                           xend = raman_x_norm,  yend = raman_y_norm),
              colour = "orange", alpha = 0.45, linewidth = 0.35,
              arrow  = ggplot2::arrow(length = ggplot2::unit(0.06, "inches"),
                                      type   = "open")
            )
        }

        ggplot2::ggsave(file.path(diag_dir, "ldir_raman_overlay_diag.png"),
                        p_overlay, width = 8, height = 8, dpi = 150)
        log_message("  Overlay diagnostic saved: ", diag_dir, "/ldir_raman_overlay_diag.png")
      }, error = function(e)
        log_message("  Could not save overlay diagnostic: ", e$message, level = "WARN"))

      # --- Persist alignment method + matrix to manifest ---
      tryCatch({
        .align_method_str <- ldir_aligned$align_method[1] %||% "ransac_icp"
        update_manifest_ldir_circle(config$output_dir, .ldir_circle_info)  # refresh
        m_path <- resolve_manifest_path(config$output_dir)
        if (file.exists(m_path) && requireNamespace("jsonlite", quietly = TRUE)) {
          m_upd <- jsonlite::fromJSON(m_path, simplifyVector = FALSE)
          m_upd$ldir_raman_alignment_method <- .align_method_str
          m_upd$ldir_raman_similarity_matrix <- lapply(
            seq_len(nrow(ldir_icp$transform)),
            function(i) as.numeric(ldir_icp$transform[i, ])
          )
          writeLines(jsonlite::toJSON(m_upd, pretty = TRUE, auto_unbox = TRUE,
                                       null = "null", na = "null"), m_path)
        }
      }, error = function(e)
        log_message("  Could not persist alignment to manifest: ", e$message,
                    level = "WARN"))

      # Step 4: Assert aligned coordinates exist before any plotting
      stopifnot(all(c("x_aligned", "y_aligned") %in% colnames(ldir_aligned)))

      # Step 5: dump traced particles + full snapshot after alignment
      if (isTRUE(config$debug)) {
        trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
        dump_particle(ldir_aligned, trace_ids, "after_alignment", config$debug_dir)
        trace_particle_snapshot(ldir_aligned, "after_alignment", config)
      }

      # 12g. LDIR↔Raman matching
      ldir_raman_match <- match_particles(
        ldir_aligned, raman_for_match, config,
        src_label = "ldir", ref_label = "raman"
      )

      # 12g-bis. TPS local refinement: a single global similarity can't overlay
      # every particle when the two coordinate systems differ by a small
      # non-rigid distortion, leaving peripheral particles unmatched despite
      # clearly corresponding. Fit a regularized thin-plate spline to the
      # residual displacement at the confident matches, warp all LDIR
      # coordinates locally, and re-match -- adopting the result only if it
      # increases matches. Displacement is capped so non-overlapping debris
      # cannot be flung onto spurious partners. (No-op under force-complete
      # matching, which already matches everything.)
      if (isTRUE(config$ldir_tps_refine) &&
          !is.null(ldir_raman_match$matched) &&
          nrow(ldir_raman_match$matched) >= (config$ldir_tps_min_controls %||% 6)) {
        mm <- ldir_raman_match$matched
        warp <- tryCatch(
          tps_fit_warp(mm$ldir_x_aligned, mm$ldir_y_aligned,
                       mm$raman_x_norm,   mm$raman_y_norm,
                       lambda = config$ldir_tps_lambda %||% 0.5),
          error = function(e) { log_message("  TPS fit failed: ",
                                            e$message, level = "WARN"); NULL })
        if (!is.null(warp)) {
          w <- tps_apply(warp, ldir_aligned$x_aligned, ldir_aligned$y_aligned)
          ldir_tps <- ldir_aligned
          ldir_tps$x_aligned <- w$x; ldir_tps$y_aligned <- w$y
          match2 <- try_or(
            match_particles(ldir_tps, raman_for_match, config,
                            src_label = "ldir", ref_label = "raman"),
            default = NULL, what = "LDIR TPS re-match")
          n0 <- ldir_raman_match$match_stats$n_matched
          n1 <- if (!is.null(match2)) match2$match_stats$n_matched else -1L
          if (!is.null(match2) && n1 > n0) {
            log_message("  LDIR TPS refinement: matches ", n0, " -> ", n1,
                        " (control residual ", round(warp$ctrl_res), " um, ",
                        warp$n, " controls)")
            ldir_aligned     <- ldir_tps
            ldir_raman_match <- match2
          } else {
            log_message("  LDIR TPS refinement: no gain (", n0, " vs ",
                        max(n1, 0L), ") — keeping global alignment")
          }
        }
      }

      # 12h. LDIR↔Raman agreement analysis
      ldir_raman_agreement <- analyze_agreement(
        ldir_raman_match, config,
        instrument_a = "LDIR", instrument_b = "Raman"
      )

      # 12i. LDIR↔FTIR matching (both already in Raman coordinate frame)
      # Set up aligned FTIR as reference: matcher expects ref with x_norm/y_norm
      ftir_as_ref <- ftir_aligned
      ftir_as_ref$x_norm <- ftir_as_ref$x_aligned
      ftir_as_ref$y_norm <- ftir_as_ref$y_aligned
      ldir_ftir_match <- match_particles(
        ldir_aligned, ftir_as_ref, config,
        src_label = "ldir", ref_label = "ftir"
      )
      log_message("  LDIR-FTIR direct matching: ",
                  ldir_ftir_match$match_stats$n_matched, " pairs")

      # 12j. Three-way triplets (FTIR↔Raman ∩ LDIR↔Raman via Raman ID)
      if (nrow(match_result$matched) > 0 &&
          nrow(ldir_raman_match$matched) > 0) {
        ftir_raman_pairs <- data.frame(
          raman_particle_id   = match_result$matched$raman_particle_id,
          ftir_particle_id    = match_result$matched$ftir_particle_id,
          ftir_material       = match_result$matched$ftir_material,
          ftir_raman_distance = match_result$matched$match_distance,
          stringsAsFactors = FALSE
        )
        ldir_raman_pairs <- data.frame(
          raman_particle_id   = ldir_raman_match$matched$raman_particle_id,
          ldir_particle_id    = ldir_raman_match$matched$ldir_particle_id,
          ldir_material       = ldir_raman_match$matched$ldir_material,
          ldir_raman_distance = ldir_raman_match$matched$match_distance,
          stringsAsFactors = FALSE
        )

        triplets <- merge(ftir_raman_pairs, ldir_raman_pairs,
                          by = "raman_particle_id")

        # Add Raman material
        raman_mat_idx <- match(triplets$raman_particle_id,
                               raman_for_match$particle_id)
        triplets$raman_material <- raman_for_match$material[raman_mat_idx]

        # Compute per-triplet material agreement quality (0–3 instruments agreeing)
        # Uses polymer family classification so minor name differences still score
        if (nrow(triplets) > 0) {
          ftir_fam  <- classify_family_vec(triplets$ftir_material)
          raman_fam <- classify_family_vec(triplets$raman_material)
          ldir_fam  <- classify_family_vec(triplets$ldir_material)

          # Count how many instrument families agree with Raman
          fr_agree   <- ftir_fam == raman_fam & ftir_fam != "Unknown"
          lr_agree   <- ldir_fam == raman_fam & ldir_fam != "Unknown"
          fl_agree   <- ftir_fam == ldir_fam  & ftir_fam != "Unknown"
          n_agree    <- as.integer(fr_agree) + as.integer(lr_agree) +
                        as.integer(fl_agree)
          # n_agree ranges 0–3 (3 = all three pairwise family comparisons agree)
          triplets$n_instrument_agreement <- n_agree
          triplets$ftir_family  <- ftir_fam
          triplets$raman_family <- raman_fam
          triplets$ldir_family  <- ldir_fam
          triplets$material_consensus <- ifelse(
            n_agree == 3, "Full agreement",
            ifelse(fr_agree & lr_agree, "All three agree",
            ifelse(fr_agree, "FTIR+Raman agree",
            ifelse(lr_agree, "LDIR+Raman agree",
            ifelse(fl_agree, "FTIR+LDIR agree",
            "No agreement"))))
          )

          n_full <- sum(n_agree == 3, na.rm = TRUE)
          n_partial <- sum(n_agree >= 2 & n_agree < 3, na.rm = TRUE)
          log_message("  Three-way triplets: ", nrow(triplets),
                      " total (", n_full, " full agreement, ",
                      n_partial, " partial, ",
                      nrow(triplets) - n_full - n_partial, " no agreement)")
        } else {
          log_message("  Three-way triplets: 0 particles matched")
        }
      }

    } else {
      log_message("  LDIR spatial alignment failed — ",
                  "material comparison only", level = "WARN")
    }

  } else {
    log_message("  LDIR: no spatial coordinates — material comparison only")
  }

  # Non-spatial material distribution (always available)
  ldir_mats <- table(ldir_clean$material)
  log_message("  LDIR material distribution (top 5):")
  for (m in names(head(sort(ldir_mats, decreasing = TRUE), 5))) {
    log_message("    ", m, ": ", ldir_mats[m])
  }

  ldir_results <- list(
    ldir_clean           = ldir_clean,
    ldir_with_coords     = ldir_with_coords,
    ldir_aligned         = ldir_aligned,
    ldir_raman_match     = ldir_raman_match,
    ldir_ftir_match      = ldir_ftir_match,
    ldir_raman_agreement = ldir_raman_agreement,
    ldir_icp             = ldir_icp,
    ldir_landmark        = if (exists("ldir_landmark_result")) ldir_landmark_result else NULL,
    ldir_alignment_method = if (exists("use_ldir_landmark") && use_ldir_landmark) "landmark" else "ransac",
    triplets             = triplets,
    has_coords           = has_ldir_coords,
    ldir_material_dist   = ldir_mats,
    ldir_image_extracted = if (exists("ldir_image_extracted")) ldir_image_extracted else NULL
  )

  # --- Debug: Branch A/B comparison + residual vectors ---
  if (isTRUE(config$debug) && has_ldir_coords && !is.null(ldir_aligned) &&
      !is.null(ldir_raman_match)) {
    debug_ldir_branches(
      ldir_aligned, raman_clean, ldir_raman_match,
      ldir_icp, config
    )
  }
}

# ---------------------------------------------------------------------------
# 12b. Cross-instrument pairwise matches (non-Raman pairs)
# ---------------------------------------------------------------------------
# Every instrument is co-registered into the Raman-normalized frame via its
# aligned coordinates, so any pair can be matched directly. We reuse the
# x_norm <- x_aligned trick (as in the LDIR<->FTIR match above) to drop the
# "reference" instrument into that shared frame. match_particles picks the
# acceptance gate automatically from the instrument labels -- a pair that
# involves LDIR uses the looser LDIR gate, an FTIR<->FTIR pair the fine gate.
bruker_perkin_match <- NULL   # FTIR (Bruker) <-> FTIR (PerkinElmer)
bruker_ldir_match   <- NULL   # FTIR (Bruker) <-> LDIR

.have_bruker_aligned <- exists("bruker_aligned") && !is.null(bruker_aligned) &&
  is.data.frame(bruker_aligned) && nrow(bruker_aligned) > 0

if (.have_bruker_aligned && exists("ftir_aligned") && !is.null(ftir_aligned) &&
    nrow(ftir_aligned) > 0) {
  perkin_as_ref <- ftir_aligned
  perkin_as_ref$x_norm <- perkin_as_ref$x_aligned
  perkin_as_ref$y_norm <- perkin_as_ref$y_aligned
  bruker_perkin_match <- match_particles(
    bruker_aligned, perkin_as_ref, config,
    src_label = "ftir_bruker", ref_label = "ftir_perkin"
  )
  log_message("FTIR (Bruker) <-> FTIR (PerkinElmer): ",
              bruker_perkin_match$match_stats$n_matched, " matched")
}

if (.have_bruker_aligned && exists("ldir_aligned") && !is.null(ldir_aligned) &&
    is.data.frame(ldir_aligned) && "x_aligned" %in% names(ldir_aligned) &&
    nrow(ldir_aligned) > 0) {
  ldir_as_ref <- ldir_aligned
  ldir_as_ref$x_norm <- ldir_as_ref$x_aligned
  ldir_as_ref$y_norm <- ldir_as_ref$y_aligned
  bruker_ldir_match <- match_particles(
    bruker_aligned, ldir_as_ref, config,
    src_label = "ftir_bruker", ref_label = "ldir"
  )
  log_message("FTIR (Bruker) <-> LDIR: ",
              bruker_ldir_match$match_stats$n_matched, " matched")
}

# ---------------------------------------------------------------------------
# 13. Diagnostics (use full datasets for overlay, matched pairs for detail)
# ---------------------------------------------------------------------------

diagnostics <- generate_diagnostics(
  ftir_aligned_all, raman_clean,
  match_result, icp_result, agreement,
  config
)

# Add LDIR diagnostics if spatial matching was performed
if (!is.null(ldir_results) && !is.null(ldir_results$ldir_aligned) &&
    !is.null(ldir_results$ldir_raman_match)) {

  # Step 4: Assert LDIR overlay uses aligned coordinates
  ldir_plot_df <- ldir_results$ldir_aligned
  stopifnot(all(c("x_aligned", "y_aligned") %in% colnames(ldir_plot_df)))

  # Step 5: dump traced particles + full snapshot for final overlay
  if (isTRUE(config$debug)) {
    trace_ids <- config$debug_trace_ids %||% c("A3", "MP_11")
    dump_particle(ldir_plot_df, trace_ids, "final_overlay", config$debug_dir)
    trace_particle_snapshot(ldir_plot_df, "final_overlay", config)

    # Log the exact plotted coordinates for traced particles
    for (pid in trace_ids) {
      ri <- which(ldir_plot_df$particle_id == pid)
      if (length(ri) > 0) {
        message(sprintf("OVERLAY CHECK '%s': x_aligned=%.2f  y_aligned=%.2f",
                        pid, ldir_plot_df$x_aligned[ri[1]],
                        ldir_plot_df$y_aligned[ri[1]]))
      } else {
        message(sprintf("OVERLAY CHECK '%s': NOT IN ldir_plot_df", pid))
      }
    }
  }

  ldir_diag <- generate_ldir_diagnostics(
    ldir_plot_df, raman_clean,
    ldir_results$ldir_raman_match, ldir_results$ldir_raman_agreement,
    debug_subtitle = if (isTRUE(config$debug)) {
      # Build per-particle subtitle for A3 etc.
      parts <- character(0)
      for (pid in (config$debug_trace_ids %||% c("A3", "MP_11"))) {
        ri <- which(ldir_plot_df$particle_id == pid)
        if (length(ri) > 0) {
          parts <- c(parts, sprintf("%s: (%.1f, %.1f)",
                                    pid,
                                    ldir_plot_df$x_aligned[ri[1]],
                                    ldir_plot_df$y_aligned[ri[1]]))
        }
      }
      if (length(parts) > 0) paste(parts, collapse = "  ") else NULL
    } else NULL
  )
  diagnostics <- c(diagnostics, ldir_diag)
  log_message("  Added ", length(ldir_diag), " LDIR diagnostic plots")
}

# ---------------------------------------------------------------------------
# 14. Export
# ---------------------------------------------------------------------------

export_results(
  match_result, agreement, diagnostics,
  icp_result, norm_result, config,
  ftir_scan_bounds = ftir_scan_bounds,
  ldir_results     = ldir_results,
  input_paths      = .run_input_paths,
  ldir_image_info  = .ldir_image_info,
  image_infos      = list(
    ftir_image = .ftir_image_info,
    raman_image = .raman_image_info,
    ldir_image = .ldir_image_info
  )
)

# Use stage-based output directories for remaining exports
.export_dirs <- get_output_dirs(config$output_dir)

# Export FTIR Bruker results
if (!is.null(bruker_match_result)) {
  if (nrow(bruker_match_result$matched) > 0) {
    write.csv(bruker_match_result$matched,
              file.path(.export_dirs$matches, "matched_ftir_bruker_raman.csv"),
              row.names = FALSE)
    log_message("  Wrote matched_ftir_bruker_raman.csv (",
                nrow(bruker_match_result$matched), " pairs)")
  }
  if (nrow(bruker_match_result$unmatched_ftir) > 0) {
    write.csv(bruker_match_result$unmatched_ftir,
              file.path(.export_dirs$matches, "unmatched_ftir_bruker_vs_raman.csv"),
              row.names = FALSE)
    log_message("  Wrote unmatched_ftir_bruker_vs_raman.csv (",
                nrow(bruker_match_result$unmatched_ftir), " particles)")
  }
  if (!is.null(bruker_icp_result) && !is.null(bruker_norm_result)) {
    bp <- bruker_icp_result$params
    bruker_param_lines <- c(
      "# FTIR_bruker-to-Raman Transform Parameters",
      paste0("# Generated: ", Sys.time()),
      "",
      paste0("scale:        ", round(bp$scale, 6)),
      paste0("rotation_deg: ", round(bp$rotation_deg, 4)),
      paste0("tx:           ", round(bp$tx, 4)),
      paste0("ty:           ", round(bp$ty, 4)),
      paste0("reflected:    ", bp$reflected),
      "",
      "# Normalization parameters (applied before transform)",
      paste0("ftir_centroid_x:  ", round(bruker_norm_result$ftir_centroid[1], 4)),
      paste0("ftir_centroid_y:  ", round(bruker_norm_result$ftir_centroid[2], 4)),
      paste0("raman_centroid_x: ", round(bruker_norm_result$raman_centroid[1], 4)),
      paste0("raman_centroid_y: ", round(bruker_norm_result$raman_centroid[2], 4)),
      "",
      "# ICP refinement info",
      paste0("icp_converged:    ", bruker_icp_result$converged),
      paste0("icp_iterations:   ", bruker_icp_result$n_iterations),
      paste0("icp_final_rms:    ",
             round(tail(bruker_icp_result$rms_history, 1), 4), " um"),
      "",
      "# 3x3 Transform matrix (homogeneous, FTIR_bruker_norm -> Raman_norm)",
      paste0("matrix_row1: ", paste(round(bruker_icp_result$transform[1, ], 8), collapse = ", ")),
      paste0("matrix_row2: ", paste(round(bruker_icp_result$transform[2, ], 8), collapse = ", ")),
      paste0("matrix_row3: ", paste(round(bruker_icp_result$transform[3, ], 8), collapse = ", "))
    )
    writeLines(bruker_param_lines,
               file.path(.export_dirs$alignment, "transform_params_ftir_bruker_raman.txt"))
    log_message("  Wrote transform_params_ftir_bruker_raman.txt")
  }
  if (!is.null(bruker_agreement) && !is.null(bruker_agreement$agreement_detail) &&
      nrow(bruker_agreement$agreement_detail) > 0) {
    write.csv(bruker_agreement$agreement_detail,
              file.path(.export_dirs$agreement, "agreement_pairwise_ftir_bruker_raman.csv"),
              row.names = FALSE)
    log_message("  Wrote agreement_pairwise_ftir_bruker_raman.csv")
  }
} else if (has_ftir_bruker && !is.null(ftir_bruker_raw) && nrow(ftir_bruker_raw) > 0) {
  # Fallback for skipped alignment (e.g. too few plastic anchors)
  write.csv(ftir_bruker_raw,
            file.path(.export_dirs$matches, "unmatched_ftir_bruker.csv"),
            row.names = FALSE)
  log_message("  Wrote unmatched_ftir_bruker.csv (", nrow(ftir_bruker_raw), " particles)")
}

# Export cross-instrument matches (Bruker<->Perkin, Bruker<->LDIR). Perkin<->LDIR
# is written by export_ldir_results() as matched_ldir_ftir_perkin.csv.
if (!is.null(bruker_perkin_match) && nrow(bruker_perkin_match$matched) > 0) {
  write.csv(bruker_perkin_match$matched,
            file.path(.export_dirs$matches, "matched_ftir_bruker_ftir_perkin.csv"),
            row.names = FALSE)
  log_message("  Wrote matched_ftir_bruker_ftir_perkin.csv (",
              nrow(bruker_perkin_match$matched), " pairs)")
}
if (!is.null(bruker_ldir_match) && nrow(bruker_ldir_match$matched) > 0) {
  write.csv(bruker_ldir_match$matched,
            file.path(.export_dirs$matches, "matched_ftir_bruker_ldir.csv"),
            row.names = FALSE)
  log_message("  Wrote matched_ftir_bruker_ldir.csv (",
              nrow(bruker_ldir_match$matched), " pairs)")
}

# Export composite matches if found
if (nrow(composites) > 0) {
  write.csv(composites,
            file.path(.export_dirs$matches, "composite_matches.csv"),
            row.names = FALSE)
  log_message("  Wrote composite_matches.csv (", nrow(composites), " composites)")
}

# Export TPS assessment
if (!is.null(tps_assessment$quadrant_residuals)) {
  write.csv(tps_assessment$quadrant_residuals,
            file.path(.export_dirs$diagnostics, "tps_quadrant_residuals.csv"),
            row.names = FALSE)
  tps_lines <- c(
    "# TPS (Thin-Plate Spline) Assessment",
    paste0("# Generated: ", Sys.time()),
    "",
    paste0("systematic_distortion: ", tps_assessment$systematic_distortion),
    paste0("recommend_tps:         ", tps_assessment$recommend_tps),
    paste0("dx_range_um:           ", tps_assessment$dx_range_um),
    paste0("dy_range_um:           ", tps_assessment$dy_range_um),
    paste0("assessment:            ", tps_assessment$message)
  )
  writeLines(tps_lines, file.path(.export_dirs$diagnostics, "tps_assessment.txt"))
}

# ---------------------------------------------------------------------------
# 15. Summary
# ---------------------------------------------------------------------------

log_message(strrep("=", 60))
log_message("Pipeline complete!")
log_message(strrep("=", 60))
log_message("  Alignment method:  ", alignment_method,
            if (alignment_method == "landmark")
              paste0(" (", landmark_result$n_inliers, " landmark inliers, ",
                     round(landmark_result$mean_residual, 1), " µm mean residual)")
            else "")
log_message("  FTIR landmarks:    ", landmark_result$n_ftir_landmarks,
            " (confident: ", landmark_result$confident, ")")
log_message("  FTIR plastic anchors (alignment): ", nrow(ftir_for_align))
log_message("  Raman targets (alignment):       ", nrow(raman_for_align))
log_message("  FTIR particles (matching):       ", nrow(ftir_for_match))
log_message("  Raman particles (matching):      ", nrow(raman_for_match))
log_message("  Matched pairs:      ", match_result$match_stats$n_matched)
log_message("  Unmatched FTIR:     ", match_result$match_stats$n_unmatched_ftir)
log_message("  Unmatched Raman:    ", match_result$match_stats$n_unmatched_raman)
log_message("  FTIR match rate:    ",
            round(match_result$match_stats$match_rate_ftir * 100, 1), "%")
if (!is.na(agreement$agreement_rate)) {
  log_message("  Material agreement: ", round(agreement$agreement_rate * 100, 1), "%")
}
if (!is.null(agreement$tiered_rates) && agreement$tiered_rates$n_total > 0) {
  tr <- agreement$tiered_rates
  log_message("  Tiered agreement:  Exact ", tr$exact_pct, "%, Family+ ",
              tr$family_or_better_pct, "%")
}
if (nrow(composites) > 0) {
  log_message("  Composite matches: ", nrow(composites))
}
log_message("  TPS assessment:    ", tps_assessment$message)
if (has_ldir && !is.null(ldir_results)) {
  log_message("  --- LDIR ---")
  log_message("  LDIR particles:    ", nrow(ldir_results$ldir_clean))
  log_message("  LDIR coordinates:  ",
              if (ldir_results$has_coords) "spatial (from image)" else "none")
  if (!is.null(ldir_results$ldir_landmark)) {
    log_message("  LDIR alignment:    ", ldir_results$ldir_alignment_method,
                " (landmarks: ", ldir_results$ldir_landmark$n_ftir_landmarks,
                ", confident: ", ldir_results$ldir_landmark$confident, ")")
  }
  if (!is.null(ldir_results$ldir_raman_match)) {
    log_message("  LDIR-Raman matched: ",
                ldir_results$ldir_raman_match$match_stats$n_matched)
  }
  if (!is.null(ldir_results$ldir_ftir_match)) {
    log_message("  LDIR-FTIR matched:  ",
                ldir_results$ldir_ftir_match$match_stats$n_matched)
  }
  if (nrow(ldir_results$triplets) > 0) {
    log_message("  Three-way triplets: ", nrow(ldir_results$triplets))
  }
  if (!is.null(ldir_results$ldir_raman_agreement) &&
      !is.null(ldir_results$ldir_raman_agreement$tiered_rates) &&
      ldir_results$ldir_raman_agreement$tiered_rates$n_total > 0) {
    ltr <- ldir_results$ldir_raman_agreement$tiered_rates
    log_message("  LDIR-Raman agreement: Exact ", ltr$exact_pct,
                "%, Family+ ", ltr$family_or_better_pct, "%")
  }
}
log_message("  Results in: ", config$output_dir)

# ---------------------------------------------------------------------------
# 16. Generate PDF Report
# ---------------------------------------------------------------------------

tryCatch({
  log_message("Generating PDF report...")

  # Source report helpers from Shiny app
  source("shiny_app/global.R", local = FALSE)

  # --- Report-only quality gate --------------------------------------------
  # Everything above this point (ingest, image recognition, alignment, ICP,
  # matching, agreement) deliberately runs on the FULL particle set — narrowing
  # it there would change which pairs the registration can find. Only the
  # report is filtered, and it uses the same per-instrument defaults as the
  # Shiny viewer's quality sliders so the two agree on what a reported
  # particle is.
  #
  # The scales differ per instrument and are NOT interchangeable: FTIR
  # (PerkinElmer and Bruker) and LDIR carry quality on 0-1, while Raman
  # carries HQI on 0-100. Confirmed against the ingested CSVs.
  REPORT_QUALITY_RANGE <- list(
    "FTIR (PerkinElmer)" = c(0.70, 1),
    "FTIR (Bruker)"      = c(0.70, 1),
    "Raman"              = c(70,   100),   # HQI scale
    "LDIR"               = c(0.80, 1)
  )

  # Mirrors filter_instrument() in the viewer, including dropping NA quality.
  .report_quality_filter <- function(df, label) {
    rng <- REPORT_QUALITY_RANGE[[label]]
    if (is.null(df) || is.null(rng) || !("quality" %in% names(df))) return(df)
    q <- suppressWarnings(as.numeric(df$quality))
    df[!is.na(q) & q >= rng[1] & q <= rng[2], , drop = FALSE]
  }

  # Prepare data for report
  report_devices <- list()

  .add_report_device <- function(devices, label, df) {
    if (is.null(df) || nrow(df) == 0) return(devices)
    n_before <- nrow(df)
    df <- .report_quality_filter(df, label)
    rng <- REPORT_QUALITY_RANGE[[label]]
    log_message(sprintf(
      "  Report filter %-19s quality %s-%s: %d of %d particles kept",
      label, format(rng[1]), format(rng[2]), nrow(df), n_before))
    if (nrow(df) == 0) return(devices)
    devices[[label]] <- df
    devices
  }

  report_devices <- .add_report_device(report_devices, "FTIR (PerkinElmer)",
                                       ftir_clean)
  report_devices <- .add_report_device(report_devices, "FTIR (Bruker)",
                                       ftir_bruker_clean)
  report_devices <- .add_report_device(report_devices, "Raman", raman_clean)
  if (has_ldir)
    report_devices <- .add_report_device(report_devices, "LDIR",
                                         ldir_results$ldir_clean)

  # Generate report PDF
  report_file <- file.path(config$output_dir, "particle_report.pdf")

  # Build report pages
  report_pages <- list()

  # Title page. The quality gate is stated explicitly: these counts are a
  # filtered subset, while the alignment and matching upstream used every
  # particle, so a reader comparing the two needs to know why they differ.
  title_text <- paste0(
    "Particle Analysis Report\n",
    "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "Run: ", basename(config$output_dir), "\n",
    "\n",
    "Quality filter applied to this report:\n",
    paste(vapply(names(REPORT_QUALITY_RANGE), function(k) {
      r <- REPORT_QUALITY_RANGE[[k]]
      sprintf("  %-19s %s - %s%s", k, format(r[1]), format(r[2]),
              if (identical(k, "Raman")) "  (HQI)" else "")
    }, character(1)), collapse = "\n"), "\n",
    "\n",
    "Alignment and matching upstream used ALL particles; the filter\n",
    "affects only the figures and tables below."
  )
  report_pages[[1]] <- report_text_page("Report", title_text, "")

  # Summary tables
  if (length(report_devices) > 0) {
    plastics_tbl <- report_plastics_table(report_devices)
    if (nrow(plastics_tbl) > 0) {
      report_pages[[length(report_pages) + 1]] <-
        report_table_page(plastics_tbl, "Plastic Families",
                          paste0("Family counts per instrument, quality-filtered ",
                                 "(see the title page for the per-instrument ",
                                 "ranges)."))
    }

    size_tbl <- report_size_stats_table(report_devices)
    if (nrow(size_tbl) > 0) {
      report_pages[[length(report_pages) + 1]] <-
        report_table_page(size_tbl, "Size Statistics",
                          paste0("Feret Max statistics per instrument, ",
                                 "quality-filtered (see the title page for ",
                                 "the per-instrument ranges)."))
    }
  }

  # Write PDF
  n_pages <- write_report_pdf(report_pages, report_file)
  log_message("  Report written: ", report_file, " (", n_pages, " pages)")

  # Write HTML (same content + interactive plotly barplot)
  html_file <- file.path(config$output_dir, "particle_report.html")
  tryCatch({
    device_counts <- lapply(
      setNames(nm = names(report_devices)),
      function(lbl) {
        x <- report_devices[[lbl]]
        if (is.null(x) || nrow(x) == 0 || !"material" %in% names(x)) return(NULL)
        table(classify_family_vec(x$material))
      }
    )
    plotly_fig <- tryCatch(
      build_plotly_barplot(device_counts),
      error = function(e) { log_message("  WARNING: plotly build failed: ", e$message); NULL }
    )
    write_report_html(report_pages, html_file, plotly_fig = plotly_fig,
                      plotly_insert_after = 2L)
    log_message("  HTML report written: ", html_file)
  }, error = function(e) {
    log_message("  WARNING: HTML report generation failed: ", e$message)
  })

}, error = function(e) {
  log_message("  WARNING: Report generation failed: ", e$message)
})
