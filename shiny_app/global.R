# =============================================================================
# global.R -- Load pipeline output and prepare data for Shiny
# =============================================================================

library(shiny)
library(ggplot2)
library(ggrepel)
library(png)
library(rintrojs)

# %||% is used throughout this app (image offsets, placement fallbacks, etc.)
# but is only DEFINED (guarded) in R/07_match.R / R/09_diagnostics.R, which
# this file does not source. It was silently relying on main.R having been
# run earlier in the same R session (which sources those files) to define it
# globally; a Shiny session started fresh (Rscript/rsconnect/a new R process)
# never gets it, so any code path that actually uses %||% (e.g. placing an
# uploaded background image) fails with "could not find function \"%||%\"".
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# Allow large instrument images (TIFF micrographs are commonly 10-40 MB;
# Shiny's default cap is only 5 MB).  NOTE: this limit is enforced by Shiny
# at the HTTP layer BEFORE any server code runs, so an oversized upload can
# never be rescued server-side -- it is rejected outright.  Uploads that pass
# this ceiling are immediately downsized in memory for display (see
# downsample_raster() below), so storage and rendering stay small.
options(shiny.maxRequestSize = 50 * 1024^2)   # 50 MB

# Longest edge (px) to which uploaded background images are downsized.
# 2000 px keeps enough resolution to visually align particles against the
# membrane/micrograph while keeping annotation_raster() rendering fast.
BG_IMAGE_MAX_DIM <- 2000L

# Default lower bound (um) for every instrument's Feret Max filter, in the
# viewer's sliders AND in the generated report. Particles below this are
# display/reporting noise. The slider minimum stays 0 so an operator can still
# dial it down; this only sets where it starts.
#
# The PIPELINE is deliberately unaffected: ingest, image recognition,
# alignment, ICP and matching all keep running on every particle, because
# dropping small ones there would change which pairs registration can find.
DEFAULT_MIN_SIZE_UM <- 20


# Locate the pipeline's R/ directory. global.R is sourced from two very
# different working directories: the Shiny app runs with getwd() == shiny_app/
# (so R/ is at ../R), while main.R sources this file from the project root
# (so R/ is at ./R). Hard-coding "../R" broke report generation in the
# pipeline with "cannot open the connection". Resolve it instead of assuming.
.pipeline_r_dir <- local({
  cands <- c(file.path("..", "R"), "R",
             file.path("..", "..", "R"))
  hit <- Filter(function(d) file.exists(file.path(d, "08b_material_map.R")), cands)
  if (length(hit) == 0)
    stop("global.R: cannot locate the pipeline R/ directory from ", getwd())
  hit[[1]]
})

# Source canonical material classification from pipeline
# (classify_family_vec, classify_category, classify_category_vec, etc.)
source(file.path(.pipeline_r_dir, "08b_material_map.R"), local = TRUE)

# Dependency-free BMP reader (read_bmp_raster) -- lets load_image_raster()
# handle instrument BMP exports even when magick is not installed.
source(file.path(.pipeline_r_dir, "read_bmp.R"), local = TRUE)

# ---------------------------------------------------------------------------
# List ALL available runs in the output directory (newest first).
# Returns a named character vector suitable for selectInput choices:
#   names = display labels (run_id + git short + timestamp)
#   values = run directory paths
# ---------------------------------------------------------------------------
list_all_runs <- function(output_dir = file.path("..", "output")) {
  if (!dir.exists(output_dir)) return(character(0))

  # Subdirectory runs: detect both staged (05_matches/) and legacy (flat) layouts
  runs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  has_staged <- file.exists(file.path(runs, "05_matches", "matched_ftir_perkin_raman.csv"))
  has_legacy <- file.exists(file.path(runs, "matched_particles.csv"))
  runs <- runs[has_staged | has_legacy]

  if (length(runs) == 0) return(character(0))

  # Sort newest first
  runs <- runs[order(file.mtime(runs), decreasing = TRUE)]

  # Build display labels (include manifest info if available)
  labels <- vapply(runs, function(run_dir) {
    run_id <- basename(run_dir)
    m_path <- file.path(run_dir, "00_manifest", "manifest.json")
    if (!file.exists(m_path)) m_path <- file.path(run_dir, "manifest.json")
    if (file.exists(m_path) && requireNamespace("jsonlite", quietly = TRUE)) {
      m <- tryCatch(
        jsonlite::fromJSON(m_path, simplifyVector = FALSE),
        error = function(e) NULL
      )
      if (!is.null(m)) {
        ts  <- if (!is.null(m$timestamp)) substr(m$timestamp, 1, 16) else ""
        git <- if (!is.null(m$git_commit) && !is.na(m$git_commit) &&
                   nzchar(m$git_commit))
                 paste0(" [", substr(m$git_commit, 1, 7), "]") else ""
        stg <- if (!is.null(m$stage)) paste0(" (", m$stage, ")") else ""
        return(paste0(run_id, git, stg))
      }
    }
    # Fallback: run_id + mtime
    mtime <- format(file.mtime(run_dir), "%Y-%m-%d %H:%M")
    paste0(run_id, " [", mtime, "]")
  }, character(1))

  setNames(runs, labels)
}

# ---------------------------------------------------------------------------
# Load manifest.json from a run directory (graceful fallback if absent).
# ---------------------------------------------------------------------------
load_run_manifest <- function(run_dir) {
  # Check staged location first, then legacy
  m_path <- file.path(run_dir, "00_manifest", "manifest.json")
  if (!file.exists(m_path)) m_path <- file.path(run_dir, "manifest.json")
  if (!file.exists(m_path)) {
    return(list(
      is_missing = TRUE,
      run_id     = basename(run_dir),
      timestamp  = NA_character_,
      git_commit = NA_character_,
      stage      = NA_character_,
      inputs     = list(),
      ldir_image = NULL,
      config_snapshot = list()
    ))
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list(is_missing = FALSE, run_id = basename(run_dir),
                timestamp = NA_character_, git_commit = NA_character_,
                stage = "unknown", inputs = list(),
                ldir_image = NULL, config_snapshot = list()))
  }
  tryCatch({
    m <- jsonlite::fromJSON(m_path, simplifyVector = FALSE)
    m$is_missing <- FALSE
    m
  }, error = function(e) {
    list(is_missing = TRUE, run_id = basename(run_dir),
         timestamp = NA_character_, git_commit = NA_character_,
         stage = "error", error = e$message, inputs = list(),
         ldir_image = NULL, config_snapshot = list())
  })
}


# Resolve a manifest image asset path (preview/canonical/original) for an instrument.
manifest_image_path <- function(manifest, input_name, preferred = c("preview", "canonical", "original")) {
  preferred <- match.arg(preferred)
  if (is.null(manifest$image_assets) || is.null(manifest$image_assets[[input_name]])) return(NULL)
  asset <- manifest$image_assets[[input_name]]
  node <- asset[[preferred]]
  if (!is.null(node$path) && file.exists(node$path)) return(node$path)
  # Fallback order
  for (alt in c("preview", "canonical", "original")) {
    node_alt <- asset[[alt]]
    if (!is.null(node_alt$path) && file.exists(node_alt$path)) return(node_alt$path)
  }
  NULL
}

load_png_raster <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  png::readPNG(path)
}

# points_df must have numeric columns x and y (in um)
# If your columns are named differently, rename before calling.
build_single_view_plot <- function(points_df, bg_png_path = NULL) {
  stopifnot(is.data.frame(points_df))
  if (!all(c("x", "y") %in% names(points_df))) {
    stop("points_df must contain columns x and y")
  }
  
  x_rng <- range(points_df$x, na.rm = TRUE)
  y_rng <- range(points_df$y, na.rm = TRUE)
  
  # background
  bg_layer <- NULL
  if (!is.null(bg_png_path) && file.exists(bg_png_path)) {
    bg <- load_png_raster(bg_png_path)
    if (!is.null(bg)) {
      bg_layer <- annotation_raster(bg, xmin=x_rng[1], xmax=x_rng[2], ymin=y_rng[1], ymax=y_rng[2])
    }
  }
  
  # aesthetics fallbacks
  if (!("match_status" %in% names(points_df))) points_df$match_status <- "unknown"
  if (!("feret_max_um" %in% names(points_df))) points_df$feret_max_um <- 50
  
  ggplot(points_df, aes(x = x, y = y)) +
    bg_layer +
    geom_point(aes(color = match_status, size = feret_max_um), alpha = 0.9) +
    coord_fixed(expand = FALSE) +
    scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
    theme_minimal()
}

# ---------------------------------------------------------------------------
# Extract canonical image paths from a run manifest.
# Returns a named list: list(ftir = path|NULL, raman = path|NULL, ldir = path|NULL)
# Paths are verified to exist; falls back to inputs/ folder filenames if needed.
# ---------------------------------------------------------------------------
get_run_image_paths <- function(manifest, run_dir) {
  out <- list(ftir = NULL, raman = NULL, ldir = NULL, ftir_bruker = NULL)

  # Highest priority: images the user uploaded in the viewer, persisted to the
  # run directory by persist_uploaded_image() (app.R).  An explicit upload
  # overrides pipeline-generated images and survives run switches / restarts.
  for (instr in names(out)) {
    up <- file.path(run_dir, "inputs", paste0(instr, "_image_uploaded.png"))
    if (file.exists(up)) out[[instr]] <- up
  }
  if (is.null(manifest) || isTRUE(manifest$is_missing)) return(out)

  # Then manifest$image_assets.  NB: assign only non-NULL results -- writing
  # NULL into a list DROPS the element, after which out$ftir would partial-
  # match out$ftir_bruker and cross-wire the two instruments' images.
  for (nm in c("ftir", "raman", "ldir")) {
    if (is.null(out[[nm]])) {
      mp <- manifest_image_path(manifest, paste0(nm, "_image"), preferred = "canonical")
      if (!is.null(mp)) out[[nm]] <- mp
    }
  }

  # Fallback to run_dir/inputs naming convention (relative)
  for (nm in c("ftir", "raman", "ldir")) {
    if (is.null(out[[nm]])) {
      fb <- file.path(run_dir, "inputs", paste0(nm, "_image_canonical.png"))
      if (file.exists(fb)) out[[nm]] <- fb
    }
  }

  # Last resort: the original source file the user selected at run time
  # (manifest$inputs$<instr>_image$path).  Runs processed without magick have
  # no canonical/preview PNGs in inputs/, but the source image is usually
  # still on disk on the same machine -- and load_image_raster() can read it
  # directly (PNG/JPEG natively, BMP via read_bmp_raster, TIFF via magick).
  for (nm in c("ftir", "raman", "ldir", "ftir_bruker")) {
    if (is.null(out[[nm]])) {
      src <- tryCatch(manifest$inputs[[paste0(nm, "_image")]]$path,
                      error = function(e) NULL)
      if (!is.null(src) && length(src) == 1 && is.character(src) &&
          nzchar(src) && file.exists(src))
        out[[nm]] <- src
    }
  }

  out
}

# ---------------------------------------------------------------------------
# Locate pipeline output.
# Supports two layouts:
#   1. Subdirectory runs:  output/2026-02-16_6/matched_particles.csv
#   2. Flat timestamped:   output/matched_particles_20260216_172156.csv
# Returns a list(dir, format) or NULL.
# ---------------------------------------------------------------------------
find_latest_run <- function(output_dir = file.path("..", "output")) {
  if (!dir.exists(output_dir)) return(NULL)

  # --- Try subdirectory format first (staged or legacy) ---
  runs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  has_staged <- file.exists(file.path(runs, "05_matches", "matched_ftir_perkin_raman.csv"))
  has_legacy <- file.exists(file.path(runs, "matched_particles.csv"))
  runs <- runs[has_staged | has_legacy]
  if (length(runs) > 0) {
    runs <- runs[order(file.mtime(runs), decreasing = TRUE)]
    return(list(dir = runs[1], format = "subdir"))
  }

  # --- Try flat timestamped format ---
  flat <- list.files(output_dir, pattern = "^matched_particles.*\\.csv$",
                     full.names = TRUE)
  if (length(flat) > 0) {
    flat <- flat[order(file.mtime(flat), decreasing = TRUE)]
    return(list(dir = output_dir, format = "flat",
                matched_file = flat[1]))
  }

  NULL
}

# ---------------------------------------------------------------------------
# Load pipeline CSVs from either subdirectory or flat format
# ---------------------------------------------------------------------------
load_run_data <- function(run_info) {
  data <- list(run_dir = run_info$dir)

  if (run_info$format == "subdir") {
    rd <- run_info$dir
    is_staged <- dir.exists(file.path(rd, "05_matches"))

    # Helper: resolve file from staged path, fall back to legacy flat path
    resolve <- function(staged_rel, legacy_rel) {
      if (is_staged) {
        fp <- file.path(rd, staged_rel)
        if (file.exists(fp)) return(fp)
      }
      fp2 <- file.path(rd, legacy_rel)
      if (file.exists(fp2)) return(fp2)
      NULL
    }

    # Core match files (staged names <- Part D pairwise naming)
    file_map <- list(
      matched               = c("05_matches/matched_ftir_perkin_raman.csv",
                                 "matched_particles.csv"),
      unmatched_ftir        = c("05_matches/unmatched_ftir_perkin_vs_raman.csv",
                                 "unmatched_ftir.csv"),
      unmatched_raman       = c("05_matches/unmatched_raman_vs_ftir_perkin.csv",
                                 "unmatched_raman.csv"),
      agreement             = c("06_agreement/agreement_summary.csv",
                                 "agreement_summary.csv"),
      matched_ftir_bruker   = c("05_matches/matched_ftir_bruker_raman.csv",
                                 "matched_ftir_bruker_raman.csv"),
      unmatched_ftir_bruker = c("05_matches/unmatched_ftir_bruker_vs_raman.csv",
                                 "unmatched_ftir_bruker.csv"),
      # Cross-instrument (non-Raman) pairwise matches
      matched_bruker_perkin = c("05_matches/matched_ftir_bruker_ftir_perkin.csv",
                                 "matched_ftir_bruker_ftir_perkin.csv"),
      matched_bruker_ldir   = c("05_matches/matched_ftir_bruker_ldir.csv",
                                 "matched_ftir_bruker_ldir.csv"),
      matched_perkin_ldir   = c("05_matches/matched_ldir_ftir_perkin.csv",
                                 "matched_ldir_ftir_perkin.csv")
    )
    for (nm in names(file_map)) {
      fp <- resolve(file_map[[nm]][1], file_map[[nm]][2])
      if (!is.null(fp)) data[[nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }

    # Transform params
    tp <- resolve("04_alignment/transform_params_ftir_perkin_raman.txt",
                  "transform_params.txt")
    if (!is.null(tp)) data$transform <- parse_transform_params(tp)

    # LDIR files (optional)
    ldir_map <- list(
      ldir_raman_matched   = c("05_matches/matched_ldir_raman.csv",
                                "ldir_raman_matched.csv"),
      unmatched_ldir       = c("05_matches/unmatched_ldir_vs_raman.csv",
                                "unmatched_ldir.csv"),
      triplets             = c("05_matches/triplets_3way.csv",
                                "triplets_3way.csv"),
      ldir_image_extracted = c("02_joined/ldir_image_extracted.csv",
                                "ldir_image_extracted.csv")
    )
    for (nm in names(ldir_map)) {
      fp <- resolve(ldir_map[[nm]][1], ldir_map[[nm]][2])
      if (!is.null(fp)) data[[nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }

  } else {
    # Flat: files have timestamps in the name. Find each by prefix.
    find_flat <- function(prefix, ext = "csv") {
      pat <- paste0("^", prefix, ".*\\.", ext, "$")
      hits <- list.files(run_info$dir, pattern = pat, full.names = TRUE)
      if (length(hits) == 0) return(NULL)
      hits[order(file.mtime(hits), decreasing = TRUE)][1]
    }

    for (info in list(
      list(nm = "matched",         prefix = "matched_particles"),
      list(nm = "unmatched_ftir",  prefix = "unmatched_ftir"),
      list(nm = "unmatched_raman", prefix = "unmatched_raman"),
      list(nm = "agreement",       prefix = "agreement_summary")
    )) {
      fp <- find_flat(info$prefix)
      if (!is.null(fp)) data[[info$nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }
    tp <- find_flat("transform_params", "txt")
    if (!is.null(tp)) data$transform <- parse_transform_params(tp)

    # LDIR files (flat format)
    for (info in list(
      list(nm = "ldir_raman_matched",  prefix = "ldir_raman_matched"),
      list(nm = "unmatched_ldir",      prefix = "unmatched_ldir"),
      list(nm = "triplets",            prefix = "triplets_3way"),
      list(nm = "ldir_image_extracted", prefix = "ldir_image_extracted"),
      list(nm = "matched_bruker_perkin", prefix = "matched_ftir_bruker_ftir_perkin"),
      list(nm = "matched_bruker_ldir",   prefix = "matched_ftir_bruker_ldir"),
      list(nm = "matched_perkin_ldir",   prefix = "matched_ldir_ftir_perkin")
    )) {
      fp <- find_flat(info$prefix)
      if (!is.null(fp)) data[[info$nm]] <- read.csv(fp, stringsAsFactors = FALSE)
    }
  }

  data <- annotate_ldir_gate(data, run_info$dir)
  enrich_material_family(data)
}

# ---------------------------------------------------------------------------
# LDIR<->Raman acceptance gate + genuine-match classification
# ---------------------------------------------------------------------------
# ldir_force_complete_match = TRUE makes the pipeline pair EVERY LDIR particle
# with a Raman particle regardless of distance, so "matched" in the raw CSV is
# vacuous -- every LDIR is always paired. The per-pair `match_distance` (aligned
# coordinate Euclidean distance, recomputed after any TPS refinement) is the
# real signal: a pair is a genuine match only when it falls within the LDIR
# acceptance gate. Classifying here, once, keeps the summary, the overlay plot,
# and the hover/tables perfectly consistent -- the image agrees with the table.

# Return the LDIR<->Raman acceptance gate (um) for a run. Read from the run
# manifest's config snapshot; falls back to the pipeline default (250 um) for
# runs whose manifest predates the gate being recorded.
ldir_acceptance_gate <- function(run_dir) {
  default_gate <- 250
  if (is.null(run_dir) || !nzchar(run_dir)) return(default_gate)
  man <- tryCatch(load_run_manifest(run_dir), error = function(e) NULL)
  cs  <- if (!is.null(man)) man$config_snapshot else NULL
  g <- NULL
  if (!is.null(cs)) {
    g <- cs$match_dist_threshold_ldir_um
    if (is.null(g)) g <- cs$match_dist_threshold_um
  }
  g <- suppressWarnings(as.numeric(unlist(g)))
  g <- g[is.finite(g)]
  if (length(g) == 0 || g[1] <= 0) return(default_gate)
  g[1]
}

# Tag each LDIR<->Raman pair with `within_gate` (TRUE = genuine match). Records
# the resolved gate on data$ldir_match_gate_um for the summary line.
annotate_ldir_gate <- function(data, run_dir) {
  m <- data$ldir_raman_matched
  if (is.null(m) || nrow(m) == 0) return(data)
  gate <- ldir_acceptance_gate(run_dir)
  data$ldir_match_gate_um <- gate
  if ("match_distance" %in% names(m)) {
    m$within_gate <- !is.na(m$match_distance) & m$match_distance <= gate
  } else {
    # Legacy CSVs without a per-pair distance: preserve prior behaviour
    # (treat all forced pairs as matched) rather than silently dropping them.
    m$within_gate <- TRUE
  }
  data$ldir_raman_matched <- m
  data
}

# Recompute the within_gate classification against an explicit gate (um). Used
# by the viewer's live gate slider so the user can retune the acceptance
# distance in real time; annotate_ldir_gate() seeds the default from the
# manifest, this overrides it with the slider value.
regate_ldir <- function(data, gate) {
  m <- data$ldir_raman_matched
  if (is.null(m) || nrow(m) == 0) return(data)
  if (is.null(gate) || !is.finite(gate) || gate <= 0) return(data)
  data$ldir_match_gate_um <- gate
  if ("match_distance" %in% names(m)) {
    m$within_gate <- !is.na(m$match_distance) & m$match_distance <= gate
    data$ldir_raman_matched <- m
  }
  data
}

# Subset of an ldir_raman_matched frame that are genuine (within-gate) matches.
# When the frame has not been annotated (defensive: no within_gate column) every
# pair is treated as genuine, preserving legacy behaviour.
ldir_genuine_pairs <- function(m) {
  if (is.null(m) || nrow(m) == 0) return(m)
  if (!"within_gate" %in% names(m)) return(m)
  m[isTRUE_vec(m$within_gate), , drop = FALSE]
}

# Vectorised isTRUE (NA-safe): TRUE only where the value is exactly TRUE.
isTRUE_vec <- function(x) !is.na(x) & x

# Helper: add <prefix>_material_family columns for every <prefix>_material
# column in the named match tables, so the overlay can filter any instrument
# side on harmonized family names -- including the cross-instrument pair tables.
enrich_material_family <- function(data) {
  match_tables <- c("matched", "matched_ftir_bruker", "ldir_raman_matched",
                    "matched_bruker_perkin", "matched_bruker_ldir",
                    "matched_perkin_ldir")
  for (tbl in match_tables) {
    d <- data[[tbl]]
    if (is.null(d) || nrow(d) == 0) next
    mat_cols <- grep("_material$", names(d), value = TRUE)
    for (mc in mat_cols) {
      fam <- sub("_material$", "_material_family", mc)
      d[[fam]] <- classify_family_vec(d[[mc]])
    }
    data[[tbl]] <- d
  }
  data
}

# ---------------------------------------------------------------------------
# Build run_data from user-uploaded CSVs (fallback when pipeline output is
# not available).  Minimum required: matched_particles.csv.
# ---------------------------------------------------------------------------
load_uploaded_data <- function(matched_path,
                                unmatched_ftir_path = NULL,
                                unmatched_raman_path = NULL,
                                transform_path = NULL) {
  data <- list(run_dir = dirname(matched_path))
  data$matched <- read.csv(matched_path, stringsAsFactors = FALSE)

  if (!is.null(unmatched_ftir_path) && file.exists(unmatched_ftir_path))
    data$unmatched_ftir <- read.csv(unmatched_ftir_path, stringsAsFactors = FALSE)
  if (!is.null(unmatched_raman_path) && file.exists(unmatched_raman_path))
    data$unmatched_raman <- read.csv(unmatched_raman_path, stringsAsFactors = FALSE)
  if (!is.null(transform_path) && file.exists(transform_path))
    data$transform <- parse_transform_params(transform_path)

  enrich_material_family(data)
}

# ---------------------------------------------------------------------------
# Parse transform_params.txt into a structured list
# ---------------------------------------------------------------------------
parse_transform_params <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)

  get_val <- function(key) {
    pat <- paste0("^", key, ":\\s*")
    idx <- grep(pat, lines)
    if (length(idx) == 0) return(NA)
    trimws(sub(pat, "", lines[idx[1]]))
  }

  get_num <- function(key) as.numeric(get_val(key))

  get_row <- function(key) {
    v <- get_val(key)
    if (is.na(v)) return(NULL)
    as.numeric(strsplit(gsub("\\s+", "", v), ",")[[1]])
  }

  # 3x3 transform matrix (FTIR_norm -> Raman_norm)
  r1 <- get_row("matrix_row1")
  r2 <- get_row("matrix_row2")
  r3 <- get_row("matrix_row3")
  M <- NULL
  if (!is.null(r1) && !is.null(r2) && !is.null(r3)) {
    M <- rbind(r1, r2, r3)
  }

  # FTIR scan bounds (optional -- present in newer pipeline output)
  scan_xmin <- get_num("ftir_scan_xmin")
  scan_xmax <- get_num("ftir_scan_xmax")
  scan_ymin <- get_num("ftir_scan_ymin")
  scan_ymax <- get_num("ftir_scan_ymax")
  scan_bounds <- NULL
  if (!is.na(scan_xmax)) {
    scan_bounds <- list(xmin = scan_xmin, xmax = scan_xmax,
                        ymin = scan_ymin, ymax = scan_ymax)
  }

  list(
    scale         = get_num("scale"),
    rotation_deg  = get_num("rotation_deg"),
    reflected     = tolower(get_val("reflected")) == "true",
    M             = M,
    ftir_centroid  = c(get_num("ftir_centroid_x"), get_num("ftir_centroid_y")),
    raman_centroid = c(get_num("raman_centroid_x"), get_num("raman_centroid_y")),
    ftir_scan_bounds = scan_bounds
  )
}

# ---------------------------------------------------------------------------
# Build the full 3x3 transform: FTIR original coords -> Raman coords
#   x_aligned = M * (x_orig - ftir_centroid) + raman_centroid
# ---------------------------------------------------------------------------
build_full_transform <- function(transform) {
  # Both FTIR aligned and Raman are displayed in the normalized (centered) frame.
  # Pipeline: x_norm = x_orig - ftir_centroid, then x_aligned = M * x_norm.
  # So full transform: subtract ftir centroid, then apply M.  No raman centroid added.
  T1 <- diag(3)
  T1[1, 3] <- -transform$ftir_centroid[1]
  T1[2, 3] <- -transform$ftir_centroid[2]

  # Full: M %*% T1  (maps FTIR original -> normalized aligned frame)
  transform$M %*% T1
}

# ---------------------------------------------------------------------------
# Transform 2D points with a 3x3 homogeneous matrix
# ---------------------------------------------------------------------------
transform_points <- function(x, y, M_full) {
  pts <- rbind(x, y, rep(1, length(x)))
  result <- M_full %*% pts
  list(x = result[1, ], y = result[2, ])
}

# ---------------------------------------------------------------------------
# Estimate FTIR scan bounds from image dimensions and particle coordinates.
#
# The PerkinElmer Spotlight exports images at ~6 rendering pixels per 25um
# grid cell.  From a 2993x2993 image: (2993+1)/6 ~ 499 grid positions,
# giving a 499 * 25 = 12475 um scan extent.  This function computes the
# bounds robustly from the image dimensions and grid step.
# ---------------------------------------------------------------------------
estimate_ftir_scan_bounds <- function(img_raster, particle_x_um = NULL,
                                       particle_y_um = NULL,
                                       grid_step_um = 25) {
  img_w <- ncol(img_raster)
  img_h <- nrow(img_raster)

  # Estimate grid positions from image pixel count
  # Empirical: (image_px + 1) / 6 gives the grid count
  render_px_per_cell <- 6
  grid_nx <- round((img_w + 1) / render_px_per_cell)
  grid_ny <- round((img_h + 1) / render_px_per_cell)

  # Scan extent: grid_count * grid_step
  x_extent <- grid_nx * grid_step_um
  y_extent <- grid_ny * grid_step_um

  # Sanity check against particle positions if available
  if (!is.null(particle_x_um)) {
    max_px <- max(particle_x_um, na.rm = TRUE)
    if (x_extent < max_px) x_extent <- ceiling(max_px / 500) * 500
  }
  if (!is.null(particle_y_um)) {
    max_py <- max(particle_y_um, na.rm = TRUE)
    if (y_extent < max_py) y_extent <- ceiling(max_py / 500) * 500
  }

  # Center the image extent on the particle distribution midpoint.
  # The scan origin is unknown, so centering on particle positions is the
  # best estimate (particles can only appear within the scan area).
  if (!is.null(particle_x_um) && length(particle_x_um) > 0) {
    mid_x <- (min(particle_x_um, na.rm = TRUE) +
              max(particle_x_um, na.rm = TRUE)) / 2
    xmin <- mid_x - x_extent / 2
    xmax <- mid_x + x_extent / 2
  } else {
    xmin <- 0; xmax <- x_extent
  }

  if (!is.null(particle_y_um) && length(particle_y_um) > 0) {
    mid_y <- (min(particle_y_um, na.rm = TRUE) +
              max(particle_y_um, na.rm = TRUE)) / 2
    ymin <- mid_y - y_extent / 2
    ymax <- mid_y + y_extent / 2
  } else {
    ymin <- 0; ymax <- y_extent
  }

  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

# ---------------------------------------------------------------------------
# Build per-instrument data frames from pipeline output
# ---------------------------------------------------------------------------
build_instrument_dfs <- function(data) {
  result <- list()

  # --- FTIR (PerkinElmer) ---
  ftir_parts <- list()
  if (!is.null(data$matched) && nrow(data$matched) > 0) {
    m <- data$matched
    ftir_parts[[1]] <- data.frame(
      particle_id  = m$ftir_particle_id,
      x = m$ftir_x_aligned, y = m$ftir_y_aligned,
      x_orig = m$ftir_x_um, y_orig = m$ftir_y_um,
      area_um2 = m$ftir_area_um2,
      major_um = m$ftir_major_um, minor_um = m$ftir_minor_um,
      feret_max = m$ftir_feret_max_um,
      material = m$ftir_material, quality = m$ftir_quality,
      match_status = "matched", match_id = m$match_id,
      matched_to_raman = TRUE,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_ftir) && nrow(data$unmatched_ftir) > 0) {
    u <- data$unmatched_ftir
    ftir_parts[[length(ftir_parts) + 1]] <- data.frame(
      particle_id = u$particle_id,
      x = u$x_aligned, y = u$y_aligned,
      x_orig = u$x_um, y_orig = u$y_um,
      area_um2 = u$area_um2,
      major_um = u$major_um, minor_um = u$minor_um,
      feret_max = u$feret_max_um,
      material = u$material, quality = u$quality,
      match_status = "unmatched", match_id = NA_integer_,
      matched_to_raman = FALSE,
      stringsAsFactors = FALSE)
  }
  result$ftir <- do.call(rbind, ftir_parts)
  if (!is.null(result$ftir))
    result$ftir$material_family <- classify_family_vec(result$ftir$material)

  # --- Raman ---
  raman_parts <- list()
  if (!is.null(data$matched) && nrow(data$matched) > 0) {
    m <- data$matched
    raman_parts[[1]] <- data.frame(
      particle_id = m$raman_particle_id,
      x = m$raman_x_norm, y = m$raman_y_norm,
      x_orig = m$raman_x_um, y_orig = m$raman_y_um,
      area_um2 = m$raman_area_um2,
      major_um = m$raman_major_um, minor_um = m$raman_minor_um,
      feret_max = m$raman_feret_max_um,
      material = m$raman_material, quality = m$raman_quality,
      match_status = "matched", match_id = m$match_id,
      matched_to_ftir_perkin = TRUE,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_raman) && nrow(data$unmatched_raman) > 0) {
    u <- data$unmatched_raman
    raman_parts[[length(raman_parts) + 1]] <- data.frame(
      particle_id = u$particle_id,
      x = u$x_norm, y = u$y_norm,
      x_orig = u$x_um, y_orig = u$y_um,
      area_um2 = u$area_um2,
      major_um = u$major_um, minor_um = u$minor_um,
      feret_max = u$feret_max_um,
      material = u$material, quality = u$quality,
      match_status = "unmatched", match_id = NA_integer_,
      matched_to_ftir_perkin = FALSE,
      stringsAsFactors = FALSE)
  }
  result$raman <- do.call(rbind, raman_parts)
  if (!is.null(result$raman))
    result$raman$material_family <- classify_family_vec(result$raman$material)

  # Add LDIR->Raman match flag if LDIR-Raman match data available. Only
  # genuine (within-gate) pairs count as matched -- an over-gate forced pairing
  # leaves the Raman particle effectively unmatched to LDIR.
  if (!is.null(result$raman) && nrow(result$raman) > 0 &&
      !is.null(data$ldir_raman_matched) && nrow(data$ldir_raman_matched) > 0) {
    genuine <- ldir_genuine_pairs(data$ldir_raman_matched)
    result$raman$matched_to_ldir <- result$raman$particle_id %in%
      genuine$raman_particle_id
  } else if (!is.null(result$raman) && nrow(result$raman) > 0) {
    result$raman$matched_to_ldir <- FALSE
  }

  # --- LDIR ---
  ldir_parts <- list()
  if (!is.null(data$ldir_raman_matched) && nrow(data$ldir_raman_matched) > 0) {
    m <- data$ldir_raman_matched
    # Safeguard: ensure aligned columns exist; fall back to originals with warning
    has_aligned <- "ldir_x_aligned" %in% names(m)
    if (!has_aligned) {
      message("[Particle Viewer] WARNING: ldir_raman_matched.csv missing ldir_x_aligned; ",
              "falling back to ldir_x_um (original coordinates)")
    }
    # Genuine match = within the acceptance gate. Over-gate rows are forced
    # pairings (ldir_force_complete_match); report them as unmatched LDIR so the
    # summary, the overlay (lone red point), and the tables all agree.
    within <- if ("within_gate" %in% names(m)) isTRUE_vec(m$within_gate) else rep(TRUE, nrow(m))
    ldir_parts[[1]] <- data.frame(
      particle_id = m$ldir_particle_id,
      x = if (has_aligned) m$ldir_x_aligned else m$ldir_x_um,
      y = if (has_aligned) m$ldir_y_aligned else m$ldir_y_um,
      x_orig = m$ldir_x_um, y_orig = m$ldir_y_um,
      area_um2 = m$ldir_area_um2,
      major_um = m$ldir_major_um, minor_um = m$ldir_minor_um,
      feret_max = m$ldir_feret_max_um,
      material = m$ldir_material, quality = m$ldir_quality,
      match_status = ifelse(within, "matched", "unmatched"),
      match_id = ifelse(within, m$match_id, NA_integer_),
      matched_to_raman = within,
      match_score       = if ("match_score"           %in% names(m)) m$match_score           else NA_real_,
      coord_match_cost  = if ("ldir_coord_match_cost" %in% names(m)) m$ldir_coord_match_cost else NA_real_,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_ldir) && nrow(data$unmatched_ldir) > 0) {
    u <- data$unmatched_ldir
    ldir_parts[[length(ldir_parts) + 1]] <- data.frame(
      particle_id = u$particle_id,
      x = if ("x_aligned" %in% names(u)) u$x_aligned else u$x_um,
      y = if ("y_aligned" %in% names(u)) u$y_aligned else u$y_um,
      x_orig = u$x_um, y_orig = u$y_um,
      area_um2 = u$area_um2,
      major_um = u$major_um, minor_um = u$minor_um,
      feret_max = u$feret_max_um,
      material = u$material, quality = u$quality,
      match_status = "unmatched", match_id = NA_integer_,
      matched_to_raman = FALSE,
      match_score      = NA_real_,
      coord_match_cost = if ("coord_match_cost" %in% names(u)) u$coord_match_cost else NA_real_,
      stringsAsFactors = FALSE)
  }
  result$ldir <- if (length(ldir_parts) > 0) do.call(rbind, ldir_parts) else NULL
  if (!is.null(result$ldir))
    result$ldir$material_family <- classify_family_vec(result$ldir$material)

  # --- FTIR Bruker (matched + unmatched, with aligned coordinates when available) ---
  ftir_bruker_parts <- list()
  if (!is.null(data$matched_ftir_bruker) && nrow(data$matched_ftir_bruker) > 0) {
    m <- data$matched_ftir_bruker
    ftir_bruker_parts[[1]] <- data.frame(
      particle_id      = m$ftir_particle_id,
      x                = m$ftir_x_aligned,
      y                = m$ftir_y_aligned,
      x_orig           = m$ftir_x_um,
      y_orig           = m$ftir_y_um,
      area_um2         = m$ftir_area_um2,
      major_um         = m$ftir_major_um,
      minor_um         = m$ftir_minor_um,
      feret_max        = m$ftir_feret_max_um,
      material         = m$ftir_material,
      quality          = m$ftir_quality,
      match_status     = "matched",
      match_id         = m$match_id,
      matched_to_raman = TRUE,
      stringsAsFactors = FALSE)
  }
  if (!is.null(data$unmatched_ftir_bruker) && nrow(data$unmatched_ftir_bruker) > 0) {
    u <- data$unmatched_ftir_bruker
    has_aligned <- "x_aligned" %in% names(u)
    ftir_bruker_parts[[length(ftir_bruker_parts) + 1]] <- data.frame(
      particle_id      = u$particle_id,
      x                = if (has_aligned) u$x_aligned else u$x_um,
      y                = if (has_aligned) u$y_aligned else u$y_um,
      x_orig           = u$x_um,
      y_orig           = u$y_um,
      area_um2         = u$area_um2,
      major_um         = u$major_um,
      minor_um         = u$minor_um,
      feret_max        = u$feret_max_um,
      material         = u$material,
      quality          = u$quality,
      match_status     = "unmatched",
      match_id         = NA_integer_,
      matched_to_raman = FALSE,
      stringsAsFactors = FALSE)
  }
  result$ftir_bruker <- if (length(ftir_bruker_parts) > 0) do.call(rbind, ftir_bruker_parts) else NULL
  if (!is.null(result$ftir_bruker))
    result$ftir_bruker$material_family <- classify_family_vec(result$ftir_bruker$material)

  result
}


# ---------------------------------------------------------------------------
# Plastics summary using canonical classify_family_vec() from 08b_material_map.R
# ---------------------------------------------------------------------------

#' Count unique particles per polymer family and category in a device data frame.
#'
#' Uses classify_family_vec() to map raw material names to canonical families,
#' then classify_category() to group into Synthetic/Semi-synthetic/Natural.
#'
#' @param df Device data frame from build_instrument_dfs() (must have a
#'   \code{material} column).
#' @return data.frame(family, category, n) sorted by category then n desc,
#'   or an empty frame when no particles have classifiable materials.
summarise_plastics <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !"material" %in% names(df))
    return(data.frame(family = character(0), category = character(0),
                      n = integer(0)))
  families <- classify_family_vec(df$material)
  categories <- classify_category_vec(families)
  # Exclude Unknown
  keep <- families != "Unknown"
  if (!any(keep))
    return(data.frame(family = character(0), category = character(0),
                      n = integer(0)))
  tbl <- as.data.frame(table(family = families[keep]), stringsAsFactors = FALSE)
  names(tbl) <- c("family", "n")
  tbl$n <- as.integer(tbl$n)
  tbl$category <- classify_category_vec(tbl$family)
  # Sort: Synthetic first, then Semi-synthetic, then Natural, within each by n desc
  cat_order <- c("Synthetic", "Semi-synthetic", "Natural/Organic")
  tbl$cat_rank <- match(tbl$category, cat_order, nomatch = 99)
  tbl <- tbl[order(tbl$cat_rank, -tbl$n), ]
  tbl$cat_rank <- NULL
  rownames(tbl) <- NULL
  tbl
}

# ---------------------------------------------------------------------------
# Coordinate bounds from particle data (safe for NULL / empty inputs)
# ---------------------------------------------------------------------------
compute_bounds <- function(...) {
  dfs <- list(...)
  all_x <- unlist(lapply(dfs, function(d) if (!is.null(d) && nrow(d) > 0) d$x))
  all_y <- unlist(lapply(dfs, function(d) if (!is.null(d) && nrow(d) > 0) d$y))
  all_x <- all_x[is.finite(all_x)]
  all_y <- all_y[is.finite(all_y)]
  if (length(all_x) == 0 || length(all_y) == 0) {
    return(list(x = c(-1000, 1000), y = c(-1000, 1000)))
  }
  pad <- 200
  list(
    x = c(min(all_x) - pad, max(all_x) + pad),
    y = c(min(all_y) - pad, max(all_y) + pad)
  )
}



# Detect image format from file signature (magic bytes)
# Returns: PNG, JPEG, TIFF, BMP, WEBP, or unknown
sniff_image_type <- function(path) {
  if (is.null(path) || !file.exists(path)) return("unknown")
  hdr <- tryCatch(as.integer(readBin(path, "raw", n = 12)), error = function(e) integer(0))
  if (length(hdr) < 4) return("unknown")
  if (hdr[1] == 137 && hdr[2] == 80 && hdr[3] == 78 && hdr[4] == 71) return("PNG")
  if (hdr[1] == 255 && hdr[2] == 216 && hdr[3] == 255) return("JPEG")
  if ((hdr[1] == 73 && hdr[2] == 73 && hdr[3] == 42 && hdr[4] == 0) ||
      (hdr[1] == 77 && hdr[2] == 77 && hdr[3] == 0 && hdr[4] == 42)) return("TIFF")
  if (hdr[1] == 66 && hdr[2] == 77) return("BMP")
  if (length(hdr) >= 12 && hdr[1] == 82 && hdr[2] == 73 && hdr[3] == 70 && hdr[4] == 70 &&
      hdr[9] == 87 && hdr[10] == 69 && hdr[11] == 66 && hdr[12] == 80) return("WEBP")
  "unknown"
}

# ---------------------------------------------------------------------------
# Load an image (any format) as a raster array for ggplot annotation.
# Accepts PNG, JPEG, TIFF, BMP, WEBP regardless of extension.
# Uses magick for broadest format support, with pkg-specific fallbacks.
# annotation_raster places row 1 at ymax (top of plot).  Standard images
# already have row 1 = top of the visual image, which corresponds to max-Y
# in Cartesian / stage coordinates.  So NO vertical flip is needed.
# ---------------------------------------------------------------------------
load_image_raster <- function(path) {
  # fileInput(multiple = TRUE) delivers a vector of datapaths; the scalar
  # conditions below would then error ("condition has length > 1") and the
  # calling observer would die without ever setting the image.
  if (length(path) > 1) {
    message("[Particle Viewer] load_image_raster: got ", length(path),
            " paths; using the first one only")
    path <- path[1]
  }
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)

  typ <- sniff_image_type(path)

  # PNG: always use png::readPNG -- avoids magick's image_data() producing
  # tiled/colour-distorted arrays for RGBA PNGs (e.g. FTIR false-colour images).
  if (typ == "PNG") {
    return(tryCatch(png::readPNG(path), error = function(e) NULL))
  }

  # JPEG: prefer the jpeg package (correct sRGB, no intermediate conversion).
  if (typ == "JPEG" && requireNamespace("jpeg", quietly = TRUE)) {
    return(tryCatch(jpeg::readJPEG(path), error = function(e) NULL))
  }

  # Other formats (TIFF, BMP, WEBP, ...): use magick when available.
  if (requireNamespace("magick", quietly = TRUE)) {
    arr <- tryCatch({
      img_mg  <- magick::image_read(path)
      if (length(img_mg) > 1) img_mg <- img_mg[1]
      img_rgb <- magick::image_convert(img_mg, colorspace = "sRGB")
      raw_data <- magick::image_data(img_rgb, channels = "rgba")
      vals <- if (is.character(raw_data)) strtoi(raw_data, base = 16L)
              else as.integer(raw_data)
      arr <- array(vals / 255, dim = dim(raw_data))
      aperm(arr, c(3, 2, 1))   # [4, W, H] -> [H, W, 4] (RGBA raster)
    }, error = function(e) NULL)
    if (!is.null(arr)) return(arr)
  }

  # BMP without magick (or magick failed): dependency-free reader from
  # R/read_bmp.R covers the uncompressed 8/24/32-bit BMPs that instrument
  # software exports.
  if (typ == "BMP") {
    arr <- tryCatch(read_bmp_raster(path), error = function(e) NULL)
    if (!is.null(arr)) return(arr)
    message("[Particle Viewer] read_bmp_raster could not decode ",
            basename(path), " (compressed or unusual BMP variant)")
  }

  # Last-resort extension-based fallback
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tif", "tiff", "bmp", "webp"))
    warning("Install magick for TIFF/WEBP and compressed-BMP support ",
            "in the Shiny viewer.")
  raw <- tryCatch(png::readPNG(path), error = function(e) NULL)
  if (is.null(raw) && requireNamespace("jpeg", quietly = TRUE))
    raw <- tryCatch(jpeg::readJPEG(path), error = function(e) NULL)
  raw
}

# ---------------------------------------------------------------------------
# Downsample a raster array so its longest edge is <= max_dim pixels.
# Uses block-average (box filter) pooling -- antialiased and dependency-free --
# with stride subsampling as a fallback for degenerate aspect ratios.
# Accepts 2D (grayscale) or 3D (H x W x channels) arrays; returns same form.
# The original pixel dimensions are recorded as attributes so callers that
# convert pixels to um (e.g. TIFF DPI metadata, which refers to the ORIGINAL
# file) can correct their scale for the reduced raster.
# ---------------------------------------------------------------------------
downsample_raster <- function(raw, max_dim = BG_IMAGE_MAX_DIM) {
  d <- dim(raw)
  if (is.null(d) || length(d) < 2) return(raw)
  h <- d[1]; w <- d[2]
  if (max(h, w) <= max_dim) return(raw)
  k <- ceiling(max(h, w) / max_dim)

  out <- if (h >= k && w >= k) {
    pool <- function(m) {
      hh <- (nrow(m) %/% k) * k
      ww <- (ncol(m) %/% k) * k
      m <- m[seq_len(hh), seq_len(ww), drop = FALSE]
      m <- rowsum(m, rep(seq_len(hh %/% k), each = k)) / k
      t(rowsum(t(m), rep(seq_len(ww %/% k), each = k))) / k
    }
    if (length(d) == 2) pool(raw)
    else vapply(seq_len(d[3]), function(ch) pool(raw[, , ch]),
                matrix(0, h %/% k, w %/% k))
  } else {
    # Extreme aspect ratio: block pooling would collapse the short axis
    rows <- seq(1, h, by = k)
    cols <- seq(1, w, by = k)
    if (length(d) == 2) raw[rows, cols, drop = FALSE]
    else raw[rows, cols, , drop = FALSE]
  }

  attr(out, "orig_width_px")  <- w
  attr(out, "orig_height_px") <- h
  out
}

# ---------------------------------------------------------------------------
# Auto-detect um-per-pixel scale from TIFF resolution metadata.
# Returns NULL silently when: not a TIFF, magick unavailable, metadata absent,
# or the value looks like a screen-default (72/96/150 DPI) rather than a real
# instrument-calibrated resolution.
# ---------------------------------------------------------------------------
extract_tiff_um_per_px <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  if (sniff_image_type(path) != "TIFF") return(NULL)
  if (!requireNamespace("magick", quietly = TRUE)) return(NULL)
  tryCatch({
    img   <- magick::image_read(path)
    info  <- magick::image_info(img)
    units <- info$units
    if (is.na(units) || units == "Undefined") return(NULL)

    # density may be "NxN" string or numeric
    dens <- info$density
    if (is.character(dens)) dens <- as.numeric(strsplit(dens, "x")[[1]][1])
    if (is.na(dens) || dens <= 0) return(NULL)

    # Convert to um/pixel
    if (units == "PixelsPerCentimeter") {
      um_per_px <- 10000 / dens   # 1 cm = 10 000 \u00b5m
    } else {
      um_per_px <- 25400 / dens   # 1 inch = 25 400 \u00b5m
    }

    # Reject common screen defaults -- these are never real instrument values
    screen_dpis <- c(72, 96, 150, 300)
    effective_dpi <- if (units == "PixelsPerCentimeter") dens * 2.54 else dens
    if (round(effective_dpi) %in% screen_dpis) return(NULL)

    # Sanity: instrument images typically 0.5-50 um/px
    if (um_per_px > 0 && um_per_px < 200) um_per_px else NULL
  }, error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Resolve the Raman image's physical extent from the WITec Particle Scout
# values stored in the run manifest (config_snapshot$raman_image_width_um /
# height_um / center_x_um / center_y_um).
#
# WITec's panel reports the image center in its video/image frame, whose Y
# axis points DOWN, while the particle export ("Visual Center Point Y") is
# in stage coordinates with Y UP -- the stage-frame center is (cx, -cy).
# Verified on real data (PET A, 2026-07: particles Y [-68, 4877], panel
# Center Y = -4394): with Y negated 100% of particles fall inside the image;
# taken as-reported only 30% do.  Because the convention may vary across
# WITec versions/exports, BOTH interpretations are scored by the fraction of
# particles they contain and the better one wins; below min_frac the
# function returns NULL and the caller falls back to heuristic placement.
#
# cfg            : config_snapshot list from the run manifest
# x_orig, y_orig : particle stage coordinates (um) used to score candidates
# min_frac       : minimum containment fraction to accept
# Returns list(xmin, xmax, ymin, ymax, y_negated, frac_inside) or NULL.
# ---------------------------------------------------------------------------
# Coarse um-per-pixel estimate from the analysed particle areas: in a
# dark-field micrograph the bright (particle) pixels should cover the same
# physical area the instrument reported for those particles, so
# sqrt(sum(area) / bright_px) recovers the scale -- the same self-calibration
# the LDIR processed-image join uses.
#
# Accuracy is only about +/-35% (it moves with where the brightness threshold
# falls), so this is a SANITY CHECK against a grossly misplaced image, never a
# way to set the scale. Returns NULL when the image is not dark-field enough
# for the pixel count to mean anything.
image_scale_from_particle_area <- function(raster, areas_um2, thr = 0.45) {
  if (is.null(raster) || length(dim(raster)) < 2) return(NULL)
  a <- sum(areas_um2[is.finite(areas_um2)])
  if (!is.finite(a) || a <= 0) return(NULL)
  g <- if (length(dim(raster)) == 3) {
    m <- raster[, , 1]
    for (k in seq_len(dim(raster)[3])[-1]) m <- pmax(m, raster[, , k])
    m
  } else raster
  n <- sum(g > thr, na.rm = TRUE)
  if (n < 50 || n / length(g) > 0.20) return(NULL)
  sqrt(a / n)
}

raman_image_extent_from_config <- function(cfg, x_orig, y_orig, min_frac = 0.5) {
  w  <- cfg$raman_image_width_um
  h  <- cfg$raman_image_height_um
  cx <- cfg$raman_image_center_x_um
  cy <- cfg$raman_image_center_y_um
  ok <- function(v) !is.null(v) && is.numeric(v) && length(v) == 1 && is.finite(v)
  if (!ok(w) || !ok(h) || w <= 0 || h <= 0 || !ok(cx) || !ok(cy)) return(NULL)

  fin <- is.finite(x_orig) & is.finite(y_orig)
  x <- x_orig[fin]; y <- y_orig[fin]

  score <- function(cy_stage) {
    ext <- list(xmin = cx - w / 2, xmax = cx + w / 2,
                ymin = cy_stage - h / 2, ymax = cy_stage + h / 2)
    ext$frac_inside <- if (length(x) == 0) 1 else
      mean(x >= ext$xmin & x <= ext$xmax & y >= ext$ymin & y <= ext$ymax)
    ext
  }

  neg <- score(-cy)   # panel value in Y-down video frame -> negate (expected)
  raw <- score(cy)    # panel value already in stage frame
  best <- if (neg$frac_inside >= raw$frac_inside) {
    c(neg, list(y_negated = TRUE))
  } else {
    c(raw, list(y_negated = FALSE))
  }
  if (best$frac_inside < min_frac) return(NULL)
  best
}

# ---------------------------------------------------------------------------
# Multi-Run image placement -- mirror each single-instrument tab so the
# reproducibility overlay reproduces the exact image<->coordinate relationship.
# Each returns a bare extent list(xmin,xmax,ymin,ymax) in the run-1 (raw) frame,
# or NULL when the required metadata is absent.
# ---------------------------------------------------------------------------

# Single numeric field from a one-row meta data frame (NA if absent/non-finite).
.repro_meta_num <- function(meta, field) {
  if (is.null(meta) || !field %in% names(meta)) return(NA_real_)
  v <- suppressWarnings(as.numeric(meta[[field]][1]))
  if (length(v) == 0 || !is.finite(v)) NA_real_ else v
}

# FTIR / Bruker P1: physical extent recorded by tools/reproducibility.R, in the
# native FTIR scan frame (coordinates are um from the scan origin, so the image
# spans [0,w] x [0,h] unless an explicit centre is given). Resize-invariant:
# re-exporting the image at a different pixel resolution does not move it.
# Mirrors the Raman WITec tier. NULL when the metadata is absent.
place_image_ftir_meta <- function(meta, x, y) {
  w <- .repro_meta_num(meta, "ftir_image_width_um")
  h <- .repro_meta_num(meta, "ftir_image_height_um")
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) return(NULL)
  cx <- .repro_meta_num(meta, "ftir_image_center_x_um")
  cy <- .repro_meta_num(meta, "ftir_image_center_y_um")
  if (!is.finite(cx)) cx <- w / 2      # default: scan origin at (0,0)
  if (!is.finite(cy)) cy <- h / 2
  list(xmin = cx - w / 2, xmax = cx + w / 2,
       ymin = cy - h / 2, ymax = cy + h / 2)
}

# FTIR / Bruker P2: aspect-preserving fit to the particle extent.
#
# `raw` is REQUIRED to preserve the aspect ratio. Without it this returns the
# bare particle bounding box, which is what the multi-run overlay used to do --
# annotation_raster() stretches the image to whatever box it is given, so a
# non-square particle hull sheared the micrograph and its features stopped
# lining up with the points (a rectangular scan squeezed into a square hull
# renders as mismatched bands). tools/reproducibility.R documents this backdrop
# as being placed "at the point extent (aspect-preserving)"; passing `raw` is
# what actually honours that.
place_image_particle_extent <- function(x, y, raw = NULL) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NULL)
  if (is.null(raw))
    return(list(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)))
  # A single point (or perfectly coincident points) has no extent to fit to.
  if (max(x) - min(x) <= 0 && max(y) - min(y) <= 0) return(NULL)
  compute_image_bounds(raw, x, y, padding_um = 0)
}

# Raman P1: WITec width/height/center from meta; Y auto-detected against the
# particles (reuses raman_image_extent_from_config, exactly like the Raman tab).
place_image_raman_meta <- function(meta, x, y) {
  cfg <- list(
    raman_image_width_um    = .repro_meta_num(meta, "raman_image_width_um"),
    raman_image_height_um   = .repro_meta_num(meta, "raman_image_height_um"),
    raman_image_center_x_um = .repro_meta_num(meta, "raman_image_center_x_um"),
    raman_image_center_y_um = .repro_meta_num(meta, "raman_image_center_y_um"))
  ext <- raman_image_extent_from_config(cfg, x, y)
  if (is.null(ext)) return(NULL)
  list(xmin = ext$xmin, xmax = ext$xmax, ymin = ext$ymin, ymax = ext$ymax)
}

# Raman P2: known um-per-pixel scale (from meta, or read from a TIFF backdrop),
# centred on the particle mean -- mirrors raman_native_image_info Priority 2. This
# is the tier the Raman tab uses when no WITec metadata is present.
place_image_raman_umpx <- function(meta, x, y, raw, bg_path = NULL) {
  if (is.null(raw)) return(NULL)
  upp <- .repro_meta_num(meta, "raman_um_per_px")
  if (!is.finite(upp) && !is.null(bg_path))
    upp <- tryCatch(extract_tiff_um_per_px(bg_path), error = function(e) NULL)
  if (is.null(upp) || !is.finite(upp) || upp <= 0) return(NULL)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NULL)
  cx <- mean(x); cy <- mean(y)
  hw <- ncol(raw) * upp / 2; hh <- nrow(raw) * upp / 2
  list(xmin = cx - hw, xmax = cx + hw, ymin = cy - hh, ymax = cy + hh)
}

# LDIR: scan-circle calibration from meta (mirrors ldir_native_image_info).
place_image_ldir_meta <- function(meta) {
  s  <- .repro_meta_num(meta, "ldir_scale_um_per_px")
  cx <- .repro_meta_num(meta, "ldir_cx_px");  cy <- .repro_meta_num(meta, "ldir_cy_px")
  w  <- .repro_meta_num(meta, "ldir_image_width_px")
  h  <- .repro_meta_num(meta, "ldir_image_height_px")
  if (any(!is.finite(c(s, cx, cy, w, h))) || s <= 0) return(NULL)
  list(xmin = -cx * s, xmax = (w - cx) * s, ymin = (cy - h) * s, ymax = cy * s)
}

# Dispatch to the instrument-appropriate placement; NULL if unavailable. For
# Raman this runs the same cascade as raman_native_image_info: WITec extent
# (P1) then um-per-pixel scale (P2); P3 (particle-extent fit) is left to the
# caller's fallback. FTIR/Bruker run the analogous two tiers: recorded physical
# extent (P1) then an aspect-preserving fit to the particle extent (P2).
# `raw` is needed by the Raman P2 tier and by the FTIR P2 fit (which cannot
# preserve the aspect ratio without knowing the raster's pixel dimensions);
# `bg_path` only by Raman P2.
place_image_multirun <- function(instrument, meta, x, y, raw = NULL, bg_path = NULL) {
  switch(as.character(instrument),
    raman = {
      ext <- place_image_raman_meta(meta, x, y)                 # P1: WITec
      if (is.null(ext)) ext <- place_image_raman_umpx(meta, x, y, raw, bg_path)  # P2
      ext
    },
    ftir_perkin = {
      ext <- place_image_ftir_meta(meta, x, y)                  # P1: physical
      if (is.null(ext)) ext <- place_image_particle_extent(x, y, raw)  # P2: fit
      ext
    },
    ftir_bruker = {
      ext <- place_image_ftir_meta(meta, x, y)
      if (is.null(ext)) ext <- place_image_particle_extent(x, y, raw)
      ext
    },
    ldir        = place_image_ldir_meta(meta),
    NULL)
}

# ---------------------------------------------------------------------------
# LDIR view rotation -- rotate the whole native LDIR scene (image raster,
# extent, particle coordinates) by a multiple of 90 deg about the origin so
# the LDIR tab can be displayed in the Raman orientation for side-by-side
# comparison.  Display-only: no stored coordinate is modified.
# deg convention: +90 = counter-clockwise, -90 = clockwise, in {0,90,-90,180}.
# ---------------------------------------------------------------------------
rotate_xy_view <- function(x, y, deg) {
  switch(as.character(((deg %% 360) + 360) %% 360),
         "90"  = list(x = -y, y =  x),
         "180" = list(x = -x, y = -y),
         "270" = list(x =  y, y = -x),
         list(x = x, y = y))
}

rotate_extent_view <- function(ext, deg) {
  d <- ((deg %% 360) + 360) %% 360
  if (d == 90)  return(list(xmin = -ext$ymax, xmax = -ext$ymin,
                            ymin =  ext$xmin, ymax =  ext$xmax))
  if (d == 180) return(list(xmin = -ext$xmax, xmax = -ext$xmin,
                            ymin = -ext$ymax, ymax = -ext$ymin))
  if (d == 270) return(list(xmin =  ext$ymin, xmax =  ext$ymax,
                            ymin = -ext$xmax, ymax = -ext$xmin))
  ext[c("xmin", "xmax", "ymin", "ymax")]
}

rotate_raster_view <- function(r, deg) {
  d <- ((deg %% 360) + 360) %% 360
  if (d == 0 || is.null(r)) return(r)
  rot1 <- function(m) {
    if (d == 90)  return(t(m)[ncol(m):1, , drop = FALSE])   # CCW
    if (d == 180) return(m[nrow(m):1, ncol(m):1, drop = FALSE])
    t(m[nrow(m):1, , drop = FALSE])                          # 270 = CW
  }
  if (length(dim(r)) == 2) return(rot1(r))
  ch <- lapply(seq_len(dim(r)[3]), function(k) rot1(r[, , k]))
  # unlist (not array(0, ...)) so an integer raster stays integer
  array(unlist(ch, use.names = FALSE),
        dim = c(nrow(ch[[1]]), ncol(ch[[1]]), length(ch)))
}

# ---------------------------------------------------------------------------
# View mirror -- reflect the scene about the X axis (y -> -y), about the origin
# like the rotations above.  A rotation alone cannot undo a handedness
# difference between two instrument exports, which is why a Y flip is needed
# on top of the four rotations.
#
# Only ONE mirror axis is offered: an X mirror is the same scene as a Y mirror
# followed by a 180 deg rotation, so {0,90,-90,180} x {no flip, flip Y} already
# spans all eight orientations (the dihedral group).
# ---------------------------------------------------------------------------
flip_xy_view <- function(x, y, flip) {
  if (!isTRUE(flip)) return(list(x = x, y = y))
  list(x = x, y = -y)
}

flip_extent_view <- function(ext, flip) {
  if (!isTRUE(flip)) return(ext[c("xmin", "xmax", "ymin", "ymax")])
  list(xmin = ext$xmin, xmax = ext$xmax, ymin = -ext$ymax, ymax = -ext$ymin)
}

flip_raster_view <- function(r, flip) {
  if (!isTRUE(flip) || is.null(r)) return(r)
  # Raster row 1 is the top edge, so reversing rows mirrors vertically.
  f1 <- function(m) m[nrow(m):1, , drop = FALSE]
  if (length(dim(r)) == 2) return(f1(r))
  ch <- lapply(seq_len(dim(r)[3]), function(k) f1(r[, , k]))
  array(unlist(ch, use.names = FALSE),
        dim = c(nrow(ch[[1]]), ncol(ch[[1]]), length(ch)))
}

# Composite view transform: MIRROR FIRST, THEN ROTATE.  All three helpers below
# use that same order, so points, extent and raster stay consistent.
view_transform_xy <- function(x, y, deg, flip = FALSE) {
  f <- flip_xy_view(x, y, flip)
  rotate_xy_view(f$x, f$y, deg)
}

view_transform_extent <- function(ext, deg, flip = FALSE) {
  rotate_extent_view(flip_extent_view(ext, flip), deg)
}

view_transform_raster <- function(r, deg, flip = FALSE) {
  rotate_raster_view(flip_raster_view(r, flip), deg)
}

# Auto LDIR view rotation: directly MEASURE which 90 deg rotation brings the
# LDIR particle cloud into the Raman particle cloud's orientation, by scoring
# each of {0, 90, -90, 180} on how many LDIR particles land on a Raman
# particle after best translation.  Operates on the exact native display
# coordinates (x_orig/y_orig) the viewer plots, so it is immune to the
# instrument export convention, Y-flips, and pipeline transform quirks that
# made the transform-file reconstruction (ldir_total_rotation_deg) unreliable.
# Scale-free: both clouds are centered and normalized to unit RMS radius
# first, so a scale mismatch (e.g. LDIR circle-calibration inflation) does
# not affect the rotation choice.  Returns an integer in {0,90,-90,180};
# falls back to 0 when no rotation clearly beats leaving it unrotated.
#
# .auto_view_best() is the shared engine: it scores a list of candidate view
# transforms (each list(deg=, flip=)) and returns the winning index.  The FIRST
# candidate must be the identity -- it is the fallback whenever nothing scores
# well enough, or when no candidate clearly beats leaving the view as-is.
.auto_view_best <- function(src_x, src_y, ref_x, ref_y, cands) {
  fk <- is.finite(src_x) & is.finite(src_y)
  fr <- is.finite(ref_x) & is.finite(ref_y)
  lx <- src_x[fk]; ly <- src_y[fk]
  rx <- ref_x[fr]; ry <- ref_y[fr]
  if (length(lx) < 4 || length(rx) < 4) return(1L)

  # Scoring is O(n_src * n_ref) per candidate. LDIR clouds are tiny, but an
  # FTIR run can carry thousands of particles -- thin deterministically (no RNG,
  # so the result stays reproducible and cacheable) to keep the tab responsive.
  cap <- 400L
  thin <- function(v, n) if (n <= cap) v else v[round(seq(1, n, length.out = cap))]
  nl <- length(lx); nr <- length(rx)
  lx <- thin(lx, nl); ly <- thin(ly, nl)
  rx <- thin(rx, nr); ry <- thin(ry, nr)

  nrm <- function(x, y) {
    x <- x - mean(x); y <- y - mean(y)
    s <- sqrt(mean(x^2 + y^2)); if (!is.finite(s) || s <= 0) s <- 1
    list(x = x / s, y = y / s)
  }
  L <- nrm(lx, ly); R <- nrm(rx, ry)
  tol <- 0.10   # normalized units (~10% of cloud radius)

  score <- function(cand) {
    p <- view_transform_xy(L$x, L$y, cand$deg, cand$flip)
    dx <- outer(R$x, p$x, "-"); dy <- outer(R$y, p$y, "-")
    # translation voting: densest bin of pairwise offsets, then inlier count
    key <- paste(round(dx / tol), round(dy / tol))
    tb  <- sort(table(key), decreasing = TRUE)
    best <- 0L
    for (k in names(tb)[seq_len(min(6, length(tb)))]) {
      sel <- key == k
      px <- p$x + mean(dx[sel]); py <- p$y + mean(dy[sel])
      d <- sqrt(outer(R$x, px, "-")^2 + outer(R$y, py, "-")^2)
      best <- max(best, sum(apply(d, 2, min) <= tol))
    }
    best
  }

  ns <- vapply(cands, score, integer(1))
  bi <- which.max(ns)
  if (ns[bi] < 4) return(1L)
  # only transform when the winner clearly beats leaving the view as-is
  if (bi != 1L && ns[bi] <= ns[1] + 1L) return(1L)
  bi
}

ldir_auto_view_rotation <- function(ldir_x, ldir_y, raman_x, raman_y) {
  degs  <- c(0L, 90L, -90L, 180L)
  cands <- lapply(degs, function(d) list(deg = d, flip = FALSE))
  degs[.auto_view_best(ldir_x, ldir_y, raman_x, raman_y, cands)]
}

# Auto view transform including a possible mirror -- same measurement as
# ldir_auto_view_rotation() but over all eight orientations, so it can tell a
# 180 deg rotation apart from a Y flip (which look identical for a symmetric
# particle cloud but are different scenes).  Used by the FTIR tabs, where the
# instrument export can differ from Raman in handedness, not just rotation.
# Returns list(deg = <0|90|-90|180>, flip = <TRUE|FALSE>).
auto_view_dihedral <- function(src_x, src_y, ref_x, ref_y) {
  cands <- list()
  for (fl in c(FALSE, TRUE))
    for (d in c(0L, 90L, -90L, 180L))
      cands[[length(cands) + 1L]] <- list(deg = d, flip = fl)
  cands[[.auto_view_best(src_x, src_y, ref_x, ref_y, cands)]]
}

# ---------------------------------------------------------------------------
# Auto view orientation measured from PAIRED coordinates.
#
# Every instrument row already carries both frames: x_orig/y_orig (native --
# what the single-instrument tab plots) and x/y (aligned into Raman space --
# what the overlay plots). That is a per-particle correspondence, so the
# native -> Raman orientation can be measured exactly with a 2D Kabsch fit,
# with no point matching at all.
#
# This is strictly better than cloud matching (auto_view_dihedral): it cannot
# be defeated by the two instruments detecting different particles, by very
# different particle counts, or by a near-symmetric layout where several
# rotations score alike -- the failure mode that left the view unrotated even
# though aligned mode was demonstrably correct.
#
# The fitted angle is snapped to a multiple of 90, since that is all a view
# rotation offers. Mirror-first ordering matches view_transform_xy().
# Returns list(deg, flip, angle, rmse) -- rmse is the residual as a fraction of
# the cloud radius, so the caller can reject a fit that did not converge --
# or NULL when there is too little paired data to decide.
# ---------------------------------------------------------------------------
auto_view_from_pairs <- function(x_native, y_native, x_aligned, y_aligned) {
  if (is.null(x_native) || is.null(x_aligned)) return(NULL)
  ok <- is.finite(x_native) & is.finite(y_native) &
        is.finite(x_aligned) & is.finite(y_aligned)
  if (sum(ok) < 3L) return(NULL)

  ax <- x_native[ok];  ay <- y_native[ok]
  bx <- x_aligned[ok]; by <- y_aligned[ok]
  ax <- ax - mean(ax); ay <- ay - mean(ay)
  bx <- bx - mean(bx); by <- by - mean(by)
  sa <- sum(ax^2 + ay^2); sb <- sum(bx^2 + by^2)
  if (!is.finite(sa) || !is.finite(sb) || sa <= 0 || sb <= 0) return(NULL)

  rms_b <- sqrt(sb / length(bx))
  # Optimal rotation for one handedness: theta = atan2(Sxy - Syx, Sxx + Syy).
  fit <- function(px, py) {
    th <- atan2(sum(px * by) - sum(py * bx), sum(px * bx) + sum(py * by))
    s  <- sqrt(sb / sa)                      # isotropic scale between frames
    qx <- s * (px * cos(th) - py * sin(th))
    qy <- s * (px * sin(th) + py * cos(th))
    list(theta = th, rmse = sqrt(mean((qx - bx)^2 + (qy - by)^2)) / rms_b)
  }
  f_plain <- fit(ax,  ay)
  f_flip  <- fit(ax, -ay)                    # mirror first, then rotate
  best <- if (f_flip$rmse < f_plain$rmse) list(f = f_flip, flip = TRUE)
          else                             list(f = f_plain, flip = FALSE)

  deg <- ((round(best$f$theta * 180 / pi / 90) * 90) %% 360 + 360) %% 360
  if (deg == 270) deg <- -90
  list(deg = as.integer(deg), flip = best$flip,
       angle = best$f$theta * 180 / pi, rmse = best$f$rmse)
}

# Total LDIR -> Raman rotation for a run, snapped to the nearest 90 deg --
# reconstructed from the pipeline's residual rotation
# (04_alignment/transform_params_ldir_raman.txt) plus the pre-rotation
# recorded in the manifest snapshot (ldir_rotate_deg_for_alignment; -90
# historical default).  Kept as a secondary reference; the viewer's Auto
# mode uses ldir_auto_view_rotation() (direct measurement) instead, which
# does not depend on parsing pipeline internals.  Returns NULL when the
# transform is unavailable or reflected (a mirror is not a pure rotation).
ldir_total_rotation_deg <- function(run_dir, manifest = NULL) {
  if (is.null(run_dir)) return(NULL)
  tp <- file.path(run_dir, "04_alignment", "transform_params_ldir_raman.txt")
  if (!file.exists(tp)) return(NULL)
  ln <- tryCatch(readLines(tp), error = function(e) NULL)
  if (is.null(ln)) return(NULL)
  gv <- function(p) suppressWarnings(
    as.numeric(sub(".*:\\s*", "", grep(p, ln, value = TRUE)[1])))
  resid <- gv("^rotation_deg")
  refl_line <- grep("^reflected", ln, value = TRUE)
  refl <- length(refl_line) > 0 && grepl("TRUE", refl_line[1])
  if (!is.finite(resid) || isTRUE(refl)) return(NULL)
  pre <- tryCatch(manifest$config_snapshot$ldir_rotate_deg_for_alignment,
                  error = function(e) NULL)
  if (is.null(pre) || !is.numeric(pre)) pre <- -90
  total   <- ((resid + pre + 180) %% 360) - 180
  snapped <- (round(total / 90) * 90) %% 360
  if (snapped > 180) snapped <- snapped - 360
  as.integer(snapped)
}

# ---------------------------------------------------------------------------
# Compute image bounds that preserve the image's native aspect ratio while
# centering on a set of particles.  The image is expanded (never cropped)
# so that all particles fit inside, plus padding.
#
# img_raster : raster array (h x w x channels)
# x_vals, y_vals : particle coordinate vectors (determines center + extent)
# padding_um : padding in um on each side
#
# Returns: list(xmin, xmax, ymin, ymax) for annotation_raster
# ---------------------------------------------------------------------------
compute_image_bounds <- function(img_raster, x_vals, y_vals, padding_um = 300) {
  img_h <- nrow(img_raster)
  img_w <- ncol(img_raster)
  img_aspect <- img_w / img_h  # width / height

  # Particle extent
  px_min <- min(x_vals, na.rm = TRUE)
  px_max <- max(x_vals, na.rm = TRUE)
  py_min <- min(y_vals, na.rm = TRUE)
  py_max <- max(y_vals, na.rm = TRUE)

  # Center on particles
  cx <- (px_min + px_max) / 2
  cy <- (py_min + py_max) / 2

  # Required span to contain all particles + padding
  span_x <- (px_max - px_min) + 2 * padding_um
  span_y <- (py_max - py_min) + 2 * padding_um

  # Expand the shorter dimension to match the image aspect ratio
  if (span_x / span_y > img_aspect) {
    # Image needs to be taller for this width
    span_y <- span_x / img_aspect
  } else {
    # Image needs to be wider for this height
    span_x <- span_y * img_aspect
  }

  list(xmin = cx - span_x / 2, xmax = cx + span_x / 2,
       ymin = cy - span_y / 2, ymax = cy + span_y / 2)
}

# ---------------------------------------------------------------------------
# Natural/numeric sort for particle IDs with prefixes.
# Extracts the first numeric substring and sorts by it, keeping the full ID.
# E.g. MP_1, MP_2, ..., MP_10 (not MP_1, MP_10, MP_11, ..., MP_2).
# ---------------------------------------------------------------------------
natural_sort_ids <- function(ids) {
  nums <- as.numeric(sub(".*?(\\d+).*", "\\1", ids))
  ids[order(nums, na.last = TRUE)]
}

# ---------------------------------------------------------------------------
# Parse a particle selection string into matching IDs.
#
# Supported syntax:
#   - Numeric range: "1-10" or "MP_1-MP_10" (matches IDs with numeric suffix 1..10)
#   - Comma list:    "MP_1,MP_5,MP_10"  (exact IDs)
#   - Glob pattern:  "MP_*" or "A*"     (wildcard matching)
#
# all_ids: character vector of all available particle IDs
# Returns: character vector of matching IDs (subset of all_ids)
# ---------------------------------------------------------------------------
parse_particle_selection <- function(text, all_ids) {
  text <- trimws(text)
  if (nchar(text) == 0) return(character(0))

  # Try numeric range: "1-10" or "MP_1-MP_10"
  range_match <- regmatches(text, regexec("^(.*)?(\\d+)\\s*-\\s*(.*)?(\\d+)$", text))[[1]]
  if (length(range_match) == 5) {
    lo <- as.integer(range_match[3])
    hi <- as.integer(range_match[5])
    if (!is.na(lo) && !is.na(hi)) {
      nums <- as.numeric(sub(".*?(\\d+).*", "\\1", all_ids))
      return(all_ids[!is.na(nums) & nums >= lo & nums <= hi])
    }
  }

  # Try comma-separated list
  if (grepl(",", text)) {
    parts <- trimws(strsplit(text, ",")[[1]])
    return(intersect(parts, all_ids))
  }

  # Try glob pattern (contains * or ?)
  if (grepl("[*?]", text)) {
    pat <- utils::glob2rx(text)
    return(all_ids[grepl(pat, all_ids)])
  }

  # Exact match
  if (text %in% all_ids) return(text)

  character(0)
}


# ===========================================================================
# PDF report content layer
# ===========================================================================
# Pure builders shared by the Summary tab's HTML tables and the PDF report, so
# the two can never drift apart. Nothing here touches Shiny reactives or
# devices: each function takes plain data and returns a data frame or a grob.
#
# The report is assembled with grDevices::pdf(), matching the pipeline's own
# plots/all_diagnostics.pdf. That keeps it dependency-free -- no pandoc, no
# LaTeX -- so it works wherever the app runs.

# Family x device count table behind "Plastics by Instrument".
#
# `devices` is a named list of particle data frames (name = device label);
# entries that are NULL or empty are dropped. Rows are ordered by material
# category then family name, and a Total row is appended. Returns a zero-row
# data frame when nothing is classifiable, so callers can test nrow() alone.
report_plastics_table <- function(devices) {
  devices <- Filter(function(d) !is.null(d) && nrow(d) > 0, devices)
  empty <- data.frame(Family = character(0), Category = character(0),
                      stringsAsFactors = FALSE)
  if (length(devices) == 0) return(empty)

  per_dev  <- lapply(devices, summarise_plastics)
  all_fams <- unique(unlist(lapply(per_dev, `[[`, "family")))
  if (length(all_fams) == 0) return(empty)

  fam_cats  <- classify_category_vec(all_fams)
  cat_order <- c("Synthetic", "Semi-synthetic", "Natural/Organic", "Unknown")
  fam_ord   <- order(match(fam_cats, cat_order, nomatch = 99), all_fams)
  all_fams  <- all_fams[fam_ord]
  fam_cats  <- fam_cats[fam_ord]

  out <- data.frame(Family = all_fams, Category = fam_cats,
                    stringsAsFactors = FALSE)
  for (dev in names(per_dev)) {
    dt <- per_dev[[dev]]
    idx <- match(all_fams, dt$family)
    out[[dev]] <- ifelse(is.na(idx), 0L, dt$n[idx])
  }
  # Total row: per-device sum over the families shown (same as the HTML footer).
  totals <- data.frame(Family = "Total", Category = "", stringsAsFactors = FALSE)
  for (dev in names(per_dev)) totals[[dev]] <- sum(per_dev[[dev]]$n)
  rbind(out, totals)
}


# Per-instrument size statistics behind "Size Statistics".
#
# `dfs` is a named list with $ftir / $raman / $ldir (as summary_dfs() returns).
# Numbers are pre-formatted for display, matching the HTML table's rounding.
report_size_stats_table <- function(dfs) {
  labels <- c(FTIR = "ftir", Raman = "raman", LDIR = "ldir")
  rows <- list()
  for (nm in names(labels)) {
    df <- dfs[[labels[[nm]]]]
    if (is.null(df) || nrow(df) == 0) next
    fm <- df$feret_max
    rows[[nm]] <- data.frame(
      Instrument = nm,
      Total      = nrow(df),
      Matched    = sum(df$match_status == "matched",   na.rm = TRUE),
      Unmatched  = sum(df$match_status == "unmatched", na.rm = TRUE),
      Mean       = paste0(round(mean(fm,   na.rm = TRUE), 1), " \u00b5m"),
      Median     = paste0(round(stats::median(fm, na.rm = TRUE), 1), " \u00b5m"),
      `Std Dev`  = paste0(round(stats::sd(fm, na.rm = TRUE), 1), " \u00b5m"),
      Range      = paste0(round(min(fm, na.rm = TRUE), 1), "\u2013",
                          round(max(fm, na.rm = TRUE), 1), " \u00b5m"),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0)
    return(data.frame(Instrument = character(0), stringsAsFactors = FALSE))
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}


# --- PDF page primitives ---------------------------------------------------
# Each returns a ggplot/grob that can be printed to one page of the report.

# Wrap long text so it fits the page width instead of running off it.
# Lines already within `width` are passed through VERBATIM: strwrap() collapses
# runs of spaces, which would destroy the column alignment of the monospaced
# "label : value" blocks on the title page.
.report_wrap <- function(txt, width = 100) {
  if (length(txt) == 0) return(character(0))
  unlist(lapply(txt, function(s) {
    if (is.na(s) || !nzchar(s)) return("")
    if (nchar(s) <= width) return(s)
    strwrap(s, width = width)
  }), use.names = FALSE)
}

# A text-only page: bold title, then left-aligned body lines.
# ---------------------------------------------------------------------------
# Input-file provenance for the report title page.
#
# Reads the run manifest's `inputs` block, which records one entry per source
# file (basename, format, size, and pixel dimensions for images). Shared by the
# pipeline report and the viewer's download so both list the same thing.
#
# Names are middle-truncated rather than wrapped: report_text_page() wraps with
# strwrap(), which breaks on whitespace and leaves a long unbroken filename to
# run off the page edge. Truncating in the middle keeps the extension visible,
# which is the part that identifies the format.
# ---------------------------------------------------------------------------
.report_trunc_mid <- function(s, width) {
  s <- as.character(s)
  n <- nchar(s)
  if (is.na(n) || n <= width) return(s)
  keep  <- width - 3L
  left  <- ceiling(keep / 2)
  right <- keep - left
  paste0(substr(s, 1, left), "...", substr(s, n - right + 1L, n))
}

.report_fmt_bytes <- function(b) {
  b <- suppressWarnings(as.numeric(b))
  if (length(b) != 1 || is.na(b)) return("")
  if (b >= 1048576) return(sprintf("%.1f MB", b / 1048576))
  if (b >= 1024)    return(sprintf("%.1f KB", b / 1024))
  sprintf("%d B", as.integer(b))
}

report_input_file_lines <- function(manifest, name_width = 46L) {
  inp <- manifest$inputs
  if (is.null(inp) || length(inp) == 0)
    return(c("Input files:", "  (not recorded in this run's manifest)"))

  friendly <- c(ftir = "FTIR (PerkinElmer)", ftir_perkin = "FTIR (PerkinElmer)",
                ftir_bruker = "FTIR (Bruker)", raman = "Raman", ldir = "LDIR")
  rows <- list()
  for (k in names(inp)) {
    e <- inp[[k]]
    if (!is.list(e)) next
    base <- e$basename %||% basename(e$path %||% "")
    if (!nzchar(base)) next
    is_img <- grepl("_image$", k) ||
      toupper(e$format %||% "") %in% c("PNG", "BMP", "TIF", "TIFF", "JPG", "JPEG")
    inst <- sub("_image$", "", k)
    dims <- if (!is.null(e$width) && !is.null(e$height) &&
                is.finite(suppressWarnings(as.numeric(e$width))))
      sprintf("%dx%d px", as.integer(e$width), as.integer(e$height)) else ""
    rows[[length(rows) + 1]] <- list(
      img = is_img,
      label = friendly[[inst]] %||% inst,
      base = .report_trunc_mid(base, name_width),
      fmt = e$format %||% "",
      size = .report_fmt_bytes(e$size_bytes),
      dims = dims)
  }
  if (length(rows) == 0)
    return(c("Input files:", "  (not recorded in this run's manifest)"))

  fmt_row <- function(r) sprintf("    %-19s %-*s %-5s %9s %s",
                                 r$label, name_width, r$base, r$fmt, r$size, r$dims)
  out <- "Input files:"
  data_rows <- Filter(function(r) !r$img, rows)
  img_rows  <- Filter(function(r)  r$img, rows)
  if (length(data_rows)) out <- c(out, "  Data:",
                                  vapply(data_rows, fmt_row, character(1)))
  if (length(img_rows))  out <- c(out, "  Images:",
                                  vapply(img_rows, fmt_row, character(1)))
  trimws(out, which = "right")
}

report_text_page <- function(title, lines = character(0), subtitle = NULL) {
  lines <- .report_wrap(lines)
  body  <- if (length(lines)) paste(lines, collapse = "\n") else ""
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1, label = title, hjust = 0, vjust = 1,
                      size = 7, fontface = "bold") +
    { if (!is.null(subtitle))
        ggplot2::annotate("text", x = 0, y = 0.93,
                          label = paste(.report_wrap(subtitle), collapse = "\n"),
                          hjust = 0, vjust = 1, size = 4.2, colour = "grey35")
      else NULL } +
    ggplot2::annotate("text", x = 0, y = 0.86, label = body, hjust = 0, vjust = 1,
                      size = 4, family = "mono") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(24, 24, 24, 24),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# A data frame rendered as a table page, anchored top-left under the title.
# Falls back to a text page when gridExtra is unavailable, so the report
# degrades rather than failing.
report_table_page <- function(df, title, caption = NULL) {
  if (is.null(df) || nrow(df) == 0)
    return(report_text_page(title, "No data available.", caption))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    txt <- utils::capture.output(print(df, row.names = FALSE))
    return(report_text_page(title, txt, caption))
  }
  chr <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE,
                       check.names = FALSE)
  nr <- nrow(chr); nc <- ncol(chr)
  # Shrink the font as the table grows so wide/long tables stay on one page.
  fs <- if (nr > 24 || nc > 9) 7 else if (nr > 14) 8.5 else 10

  # Bold a trailing "Total" row, mirroring the HTML table's bold footer.
  total_row <- if (nr > 0 && identical(as.character(chr[[1]][nr]), "Total")) nr else NA_integer_
  faces <- matrix(1L, nrow = nr, ncol = nc)
  if (!is.na(total_row)) faces[total_row, ] <- 2L

  tg <- gridExtra::tableGrob(
    chr, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = fs,
      core    = list(fg_params = list(hjust = 0, x = 0.03,
                                      fontface = as.vector(faces))),
      colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.03))))

  # Rule above the Total row (the HTML uses a 2px top border).
  if (!is.na(total_row) && requireNamespace("gtable", quietly = TRUE)) {
    tg <- gtable::gtable_add_grob(
      tg, grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1,
                             gp = grid::gpar(lwd = 1.6, col = "grey35")),
      t = total_row + 1L, b = total_row + 1L, l = 1L, r = nc, name = "total-rule")
  }
  # Rule under the header, so the column names read as a header.
  if (requireNamespace("gtable", quietly = TRUE)) {
    tg <- gtable::gtable_add_grob(
      tg, grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                             gp = grid::gpar(lwd = 1, col = "grey60")),
      t = 1L, b = 1L, l = 1L, r = nc, name = "head-rule")
  }

  # Pin the table to the top-left: pad to the right and below with null space
  # instead of letting arrangeGrob centre it in the page. The caption sits
  # directly under the table, where it still reads as belonging to it.
  row <- gridExtra::arrangeGrob(
    tg, grid::nullGrob(), ncol = 2,
    widths = grid::unit.c(sum(tg$widths), grid::unit(1, "null")))
  cap_g <- if (!is.null(caption) && any(nzchar(caption)))
    grid::textGrob(paste(.report_wrap(caption, 120), collapse = "\n"),
                   x = 0, hjust = 0, vjust = 1,
                   gp = grid::gpar(fontsize = 9, col = "grey35")) else NULL
  parts   <- list(row)
  heights <- list(sum(tg$heights))
  if (!is.null(cap_g)) {
    parts   <- c(parts, list(cap_g))
    heights <- c(heights, list(grid::unit(2.2, "lines")))
  }
  parts   <- c(parts, list(grid::nullGrob()))
  heights <- c(heights, list(grid::unit(1, "null")))
  gridExtra::arrangeGrob(
    grobs = parts, ncol = 1,
    heights = do.call(grid::unit.c, heights),
    top = grid::textGrob(title, x = 0, hjust = 0,
                         gp = grid::gpar(fontsize = 16, fontface = "bold")),
    # Page margins -- without these the title sits flush against the paper edge.
    vp = grid::viewport(width  = grid::unit(1, "npc") - grid::unit(1, "cm"),
                        height = grid::unit(1, "npc") - grid::unit(1, "cm")))
}

# A figure page: an existing ggplot plus a caption recording the filters that
# produced it, so a page lifted out of the PDF still says what it is showing.
report_figure_page <- function(plot, title, caption = NULL) {
  if (is.null(plot)) return(NULL)
  cap <- if (!is.null(caption) && length(caption) && any(nzchar(caption)))
    paste(.report_wrap(caption, 140), collapse = "\n") else NULL
  # Demote the plot's OWN title to a subtitle rather than discarding it: the
  # viewer puts the particle count and the coordinate frame there, and on an
  # empty view it carries the "no particles match the current filters" message.
  # Overwriting it left an unexplained blank figure in the report.
  own <- plot$labels$title
  sub <- if (!is.null(own) && length(own) == 1L && !is.na(own) && nzchar(own) &&
             !identical(own, title)) own else NULL
  plot +
    ggplot2::labs(title = title, subtitle = sub, caption = cap) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10.5, colour = "grey25"),
      plot.caption  = ggplot2::element_text(size = 8.5, colour = "grey35",
                                            hjust = 0),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# A full-bleed figure page: title and optional caption are grobs outside the
# plot area so the ggplot fills the remaining space without title/margin waste.
# Use this for the barplot and any chart that is clipped by report_figure_page.
report_full_figure_page <- function(plot, title, caption = NULL) {
  if (is.null(plot)) return(NULL)
  # Strip the plot's own title so we don't double-render it
  plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.margin = ggplot2::margin(4, 6, 4, 6))

  title_g <- grid::textGrob(title, x = 0, hjust = 0,
                             gp = grid::gpar(fontsize = 15, fontface = "bold"))
  cap_g   <- if (!is.null(caption) && any(nzchar(caption)))
    grid::textGrob(paste(.report_wrap(caption, 140), collapse = "\n"),
                   x = 0, hjust = 0,
                   gp = grid::gpar(fontsize = 8.5, col = "grey35")) else NULL

  parts   <- list(title_g, plot)
  heights <- list(grid::unit(1.6, "lines"), grid::unit(1, "null"))
  if (!is.null(cap_g)) {
    parts   <- c(parts, list(cap_g))
    heights <- c(heights, list(grid::unit(1.8, "lines")))
  }
  gridExtra::arrangeGrob(
    grobs = parts, ncol = 1,
    heights = do.call(grid::unit.c, heights),
    vp = grid::viewport(width  = grid::unit(1, "npc") - grid::unit(1.2, "cm"),
                        height = grid::unit(1, "npc") - grid::unit(1.2, "cm")))
}

# Arrange several plots on one page (used for the pies and size histograms).
report_grid_page <- function(plots, title, caption = NULL, ncol = 2) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(NULL)
  if (!requireNamespace("gridExtra", quietly = TRUE)) return(plots[[1]])
  top <- grid::textGrob(title, x = 0.02, hjust = 0,
                        gp = grid::gpar(fontsize = 15, fontface = "bold"))
  bottom <- if (!is.null(caption) && any(nzchar(caption)))
    grid::textGrob(paste(.report_wrap(caption, 140), collapse = "\n"),
                   x = 0.02, hjust = 0,
                   gp = grid::gpar(fontsize = 8.5, col = "grey35")) else NULL
  gridExtra::arrangeGrob(
    grobs = plots, ncol = min(ncol, length(plots)),
    top = top, bottom = bottom,
    # Page margins, so the title is not flush against the paper edge.
    vp = grid::viewport(width  = grid::unit(1, "npc") - grid::unit(1, "cm"),
                        height = grid::unit(1, "npc") - grid::unit(1, "cm")))
}



# ===========================================================================
# Shared plot renderers (hoisted out of app.R's server scope)
# ===========================================================================
# These were defined inside server(), which made them unreachable from
# main.R -- so the pipeline's report could only ever contain tables while the
# viewer's download contained figures. They reference no reactives and no
# input$, so they live here unchanged and BOTH callers now use them. That is
# the point: one renderer, so the two reports cannot drift apart.
# ---------------------------------------------------------------------------

# Degenerate-extent guard, used by sanitize_bounds() below. Lives at top
# level so both the viewer and the pipeline report can reach it.
.MIN_SPAN <- 1e-6

plot_size_distribution <- function(df, inst_name, color_matched = "#d62728", color_unmatched = "#bcbd22") {

  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No data"),
                                 size = 5, colour = "grey50") +
           theme_void())
  }
  # Show all particles in one histogram with smoothed density overlay (normalized to %)
  ggplot(df, aes(x = feret_max)) +
    geom_histogram(aes(y = after_stat(density) * 100), alpha = 0.6, bins = 20, fill = "#4472C4", color = "white") +
    geom_density(aes(y = after_stat(density) * 100), alpha = 0.5, fill = "#70AD47", color = "#70AD47", linewidth = 1.2) +
    labs(title = inst_name, x = "Feret Max (\u00b5m)", y = "Relative Frequency (%)") +
    theme_minimal() + theme(plot.title = element_text(size = 11, face = "bold"))
}

add_image_bg <- function(p, img_info, alpha = 0.4) {
  if (is.null(img_info)) return(p)
  p + annotation_raster(img_info$raster,
        xmin = img_info$xmin, xmax = img_info$xmax,
        ymin = img_info$ymin, ymax = img_info$ymax,
        interpolate = TRUE)
}

breaks_adaptive <- function(rng) {
  if (is.null(rng) || length(rng) < 2 || rng[1] >= rng[2]) return(NULL)

  span <- rng[2] - rng[1]

  # Choose interval to get ~4-8 breaks (target: 6)
  intervals <- c(1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)
  best_int <- 1000
  for (int in intervals) {
    n_breaks <- span / int
    if (n_breaks >= 4 && n_breaks <= 8) {
      best_int <- int
      break
    }
    if (n_breaks < 4) {
      best_int <- int
      break
    }
  }

  seq(floor(rng[1] / best_int) * best_int, ceiling(rng[2] / best_int) * best_int, by = best_int)
}

safe_size_limits <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NULL)
  r <- range(v)
  if (r[1] == r[2]) r <- c(0, r[2] + 1)
  r
}

add_particle_labels <- function(p, df, bounds, id_col = "particle_id",
                                size = 3) {
  if (is.null(df) || nrow(df) == 0 || !(id_col %in% names(df))) return(p)
  lab <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (nrow(lab) == 0) return(p)
  lab$.lbl <- as.character(lab[[id_col]])
  lab <- lab[!is.na(lab$.lbl) & nzchar(lab$.lbl), , drop = FALSE]
  if (nrow(lab) == 0) return(p)
  size <- if (is.numeric(size) && length(size) == 1 && is.finite(size)) size else 3
  xr <- diff(bounds$x); yr <- diff(bounds$y)
  # Nudge/shadow scale gently with text size so bigger labels stay clear of
  # the marker and keep a proportional shadow.
  lab$.ly  <- lab$y + yr * 0.006 * size     # nudge above the marker
  lab$.sx  <- lab$x + xr * 0.0006 * size    # shadow offset
  lab$.sy  <- lab$.ly - yr * 0.0006 * size
  p +
    geom_text(data = lab, aes(x = .sx, y = .sy, label = .lbl),
              vjust = 0, size = size, colour = "black", alpha = 0.85,
              inherit.aes = FALSE) +
    geom_text(data = lab, aes(x = x, y = .ly, label = .lbl),
              vjust = 0, size = size, colour = "white", fontface = "bold",
              inherit.aes = FALSE)
}

make_scatter <- function(df, img_info, bounds, title,
                          match_colours = NULL, highlight_id = NULL,
                          full_df = NULL, match_labels = NULL,
                          plain = FALSE, show_labels = FALSE, label_size = 3,
                          subtitle = NULL) {

  p <- ggplot(df, aes(x = x, y = y))

  # Background image (with per-image bounds)
  p <- add_image_bg(p, img_info)

  # Points. In "plain" mode (Show all detected) every particle is drawn in a
  # single colour with no matched/unmatched distinction or legend.
  if (isTRUE(plain)) {
    p <- p + geom_point(aes(size = feret_max), colour = "#1f77b4",
                        alpha = 0.7)
  } else {
    p <- p + geom_point(aes(colour = match_status, size = feret_max),
                        alpha = 0.7)
    if (!is.null(match_colours)) {
      if (!is.null(match_labels))
        p <- p + scale_colour_manual(values = match_colours, labels = match_labels)
      else
        p <- p + scale_colour_manual(values = match_colours)
    }
  }

  p <- p +
    scale_size_continuous(name = "Feret Max (\u00b5m)", range = c(2, 12),
                          limits = safe_size_limits(df$feret_max)) +
    scale_x_continuous(breaks = breaks_adaptive(bounds$x)) +
    scale_y_continuous(breaks = breaks_adaptive(bounds$y)) +
    coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = "X (\u00b5m)", y = "Y (\u00b5m)") +
    theme_minimal(base_size = 15) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "grey98", colour = NA),
      panel.grid       = element_line(colour = "grey90"),
      plot.subtitle    = element_text(size = 11.5, colour = "#b02a37"),
      legend.position  = "right",
      legend.title     = element_text(size = 13),
      legend.text      = element_text(size = 11)
    )

  # Highlight selected particle(s) -- ALWAYS shown even if filtered out.
  # First try the filtered df, then fall back to full_df (unfiltered).
  # highlight_id can be a character vector (multiple IDs from pattern select)
  # or a single ID (from the selectInput dropdown).
  if (!is.null(highlight_id) && length(highlight_id) > 0 &&
      !identical(highlight_id, "None") && !identical(highlight_id, character(0))) {
    hl <- NULL
    if ("particle_id" %in% names(df))
      hl <- df[df$particle_id %in% highlight_id, ]
    if ((is.null(hl) || nrow(hl) == 0) && !is.null(full_df) &&
        "particle_id" %in% names(full_df))
      hl <- full_df[full_df$particle_id %in% highlight_id, ]
    if (!is.null(hl) && nrow(hl) > 0) {
      # Y-offset scales with plot extent so label doesn't overlap the circle
      y_span <- diff(bounds$y)
      y_nudge <- y_span * 0.03   # 3% of visible y-range
      hl$label_y <- hl$y + y_nudge
      p <- p + geom_point(data = hl, aes(x = x, y = y),
                           shape = 21, size = 10, stroke = 2,
                           fill = NA, colour = "#FFD700") +
               geom_text(data = hl, aes(x = x, y = label_y, label = particle_id),
                          vjust = 0, size = 4.0, fontface = "bold",
                          colour = "#FFD700")
    }
  }

  # Number every displayed particle (opt-in)
  if (isTRUE(show_labels)) p <- add_particle_labels(p, df, bounds, size = label_size)

  p
}

sanitize_bounds <- function(b, fallback = list(x = c(-1000, 1000),
                                               y = c(-1000, 1000))) {
  ok <- function(v) is.numeric(v) && length(v) == 2L && all(is.finite(v))
  if (is.null(b) || !ok(b$x) || !ok(b$y)) return(fallback)
  widen <- function(v) {
    v <- sort(v)
    if (diff(v) > .MIN_SPAN) v else c(mean(v) - 0.5, mean(v) + 0.5)
  }
  list(x = widen(b$x), y = widen(b$y))
}

# --- Report assembly ------------------------------------------------------

# Write `pages` (ggplots / grobs, NULLs skipped) to a multi-page PDF at `path`.
# Returns the number of pages written. A page that fails to draw is replaced by
# an error page rather than aborting the whole report -- a single bad figure
# should not cost the user the other twenty.
write_report_pdf <- function(pages, path, width = 11, height = 8.5) {
  pages <- Filter(Negate(is.null), pages)
  if (length(pages) == 0) pages <- list(report_text_page(
    "Report", "Nothing to report \u2014 no data is loaded in the viewer."))

  # cairo_pdf handles UTF-8 text properly regardless of the R session's locale.
  # The report is full of \u00b5m and en-dashes, and the base pdf() device encodes
  # text using the locale's charset -- under a C locale that mangles them (um
  # became garbage, en-dashes became "..."). Fall back to pdf() with a Latin-1
  # encoding, which still covers u, if cairo is not compiled in.
  if (isTRUE(unname(capabilities("cairo")))) {
    grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::pdf(path, width = width, height = height, onefile = TRUE,
                   encoding = "ISOLatin1.enc")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  n <- 0L
  for (pg in pages) {
    ok <- tryCatch({
      if (inherits(pg, "ggplot")) print(pg) else { grid::grid.newpage(); grid::grid.draw(pg) }
      TRUE
    }, error = function(e) { message("[report] page failed: ", conditionMessage(e)); FALSE })
    if (!ok) {
      tryCatch(print(report_text_page("Figure unavailable",
        "This page could not be rendered. The rest of the report is unaffected.")),
        error = function(e) NULL)
    }
    n <- n + 1L
  }
  n
}

# ---------------------------------------------------------------------------
# Standalone plotly barplot builder (non-reactive; used by write_report_html
# and by main.R so the interactive chart is reproducible outside Shiny).
#
# device_counts : named list, device_label -> named integer table of families
#                 (NULL entries are skipped). Mirrors active_material_counts().
# sel_fam       : "__all_plastics__" or a single family name.
# rel_mode      : if TRUE, show share (%) instead of counts.
# show_all      : if TRUE, include non-plastic families; otherwise only plastics.
# cat_mode      : "synthetic", "semi", or "both" (only used when !show_all).
# ---------------------------------------------------------------------------
build_plotly_barplot <- function(device_counts,
                                 sel_fam  = "__all_plastics__",
                                 rel_mode = FALSE,
                                 show_all = TRUE,
                                 cat_mode = "both") {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("plotly is required for build_plotly_barplot()")

  device_colors <- c("FTIR (PerkinElmer)" = "#2ca02c", "FTIR (Bruker)" = "#9467bd",
                     "Raman" = "#1f77b4", "LDIR" = "#d62728")
  .fam_palette  <- c(
    PE = "#e41a1c", PP = "#377eb8", PS = "#4daf4a", PET = "#984ea3",
    PVC = "#ff7f00", PA = "#a65628", PU = "#f781bf", PC = "#999999",
    PMMA = "#66c2a5", PTFE = "#fc8d62", ABS = "#e78ac3", Rubber = "#7570b3",
    Cellulose = "#bcbd22", Acrylate = "#17becf", Other = "#e5c494"
  )

  cts <- device_counts
  has_data    <- vapply(names(cts), function(l) !is.null(cts[[l]]), logical(1))
  cts_present <- cts[has_data]
  inst_levels <- names(cts_present)

  if (length(cts_present) == 0)
    return(plotly::plot_ly() |>
             plotly::layout(title = "No data available",
                            xaxis = list(visible = FALSE),
                            yaxis = list(visible = FALSE)))

  plastic_fams       <- c(synthetic_families, semi_synthetic_families)
  if (sel_fam == "__all_plastics__") {
    keep_fams <- if (show_all) NULL
                 else if (identical(cat_mode, "synthetic")) synthetic_families
                 else plastic_fams
  } else {
    keep_fams <- sel_fam
  }

  rows <- do.call(rbind, lapply(inst_levels, function(dev_label) {
    tbl       <- cts_present[[dev_label]]
    fams_here <- names(tbl)
    if (!is.null(keep_fams)) fams_here <- fams_here[fams_here %in% keep_fams]
    if (length(fams_here) == 0) return(NULL)
    data.frame(instrument = dev_label, family = fams_here,
               count = as.integer(tbl[fams_here]), stringsAsFactors = FALSE)
  }))

  if (is.null(rows) || nrow(rows) == 0)
    return(plotly::plot_ly() |>
             plotly::layout(title = "No matching data",
                            xaxis = list(visible = FALSE),
                            yaxis = list(visible = FALSE)))

  inst_totals  <- tapply(rows$count, rows$instrument, sum)
  rows$value   <- if (rel_mode)
    round(rows$count / inst_totals[rows$instrument] * 100, 1)
  else
    rows$count
  rows$tooltip <- if (rel_mode)
    paste0(rows$family, ": ", rows$value, "% (n=", rows$count, ")")
  else
    paste0(rows$family, ": ", rows$value)

  y_label    <- if (rel_mode) "Share (%)" else "Particle Count"
  title_str  <- paste0(if (sel_fam == "__all_plastics__") "All Plastics" else sel_fam,
                       " across instruments")

  fam_totals <- if (sel_fam == "__all_plastics__") {
    fams_present <- unique(rows$family)
    totals       <- vapply(fams_present, function(f) sum(rows$value[rows$family == f]), numeric(1))
    names(totals) <- fams_present
    sort(totals, decreasing = TRUE)
  } else {
    setNames(sum(rows$value), sel_fam)
  }
  fam_order <- names(fam_totals)

  p <- plotly::plot_ly()
  for (fam in fam_order) {
    sub         <- rows[rows$family == fam, ]
    sub_aligned <- data.frame(instrument = inst_levels, stringsAsFactors = FALSE)
    sub_aligned <- merge(sub_aligned, sub, by = "instrument", all.x = TRUE)
    sub_aligned$value[is.na(sub_aligned$value)]     <- 0
    sub_aligned$tooltip[is.na(sub_aligned$tooltip)] <- paste0(fam, ": 0")
    bar_label <- ifelse(sub_aligned$value > 0, as.character(sub_aligned$value), "")

    p <- plotly::add_trace(p,
      x                = sub_aligned$instrument,
      y                = sub_aligned$value,
      type             = "bar",
      name             = fam,
      text             = bar_label,
      textposition     = "inside",
      insidetextanchor = "middle",
      textfont         = list(size = 13, color = "white"),
      hovertext        = sub_aligned$tooltip,
      hoverinfo        = "text",
      marker           = list(color = if (sel_fam == "__all_plastics__")
                                (.fam_palette[fam] %||% "#cccccc")
                              else unname(device_colors[sub_aligned$instrument]))
    )
  }

  plotly::layout(p,
    barmode = if (sel_fam == "__all_plastics__") "stack" else "group",
    title   = list(text = title_str, x = 0.5, xanchor = "center",
                   font = list(size = 17)),
    xaxis   = list(title = "", tickfont = list(size = 14)),
    yaxis   = list(title = y_label,
                   titlefont = list(size = 14), tickfont = list(size = 13)),
    legend  = list(title = list(text = "<b>Family</b>", font = list(size = 14)),
                   font  = list(size = 13)),
    margin      = list(t = 65, r = 20, b = 50, l = 65),
    uniformtext = list(minsize = 10, mode = "hide")
  )
}

# ---------------------------------------------------------------------------
# Write a self-contained HTML report from the same `pages` list used by
# write_report_pdf().  Each page is rendered to a PNG and embedded as a
# base64 data URI so the file has no external dependencies except for the
# plotly.js CDN link used by the optional interactive chart.
#
# pages              : list of ggplot / grob objects (NULLs are skipped).
# path               : destination .html file path.
# plotly_fig         : optional plotly object inserted after page
#                      `plotly_insert_after`.
# plotly_insert_after: page index after which to inject the plotly widget
#                      (default 2 = after the static barplot page).
# width / height     : PNG dimensions in inches (matches PDF defaults).
# ---------------------------------------------------------------------------
write_report_html <- function(pages, path,
                              plotly_fig         = NULL,
                              plotly_insert_after = 2L,
                              plotly_figs        = NULL,
                              width = 11, height = 8.5) {

  pages <- Filter(Negate(is.null), pages)
  if (length(pages) == 0)
    pages <- list(report_text_page("Report",
      "Nothing to report \u2014 no data is loaded in the viewer."))

  # --- Render each page to a base64-encoded PNG ----------------------------
  encode_page <- function(pg, w = width, h = height) {
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    grDevices::png(tmp, width = w, height = h, units = "in", res = 144)
    ok <- tryCatch({
      if (inherits(pg, "ggplot")) print(pg)
      else { grid::grid.newpage(); grid::grid.draw(pg) }
      TRUE
    }, error = function(e) FALSE)
    grDevices::dev.off()
    if (!ok) return(NULL)
    paste0("data:image/png;base64,", base64enc::base64encode(tmp))
  }

  page_uris <- lapply(pages, encode_page)

  # --- Build plotly HTML snippets ------------------------------------------
  # `figs` is a list of list(fig=, after=, title=, caption=) so the report can
  # carry more than one interactive chart (absolute and relative share).
  .plotly_snippet <- function(fig, title, caption, idx) {
    if (is.null(fig) || !requireNamespace("plotly", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)) return("")
    built    <- plotly::plotly_build(fig)
    spec     <- built$x[c("data", "layout")]
    fig_json <- jsonlite::toJSON(spec, auto_unbox = TRUE, null = "null",
                                 na = "null", digits = 6)
    uid <- paste0("plotly-", format(Sys.time(), "%Y%m%d%H%M%S"), "-", idx)
    paste0(
      '<div class="report-section plotly-section">',
      '<h2>', title, '</h2>',
      '<p class="caption">', caption, '</p>',
      '<div id="', uid, '" style="width:100%;height:520px;"></div>',
      '<script>',
      '(function(){',
      'var spec=', fig_json, ';',
      'Plotly.newPlot("', uid, '",spec.data,spec.layout,',
      '{responsive:true,displayModeBar:true});',
      '})();',
      '</script>',
      '</div>'
    )
  }

  # Back-compat: a single plotly_fig/plotly_insert_after still works.
  if (is.null(plotly_figs) && !is.null(plotly_fig))
    plotly_figs <- list(list(fig = plotly_fig, after = plotly_insert_after))
  if (is.null(plotly_figs)) plotly_figs <- list()

  snippets <- list()
  for (i in seq_along(plotly_figs)) {
    spec <- plotly_figs[[i]]
    s <- .plotly_snippet(
      spec$fig,
      spec$title   %||% "Material Comparison \u2014 Interactive Chart",
      spec$caption %||% paste0("Hover over bars for exact values. ",
                               "Use the legend to show/hide families."),
      i)
    if (nzchar(s))
      snippets[[length(snippets) + 1]] <- list(after = spec$after %||% 2L,
                                               html = s)
  }
  has_plotly <- length(snippets) > 0

  # --- Assemble HTML -------------------------------------------------------
  css <- paste0(
    'body{font-family:Arial,sans-serif;background:#f0f2f5;margin:0;padding:20px;}',
    '.report-header{max-width:1200px;margin:0 auto 24px;padding:16px 24px;',
    'background:#1a3a5c;color:#fff;border-radius:6px;}',
    '.report-header h1{margin:0;font-size:1.5em;font-weight:600;}',
    '.report-header p{margin:4px 0 0;opacity:.8;font-size:.9em;}',
    '.report-section{max-width:1200px;margin:0 auto 20px;background:#fff;',
    'border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.12);overflow:hidden;}',
    '.report-section img{width:100%;height:auto;display:block;}',
    '.plotly-section{padding:20px;}',
    '.plotly-section h2{margin:0 0 8px;font-size:1.15em;color:#1a3a5c;}',
    '.caption{font-size:.8em;color:#666;margin:0 0 12px;}'
  )

  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  hdr <- paste0(
    '<div class="report-header">',
    '<h1>Multi-Instrument Particle Matching \u2014 Report</h1>',
    '<p>Generated: ', ts, '</p>',
    '</div>'
  )

  body_parts <- character(0)
  for (i in seq_along(page_uris)) {
    uri <- page_uris[[i]]
    if (!is.null(uri))
      body_parts <- c(body_parts, paste0(
        '<div class="report-section">',
        '<img src="', uri, '" alt="Report page ', i, '">',
        '</div>'
      ))
    for (sn in snippets)
      if (i == sn$after) body_parts <- c(body_parts, sn$html)
  }
  # Any snippet anchored beyond the last page is appended at the end.
  for (sn in snippets)
    if (sn$after > length(page_uris)) body_parts <- c(body_parts, sn$html)

  html <- paste0(
    '<!DOCTYPE html>\n<html lang="en">\n<head>\n',
    '<meta charset="UTF-8">\n',
    '<meta name="viewport" content="width=device-width,initial-scale=1">\n',
    '<title>Particle Analysis Report</title>\n',
    if (has_plotly)
      '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>\n'
    else "",
    '<style>', css, '</style>\n',
    '</head>\n<body>\n',
    hdr, '\n',
    paste(body_parts, collapse = "\n"),
    '\n</body>\n</html>\n'
  )

  writeLines(html, path, useBytes = FALSE)
  invisible(length(pages))
}

# ===========================================================================
# Shared report assembly (used by the pipeline's automatic report)
# ===========================================================================
# The viewer's download handler builds its pages from Shiny reactives, which
# main.R cannot call -- which is why the automatic report used to contain only
# tables. Everything below takes plain data frames and therefore runs equally
# well inside or outside a Shiny session. It reuses the SAME renderers the
# viewer uses (make_scatter, plot_size_distribution, report_*_page), so the
# two reports show the same figures for the same run.
#
# Fidelity note: the viewer's report reflects whatever filters, zoom, view
# rotation and highlight selections the operator had set. The automatic
# report has no operator, so it renders the DEFAULT view state -- which is
# what the viewer's report shows on a freshly opened session.
# ---------------------------------------------------------------------------

# ggplot twin of build_plotly_barplot(), for the PDF (plotly cannot be printed
# to a PDF device). Same families, same ordering, same palette.
build_material_barplot_gg <- function(device_counts,
                                      sel_fam  = "__all_plastics__",
                                      rel_mode = FALSE,
                                      show_all = TRUE,
                                      cat_mode = "both") {
  cts <- device_counts[!vapply(device_counts, is.null, logical(1))]
  if (length(cts) == 0) return(NULL)

  .fam_palette <- c(
    PE = "#e41a1c", PP = "#377eb8", PS = "#4daf4a", PET = "#984ea3",
    PVC = "#ff7f00", PA = "#a65628", PU = "#f781bf", PC = "#999999",
    PMMA = "#66c2a5", PTFE = "#fc8d62", ABS = "#e78ac3", Rubber = "#7570b3",
    Cellulose = "#bcbd22", Acrylate = "#17becf", Other = "#e5c494")

  plastic_fams <- c(synthetic_families, semi_synthetic_families)
  keep_fams <- if (sel_fam == "__all_plastics__") {
    if (show_all) NULL
    else if (identical(cat_mode, "synthetic")) synthetic_families
    else plastic_fams
  } else sel_fam

  inst_levels <- names(cts)
  rows <- do.call(rbind, lapply(inst_levels, function(lbl) {
    tbl <- cts[[lbl]]
    fams <- names(tbl)
    if (!is.null(keep_fams)) fams <- fams[fams %in% keep_fams]
    if (length(fams) == 0) return(NULL)
    data.frame(instrument = lbl, family = fams,
               count = as.integer(tbl[fams]), stringsAsFactors = FALSE)
  }))
  if (is.null(rows) || nrow(rows) == 0) return(NULL)

  inst_totals <- tapply(rows$count, rows$instrument, sum)
  rows$value <- if (rel_mode)
    round(rows$count / inst_totals[rows$instrument] * 100, 1) else rows$count

  # Largest total at the bottom of the stack, mirroring the plotly version.
  fams_present <- unique(rows$family)
  totals <- vapply(fams_present, function(f) sum(rows$value[rows$family == f]),
                   numeric(1))
  fam_order <- names(sort(setNames(totals, fams_present), decreasing = TRUE))
  rows$family     <- factor(rows$family, levels = rev(fam_order))
  rows$instrument <- factor(rows$instrument, levels = inst_levels)

  pal <- vapply(levels(rows$family),
                function(f) unname(.fam_palette[f]) %||% "#cccccc", character(1))
  names(pal) <- levels(rows$family)

  ggplot2::ggplot(rows, ggplot2::aes(x = instrument, y = value, fill = family)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(value > 0, value, "")),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.2, colour = "white") +
    ggplot2::scale_fill_manual(values = pal, name = "Family",
                               breaks = fam_order) +
    ggplot2::labs(x = NULL,
                  y = if (rel_mode) "Share (%)" else "Particle Count") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# Viewport for an instrument page: frame on the image when there is one
# (pad 200), else on the particles (pad 300). Mirrors the viewer's tabs.
.report_view_bounds <- function(img, df) {
  if (!is.null(img)) {
    pad <- 200
    return(sanitize_bounds(list(x = c(img$xmin - pad, img$xmax + pad),
                                y = c(img$ymin - pad, img$ymax + pad))))
  }
  if (!is.null(df) && nrow(df) > 0 && any(is.finite(df$x))) {
    pad <- 300
    return(sanitize_bounds(list(
      x = c(min(df$x, na.rm = TRUE) - pad, max(df$x, na.rm = TRUE) + pad),
      y = c(min(df$y, na.rm = TRUE) - pad, max(df$y, na.rm = TRUE) + pad))))
  }
  sanitize_bounds(list(x = c(-1000, 1000), y = c(-1000, 1000)))
}

# Resolve an instrument's background image + physical extent, running the same
# placement cascade the viewer runs (place_image_multirun), with the viewer's
# P3 fallback (aspect-preserving fit to the particle extent, 300um padding).
report_instrument_image <- function(key, df, meta, bg_path, native = TRUE) {
  if (is.null(bg_path) || !nzchar(bg_path) || !file.exists(bg_path)) return(NULL)
  raw <- tryCatch(load_image_raster(bg_path), error = function(e) NULL)
  if (is.null(raw)) return(NULL)
  xs <- if (native && "x_orig" %in% names(df)) df$x_orig else df$x
  ys <- if (native && "y_orig" %in% names(df)) df$y_orig else df$y
  ext <- tryCatch(place_image_multirun(key, meta, xs, ys, raw, bg_path),
                  error = function(e) NULL)
  if (is.null(ext))
    ext <- tryCatch(compute_image_bounds(raw, xs, ys, padding_um = 300),
                    error = function(e) NULL)
  if (is.null(ext)) return(NULL)
  c(list(raster = raw), ext[c("xmin", "xmax", "ymin", "ymax")])
}

# Assemble the full report page list, matching the viewer's download page for
# page: title, material barplot, plastics table, size distributions, size
# statistics, one page per instrument over its image, and the overlay.
#
# dfs        keyed list (ftir / ftir_bruker / raman / ldir) as returned by
#            build_instrument_dfs(), already quality-filtered by the caller
# meta       manifest config_snapshot (drives image placement)
# img_paths  keyed list of background image paths
build_report_pages <- function(dfs, meta = list(), img_paths = list(),
                               run_label = "unknown", run_id = "unknown",
                               quality_note = NULL, manifest = list()) {

  DEV <- c("FTIR (PerkinElmer)" = "ftir", "FTIR (Bruker)" = "ftir_bruker",
           "Raman" = "raman", "LDIR" = "ldir")
  n_of <- function(k) { d <- dfs[[k]]; if (is.null(d)) 0L else nrow(d) }
  pages <- list()

  # --- 1. Title / provenance ------------------------------------------------
  pages <- c(pages, list(report_text_page(
    "Multi-Instrument Particle Matching \u2014 Report",
    c(paste0("Run                 : ", run_label),
      paste0("Run ID              : ", run_id),
      paste0("Generated           : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      "Particles included in this report:",
      paste0("  FTIR (PerkinElmer): ", n_of("ftir")),
      paste0("  FTIR (Bruker)     : ", n_of("ftir_bruker")),
      paste0("  Raman             : ", n_of("raman")),
      paste0("  LDIR              : ", n_of("ldir")),
      "",
      report_input_file_lines(manifest),
      if (!is.null(quality_note)) "" else NULL,
      quality_note),
    subtitle = "Generated automatically by the pipeline")))

  # --- 2. Material comparison ----------------------------------------------
  counts <- lapply(setNames(nm = names(DEV)), function(lbl) {
    d <- dfs[[DEV[[lbl]]]]
    if (is.null(d) || nrow(d) == 0 || !"material" %in% names(d)) return(NULL)
    table(classify_family_vec(d$material))
  })
  # Two views of the same data: absolute counts show how much each instrument
  # found, relative share shows composition independent of that. Reporting only
  # one hides half the picture, so the report carries both.
  pages <- c(pages, list(report_full_figure_page(
    build_material_barplot_gg(counts, rel_mode = FALSE),
    "Material Comparison Across Instruments \u2014 All Plastics (stacked, absolute)",
    "Family counts per instrument, stacked. Absolute particle counts.")))
  pages <- c(pages, list(report_full_figure_page(
    build_material_barplot_gg(counts, rel_mode = TRUE),
    "Material Comparison Across Instruments \u2014 All Plastics (stacked, relative)",
    paste0("The same families as a share of each instrument's own total (%), ",
           "so composition can be compared across instruments that found ",
           "different numbers of particles."))))

  # --- 3. Plastics table ----------------------------------------------------
  devices <- Filter(Negate(is.null), lapply(setNames(nm = names(DEV)),
                                            function(lbl) dfs[[DEV[[lbl]]]]))
  pages <- c(pages, list(report_table_page(
    report_plastics_table(devices), "Plastics by Instrument",
    paste0("Material family counts per device. 'Unknown' families are ",
           "excluded; the Total row sums the families shown."))))

  # --- 4. Size distributions ------------------------------------------------
  size_plots <- list()
  if (n_of("ftir") > 0)
    size_plots <- c(size_plots,
                    list(plot_size_distribution(dfs$ftir, "FTIR (PerkinElmer)")))
  if (n_of("raman") > 0)
    size_plots <- c(size_plots,
                    list(plot_size_distribution(dfs$raman, "Raman",
                                                color_matched = "#1f77b4")))
  if (n_of("ldir") > 0)
    size_plots <- c(size_plots,
                    list(plot_size_distribution(dfs$ldir, "LDIR",
                                                color_matched = "#ff7f0e")))
  pages <- c(pages, list(report_grid_page(
    size_plots, "Size Distribution by Instrument",
    paste0("Feret Max (\u00b5m) per instrument. Solid bars: matched. ",
           "Outline bars: unmatched."),
    ncol = min(3, max(1, length(size_plots))))))

  # --- 5. Size statistics ---------------------------------------------------
  pages <- c(pages, list(report_table_page(
    report_size_stats_table(dfs), "Size Statistics",
    "Feret Max summary statistics per instrument.")))

  # --- 6-9. Per-instrument views over the instrument image ------------------
  # Native instrument frame (x_orig/y_orig), as each viewer tab shows it.
  inst_spec <- list(
    list(key = "ftir",        pkey = "ftir_perkin",
         title = "FTIR (PerkinElmer) \u2014 particles over instrument image",
         cols = c(matched = "#2ca02c", unmatched = "#d62728"),
         labs = c(matched = "matched to Raman", unmatched = "unmatched")),
    list(key = "ftir_bruker", pkey = "ftir_bruker",
         title = "FTIR (Bruker) \u2014 particles over instrument image",
         cols = c(matched = "#9467bd", unmatched = "#d62728"),
         labs = c(matched = "matched to Raman", unmatched = "unmatched")),
    list(key = "raman",       pkey = "raman",
         title = "Raman \u2014 particles over instrument image",
         cols = c(matched = "#1f77b4", unmatched = "#ff7f0e"),
         labs = c(matched = "matched to FTIR", unmatched = "unmatched")),
    list(key = "ldir",        pkey = "ldir",
         title = "LDIR \u2014 particles over instrument image",
         cols = c(matched = "#d62728", unmatched = "#ff7f0e"),
         labs = c(matched = "matched to Raman", unmatched = "unmatched")))

  for (sp in inst_spec) {
    d <- dfs[[sp$key]]
    if (is.null(d) || nrow(d) == 0) next
    dd <- d
    if (all(c("x_orig", "y_orig") %in% names(dd))) {
      dd$x <- dd$x_orig; dd$y <- dd$y_orig
    }
    img <- report_instrument_image(sp$pkey, dd, meta,
                                   img_paths[[sp$key]], native = FALSE)
    bounds <- .report_view_bounds(img, dd)
    pg <- tryCatch(
      make_scatter(dd, img, bounds,
                   paste0(sp$title, "  (", nrow(dd), " shown)"),
                   match_colours = sp$cols, match_labels = sp$labs),
      error = function(e) NULL)
    pages <- c(pages, list(report_figure_page(
      pg, sp$title,
      paste0("Native instrument frame. ", nrow(dd), " particles after the ",
             "report quality filter."))))
  }

  # --- 10. Overlay ----------------------------------------------------------
  pages <- c(pages, list(report_figure_page(
    build_overlay_plot(dfs),
    "Overlay \u2014 all instruments in the shared Raman frame",
    paste0("All instruments in aligned (Raman) coordinates. ",
           "One colour per instrument."))))

  Filter(Negate(is.null), pages)
}

# Overlay: every instrument's aligned coordinates in the shared Raman frame.
build_overlay_plot <- function(dfs) {
  spec <- list(
    list(key = "ftir",        lbl = "FTIR (PerkinElmer)", col = "#2ca02c"),
    list(key = "ftir_bruker", lbl = "FTIR (Bruker)",      col = "#9467bd"),
    list(key = "raman",       lbl = "Raman",              col = "#1f77b4"),
    list(key = "ldir",        lbl = "LDIR",               col = "#d62728"))
  parts <- list()
  for (s in spec) {
    d <- dfs[[s$key]]
    if (is.null(d) || nrow(d) == 0) next
    if (!all(c("x", "y") %in% names(d))) next
    ok <- is.finite(d$x) & is.finite(d$y)
    if (!any(ok)) next
    parts[[length(parts) + 1]] <- data.frame(
      x = d$x[ok], y = d$y[ok],
      feret_max = if ("feret_max" %in% names(d)) d$feret_max[ok] else 50,
      instrument = s$lbl, stringsAsFactors = FALSE)
  }
  if (length(parts) == 0) return(NULL)
  all_pts <- do.call(rbind, parts)
  present <- vapply(spec, function(s) s$lbl, character(1))
  present <- present[present %in% unique(all_pts$instrument)]
  pal <- vapply(spec, function(s) s$col, character(1))
  names(pal) <- vapply(spec, function(s) s$lbl, character(1))
  all_pts$instrument <- factor(all_pts$instrument, levels = present)

  ggplot2::ggplot(all_pts, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = instrument, size = feret_max),
                        alpha = 0.6) +
    ggplot2::scale_colour_manual(values = pal[present], name = "Instrument") +
    ggplot2::scale_size_continuous(name = "Feret Max (\u00b5m)",
                                   range = c(2, 10),
                                   limits = safe_size_limits(all_pts$feret_max)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "X (\u00b5m)", y = "Y (\u00b5m)") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "white", colour = NA))
}
