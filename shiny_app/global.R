# =============================================================================
# global.R — Load pipeline output and prepare data for Shiny
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
# never be rescued server-side — it is rejected outright.  Uploads that pass
# this ceiling are immediately downsized in memory for display (see
# downsample_raster() below), so storage and rendering stay small.
options(shiny.maxRequestSize = 50 * 1024^2)   # 50 MB

# Longest edge (px) to which uploaded background images are downsized.
# 2000 px keeps enough resolution to visually align particles against the
# membrane/micrograph while keeping annotation_raster() rendering fast.
BG_IMAGE_MAX_DIM <- 2000L

# Source canonical material classification from pipeline
# (classify_family_vec, classify_category, classify_category_vec, etc.)
source(file.path("..", "R", "08b_material_map.R"), local = TRUE)

# Dependency-free BMP reader (read_bmp_raster) — lets load_image_raster()
# handle instrument BMP exports even when magick is not installed.
source(file.path("..", "R", "read_bmp.R"), local = TRUE)

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

# points_df must have numeric columns x and y (in µm)
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

  # Then manifest$image_assets.  NB: assign only non-NULL results — writing
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
  # still on disk on the same machine — and load_image_raster() can read it
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

    # Core match files (staged names ← Part D pairwise naming)
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
# vacuous — every LDIR is always paired. The per-pair `match_distance` (aligned
# coordinate Euclidean distance, recomputed after any TPS refinement) is the
# real signal: a pair is a genuine match only when it falls within the LDIR
# acceptance gate. Classifying here, once, keeps the summary, the overlay plot,
# and the hover/tables perfectly consistent — the image agrees with the table.

# Return the LDIR<->Raman acceptance gate (µm) for a run. Read from the run
# manifest's config snapshot; falls back to the pipeline default (250 µm) for
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

# Recompute the within_gate classification against an explicit gate (µm). Used
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
# side on harmonized family names — including the cross-instrument pair tables.
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

  # FTIR scan bounds (optional — present in newer pipeline output)
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
# The PerkinElmer Spotlight exports images at ~6 rendering pixels per 25µm
# grid cell.  From a 2993×2993 image: (2993+1)/6 ≈ 499 grid positions,
# giving a 499 * 25 = 12475 µm scan extent.  This function computes the
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

  # Add LDIR→Raman match flag if LDIR-Raman match data available. Only
  # genuine (within-gate) pairs count as matched — an over-gate forced pairing
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

  # PNG: always use png::readPNG — avoids magick's image_data() producing
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
# Uses block-average (box filter) pooling — antialiased and dependency-free —
# with stride subsampling as a fallback for degenerate aspect ratios.
# Accepts 2D (grayscale) or 3D (H x W x channels) arrays; returns same form.
# The original pixel dimensions are recorded as attributes so callers that
# convert pixels to µm (e.g. TIFF DPI metadata, which refers to the ORIGINAL
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
# Auto-detect µm-per-pixel scale from TIFF resolution metadata.
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

    # Convert to µm/pixel
    if (units == "PixelsPerCentimeter") {
      um_per_px <- 10000 / dens   # 1 cm = 10 000 µm
    } else {
      um_per_px <- 25400 / dens   # 1 inch = 25 400 µm
    }

    # Reject common screen defaults — these are never real instrument values
    screen_dpis <- c(72, 96, 150, 300)
    effective_dpi <- if (units == "PixelsPerCentimeter") dens * 2.54 else dens
    if (round(effective_dpi) %in% screen_dpis) return(NULL)

    # Sanity: instrument images typically 0.5–50 µm/px
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
# in stage coordinates with Y UP — the stage-frame center is (cx, -cy).
# Verified on real data (PET A, 2026-07: particles Y [-68, 4877], panel
# Center Y = -4394): with Y negated 100% of particles fall inside the image;
# taken as-reported only 30% do.  Because the convention may vary across
# WITec versions/exports, BOTH interpretations are scored by the fraction of
# particles they contain and the better one wins; below min_frac the
# function returns NULL and the caller falls back to heuristic placement.
#
# cfg            : config_snapshot list from the run manifest
# x_orig, y_orig : particle stage coordinates (µm) used to score candidates
# min_frac       : minimum containment fraction to accept
# Returns list(xmin, xmax, ymin, ymax, y_negated, frac_inside) or NULL.
# ---------------------------------------------------------------------------
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
# Multi-Run image placement — mirror each single-instrument tab so the
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

# FTIR / Bruker: image spans the raw particle extent (ftir_native_image_info).
place_image_particle_extent <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NULL)
  list(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y))
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

# Raman P2: known µm-per-pixel scale (from meta, or read from a TIFF backdrop),
# centred on the particle mean — mirrors raman_native_image_info Priority 2. This
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
# (P1) then µm-per-pixel scale (P2); P3 (particle-extent fit) is left to the
# caller's fallback. `raw`/`bg_path` are needed only for the Raman P2 tier.
place_image_multirun <- function(instrument, meta, x, y, raw = NULL, bg_path = NULL) {
  switch(as.character(instrument),
    raman = {
      ext <- place_image_raman_meta(meta, x, y)                 # P1: WITec
      if (is.null(ext)) ext <- place_image_raman_umpx(meta, x, y, raw, bg_path)  # P2
      ext
    },
    ftir_perkin = place_image_particle_extent(x, y),
    ftir_bruker = place_image_particle_extent(x, y),
    ldir        = place_image_ldir_meta(meta),
    NULL)
}

# ---------------------------------------------------------------------------
# LDIR view rotation — rotate the whole native LDIR scene (image raster,
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
  out <- array(0, dim = c(nrow(ch[[1]]), ncol(ch[[1]]), length(ch)))
  for (k in seq_along(ch)) out[, , k] <- ch[[k]]
  out
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
ldir_auto_view_rotation <- function(ldir_x, ldir_y, raman_x, raman_y) {
  fk <- is.finite(ldir_x) & is.finite(ldir_y)
  fr <- is.finite(raman_x) & is.finite(raman_y)
  lx <- ldir_x[fk]; ly <- ldir_y[fk]
  rx <- raman_x[fr]; ry <- raman_y[fr]
  if (length(lx) < 4 || length(rx) < 4) return(0L)

  nrm <- function(x, y) {
    x <- x - mean(x); y <- y - mean(y)
    s <- sqrt(mean(x^2 + y^2)); if (!is.finite(s) || s <= 0) s <- 1
    list(x = x / s, y = y / s)
  }
  L <- nrm(lx, ly); R <- nrm(rx, ry)
  tol <- 0.10   # normalized units (~10% of cloud radius)

  score <- function(deg) {
    p <- rotate_xy_view(L$x, L$y, deg)
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

  degs <- c(0L, 90L, -90L, 180L)
  ns   <- vapply(degs, score, integer(1))
  bi   <- which.max(ns)
  if (ns[bi] < 4) return(0L)
  # only rotate when a rotation clearly beats leaving the view upright
  if (degs[bi] != 0L && ns[bi] <= ns[1] + 1L) return(0L)
  degs[bi]
}

# Total LDIR -> Raman rotation for a run, snapped to the nearest 90 deg —
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
# padding_um : padding in µm on each side
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

