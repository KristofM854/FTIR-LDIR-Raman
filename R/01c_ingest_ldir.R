# =============================================================================
# 01c_ingest_ldir.R — LDIR data ingestion
# =============================================================================
#
# The Agilent 8700 LDIR exports particle data (size, shape, material) in an
# Excel file with two sheets:
#   - "Info"      — sample name, instrument serial, etc.
#   - "Particles" — one row per particle with morphological + spectral data
#
# Critically, LDIR does NOT export particle X/Y centroid coordinates.
# This module ingests the tabular data and leaves coordinates at NA,
# to be populated later by image-based extraction (01b_ingest_image.R)
# or by non-spatial comparison methods.
# =============================================================================


#' Ingest LDIR particle data
#'
#' Reads the Agilent 8700 LDIR export and standardizes columns.
#'
#' @param filepath Path to the LDIR Excel file (.xlsx)
#' @param sheet    Sheet name to read (default "Particles")
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um (NA), area_um2, major_um, minor_um,
#'   feret_min_um, feret_max_um, material, quality, source_file
ingest_ldir <- function(filepath, sheet = "Particles") {
  log_message("Reading LDIR data from: ", filepath)

  raw <- readxl::read_excel(filepath, sheet = sheet)
  log_message("  Raw LDIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns: ", paste(names(raw), collapse = ", "))

  # Try to read sample info from "Info" sheet
  info_sheet <- tryCatch({
    readxl::read_excel(filepath, sheet = "Info")
  }, error = function(e) NULL)

  if (!is.null(info_sheet)) {
    sample_name <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Sample Name"]
    }, error = function(e) NA_character_)
    instrument <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Instrument Name"]
    }, error = function(e) NA_character_)
    log_message("  Sample: ", sample_name)
    log_message("  Instrument: ", instrument)
  }

  # Find columns by fuzzy matching
  id_col     <- find_column(raw, c("Id", "Identifier", "#"))
  width_col  <- find_column(raw, c("Width"))
  height_col <- find_column(raw, c("Height"))
  diam_col   <- find_column(raw, c("Diameter"))
  area_col   <- find_column(raw, c("Area"))
  perim_col  <- find_column(raw, c("Perimeter"))
  mat_col    <- find_column(raw, c("Identification", "Material"))
  qual_col   <- find_column(raw, c("Quality"))
  valid_col  <- find_column(raw, c("Is Valid"))
  aspect_col <- find_column(raw, c("Aspect Ratio"))

  # Particle IDs
  particle_ids <- if (!is.null(id_col)) {
    as.character(raw[[id_col]])
  } else {
    paste0("LDIR_", seq_len(nrow(raw)))
  }

  # Sizes
  width_um  <- safe_col_numeric(raw, width_col)
  height_um <- safe_col_numeric(raw, height_col)
  diam_um   <- safe_col_numeric(raw, diam_col)

  # For feret max/min: use width and height as proxies
  feret_max <- pmax(width_um, height_um, na.rm = TRUE)
  feret_min <- pmin(width_um, height_um, na.rm = TRUE)

  # Material
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Quality score
  quality <- safe_col_numeric(raw, qual_col)

  # Is Valid filter
  if (!is.null(valid_col)) {
    valid <- raw[[valid_col]]
    n_invalid <- sum(tolower(valid) != "true", na.rm = TRUE)
    if (n_invalid > 0) {
      log_message("  Flagged ", n_invalid, " invalid particles (keeping all for now)")
    }
  }

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = NA_real_,    # LDIR does not export coordinates
    y_um         = NA_real_,    # will be populated by image extraction
    area_um2     = safe_col_numeric(raw, area_col),
    major_um     = feret_max,
    minor_um     = feret_min,
    feret_min_um = feret_min,
    feret_max_um = feret_max,
    material     = material,
    quality      = quality,
    source_file  = basename(filepath),
    # Keep LDIR-specific columns
    diameter_um  = diam_um,
    aspect_ratio = safe_col_numeric(raw, aspect_col),
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " LDIR particles")
  log_message("  Coordinates: NOT available (LDIR does not export X/Y)")
  log_message("  Size range (feret max): ",
              round(min(df$feret_max_um, na.rm = TRUE), 1), " – ",
              round(max(df$feret_max_um, na.rm = TRUE), 1), " µm")
  log_message("  Materials: ", paste(names(head(sort(table(df$material), decreasing = TRUE), 10)),
                                     collapse = ", "))

  df
}


#' Pre-filter LDIR particle data
#'
#' @param df Data frame from ingest_ldir()
#' @param min_quality Minimum quality score (0 = keep all)
#' @param min_size_um Minimum particle size in µm (applied to feret_max_um)
#' @param remove_invalid Logical, remove "Is Valid" = FALSE particles
#' @return Filtered data frame
prefilter_ldir <- function(df, min_quality = 0, min_size_um = 0,
                           remove_invalid = FALSE) {
  n_start <- nrow(df)
  log_message("Pre-filtering LDIR: ", n_start, " particles")

  # Quality filter
  if (min_quality > 0 && any(!is.na(df$quality))) {
    keep <- is.na(df$quality) | df$quality >= min_quality
    df <- df[keep, ]
    log_message("  Quality filter (>= ", min_quality, "): kept ", nrow(df),
                " of ", n_start, " particles")
  }

  # Size filter
  if (min_size_um > 0 && any(!is.na(df$feret_max_um))) {
    keep <- is.na(df$feret_max_um) | df$feret_max_um >= min_size_um
    df <- df[keep, ]
    log_message("  Size filter (>= ", min_size_um, " µm): kept ", nrow(df), " particles")
  }

  log_message("  LDIR after filtering: ", nrow(df), " particles")
  df
}
