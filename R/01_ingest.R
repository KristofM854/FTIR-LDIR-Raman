# =============================================================================
# 01_ingest.R — Data ingestion and column standardization
# =============================================================================

#' Ingest FTIR particle data from Excel export
#'
#' Reads the Long_Table sheet, parses coordinates from "[x;y]" format,
#' and standardizes column names.
#'
#' @param filepath Path to the FTIR Excel file
#' @param sheet Sheet name to read (default: "Long_Table")
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um, area_um2, major_um, minor_um,
#'   feret_min_um, feret_max_um, material, quality, source_file
ingest_ftir <- function(filepath, sheet = "Long_Table") {
  log_message("Reading FTIR data from: ", filepath)

  raw <- readxl::read_excel(filepath, sheet = sheet)
  log_message("  Raw FTIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")

  # Parse "[x;y]" coordinates
  coord_col <- grep("Coord.*\\[.*m\\]", names(raw), value = TRUE)
  if (length(coord_col) == 0) {
    stop("Could not find coordinate column in FTIR data. ",
         "Expected a column matching 'Coord.*[µm]'. ",
         "Available columns: ", paste(names(raw), collapse = ", "))
  }
  # Prefer the µm column over pixels
  coord_col_um <- grep("µm|um", coord_col, value = TRUE)
  if (length(coord_col_um) > 0) coord_col <- coord_col_um[1] else coord_col <- coord_col[1]

  coords <- parse_ftir_coordinates(raw[[coord_col]])

  # Build standardized data frame
  df <- data.frame(
    particle_id  = as.character(raw[["Identifier"]]),
    x_um         = coords$x_um,
    y_um         = coords$y_um,
    area_um2     = safe_numeric(raw, "Area on map [µm²]", alt = "Area on map"),
    major_um     = safe_numeric(raw, "Major dim [µm]",    alt = "Major dim"),
    minor_um     = safe_numeric(raw, "Minor dim [µm]",    alt = "Minor dim"),
    feret_min_um = safe_numeric(raw, "Feret min [µm]",    alt = "Feret min"),
    feret_max_um = safe_numeric(raw, "Major dim [µm]",    alt = "Major dim"),
    material     = as.character(raw[["Group"]]),
    quality      = safe_numeric(raw, "Max AAU score"),
    source_file  = if ("source_file" %in% names(raw)) as.character(raw[["source_file"]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " FTIR particles")
  df
}


#' Ingest Raman particle data from Excel export
#'
#' Reads the Long_Table sheet and standardizes column names.
#'
#' @param filepath Path to the Raman Excel file
#' @param sheet Sheet name to read (default: "Long_Table")
#' @return Data frame with standardized columns (same as ingest_ftir)
ingest_raman <- function(filepath, sheet = "Long_Table") {
  log_message("Reading Raman data from: ", filepath)

  raw <- readxl::read_excel(filepath, sheet = sheet)
  log_message("  Raw Raman data: ", nrow(raw), " rows, ", ncol(raw), " columns")

  # Extract coordinates from dedicated columns
  x_col <- find_column(raw, c("Visual Center Point X [µm]",
                               "Visual Center Point X [um]",
                               "Visual Center Point X"))
  y_col <- find_column(raw, c("Visual Center Point Y [µm]",
                               "Visual Center Point Y [um]",
                               "Visual Center Point Y"))

  # Size metrics
  area_col     <- find_column(raw, c("Area [µm²]", "Area [um²]", "Area"))
  feret_max_col <- find_column(raw, c("Feret Max [µm]", "Feret Max [um]", "Feret Max",
                                       "feret_max_um"))
  feret_min_col <- find_column(raw, c("Feret Min [µm]", "Feret Min [um]", "Feret Min"))
  length_col   <- find_column(raw, c("Length [µm]", "Length [um]", "Length"))
  width_col    <- find_column(raw, c("Width [µm]", "Width [um]", "Width"))

  # Material and quality
  material_col <- find_column(raw, c("Material", "material_raw", "material_mapped"))
  hqi_col      <- find_column(raw, c("HQI", "hqi_value"))

  df <- data.frame(
    particle_id  = as.character(raw[["Particle Name"]]),
    x_um         = as.numeric(raw[[x_col]]),
    y_um         = as.numeric(raw[[y_col]]),
    area_um2     = safe_col_numeric(raw, area_col),
    major_um     = safe_col_numeric(raw, length_col),
    minor_um     = safe_col_numeric(raw, width_col),
    feret_min_um = safe_col_numeric(raw, feret_min_col),
    feret_max_um = safe_col_numeric(raw, feret_max_col),
    material     = if (!is.null(material_col)) as.character(raw[[material_col]]) else NA_character_,
    quality      = safe_col_numeric(raw, hqi_col),
    source_file  = if ("source_file" %in% names(raw)) as.character(raw[["source_file"]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " Raman particles")
  df
}


# ---------------------------------------------------------------------------
# Column-finding helpers
# ---------------------------------------------------------------------------

#' Safely extract a numeric column by name, with fallback alternatives
safe_numeric <- function(df, col_name, alt = NULL) {
  candidates <- c(col_name, alt)
  for (cn in candidates) {
    matches <- grep(cn, names(df), fixed = TRUE, value = TRUE)
    if (length(matches) > 0) return(as.numeric(df[[matches[1]]]))
  }
  rep(NA_real_, nrow(df))
}

#' Find a column by trying multiple candidate names
find_column <- function(df, candidates) {
  for (cn in candidates) {
    if (cn %in% names(df)) return(cn)
    # Try partial match
    matches <- grep(cn, names(df), fixed = TRUE, value = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  NULL
}

#' Safely extract numeric values from a column that may be NULL
safe_col_numeric <- function(df, col_name) {
  if (is.null(col_name)) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col_name]])
}
