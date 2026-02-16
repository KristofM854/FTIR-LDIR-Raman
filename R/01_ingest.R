# =============================================================================
# 01_ingest.R — Data ingestion and column standardization
# =============================================================================

#' Read a data file (CSV or Excel) with automatic format detection
#'
#' @param filepath Path to data file (.csv, .xlsx, .xls)
#' @param sheet Sheet name for Excel files (ignored for CSV)
#' @return Data frame
read_data_file <- function(filepath, sheet = "Long_Table") {
  ext <- tolower(tools::file_ext(filepath))

  if (ext %in% c("xlsx", "xls")) {
    return(readxl::read_excel(filepath, sheet = sheet))
  }

  if (ext == "csv") {
    # Read lines with latin1 encoding (handles ISO-8859 FTIR exports)
    lines <- readLines(filepath, encoding = "latin1", warn = FALSE)
    header <- lines[1]

    # Count potential delimiters in header
    n_semicolons <- nchar(gsub("[^;]", "", header))
    n_commas     <- nchar(gsub("[^,]", "", header))
    sep <- if (n_semicolons > n_commas) ";" else ","

    # Write to temp file to avoid encoding conversion issues in read.csv
    tmp <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(lines, tmp, useBytes = TRUE)

    raw <- read.csv(tmp, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
    return(raw)
  }

  stop("Unsupported file format: ", ext, ". Use .csv, .xlsx, or .xls")
}


#' Normalize column names for fuzzy matching
#'
#' Strips encoding artifacts, normalizes unicode, lowercases, and
#' reduces whitespace for robust column matching across file formats.
#' @param x Character vector of column names
#' @return Normalized character vector
normalize_colname <- function(x) {
  x <- tolower(x)
  # Replace common encoding artifacts for µ (micro sign / mu)
  x <- gsub("\xb5|\xc2\xb5|\xc2|\xb2|\xb3", "", x, useBytes = TRUE)
  x <- gsub("[µμ²³]", "", x)
  x <- gsub("\\s+", " ", trimws(x))
  x
}


#' Find a column by fuzzy matching against normalized names
#'
#' @param df Data frame
#' @param patterns Character vector of patterns to try (in priority order)
#' @return Column name (original) or NULL
find_column <- function(df, patterns) {
  col_names <- names(df)
  col_norm  <- normalize_colname(col_names)

  for (pat in patterns) {
    pat_norm <- normalize_colname(pat)

    # Exact normalized match
    idx <- which(col_norm == pat_norm)
    if (length(idx) > 0) return(col_names[idx[1]])

    # Substring match
    idx <- which(grepl(pat_norm, col_norm, fixed = TRUE))
    if (length(idx) > 0) return(col_names[idx[1]])
  }

  # Fallback: try regex partial match on original names
  for (pat in patterns) {
    # Extract the core word (e.g., "Feret Max" from "Feret Max [µm]")
    core <- gsub("\\s*\\[.*\\]", "", pat)
    idx <- grep(core, col_names, ignore.case = TRUE)
    if (length(idx) > 0) return(col_names[idx[1]])
  }

  NULL
}


#' Safely extract numeric values from a column found by fuzzy matching
safe_numeric <- function(df, ...) {
  col <- find_column(df, c(...))
  if (is.null(col)) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col]])
}


#' Safely extract numeric values from a named column (may be NULL)
safe_col_numeric <- function(df, col_name) {
  if (is.null(col_name)) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col_name]])
}


# ---------------------------------------------------------------------------
# FTIR ingestion
# ---------------------------------------------------------------------------

#' Ingest FTIR particle data
#'
#' Reads CSV or Excel, parses "[x;y]" coordinates, standardizes columns.
#'
#' @param filepath Path to the FTIR data file
#' @param sheet Sheet name for Excel files
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um, area_um2, major_um, minor_um,
#'   feret_min_um, feret_max_um, material, quality, source_file
ingest_ftir <- function(filepath, sheet = "Long_Table") {
  log_message("Reading FTIR data from: ", filepath)

  raw <- read_data_file(filepath, sheet = sheet)
  log_message("  Raw FTIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns: ", paste(names(raw), collapse = ", "))

  # Find the µm coordinate column (not the pixel one)
  col_names <- names(raw)
  coord_candidates <- grep("Coord", col_names, ignore.case = TRUE, value = TRUE)

  if (length(coord_candidates) == 0) {
    stop("Could not find any 'Coord' column in FTIR data. ",
         "Available columns: ", paste(col_names, collapse = ", "))
  }

  # Prefer the µm column — it's the one that does NOT contain "pixel"
  coord_col <- NULL
  for (cc in coord_candidates) {
    if (!grepl("pixel", cc, ignore.case = TRUE)) {
      coord_col <- cc
      break
    }
  }
  if (is.null(coord_col)) coord_col <- coord_candidates[1]

  log_message("  Using coordinate column: '", coord_col, "'")

  coords <- parse_ftir_coordinates(raw[[coord_col]])

  # Find identifier column
  id_col <- find_column(raw, c("Identifier", "ID", "Particle"))
  particle_ids <- if (!is.null(id_col)) as.character(raw[[id_col]]) else paste0("P_", seq_len(nrow(raw)))

  # Find material/group column
  mat_col <- find_column(raw, c("Group", "Material", "Polymer"))
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = coords$x_um,
    y_um         = coords$y_um,
    area_um2     = safe_numeric(raw, "Area on map"),
    major_um     = safe_numeric(raw, "Major dim"),
    minor_um     = safe_numeric(raw, "Minor dim"),
    feret_min_um = safe_numeric(raw, "Feret min"),
    feret_max_um = safe_numeric(raw, "Major dim"),
    material     = material,
    quality      = safe_numeric(raw, "Max AAU score", "AAU"),
    source_file  = if ("source_file" %in% names(raw)) as.character(raw[["source_file"]]) else basename(filepath),
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " FTIR particles")
  log_message("  Coordinate range: X [", round(min(df$x_um, na.rm = TRUE)),
              ", ", round(max(df$x_um, na.rm = TRUE)), "], Y [",
              round(min(df$y_um, na.rm = TRUE)), ", ",
              round(max(df$y_um, na.rm = TRUE)), "]")
  log_message("  Materials: ", paste(names(table(df$material)), collapse = ", "))

  df
}


# ---------------------------------------------------------------------------
# Raman ingestion
# ---------------------------------------------------------------------------

#' Ingest Raman particle data
#'
#' Reads CSV or Excel, standardizes column names.
#'
#' @param filepath Path to the Raman data file
#' @param sheet Sheet name for Excel files
#' @return Data frame with standardized columns (same as ingest_ftir)
ingest_raman <- function(filepath, sheet = "Long_Table") {
  log_message("Reading Raman data from: ", filepath)

  raw <- read_data_file(filepath, sheet = sheet)
  log_message("  Raw Raman data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns: ", paste(names(raw), collapse = ", "))

  # Find coordinate columns
  x_col <- find_column(raw, c("Visual Center Point X [µm]",
                               "Visual Center Point X"))
  y_col <- find_column(raw, c("Visual Center Point Y [µm]",
                               "Visual Center Point Y"))

  if (is.null(x_col) || is.null(y_col)) {
    stop("Could not find X/Y coordinate columns in Raman data. ",
         "Available columns: ", paste(names(raw), collapse = ", "))
  }

  log_message("  Using coordinate columns: '", x_col, "', '", y_col, "'")

  # Find particle ID column
  id_col <- find_column(raw, c("Particle Name", "Particle", "ID"))
  particle_ids <- if (!is.null(id_col)) as.character(raw[[id_col]]) else paste0("R_", seq_len(nrow(raw)))

  # Find material column
  mat_col <- find_column(raw, c("Material", "material_mapped", "material_raw",
                                 "Group", "Polymer"))
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = as.numeric(raw[[x_col]]),
    y_um         = as.numeric(raw[[y_col]]),
    area_um2     = safe_numeric(raw, "Area"),
    major_um     = safe_numeric(raw, "Length"),
    minor_um     = safe_numeric(raw, "Width"),
    feret_min_um = safe_numeric(raw, "Feret Min"),
    feret_max_um = safe_numeric(raw, "Feret Max"),
    material     = material,
    quality      = safe_numeric(raw, "HQI", "hqi_value"),
    source_file  = if ("source_file" %in% names(raw)) as.character(raw[["source_file"]]) else basename(filepath),
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " Raman particles")
  log_message("  Coordinate range: X [", round(min(df$x_um, na.rm = TRUE)),
              ", ", round(max(df$x_um, na.rm = TRUE)), "], Y [",
              round(min(df$y_um, na.rm = TRUE)), ", ",
              round(max(df$y_um, na.rm = TRUE)), "]")
  log_message("  Materials: ", paste(names(table(df$material)), collapse = ", "))

  df
}
