# =============================================================================
# 01_ingest.R — Data ingestion and column standardization
# =============================================================================

#' Read a data file (CSV or Excel) with automatic format detection
#'
#' @param filepath Path to data file (.csv, .xlsx, .xls)
#' @param sheet Sheet name for Excel files (ignored for CSV). A workbook
#'   without that sheet but with a single sheet is read from that one.
#' @return Data frame
read_data_file <- function(filepath, sheet = "Long_Table") {
  ext <- tolower(tools::file_ext(filepath))

  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(filepath)
    if (!is.null(sheet) && !sheet %in% sheets) {
      if (length(sheets) == 1) {
        log_message("  Sheet '", sheet, "' not found in ", basename(filepath),
                    "; reading its only sheet '", sheets, "'")
        sheet <- sheets
      } else {
        stop("Sheet '", sheet, "' not found in ", basename(filepath), ". Its sheets: ",
             paste0("'", sheets, "'", collapse = ", "), ". Set the sheet in the config ",
             "(ftir_sheet or raman_sheet).", call. = FALSE)
      }
    }
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
  # Sanitize encoding FIRST (prevents "invalid multibyte string" in tolower)
  x <- safe_colnames(x)
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
  as_number(df[[col]])
}


#' Safely extract numeric values from a named column (may be NULL)
safe_col_numeric <- function(df, col_name) {
  if (is.null(col_name)) return(rep(NA_real_, nrow(df)))
  as_number(df[[col_name]])
}


#' Stop unless an export has particle rows.
require_rows <- function(raw, filepath, instrument) {
  if (is.null(raw) || nrow(raw) == 0)
    stop(instrument, " file contains no particle rows: ", filepath, call. = FALSE)
  invisible(raw)
}


#' Stop if particle IDs are duplicated.
#'
#' IDs are join keys downstream (LDIR coordinate swaps, three-way triplets,
#' the viewer's selection), so a duplicate silently attaches data to the wrong
#' particle. Fail with the offending IDs instead.
require_unique_ids <- function(ids, filepath, instrument) {
  dup <- unique(ids[duplicated(ids) & !is.na(ids)])
  if (length(dup) > 0)
    stop(instrument, " file has ", length(dup), " duplicated particle ID(s): ",
         paste(utils::head(dup, 10), collapse = ", "),
         if (length(dup) > 10) ", ..." else "", " (", filepath, ")", call. = FALSE)
  invisible(ids)
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
  require_rows(raw, filepath, "FTIR")

  # Sanitize column names to valid UTF-8 immediately
  # (prevents "invalid multibyte string" from µ/²/³ in tolower/grep)
  names(raw) <- safe_colnames(names(raw))
  col_names  <- names(raw)

  log_message("  Raw FTIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns (sanitized): ", paste(col_names, collapse = ", "))

  # Find coordinate columns — ASCII-safe matching (no reliance on µ literal)
  coord_candidates <- grep("^Coord", col_names, ignore.case = TRUE, value = TRUE)

  if (length(coord_candidates) == 0) {
    # Fallback: any column containing "Coord" anywhere
    coord_candidates <- grep("Coord", col_names, ignore.case = TRUE, value = TRUE)
  }

  if (length(coord_candidates) == 0) {
    stop("Could not find any 'Coord' column in FTIR data. ",
         "Available columns: ", paste(col_names, collapse = ", "))
  }

  # Split into "pixels" column vs "µm" column (the one without "pixel")
  coord_px <- coord_candidates[grepl("pixel", coord_candidates, ignore.case = TRUE)]
  coord_um <- setdiff(coord_candidates, coord_px)

  # The µm column is required. Falling back to the pixel column would carry
  # map-pixel indices forward as micrometres (off by the pixel size, 6.25x on
  # a Spotlight map) with nothing downstream able to notice.
  if (length(coord_um) == 0)
    stop("FTIR file has no micrometre coordinate column (found only: ",
         paste(coord_candidates, collapse = ", "), "). Export 'Coord. [µm]' ",
         "or rename the micrometre column so it starts with 'Coord'. ", filepath,
         call. = FALSE)
  coord_col <- coord_um[1]

  log_message("  Coordinate candidates: ", paste(coord_candidates, collapse = ", "))
  log_message("  Using coordinate column: '", coord_col, "'")

  coords <- parse_ftir_coordinates(raw[[coord_col]])
  if (all(is.na(coords$x_um)) || all(is.na(coords$y_um)))
    stop("FTIR coordinate column '", coord_col, "' could not be parsed as ",
         "\"[x;y]\" (first value: ", raw[[coord_col]][1], "). ", filepath, call. = FALSE)

  # Find identifier column
  id_col <- find_column(raw, c("Identifier", "ID", "Particle"))
  particle_ids <- if (!is.null(id_col)) as.character(raw[[id_col]]) else paste0("P_", seq_len(nrow(raw)))
  require_unique_ids(particle_ids, filepath, "FTIR")

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
  require_rows(raw, filepath, "Raman")

  # Sanitize column names to valid UTF-8 immediately
  names(raw) <- safe_colnames(names(raw))

  log_message("  Raw Raman data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns (sanitized): ", paste(names(raw), collapse = ", "))

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
  require_unique_ids(particle_ids, filepath, "Raman")

  # Find material column
  mat_col <- find_column(raw, c("Material", "material_mapped", "material_raw",
                                 "Group", "Polymer"))
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = as_number(raw[[x_col]]),
    y_um         = as_number(raw[[y_col]]),
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


# ---------------------------------------------------------------------------
# FTIR Bruker ingestion (thin alias — same tabular format as PerkinElmer for now)
# ---------------------------------------------------------------------------

#' Ingest FTIR Bruker particle data
#'
#' Reads a Bruker OPUS / ALPHA export using the same column-mapping logic
#' as ingest_ftir().  The source_instrument column is set to "FTIR_bruker"
#' so the viewer can distinguish the two FTIR devices.
#'
#' @param filepath Path to the FTIR Bruker data file
#' @param sheet Sheet name for Excel files
#' @return Data frame with the same standardized columns as ingest_ftir()
ingest_ftir_bruker <- function(filepath, sheet = "Long_Table") {
  log_message("Reading FTIR (Bruker) data from: ", filepath)
  df <- ingest_ftir(filepath, sheet = sheet)
  df$source_instrument <- "FTIR_bruker"
  log_message("  FTIR (Bruker): ", nrow(df), " particles")
  df
}
