# =============================================================================
# 01_readers.R — one entry point per instrument role
# =============================================================================
#
# read_instrument(role, config) reads the export configured for a role and
# returns a canonical particle table (R/particle_contract.R). Per role it uses
#
#   * the generic column-mapping reader, when config[["<role>_column_map"]] is
#     set (an R list or the path of a JSON file; see docs/input_format.md), or
#   * the built-in reader for the instrument's own export format:
#       ftir         PerkinElmer Spotlight   ingest_ftir()        R/01_ingest.R
#       ftir_bruker  Bruker LUMOS (same table layout as the PerkinElmer export)
#                                            ingest_ftir_bruker() R/01_ingest.R
#       raman        WITec ParticleScout     ingest_raman()       R/01_ingest.R
#       ldir         Agilent 8700 LDIR       ingest_ldir()        R/01c_ingest_ldir.R
#
# and then, the same for every reader:
#   1. the role's material table (inst/materials/<role>.csv, or
#      config[["<role>_material_table"]]) translates instrument labels;
#      the label as exported is kept in identification_raw
#   2. LDIR only: identifications below ldir_hqi_unknown_threshold become
#      "unknown" (relabel_ldir_low_hqi)
#   3. the contract columns are completed (as_canonical_particles) and
#      material_polymer / material_family added (annotate_materials)
#   4. validate_particles() stops on a table that breaks the contract.
#
# Vendor column names live only in the built-in readers. Everything after
# read_instrument() uses the canonical columns.
# =============================================================================

instrument_roles <- list(
  ftir = list(
    label = "FTIR (PerkinElmer)", path_key = "ftir_path", quality_metric = "AAU",
    read = function(path, config) ingest_ftir(path, sheet = config$ftir_sheet %||% "Long_Table")),
  ftir_bruker = list(
    label = "FTIR (Bruker)", path_key = "ftir_bruker_path", quality_metric = "AAU",
    read = function(path, config) ingest_ftir_bruker(path)),
  raman = list(
    label = "Raman", path_key = "raman_path", quality_metric = "HQI",
    read = function(path, config) ingest_raman(path, sheet = config$raman_sheet %||% "Long_Table")),
  ldir = list(
    label = "LDIR", path_key = "ldir_path", quality_metric = "LDIR quality",
    read = function(path, config) ingest_ldir(path))
)


#' Read one instrument's export into a canonical particle table.
#'
#' @param role One of names(instrument_roles)
#' @param config Pipeline config
#' @param path Export to read (default: the role's path in config)
#' @return Canonical particle table (validated)
read_instrument <- function(role, config, path = config[[instrument_roles[[role]]$path_key]]) {
  spec <- instrument_roles[[role]]
  if (is.null(spec)) stop("Unknown instrument role '", role, "'", call. = FALSE)
  if (is.null(path) || !nzchar(path)) stop(spec$label, ": no input file configured", call. = FALSE)

  column_map <- config[[paste0(role, "_column_map")]]
  if (!is.null(column_map)) {
    column_map <- load_column_map(column_map)
    log_message("Reading ", spec$label, " with the column mapping reader: ", basename(path))
    df <- read_mapped_particles(path, column_map)
    quality_metric <- column_map$quality_metric %||% spec$quality_metric
  } else {
    df <- spec$read(path, config)
    quality_metric <- spec$quality_metric
  }

  if (!"identification_raw" %in% names(df)) df$identification_raw <- df$material
  df$material <- translate_materials(df$material, load_material_table(role, config),
                                     what = spec$label)
  if (role == "ldir") df <- relabel_ldir_low_hqi(df, config)

  df <- as_canonical_particles(df, role, quality_metric,
                               sample_id = config$sample_id %||% NA_character_)
  df <- annotate_materials(df)
  # LDIR exports carry no coordinates: they come from the image join later.
  validate_particles(df, require_coords = role != "ldir", what = paste(spec$label, "export"))
}


# ---------------------------------------------------------------------------
# Generic column-mapping reader
# ---------------------------------------------------------------------------

column_map_fields <- c("particle_id", "x_um", "y_um", "feret_max_um", "feret_min_um",
                       "major_um", "minor_um", "area_um2", "material", "quality")

#' A column mapping as a list: pass-through for a list, parsed for a JSON path.
load_column_map <- function(map) {
  if (is.character(map) && length(map) == 1) {
    if (!file.exists(map)) stop("Column mapping file not found: ", map, call. = FALSE)
    map <- jsonlite::fromJSON(map, simplifyVector = TRUE)
  }
  if (!is.list(map) || is.null(map$columns))
    stop("A column mapping needs a 'columns' entry (see docs/input_format.md)", call. = FALSE)
  unknown <- setdiff(names(map$columns), column_map_fields)
  if (length(unknown))
    stop("Column mapping: unknown field(s) ", paste(unknown, collapse = ", "),
         "; allowed: ", paste(column_map_fields, collapse = ", "), call. = FALSE)
  if (is.null(map$columns$material))
    stop("Column mapping: 'material' is required", call. = FALSE)
  map
}


#' Read any CSV / Excel particle table through a column mapping.
#'
#' The mapping names, for each canonical field, the column holding it, and
#' how to convert the values:
#'   columns            named list: canonical field -> column name in the file
#'   sheet              Excel sheet (default: the first)
#'   coordinate_unit_um micrometres per coordinate unit (1 = um, 1000 = mm)
#'   size_unit_um       micrometres per size unit (areas use its square)
#'   y_axis             "up" (default) or "down" (image-style: y grows downward;
#'                      y is negated so the table is in the contract's frame)
#'   quality_metric     label for the quality column (e.g. "HQI")
#'
#' Missing ID column: IDs are the row numbers. Missing feret_max_um / major_um
#' (or feret_min_um / minor_um): each is filled from the other.
#'
#' @param path CSV or Excel file
#' @param map Column mapping (see load_column_map)
#' @return Particle table with the canonical columns the mapping provides
read_mapped_particles <- function(path, map) {
  map <- load_column_map(map)
  raw <- if (tolower(tools::file_ext(path)) %in% c("xlsx", "xls"))
    as.data.frame(readxl::read_excel(path, sheet = map$sheet %||% 1))
  else read_data_file(path)
  require_rows(raw, path, "Particle table")

  cols <- map$columns
  missing <- setdiff(unlist(cols), names(raw))
  if (length(missing))
    stop("Column mapping: column(s) not in ", basename(path), ": ",
         paste(missing, collapse = ", "), "\n  Columns in the file: ",
         paste(names(raw), collapse = ", "), call. = FALSE)

  get_num <- function(field, unit = 1) {
    if (is.null(cols[[field]])) return(rep(NA_real_, nrow(raw)))
    as_number(raw[[cols[[field]]]]) * unit
  }
  cu <- as.numeric(map$coordinate_unit_um %||% 1)
  su <- as.numeric(map$size_unit_um %||% 1)
  y_sign <- switch(map$y_axis %||% "up", up = 1, down = -1,
                   stop("Column mapping: y_axis must be 'up' or 'down'", call. = FALSE))

  ids <- if (!is.null(cols$particle_id)) trimws(as.character(raw[[cols$particle_id]]))
         else as.character(seq_len(nrow(raw)))
  require_unique_ids(ids, path, "Particle table")

  feret_max <- get_num("feret_max_um", su); major <- get_num("major_um", su)
  feret_min <- get_num("feret_min_um", su); minor <- get_num("minor_um", su)
  if (is.null(cols$feret_max_um)) feret_max <- major
  if (is.null(cols$major_um))     major     <- feret_max
  if (is.null(cols$feret_min_um)) feret_min <- minor
  if (is.null(cols$minor_um))     minor     <- feret_min

  data.frame(
    particle_id  = ids,
    x_um         = get_num("x_um", cu),
    y_um         = y_sign * get_num("y_um", cu),
    area_um2     = get_num("area_um2", su^2),
    major_um     = major,
    minor_um     = minor,
    feret_min_um = feret_min,
    feret_max_um = feret_max,
    material     = trimws(as.character(raw[[cols$material]])),
    quality      = get_num("quality"),
    source_file  = basename(path),
    stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Per-instrument material tables
# ---------------------------------------------------------------------------

#' The material table for a role: config[["<role>_material_table"]] (a path or
#' a data frame), else inst/materials/<role>.csv. Columns: label, material.
#' @return data.frame(label, material); zero rows when there is none
load_material_table <- function(role, config) {
  src <- config[[paste0(role, "_material_table")]]
  if (is.null(src)) {
    dir <- config$material_table_dir %||% .default_material_dir()
    src <- if (!is.null(dir)) file.path(dir, paste0(role, ".csv"))
    if (is.null(src) || !file.exists(src))
      return(data.frame(label = character(), material = character()))
  }
  tab <- if (is.data.frame(src)) src
         else utils::read.csv(src, stringsAsFactors = FALSE, comment.char = "#",
                              strip.white = TRUE, encoding = "UTF-8")
  if (!all(c("label", "material") %in% names(tab)))
    stop("Material table ", if (is.character(src)) src else "", " needs the columns ",
         "'label' and 'material'", call. = FALSE)
  key <- toupper(trimws(tab$label))
  dup <- unique(key[duplicated(key)])
  if (length(dup))
    stop("Material table for ", role, " lists a label twice: ",
         paste(dup, collapse = ", "), call. = FALSE)
  tab[, c("label", "material")]
}


# inst/materials of the repository (run from the root, or from tests/testthat).
.default_material_dir <- function() {
  cands <- file.path(c(".", "..", file.path("..", "..")), "inst", "materials")
  hit <- cands[dir.exists(cands)]
  if (length(hit)) hit[1] else NULL
}


#' Replace instrument labels listed in a material table (case-insensitive,
#' surrounding spaces ignored); other labels are returned unchanged.
translate_materials <- function(material, table, what = "") {
  if (is.null(table) || nrow(table) == 0) return(material)
  hit <- match(toupper(trimws(material)), toupper(trimws(table$label)))
  n <- sum(!is.na(hit))
  if (n > 0) {
    material[!is.na(hit)] <- table$material[hit[!is.na(hit)]]
    log_message("  Material table: ", n, " ", what, " label(s) translated")
  }
  material
}
