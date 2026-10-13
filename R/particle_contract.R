# =============================================================================
# particle_contract.R — the internal data contract between readers and the
# rest of the pipeline
# =============================================================================
#
# Every reader (one per instrument export, plus the generic column-mapping
# reader) returns a CANONICAL PARTICLE TABLE and, where an image is involved, a
# CANONICAL IMAGE OBJECT. Everything downstream (alignment, matching,
# agreement, export, viewer) reads only these, never vendor column names.
#
# Canonical particle table -- one row per particle:
#   instrument      chr  role in the run: "ftir", "ftir_bruker", "raman", "ldir"
#   sample_id       chr  sample / filter label (NA when unknown)
#   particle_id     chr  unique within the table
#   x_um, y_um      num  particle centre in micrometres, in the instrument's own
#                        frame: any origin, X to the right, Y UP (Cartesian).
#                        NA only for LDIR before its image join.
#   feret_max_um    num  longest dimension (instrument's definition; FTIR
#                        exports the major axis, Raman a true Feret maximum)
#   feret_min_um    num  shortest dimension
#   major_um, minor_um num  ellipse axes (or best equivalent)
#   area_um2        num  area
#   material        chr  identification exactly as the instrument reports it
#   quality         num  the instrument's identification score (see quality_metric)
#   quality_metric  chr  what `quality` is: "AAU" (FTIR, 0-1), "HQI" (Raman,
#                        0-100), "LDIR quality" (0-1), or a user label
#   source_file     chr  basename of the export
# plus, added by annotate_materials(): material_polymer, material_family.
# Readers may add further columns (LDIR shape descriptors, identification_raw).
#
# Canonical image object -- an image placed in the particles' frame:
#   instrument, path         as above / file on disk
#   width_px, height_px      raster size
#   pixel_size_um            c(x, y) micrometres per pixel
#   origin_px                c(col, row): the pixel position that maps to (0, 0) um
#   y_axis                   "up": image rows grow downward while y_um grows
#                            upward (y_um = (origin_row - row) * pixel_size)
#   extent_um                c(xmin, xmax, ymin, ymax) of the whole raster
#   data_region_px           c(col_min, col_max, row_min, row_max) or NULL
#   provenance               how the geometry was obtained (see image_provenances)
# =============================================================================

canonical_particle_columns <- c(
  "instrument", "sample_id", "particle_id", "x_um", "y_um",
  "feret_max_um", "feret_min_um", "major_um", "minor_um", "area_um2",
  "material", "quality", "quality_metric", "source_file")

image_provenances <- c(
  "metadata",              # instrument / export metadata
  "config",                # extent entered in the config
  "registration",          # fitted to the particles visible in the image
  "scan_circle",           # LDIR scan circle detected in the image
  "processed_full_image",  # LDIR processed overlay: whole image = scan field
  "excel_area_scale",      # LDIR overlay scaled so blob areas match the table
  "particle_extent_approx" # approximate: fitted to the particle extent
)


#' Complete a reader's output to the canonical particle table.
#'
#' Adds the contract columns a reader does not produce itself (instrument,
#' sample_id, quality_metric, missing size columns as NA) and moves the
#' Columns and values the reader produced are left as they are; the new
#' columns are appended.
#'
#' @param df Reader output
#' @param instrument Role in the run ("ftir", "ftir_bruker", "raman", "ldir")
#' @param quality_metric What the quality column holds
#' @param sample_id Sample label (NA when unknown)
as_canonical_particles <- function(df, instrument, quality_metric = NA_character_,
                                   sample_id = NA_character_) {
  n <- nrow(df)
  if (!"instrument" %in% names(df))     df$instrument     <- rep(instrument, n)
  if (!"sample_id" %in% names(df))      df$sample_id      <- rep(as.character(sample_id), n)
  if (!"quality_metric" %in% names(df)) df$quality_metric <- rep(as.character(quality_metric), n)
  for (k in c("feret_max_um", "feret_min_um", "major_um", "minor_um", "area_um2", "quality"))
    if (!k %in% names(df)) df[[k]] <- rep(NA_real_, n)
  for (k in c("material", "source_file"))
    if (!k %in% names(df)) df[[k]] <- rep(NA_character_, n)
  df
}


#' Check a particle table against the contract; stop with every problem found.
#'
#' @param df Particle table
#' @param require_coords TRUE: every particle needs finite x_um/y_um.
#'   FALSE: coordinates may be NA (LDIR before the image join), but must be
#'   numeric.
#' @param what Label for messages (e.g. "FTIR (PerkinElmer) export")
#' @return df, invisibly
validate_particles <- function(df, require_coords = TRUE, what = "particle table") {
  problems <- character(0)
  if (!is.data.frame(df)) stop(what, ": not a data frame", call. = FALSE)
  if (nrow(df) == 0) problems <- c(problems, "has no rows")
  missing <- setdiff(canonical_particle_columns, names(df))
  if (length(missing))
    problems <- c(problems, paste0("missing column(s): ", paste(missing, collapse = ", ")))
  num_cols <- intersect(c("x_um", "y_um", "feret_max_um", "feret_min_um", "major_um",
                          "minor_um", "area_um2", "quality"), names(df))
  for (k in num_cols)
    if (!is.numeric(df[[k]]))
      problems <- c(problems, paste0("column '", k, "' is not numeric"))
  if ("particle_id" %in% names(df)) {
    ids <- as.character(df$particle_id)
    if (any(is.na(ids) | !nzchar(ids))) problems <- c(problems, "empty particle_id")
    dup <- unique(ids[duplicated(ids)])
    if (length(dup))
      problems <- c(problems, paste0("duplicated particle_id: ",
                                     paste(utils::head(dup, 10), collapse = ", ")))
  }
  if (all(c("x_um", "y_um") %in% names(df)) && is.numeric(df$x_um) && is.numeric(df$y_um)) {
    bad_inf <- is.infinite(df$x_um) | is.infinite(df$y_um)
    if (any(bad_inf)) problems <- c(problems, paste0(sum(bad_inf), " infinite coordinate(s)"))
    if (require_coords) {
      na <- is.na(df$x_um) | is.na(df$y_um)
      if (all(na)) problems <- c(problems, "no particle has coordinates")
      else if (any(na))
        message(what, ": ", sum(na), " of ", nrow(df),
                " particles have no coordinates and will be dropped by the prefilter")
    }
  }
  for (k in intersect(c("feret_max_um", "feret_min_um", "major_um", "minor_um", "area_um2"),
                      num_cols)) {
    if (is.numeric(df[[k]]) && any(df[[k]] < 0, na.rm = TRUE))
      problems <- c(problems, paste0("negative values in '", k, "'"))
  }
  if (length(problems))
    stop(what, " does not meet the particle contract:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  invisible(df)
}


#' Build a canonical image object from its pixel geometry.
#'
#' @param instrument,path,provenance See the header
#' @param width_px,height_px Raster size
#' @param pixel_size_um Scalar or c(x, y), micrometres per pixel
#' @param origin_px c(col, row): pixel position of the (0, 0) um point
#' @param y_axis "up" (default; y_um = (origin_row - row) * size) or "down"
#' @param data_region_px Optional c(col_min, col_max, row_min, row_max)
new_image_object <- function(instrument, path, width_px, height_px,
                             pixel_size_um, origin_px, provenance,
                             y_axis = "up", data_region_px = NULL) {
  s <- rep_len(as.numeric(pixel_size_um), 2)
  ox <- origin_px[1]; oy <- origin_px[2]
  xs <- (c(0, width_px) - ox) * s[1]
  ys <- if (identical(y_axis, "up")) (oy - c(height_px, 0)) * s[2]
        else (c(0, height_px) - oy) * s[2]
  img <- list(instrument = instrument, path = path,
              width_px = width_px, height_px = height_px,
              pixel_size_um = s, origin_px = c(ox, oy), y_axis = y_axis,
              extent_um = c(xmin = min(xs), xmax = max(xs), ymin = min(ys), ymax = max(ys)),
              data_region_px = data_region_px, provenance = provenance)
  validate_image(img)
}


#' Build a canonical image object from a placement extent (xmin..ymax), as the
#' FTIR and Raman placements produce: the raster spans the extent, row 0 at
#' ymax (y_axis "up"), or at ymin for a row-flipped image (y_axis "down").
image_object_from_extent <- function(instrument, path, width_px, height_px,
                                     extent, provenance, y_axis = "up") {
  sx <- (extent[["xmax"]] - extent[["xmin"]]) / width_px
  sy <- (extent[["ymax"]] - extent[["ymin"]]) / height_px
  oy <- if (identical(y_axis, "up")) extent[["ymax"]] / sy else -extent[["ymin"]] / sy
  new_image_object(instrument, path, width_px, height_px, c(sx, sy),
                   origin_px = c(-extent[["xmin"]] / sx, oy),
                   provenance = provenance, y_axis = y_axis)
}


#' Check an image object against the contract; stop with every problem found.
validate_image <- function(img, what = "image") {
  problems <- character(0)
  need <- c("instrument", "path", "width_px", "height_px", "pixel_size_um",
            "origin_px", "y_axis", "extent_um", "provenance")
  missing <- setdiff(need, names(img))
  if (length(missing))
    problems <- c(problems, paste0("missing field(s): ", paste(missing, collapse = ", ")))
  pos <- function(v) is.numeric(v) && length(v) >= 1 && all(is.finite(v)) && all(v > 0)
  if (!is.null(img$width_px) && !pos(img$width_px)) problems <- c(problems, "width_px must be > 0")
  if (!is.null(img$height_px) && !pos(img$height_px)) problems <- c(problems, "height_px must be > 0")
  if (!is.null(img$pixel_size_um) && !(pos(img$pixel_size_um) && length(img$pixel_size_um) == 2))
    problems <- c(problems, "pixel_size_um must be two positive numbers")
  if (!is.null(img$origin_px) && !(is.numeric(img$origin_px) && length(img$origin_px) == 2 &&
                                   all(is.finite(img$origin_px))))
    problems <- c(problems, "origin_px must be two finite numbers")
  if (!is.null(img$y_axis) && !img$y_axis %in% c("up", "down"))
    problems <- c(problems, "y_axis must be 'up' or 'down'")
  if (!is.null(img$extent_um)) {
    e <- img$extent_um
    if (!(is.numeric(e) && length(e) == 4 && all(is.finite(e)) && e[2] > e[1] && e[4] > e[3]))
      problems <- c(problems, "extent_um must be finite with xmax > xmin and ymax > ymin")
  }
  if (!is.null(img$provenance) && !img$provenance %in% image_provenances)
    problems <- c(problems, paste0("unknown provenance '", img$provenance, "'"))
  if (length(problems))
    stop(what, " does not meet the image contract:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  invisible(img)
}


#' Pixel positions in an image -> particle-frame micrometres.
#' @return data.frame(x_um, y_um)
image_pixels_to_um <- function(col_px, row_px, img) {
  x <- (col_px - img$origin_px[1]) * img$pixel_size_um[1]
  y <- if (identical(img$y_axis, "up")) (img$origin_px[2] - row_px) * img$pixel_size_um[2]
       else (row_px - img$origin_px[2]) * img$pixel_size_um[2]
  data.frame(x_um = x, y_um = y)
}

#' Particle-frame micrometres -> pixel positions in an image (inverse).
#' @return data.frame(col_px, row_px)
image_um_to_pixels <- function(x_um, y_um, img) {
  col <- x_um / img$pixel_size_um[1] + img$origin_px[1]
  row <- if (identical(img$y_axis, "up")) img$origin_px[2] - y_um / img$pixel_size_um[2]
         else img$origin_px[2] + y_um / img$pixel_size_um[2]
  data.frame(col_px = col, row_px = row)
}

#' Image object as a plain list for JSON (00_manifest/images.json).
image_object_record <- function(img) {
  img[c("instrument", "path", "width_px", "height_px", "pixel_size_um", "origin_px",
        "y_axis", "extent_um", "data_region_px", "provenance")]
}


#' Write the run's image objects to <run_dir>/00_manifest/images.json.
#'
#' @param run_dir Run output directory
#' @param images Named list of image objects (NULL entries are skipped)
write_image_objects <- function(run_dir, images) {
  images <- Filter(Negate(is.null), images)
  dir <- file.path(run_dir, "00_manifest")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, "images.json")
  writeLines(jsonlite::toJSON(lapply(images, image_object_record), pretty = TRUE,
                              auto_unbox = TRUE, digits = NA, null = "null"), path)
  invisible(path)
}
