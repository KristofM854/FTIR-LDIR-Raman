# =============================================================================
# ftir_segments.R -- segmented PerkinElmer exports (temporary)
# =============================================================================
# The current Spotlight software exports a whole-filter scan as six segments:
# one particle table (<name>M1.csv ... M6.csv) and one "Average Abs." image
# (<name>_M1.bmp ... _M6.bmp) per segment, each in its own coordinate frame.
# stitch_ftir_segments() joins them into ONE table and ONE image that the
# pipeline reads as a normal PerkinElmer export; run_pipeline() does this by
# itself when the selected FTIR table is a segment (resolve_ftir_segments()).
# Remove this file, its call in run_pipeline() and tools/stitch_ftir_segments.R
# once the software exports a single file per scan.
#
# Geometry (field sample Faanui SW1L R4 A, checked against the images):
#   - segments M1..M6 form a 2 x 3 grid filled column by column, M1 top left:
#         M1 M4
#         M2 M5
#         M3 M6
#   - image: 2.00 um per pixel; table coordinates in um with y pointing up
#     from the bottom edge of each segment; no overlap between segments.
# A segment without particles has no table; its image is still needed.
# =============================================================================

FTIR_SEGMENT_TABLE <- "M([0-9]+)\\.(csv|xlsx)$"
FTIR_SEGMENT_IMAGE <- "M([0-9]+)\\.(bmp|png|tif|tiff)$"

#' The files of one segmented scan: every file in `path`'s folder with the same
#' name stem (the name before M<k>) and an extension matching `pattern`,
#' ordered by segment number. NULL when `path` is not a segment file.
ftir_segment_files <- function(path, pattern) {
  if (is.null(path) || !nzchar(path) || !grepl(pattern, basename(path), ignore.case = TRUE))
    return(NULL)
  stem <- sub(pattern, "", basename(path), ignore.case = TRUE)
  f <- list.files(dirname(path), full.names = TRUE)
  f <- f[grepl(pattern, basename(f), ignore.case = TRUE) &
         sub(pattern, "", basename(f), ignore.case = TRUE) == stem]
  f[order(.segment_number(f, pattern))]
}

.segment_number <- function(f, pattern)
  as.integer(sub(paste0(".*", pattern), "\\1", basename(f), ignore.case = TRUE))

#' Join the segments of one PerkinElmer scan.
#'
#' @param tables Segment tables (CSV or XLSX; segments without particles may
#'   be missing)
#' @param images The cols x rows segment images
#' @param out_dir Folder for <name>_stitched.csv / .png / .txt
#' @param cols,rows Segment grid, filled column by column from the top left
#' @param um_per_px Image scale
#' @param say Logger (cat by default; log_message in the pipeline)
#' @return list(csv, png, placement, width_um, height_um, n_particles,
#'   on_particle = share of particles on a bright image spot)
stitch_ftir_segments <- function(tables, images, out_dir, cols = 2, rows = 3,
                                 um_per_px = 2,
                                 say = function(...) cat(..., "\n", sep = "")) {
  n <- cols * rows
  img_k <- .segment_number(images, FTIR_SEGMENT_IMAGE)
  if (length(images) != n || !setequal(img_k, seq_len(n)))
    stop("Expected ", n, " FTIR segment images (M1..M", n, "); found: ",
         paste(basename(images), collapse = ", "), call. = FALSE)
  images <- images[order(img_k)]
  if (length(tables) == 0) stop("No FTIR segment tables", call. = FALSE)
  tables <- tables[order(.segment_number(tables, FTIR_SEGMENT_TABLE))]

  # --- Image: tiles column by column ----------------------------------------
  tiles <- lapply(images, magick::image_read)
  info <- do.call(rbind, lapply(tiles, magick::image_info))
  if (length(unique(info$width)) > 1 || length(unique(info$height)) > 1)
    stop("FTIR segment images differ in size", call. = FALSE)
  W <- info$width[1]; H <- info$height[1]
  column <- function(j) magick::image_append(do.call(c, tiles[(j - 1) * rows + seq_len(rows)]),
                                             stack = TRUE)
  mosaic <- magick::image_append(do.call(c, lapply(seq_len(cols), column)))
  W_um <- W * um_per_px; H_um <- H * um_per_px

  # --- Tables: shift each segment into the mosaic frame (y up, origin at the
  # mosaic's bottom-left corner) -------------------------------------------
  read_seg <- function(f) {
    if (grepl("\\.xlsx$", f, ignore.case = TRUE))
      return(as.data.frame(readxl::read_excel(f), check.names = FALSE))
    txt <- readLines(f, warn = FALSE)
    if (!all(validUTF8(txt))) txt <- iconv(txt, "latin1", "UTF-8")
    utils::read.csv(text = txt, check.names = FALSE, stringsAsFactors = FALSE)
  }
  um_column <- function(d) {
    cand <- grep("^Coord\\..*m\\]$", names(d), value = TRUE)
    cand[!grepl("pixel", cand)][1]
  }
  pair <- function(s) {
    m <- do.call(rbind, lapply(strsplit(gsub("\\[|\\]", "", as.character(s)), ";"), as.numeric))
    if (is.null(dim(m))) m <- matrix(m, ncol = 2)
    m
  }
  out <- list()
  for (f in tables) {
    k <- .segment_number(f, FTIR_SEGMENT_TABLE)
    d <- read_seg(f)
    if (nrow(d) == 0) next
    um_col <- um_column(d)
    if (is.na(um_col)) stop("No 'Coord. [um]' column in ", basename(f), call. = FALSE)
    px_col <- grep("^Coord\\..*pixel", names(d), value = TRUE)[1]
    c0 <- (k - 1) %/% rows; r0 <- (k - 1) %% rows          # column, row (row 0 = top)
    dx <- c0 * W_um; dy <- (rows - 1 - r0) * H_um
    um <- pair(d[[um_col]])
    if (!is.na(px_col)) {                                 # the chemical map's own pixel grid
      map_um <- stats::median(um[, 1] / pair(d[[px_col]])[, 1], na.rm = TRUE)
      d[[px_col]] <- sprintf("[%d;%d]", round((um[, 1] + dx) / map_um), round((um[, 2] + dy) / map_um))
    }
    d[[um_col]] <- sprintf("[%d;%d]", round(um[, 1] + dx), round(um[, 2] + dy))
    d[[1]] <- paste0("M", k, "_", d[[1]])                    # identifiers unique across segments
    out[[length(out) + 1]] <- d
  }
  tab <- do.call(rbind, out)
  if (is.null(tab)) stop("The FTIR segment tables hold no particles", call. = FALSE)

  # --- Check: do the shifted particles land on particles in the mosaic? -----
  # A particle is any spot whose colour stands out from the background: the
  # absorbance map draws strong absorbers green/yellow but weak ones (much
  # cellulose and protein) blue, so brightness in one channel is not enough.
  bm <- magick::image_data(mosaic, "rgb")                   # raw [3, width, height]
  bg <- vapply(1:3, function(ch) stats::median(as.integer(bm[ch, seq(1, cols * W, 7),
                                                             seq(1, rows * H, 7)])), 0)
  um <- pair(tab[[um_column(tab)]])
  col_px <- round(um[, 1] / um_per_px); row_px <- round(rows * H - um[, 2] / um_per_px)
  hit <- mapply(function(cx, ry) {
    xs <- max(1, cx - 4):min(cols * W, cx + 4); ys <- max(1, ry - 4):min(rows * H, ry + 4)
    v <- as.integer(bm[, xs, ys]); dim(v) <- c(3, length(v) / 3)
    max(sqrt(colSums((v - bg)^2))) > 100 }, col_px, row_px)
  say(sprintf("  Stitched %d FTIR segments (%d x %d px, %.0f x %.0f um): %d particles from %d tables, %d of them on a particle in the image",
              n, cols * W, rows * H, cols * W_um, rows * H_um, nrow(tab), length(tables), sum(hit)))
  if (mean(hit) < 0.5)
    warning("Fewer than half of the FTIR particles fall on a particle in the stitched image: ",
            "check the segment order and the image scale", call. = FALSE)

  base <- sub(FTIR_SEGMENT_TABLE, "", basename(tables[1]), ignore.case = TRUE)
  base <- sub("[ _]+$", "", base)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  csv_out <- file.path(out_dir, paste0(base, "_stitched.csv"))
  png_out <- file.path(out_dir, paste0(base, "_stitched.png"))
  # write.csv escapes the "µ" of the column names in a non-UTF-8 locale
  # ("<U+00B5>"), which the reader would no longer recognise: write the bytes
  body <- utils::capture.output(utils::write.table(tab, sep = ",", row.names = FALSE,
                                                   col.names = FALSE, qmethod = "double"))
  header <- paste0('"', enc2utf8(names(tab)), '"', collapse = ",")
  con <- file(csv_out, "wb"); writeLines(enc2utf8(c(header, body)), con, useBytes = TRUE); close(con)
  magick::image_write(mosaic, png_out, format = "png")
  # Placement next to the image (read_ftir_image_metadata()): the image spans
  # [0, width] x [0, height] um of the table's frame
  txt_out <- paste0(tools::file_path_sans_ext(png_out), ".txt")
  writeLines(c("FTIR image placement (stitched FTIR segments), um in the table's frame",
               sprintf("Width: %.0f", cols * W_um), sprintf("Height: %.0f", rows * H_um),
               sprintf("Center X: %.1f", cols * W_um / 2), sprintf("Center Y: %.1f", rows * H_um / 2)),
             txt_out)
  invisible(list(csv = csv_out, png = png_out, placement = txt_out,
                 width_um = cols * W_um, height_um = rows * H_um,
                 n_particles = nrow(tab), on_particle = mean(hit)))
}

#' When the selected PerkinElmer table is one segment of a segmented export,
#' stitch the scan into the run's inputs/ folder and point the config at the
#' stitched table and image. The images are taken from the selected FTIR image
#' (any one segment) or, without one, from the table's folder.
#'
#' @return config, with ftir_path / ftir_image replaced when stitched
resolve_ftir_segments <- function(config) {
  tables <- ftir_segment_files(config$ftir_path, FTIR_SEGMENT_TABLE)
  if (length(tables) < 2) return(config)       # a single export (or one segment only)
  images <- ftir_segment_files(config$ftir_image, FTIR_SEGMENT_IMAGE)
  if (is.null(images)) {
    cand <- list.files(dirname(config$ftir_path), FTIR_SEGMENT_IMAGE, full.names = TRUE,
                       ignore.case = TRUE)
    stems <- unique(sub(FTIR_SEGMENT_IMAGE, "", basename(cand), ignore.case = TRUE))
    if (length(stems) == 1) images <- cand
  }
  if (length(images) == 0)
    stop("The FTIR table is one of ", length(tables), " segments, but the segment images ",
         "(..._M1.bmp ... _M6.bmp) were not found: select one of them as the FTIR image",
         call. = FALSE)
  log_message("FTIR (PerkinElmer): segmented export (", length(tables), " tables, ",
              length(images), " images) — stitching")
  st <- stitch_ftir_segments(tables, images, file.path(config$output_dir, "inputs", "ftir_stitched"),
                             say = function(...) log_message(...))
  config$ftir_segment_tables <- tables
  config$ftir_segment_images <- images
  config$ftir_path  <- st$csv
  config$ftir_image <- st$png
  config
}
