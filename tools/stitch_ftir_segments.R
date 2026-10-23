# =============================================================================
# stitch_ftir_segments.R -- temporary workaround for segmented PerkinElmer exports
# =============================================================================
# The current Spotlight software exports a whole-filter scan as six segments:
# one particle table (<name>M1.csv ... M6.csv) and one "Average Abs." image
# (<name>_M1.bmp ... _M6.bmp) per segment, each in its own coordinate frame.
# This script joins them into ONE table and ONE image that the pipeline reads
# as a normal PerkinElmer export. Remove it once the software exports a single
# file (announced for the next software update).
#
# Geometry (field sample Faanui SW1L R4 A, checked against the images):
#   - segments M1..M6 form a 2 x 3 grid filled column by column, M1 top left:
#         M1 M4
#         M2 M5
#         M3 M6
#   - image: 2.00 um per pixel; table coordinates in um with y pointing up
#     from the bottom edge of each segment; no overlap between segments.
# A segment without plastics has no table; its image is still needed.
#
#   Rscript tools/stitch_ftir_segments.R <folder> [output folder]
#
# Writes <name>_stitched.csv and <name>_stitched.png, checks that the particles
# land on particles in the stitched image, and prints the two config values
# that place the image exactly (ftir_image_width_um / _height_um).
# =============================================================================

stitch_ftir_segments <- function(dir, out_dir = dir, cols = 2, rows = 3,
                                 um_per_px = 2, quiet = FALSE) {
  say <- function(...) if (!quiet) cat(..., "\n", sep = "")
  seg_of <- function(f) as.integer(sub(".*M([0-9]+)\\.[A-Za-z]+$", "\\1", f))
  imgs <- list.files(dir, "M[0-9]+\\.(bmp|png|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  tabs <- list.files(dir, "M[0-9]+\\.(csv|xlsx)$", full.names = TRUE, ignore.case = TRUE)
  n <- cols * rows
  if (length(imgs) != n || !setequal(seg_of(imgs), seq_len(n)))
    stop("Expected ", n, " segment images ending in M1..M", n, " in ", dir,
         "; found: ", paste(basename(imgs), collapse = ", "), call. = FALSE)
  imgs <- imgs[order(seg_of(imgs))]
  if (length(tabs) == 0) stop("No segment tables (M<k>.csv / .xlsx) in ", dir, call. = FALSE)
  tabs <- tabs[order(seg_of(tabs))]

  # --- Image: tiles column by column ----------------------------------------
  tiles <- lapply(imgs, magick::image_read)
  info <- do.call(rbind, lapply(tiles, magick::image_info))
  if (length(unique(info$width)) > 1 || length(unique(info$height)) > 1)
    stop("Segment images differ in size", call. = FALSE)
  W <- info$width[1]; H <- info$height[1]
  column <- function(j) magick::image_append(do.call(c, tiles[(j - 1) * rows + seq_len(rows)]),
                                             stack = TRUE)
  mosaic <- magick::image_append(do.call(c, lapply(seq_len(cols), column)))
  W_um <- W * um_per_px; H_um <- H * um_per_px

  # --- Tables: shift each segment into the mosaic frame (y up, origin at the
  # mosaic's bottom-left corner) -------------------------------------------
  read_seg <- function(f) {
    if (grepl("\\.xlsx$", f, ignore.case = TRUE)) {
      d <- as.data.frame(readxl::read_excel(f), check.names = FALSE)
    } else {
      txt <- readLines(f, warn = FALSE)
      if (!all(validUTF8(txt))) txt <- iconv(txt, "latin1", "UTF-8")
      d <- utils::read.csv(text = txt, check.names = FALSE, stringsAsFactors = FALSE)
    }
    d
  }
  um_column <- function(d) {
    cand <- grep("^Coord\\..*m\\]$", names(d), value = TRUE)
    cand[!grepl("pixel", cand)][1]
  }
  pair <- function(s) {
    v <- strsplit(gsub("\\[|\\]", "", as.character(s)), ";")
    m <- do.call(rbind, lapply(v, as.numeric))
    if (is.null(dim(m))) m <- matrix(m, ncol = 2)
    m
  }
  out <- list()
  for (f in tabs) {
    k <- seg_of(f)
    d <- read_seg(f)
    if (nrow(d) == 0) next
    um_col <- um_column(d)
    px_col <- grep("^Coord\\..*pixel", names(d), value = TRUE)[1]
    if (is.na(um_col)) stop("No 'Coord. [um]' column in ", basename(f), call. = FALSE)
    c0 <- (k - 1) %/% rows; r0 <- (k - 1) %% rows          # column, row (row 0 = top)
    dx <- c0 * W_um; dy <- (rows - 1 - r0) * H_um
    um <- pair(d[[um_col]])
    # map pixel size of the chemical map (the table's own pixel grid)
    if (!is.na(px_col)) {
      px <- pair(d[[px_col]])
      map_um <- stats::median(um[, 1] / px[, 1], na.rm = TRUE)
      d[[px_col]] <- sprintf("[%d;%d]", round((um[, 1] + dx) / map_um), round((um[, 2] + dy) / map_um))
    }
    d[[um_col]] <- sprintf("[%d;%d]", round(um[, 1] + dx), round(um[, 2] + dy))
    d[[1]] <- paste0("M", k, "_", d[[1]])                    # identifiers unique across segments
    out[[length(out) + 1]] <- d
  }
  tab <- do.call(rbind, out)

  # --- Check: do the shifted particles land on particles in the mosaic? -----
  g <- as.integer(magick::image_data(mosaic, "rgb"))[, , 2]  # particles are green/yellow
  um <- pair(tab[[um_column(tab)]])
  col_px <- round(um[, 1] / um_per_px); row_px <- round(rows * H - um[, 2] / um_per_px)
  hit <- mapply(function(cx, ry) {
    xs <- max(1, cx - 4):min(cols * W, cx + 4); ys <- max(1, ry - 4):min(rows * H, ry + 4)
    max(g[ys, xs]) > 100 }, col_px, row_px)
  say(sprintf("Stitched %d segments (%d x %d px, %.0f x %.0f um); %d particles from %d tables",
              n, cols * W, rows * H, cols * W_um, rows * H_um, nrow(tab), length(tabs)))
  say(sprintf("Particles on a bright spot of the stitched image: %d of %d", sum(hit), length(hit)))
  if (mean(hit) < 0.5)
    warning("Fewer than half of the particles fall on a particle in the image: check the ",
            "segment order and um_per_px before using the output", call. = FALSE)

  base <- sub("_?M[0-9]+\\.[A-Za-z]+$", "", basename(tabs[1]))
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
  say("Wrote ", csv_out, "\nWrote ", png_out)
  say(sprintf(paste0("To place the image exactly, set in the config (00_config.R or ",
                     "config_overrides):\n  ftir_image_width_um = %.0f, ftir_image_height_um = %.0f"),
              cols * W_um, rows * H_um))
  invisible(list(table = tab, csv = csv_out, png = png_out,
                 width_um = cols * W_um, height_um = rows * H_um, on_particle = mean(hit)))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 1) stop("Usage: Rscript tools/stitch_ftir_segments.R <folder> [output folder]")
  stitch_ftir_segments(a[1], if (length(a) >= 2) a[2] else a[1])
}
