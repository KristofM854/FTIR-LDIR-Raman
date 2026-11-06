# =============================================================================
# ldir_sam.R -- LDIR particle positions from the Clarity project file (.sam)
# =============================================================================
# The LDIR Excel export has no coordinates, so join_ldir_coords() attaches the
# table rows to the particles of the analyzed image by size, shape and colour
# class. The Clarity project file (.sam, Agilent Clarity 1.7) stores each
# particle with its ID and outline, so the join can be exact: each row takes
# the analyzed-image particle its own outline covers.
#
# Format (decoded from Faanui SW1L R4 A, Clarity 1.7.27): a sequence of
# serialized objects, each introduced by its assembly and class name as
# length-prefixed strings. A particle record, after the class name
# "Agilent.Clarity.Models.Data.Analysis.LDIRParticle":
#   int32   field count
#   16 B    GUID
#   string  particle ID ("A1"; 7-bit length prefix)
#   int32, double, 8 B   (not used)
#   int64   number of outline points n
#   n x (float32 x, float32 y)  outline, in pixels of the IR scan image (the
#           grid of the analyzed image; 0-based column, y counted from the
#           BOTTOM row)
#   double  area (um^2), double width (um), double height (um)
# A1: outline x 337-467, y 1884-1920 px (130 x 36 px = 650 x 180 um at 5
# um/px), area 39300 -- the Excel values; the analyzed image has A1 at the
# same pixels. The file also holds the optical and IR images and each
# particle's image crop and spectrum, which are not read here.
# =============================================================================

LDIR_SAM_PARTICLE_CLASS <- "Agilent.Clarity.Models.Data.Analysis.LDIRParticle"

#' Read the particle records of a Clarity .sam file.
#'
#' Scans the file in chunks for the particle class name and decodes each
#' record. A file may hold more than one particle analysis; a new group starts
#' where the IDs restart (an ID seen before in the current group).
#'
#' @param path .sam file
#' @param chunk_mb Read size; records are re-read in full where they cross a
#'   chunk boundary
#' @return data.frame(group, particle_id, n_points, x_px, y_px (outline
#'   centroid, 0-based, y from the bottom), x_min, x_max, y_min, y_max,
#'   area_um2, width_um, height_um), or NULL when no record is found
read_ldir_sam_particles <- function(path, chunk_mb = 64) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  size <- file.info(path)$size
  pat <- charToRaw(LDIR_SAM_PARTICLE_CLASS)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  chunk <- as.integer(chunk_mb * 1024^2); overlap <- length(pat)
  hits <- numeric(0); pos <- 0
  while (pos < size) {
    seek(con, pos)
    b <- readBin(con, "raw", n = min(chunk + overlap, size - pos))
    h <- grepRaw(pat, b, fixed = TRUE, all = TRUE)
    if (length(h)) hits <- c(hits, pos + h - 1)
    pos <- pos + chunk
  }
  hits <- unique(hits)
  if (!length(hits)) return(NULL)
  rd <- function(at, n) { seek(con, at); readBin(con, "raw", n = n) }
  int32 <- function(r, i) readBin(r[i:(i + 3)], "integer", size = 4, endian = "little")
  dbl <- function(r, i) readBin(r[i:(i + 7)], "double", size = 8, endian = "little")
  recs <- vector("list", length(hits))
  for (k in seq_along(hits)) {
    p0 <- hits[k] + length(pat)
    r <- rd(p0, 512)
    i <- 1 + 4 + 16                                   # field count, GUID
    len <- 0L; shift <- 0L                            # 7-bit encoded string length
    repeat { byte <- as.integer(r[i]); i <- i + 1L
             len <- len + bitwAnd(byte, 127L) * 2^shift; shift <- shift + 7L
             if (byte < 128L) break }
    if (len <= 0 || len > 64) next
    pid <- rawToChar(r[i:(i + len - 1)]); i <- i + len
    i <- i + 4 + 8 + 8                                # int32, double, 8 bytes
    n <- int32(r, i) + 2^32 * int32(r, i + 4); i <- i + 8
    if (!is.finite(n) || n < 1 || n > 1e6) next
    rr <- rd(p0 + i - 1, 8 * n + 24)
    if (length(rr) < 8 * n + 24) next
    xy <- readBin(rr[seq_len(8 * n)], "double", size = 4, n = 2 * n, endian = "little")
    x <- xy[c(TRUE, FALSE)]; y <- xy[c(FALSE, TRUE)]
    tail <- rr[(8 * n + 1):(8 * n + 24)]
    # Area-weighted centroid of the outline (shoelace); vertex mean for a line
    x2 <- c(x[-1], x[1]); y2 <- c(y[-1], y[1]); cr <- x * y2 - x2 * y; A <- sum(cr) / 2
    cx <- if (abs(A) > 0.5) sum((x + x2) * cr) / (6 * A) else mean(x)
    cy <- if (abs(A) > 0.5) sum((y + y2) * cr) / (6 * A) else mean(y)
    recs[[k]] <- data.frame(particle_id = pid, n_points = n, x_px = cx, y_px = cy,
                            x_min = min(x), x_max = max(x), y_min = min(y), y_max = max(y),
                            area_um2 = dbl(tail, 1), width_um = dbl(tail, 9),
                            height_um = dbl(tail, 17), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, recs)
  if (is.null(out) || !nrow(out)) return(NULL)
  # Groups: a new analysis starts where an ID repeats
  g <- 1L; seen <- character(0); grp <- integer(nrow(out))
  for (k in seq_len(nrow(out))) {
    if (out$particle_id[k] %in% seen) { g <- g + 1L; seen <- character(0) }
    seen <- c(seen, out$particle_id[k]); grp[k] <- g
  }
  out$group <- grp
  out
}


#' The .sam file of an LDIR run: config$ldir_sam_path, else the only .sam
#' file next to the LDIR table. NULL when there is none or several.
find_ldir_sam <- function(config) {
  p <- config$ldir_sam_path
  if (!is.null(p) && nzchar(p)) return(if (file.exists(p)) p else NULL)
  if (is.null(config$ldir_path)) return(NULL)
  f <- list.files(dirname(config$ldir_path), pattern = "\\.sam$", full.names = TRUE,
                  ignore.case = TRUE)
  if (length(f) == 1) f else NULL
}


#' Exact LDIR join from the .sam outlines.
#'
#' Picks the particle group of the .sam whose IDs and areas match the Excel
#' table, then gives each Excel row the analyzed-image blob its outline
#' centroid falls on (nearest blob centroid within `max_px`, area within a
#' factor 2). Returns the table with the columns join_ldir_coords() writes
#' (coord_source "sam"), or NULL when the file does not fit the table.
#'
#' @param excel_df LDIR rows (particle_id, area_um2)
#' @param sam read_ldir_sam_particles() result
#' @param image_df analyzed-image blobs (centroid_px_x / centroid_px_y in
#'   1-based image columns / rows from the top, x_um, y_um, area_um2)
#' @param img_h analyzed image height (px)
#' @param min_share rows that must be positioned for the result to be used
join_ldir_by_sam <- function(excel_df, sam, image_df, img_h, max_px = 4, min_share = 0.95) {
  if (is.null(sam) || !nrow(sam) || is.null(image_df$centroid_px_x) || !is.finite(img_h))
    return(NULL)
  # The group that describes this table: most IDs whose area agrees
  fit <- vapply(split(sam, sam$group), function(s) {
    k <- match(excel_df$particle_id, s$particle_id)
    sum(!is.na(k) & abs(s$area_um2[k] - excel_df$area_um2) <= 0.01 * pmax(excel_df$area_um2, 1),
        na.rm = TRUE)
  }, numeric(1))
  best <- as.integer(names(fit)[which.max(fit)])
  if (max(fit) < min_share * nrow(excel_df)) {
    log_message(sprintf("  .sam particles do not fit the LDIR table (best group: %d of %d rows with ID and area)",
                        max(fit), nrow(excel_df)), level = "WARN")
    return(NULL)
  }
  s <- sam[sam$group == best, , drop = FALSE]
  k <- match(excel_df$particle_id, s$particle_id)
  col <- s$x_px[k] + 1; row <- img_h - s$y_px[k]      # to 1-based column / row from the top
  ok <- is.finite(col) & is.finite(row)
  nn <- .nn1(cbind(image_df$centroid_px_x, image_df$centroid_px_y), cbind(col[ok], row[ok]))
  j <- rep(NA_integer_, nrow(excel_df)); d <- rep(NA_real_, nrow(excel_df))
  j[ok] <- nn$idx; d[ok] <- nn$dist
  ar <- image_df$area_um2[j] / pmax(excel_df$area_um2, 1e-6)
  good <- !is.na(j) & d <= max_px & is.finite(ar) & ar > 0.5 & ar < 2
  good[good] <- !duplicated(j[good])
  if (sum(good) < min_share * nrow(excel_df)) {
    log_message(sprintf("  .sam outlines found a blob for only %d of %d LDIR rows -- not used",
                        sum(good), nrow(excel_df)), level = "WARN")
    return(NULL)
  }
  out <- excel_df
  out$x_um <- ifelse(good, image_df$x_um[j], NA_real_)
  out$y_um <- ifelse(good, image_df$y_um[j], NA_real_)
  out$coord_match_cost <- ifelse(good, 0, NA_real_)
  out$coord_source <- ifelse(good, "sam", "none")
  out$image_area_um2 <- ifelse(good, image_df$area_um2[j], NA_real_)
  out$image_feret_um <- ifelse(good, image_df$feret_max_um[j], NA_real_)
  out$coord_join_warning <- FALSE
  log_message(sprintf(paste0("  Exact LDIR join from the .sam file: %d of %d rows placed by their ",
                             "own outline (group %d of %d; median outline-to-blob distance %.1f px)"),
                      sum(good), nrow(excel_df), best, length(fit), stats::median(d[good])))
  out
}
