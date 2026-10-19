# =============================================================================
# read_bmp.R — Dependency-free BMP reader (base R only)
# =============================================================================
# Instrument PCs frequently export microscope images as Windows bitmaps, and
# the 'magick' package is not always installable on the machine that runs the
# pipeline or the Shiny viewer.  This reader covers the variants those exports
# actually use: uncompressed 8-bit palette, 24-bit and 32-bit BMPs (BI_RGB,
# plus BI_BITFIELDS with the standard BGRA channel masks), in both bottom-up
# and top-down row order.
#
# Not supported (returns NULL): RLE-compressed BMPs, 1/4/16-bit depths,
# BITMAPCOREHEADER files, and non-standard bitfield masks.  Callers should
# fall back to magick for those.
#
# Sourced by main.R (pipeline: read_image_any fallback) and by
# shiny_app/global.R (viewer: load_image_raster fallback). Also holds
# magick_bitmap_array(), the one place magick pixel data is unpacked.

#' Pixel array from a magick bitmap, as [height, width, channels] integers 0-255.
#'
#' magick::image_data() returns a bitmap stored channel-first ([C, W, H]).
#' Current magick versions define as.integer() on it to return the array
#' already transposed to [H, W, C]; reshaping that result as [C, W, H] again
#' (what the readers here used to do) scrambles every pixel -- the "tiled /
#' colour-distorted arrays" the viewer worked around by reading PNGs with png.
#' Older magick versions return a flat vector in bitmap order; both are handled.
#'
#' @param bm magick bitmap (raw or character array, dim = c(C, W, H))
#' @return Integer array, dim = c(H, W, C)
magick_bitmap_array <- function(bm) {
  d <- dim(bm)
  ints <- if (is.character(bm)) strtoi(bm, base = 16L) else as.integer(bm)
  if (length(dim(ints)) == 3L && all(dim(ints) == d[c(3L, 2L, 1L)]))
    return(ints)
  aperm(array(as.integer(ints), dim = d), c(3L, 2L, 1L))
}


#' Read an uncompressed BMP file as a raster array
#'
#' @param path Path to a .bmp file
#' @return Numeric array in [0,1] with dim = c(height, width, channels),
#'   channels = 3 (RGB) or 4 (RGBA) — the same layout as png::readPNG —
#'   or NULL if the file is not a BMP variant this reader supports.
read_bmp_raster <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !file.exists(path))
    return(NULL)
  n <- file.info(path)$size
  if (is.na(n) || n < 54) return(NULL)
  bytes <- readBin(path, "raw", n = n)
  if (bytes[1] != as.raw(0x42) || bytes[2] != as.raw(0x4D)) return(NULL) # "BM"

  # Little-endian readers; i is the 1-based index of the first byte.
  u16 <- function(i) as.integer(bytes[i]) + 256L * as.integer(bytes[i + 1L])
  u32 <- function(i) sum(as.numeric(bytes[i + 0:3]) * 256^(0:3))
  s32 <- function(i) { v <- u32(i); if (v >= 2^31) v - 2^32 else v }

  data_offset <- u32(11L)   # file offset of the pixel array
  hdr_size    <- u32(15L)   # DIB header size (40 = BITMAPINFOHEADER, ...)
  if (hdr_size < 40) return(NULL)

  width       <- s32(19L)
  height      <- s32(23L)   # negative = top-down row order
  bpp         <- u16(29L)
  compression <- u32(31L)
  colors_used <- u32(47L)

  top_down <- height < 0
  habs     <- abs(height)
  if (width <= 0 || habs == 0 || !(bpp %in% c(8L, 24L, 32L))) return(NULL)

  # BI_RGB (0) is plain; BI_BITFIELDS (3) is accepted only with the standard
  # BGRA masks (the layout GDI and most exporters write for 32-bit).
  # In every header version the masks sit at file offset 54.
  extra_masks <- 0L
  if (compression == 3L) {
    if (n < 66 || bpp != 32L) return(NULL)
    if (u32(55L) != 0xFF0000 || u32(59L) != 0xFF00 || u32(63L) != 0xFF)
      return(NULL)
    if (hdr_size == 40) extra_masks <- 12L  # masks stored after the header
  } else if (compression != 0L) {
    return(NULL)
  }

  bytes_pp <- bpp %/% 8L
  stride   <- ((width * bytes_pp + 3L) %/% 4L) * 4L  # rows pad to 4 bytes
  if (n < data_offset + stride * habs) return(NULL)

  px <- as.integer(bytes[(data_offset + 1L):(data_offset + stride * habs)])
  # 0-based offset of each pixel's first byte within px: [row, col]
  idx <- outer((seq_len(habs) - 1L) * stride,
               (seq_len(width) - 1L) * bytes_pp, "+")

  a <- NULL
  if (bpp == 8L) {
    pal_n     <- if (colors_used > 0) colors_used else 256L
    pal_start <- 14L + hdr_size + extra_masks
    if (n < pal_start + pal_n * 4L) return(NULL)
    pal <- matrix(as.integer(bytes[pal_start + seq_len(pal_n * 4L)]),
                  nrow = 4L)                      # rows: B, G, R, reserved
    ind <- px[idx + 1L] + 1L
    if (any(ind > pal_n)) return(NULL)
    r <- matrix(pal[3L, ind], habs, width)
    g <- matrix(pal[2L, ind], habs, width)
    b <- matrix(pal[1L, ind], habs, width)
  } else {
    b <- matrix(px[idx + 1L], habs, width)
    g <- matrix(px[idx + 2L], habs, width)
    r <- matrix(px[idx + 3L], habs, width)
    if (bpp == 32L) a <- matrix(px[idx + 4L], habs, width)
  }

  if (!top_down) {
    flip <- rev(seq_len(habs))
    r <- r[flip, , drop = FALSE]
    g <- g[flip, , drop = FALSE]
    b <- b[flip, , drop = FALSE]
    if (!is.null(a)) a <- a[flip, , drop = FALSE]
  }

  # 32-bit BI_RGB files routinely leave the 4th byte at 0 ("unused"); an
  # all-zero alpha channel means opaque, not fully transparent.
  if (!is.null(a) && any(a > 0L)) {
    array(c(r, g, b, a), dim = c(habs, width, 4L)) / 255
  } else {
    array(c(r, g, b), dim = c(habs, width, 3L)) / 255
  }
}


#' Read an uncompressed BMP downsized, strip by strip from disk.
#'
#' Block-averages k x k pixels (k = the smallest integer that brings the long
#' side to <= max_dim) while reading a few MB of rows at a time, so a
#' field-sample micrograph of 1-2 GB never has to fit in memory -- neither in
#' R nor in ImageMagick, whose resource limits refuse such files ("cache
#' resources exhausted"). Edge blocks average the pixels they have, so the
#' result covers the whole image.
#'
#' @param path Path to a 24- or 32-bit uncompressed (BI_RGB) BMP
#' @param max_dim Longest side of the result (px)
#' @param chunk_mb Rows read per step (MB)
#' @return Numeric RGB array in [0,1], dim c(height, width, 3), with
#'   attr(, "resized_from_px") = c(width, height) of the file; NULL for a BMP
#'   variant this reader does not support
read_bmp_downsized <- function(path, max_dim = 2000L, chunk_mb = 32) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  h <- readBin(con, "raw", n = 54L)
  if (length(h) < 54L || h[1] != as.raw(0x42) || h[2] != as.raw(0x4D)) return(NULL)
  u16 <- function(i) as.integer(h[i]) + 256L * as.integer(h[i + 1L])
  u32 <- function(i) sum(as.numeric(h[i + 0:3]) * 256^(0:3))
  s32 <- function(i) { v <- u32(i); if (v >= 2^31) v - 2^32 else v }
  data_offset <- u32(11L); width <- s32(19L); height <- s32(23L)
  bpp <- u16(29L); compression <- u32(31L)
  if (width <= 0 || height == 0 || !(bpp %in% c(24L, 32L)) || compression != 0) return(NULL)
  top_down <- height < 0; habs <- abs(height)
  bytes_pp <- bpp %/% 8L
  stride <- ((width * bytes_pp + 3L) %/% 4L) * 4L
  if (file.info(path)$size < data_offset + stride * habs) return(NULL)

  k <- max(1L, as.integer(ceiling(max(width, habs) / max_dim)))
  out_w <- as.integer(ceiling(width / k)); out_h <- as.integer(ceiling(habs / k))
  col_block <- (seq_len(width) - 1L) %/% k + 1L
  acc <- array(0, dim = c(out_w, out_h, 3L))
  cnt <- matrix(0, out_w, out_h)
  col_n <- tabulate(col_block, out_w)

  seek(con, data_offset)
  rows_per_chunk <- max(1L, as.integer(chunk_mb * 2^20 %/% stride))
  j <- 0L                                   # 0-based file row
  while (j < habs) {
    nr <- min(rows_per_chunk, habs - j)
    raw <- readBin(con, "raw", n = stride * nr)
    m <- matrix(as.integer(raw), nrow = stride)          # one column per file row
    img_row <- if (top_down) j + seq_len(nr) - 1L else habs - 1L - (j + seq_len(nr) - 1L)
    blk <- img_row %/% k + 1L
    for (ch in 1:3) {                                     # file order B, G, R
      v <- m[seq(4L - ch, width * bytes_pp, by = bytes_pp), , drop = FALSE]
      s <- rowsum(v, col_block)                           # out_w x nr
      s <- t(rowsum(t(s), blk))                           # out_w x blocks in chunk
      b <- as.integer(colnames(s))
      acc[, b, ch] <- acc[, b, ch] + s
    }
    nb <- tabulate(blk, out_h)
    cnt <- cnt + outer(col_n, nb)
    j <- j + nr
  }
  out <- array(0, dim = c(out_h, out_w, 3L))
  for (ch in 1:3) out[, , ch] <- t(acc[, , ch] / cnt) / 255
  attr(out, "resized_from_px") <- c(width, habs)
  out
}


#' Width and height of a BMP from its header (NULL if not a BMP).
bmp_dimensions <- function(path) {
  h <- tryCatch(readBin(path, "raw", n = 26L), error = function(e) raw(0))
  if (length(h) < 26L || h[1] != as.raw(0x42) || h[2] != as.raw(0x4D)) return(NULL)
  s32 <- function(i) { v <- sum(as.numeric(h[i + 0:3]) * 256^(0:3)); if (v >= 2^31) v - 2^32 else v }
  c(width = s32(19L), height = abs(s32(23L)))
}
