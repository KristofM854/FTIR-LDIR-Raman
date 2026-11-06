# LDIR particle outlines from the Clarity project file (.sam): the reader
# (R/ldir_sam.R) and the exact join. The .sam here is synthetic, written in
# the record layout decoded from a real file (header comment of ldir_sam.R).

local({
  for (m in c("01c_ingest_ldir.R", "ftir_image_placement.R", "ldir_sam.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
})

.sam_record <- function(id, x, y, area, w, h) {
  i32 <- function(v) writeBin(as.integer(v), raw(), size = 4, endian = "little")
  f64 <- function(v) writeBin(as.double(v), raw(), size = 8, endian = "little")
  f32 <- function(v) writeBin(as.double(v), raw(), size = 4, endian = "little")
  c(as.raw(0x31), charToRaw("Clarity.Models"), charToRaw(LDIR_SAM_PARTICLE_CLASS),
    i32(8), as.raw(sample(0:255, 16, TRUE)), as.raw(nchar(id)), charToRaw(id),
    i32(1234), f64(41.6), as.raw(rep(0, 8)), i32(length(x)), i32(0),
    f32(as.vector(rbind(x, y))), f64(area), f64(w), f64(h), as.raw(sample(0:255, 40, TRUE)))
}
# Rectangle outline (boundary pixels), x0..x1 by y0..y1, y counted from the bottom
.rect <- function(x0, x1, y0, y1) list(x = c(x0:x1, rep(x1, y1 - y0), x1:x0, rep(x0, y1 - y0)),
                                       y = c(rep(y0, x1 - x0 + 1), (y0 + 1):y1, rep(y1, x1 - x0 + 1), y1:(y0 + 1)))

test_that("the .sam reader decodes ID, outline and sizes; repeated IDs start a new group", {
  set.seed(1)
  r1 <- .rect(337, 467, 1884, 1920); r2 <- .rect(1580, 1632, 2108, 2156)
  b <- c(as.raw(sample(0:255, 500, TRUE)),
         .sam_record("A1", r1$x, r1$y, 39300, 650, 180), as.raw(sample(0:255, 300, TRUE)),
         .sam_record("A2", r2$x, r2$y, 30437.5, 178, 291),
         .sam_record("A1", r1$x, r1$y, 39300, 650, 180))           # a second analysis
  f <- withr::local_tempfile(fileext = ".sam"); writeBin(b, f)
  s <- read_ldir_sam_particles(f, chunk_mb = 0.0005)               # records across chunk borders
  expect_equal(s$particle_id, c("A1", "A2", "A1"))
  expect_equal(s$group, c(1L, 1L, 2L))
  expect_equal(c(s$x_min[1], s$x_max[1], s$y_min[1], s$y_max[1]), c(337, 467, 1884, 1920))
  expect_equal(s$area_um2[1:2], c(39300, 30437.5))
  expect_equal(c(s$width_um[1], s$height_um[1]), c(650, 180))
  expect_equal(s$x_px[1], 402, tolerance = 0.01)
  expect_equal(s$y_px[1], 1902, tolerance = 0.01)
})

test_that("the exact join places each row on the blob of its own outline", {
  set.seed(2)
  H <- 2605
  sam <- data.frame(group = 1L, particle_id = c("A1", "A2", "A3"), x_px = c(402, 1606, 900),
                    y_px = c(1902, 2132, 400), area_um2 = c(39300, 30437, 2000))
  ex <- data.frame(particle_id = c("A1", "A2", "A3"), area_um2 = c(39300, 30437, 2000))
  # Blobs in 1-based column / row from the top, shuffled, plus a decoy of A2's size
  img <- data.frame(centroid_px_x = c(1607, 403, 901, 1200), centroid_px_y = H - c(2132, 1902, 400, 1000),
                    area_um2 = c(31000, 41000, 2100, 30000), feret_max_um = 1)
  img$x_um <- img$centroid_px_x * 5; img$y_um <- -img$centroid_px_y * 5
  j <- suppressMessages(join_ldir_by_sam(ex, sam, img, H))
  expect_equal(j$x_um, c(403, 1607, 901) * 5)
  expect_true(all(j$coord_source == "sam"))
  # A .sam of another scan (areas do not match): not used
  bad <- sam; bad$area_um2 <- bad$area_um2 * 3
  expect_null(suppressMessages(join_ldir_by_sam(ex, bad, img, H)))
})

test_that("the exact join returns the outlines in the frame of x_um / y_um", {
  set.seed(3)
  H <- 2605
  r1 <- .rect(337, 467, 1884, 1920)
  b <- c(as.raw(sample(0:255, 200, TRUE)), .sam_record("A1", r1$x, r1$y, 39300, 650, 180))
  f <- withr::local_tempfile(fileext = ".sam"); writeBin(b, f)
  s <- read_ldir_sam_particles(f)
  expect_equal(nrow(attr(s, "outlines")), length(r1$x))
  # Blobs with a known pixel -> um map: x = 5 * col - 100, y = -5 * row + 7000
  img <- data.frame(centroid_px_x = c(403, 1000, 2000), centroid_px_y = H - c(1902, 300, 1500),
                    area_um2 = c(41000, 500, 500), feret_max_um = 1)
  img$x_um <- 5 * img$centroid_px_x - 100; img$y_um <- -5 * img$centroid_px_y + 7000
  j <- suppressMessages(join_ldir_by_sam(data.frame(particle_id = "A1", area_um2 = 39300), s, img, H))
  ol <- attr(j, "outlines")
  expect_equal(range(ol$x_um), 5 * (c(337, 467) + 1) - 100)
  expect_equal(range(ol$y_um), sort(-5 * (H - c(1920, 1884)) + 7000))
})
