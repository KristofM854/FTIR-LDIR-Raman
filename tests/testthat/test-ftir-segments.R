# Segmented PerkinElmer exports (R/ftir_segments.R, temporary): six segment
# tables and images in a 2 x 3 grid filled column by column, each segment in
# its own frame (y up from its bottom edge), stitched into one table and image.

local({
  sys.source(file.path(REPO_ROOT, "R", "ftir_image_placement.R"), envir = globalenv())
  sys.source(file.path(REPO_ROOT, "R", "ftir_segments.R"), envir = globalenv())
})

.write_segments <- function(dir, um_per_px = 2, W = 40, H = 30) {
  # One particle per segment (none in M3), drawn as a bright green 3 x 3 spot
  # at its local position; local coordinates in um, y up from the bottom edge.
  loc <- data.frame(k = c(1, 2, 4, 5, 6), x = c(10, 20, 30, 50, 70), y = c(10, 40, 20, 30, 50))
  for (k in 1:6) {
    a <- array(0L, c(H, W, 3)); a[, , 3] <- 80L
    p <- loc[loc$k == k, ]
    if (nrow(p)) {
      cx <- round(p$x / um_per_px); ry <- round(H - p$y / um_per_px)
      a[(ry - 1):(ry + 1), (cx - 1):(cx + 1), 2] <- 255L
    }
    magick::image_write(magick::image_read(a / 255), file.path(dir, sprintf("Average Abs.(scan)_M%d.bmp", k)),
                        format = "bmp")
  }
  for (i in seq_len(nrow(loc))) {
    p <- loc[i, ]
    writeLines(c("Identifier,Coord. [pixels],Coord. [µm],Max AAU score,Group,Major dim [µm]",
                 sprintf("MP_0,[%d;%d],[%d;%d],0.8,PET,30", round(p$x / 6.25), round(p$y / 6.25), p$x, p$y)),
               file.path(dir, sprintf("scan M%d.csv", p$k)))
  }
  loc
}

test_that("segments are stitched column by column and shifted into one frame", {
  skip_if_not_installed("magick")
  d <- withr::local_tempdir(); out <- withr::local_tempdir()
  loc <- .write_segments(d)
  tabs <- ftir_segment_files(file.path(d, "scan M1.csv"), FTIR_SEGMENT_TABLE)
  imgs <- ftir_segment_files(file.path(d, "Average Abs.(scan)_M4.bmp"), FTIR_SEGMENT_IMAGE)
  expect_length(tabs, 5); expect_length(imgs, 6)
  st <- suppressWarnings(stitch_ftir_segments(tabs, imgs, out, say = function(...) NULL))
  expect_equal(c(st$width_um, st$height_um), c(160, 180))     # 2 x 40 px, 3 x 30 px at 2 um
  expect_equal(st$on_particle, 1)                             # every particle on its spot
  t <- utils::read.csv(st$csv, check.names = FALSE, encoding = "UTF-8")
  um <- grep("^Coord\\..*m\\]$", names(t), value = TRUE); um <- um[!grepl("pixel", um)]
  xy <- do.call(rbind, lapply(strsplit(gsub("\\[|\\]", "", t[[um]]), ";"), as.numeric))
  # column c = (k-1) %/% 3, row r = (k-1) %% 3 from the top; y up from the bottom
  c0 <- (loc$k - 1) %/% 3; r0 <- (loc$k - 1) %% 3
  expect_equal(xy[, 1], loc$x + c0 * 80)
  expect_equal(xy[, 2], loc$y + (2 - r0) * 60)
  expect_equal(t[[1]], paste0("P_", seq_len(nrow(loc))))     # unique across segments
  expect_equal(t[["Identifier (export)"]], rep("MP_0", nrow(loc)))
  expect_equal(t$Segment, paste0("M", loc$k))
  expect_equal(read_ftir_image_metadata(st$png)$width_um, 160)
})

test_that("run_pipeline's resolver replaces the selected segment by the stitched files", {
  skip_if_not_installed("magick")
  d <- withr::local_tempdir(); run <- withr::local_tempdir()
  .write_segments(d)
  cfg <- list(ftir_path = file.path(d, "scan M2.csv"), ftir_image = NULL, output_dir = run)
  cfg2 <- suppressWarnings(suppressMessages(resolve_ftir_segments(cfg)))
  expect_match(cfg2$ftir_path, "_stitched\\.csv$")
  expect_match(cfg2$ftir_image, "_stitched\\.png$")
  expect_length(cfg2$ftir_segment_images, 6)
  # A single (non-segment) export is left alone
  writeLines("Identifier", file.path(d, "single_export.csv"))
  one <- list(ftir_path = file.path(d, "single_export.csv"), output_dir = run)
  expect_identical(resolve_ftir_segments(one), one)
})
