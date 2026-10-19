# Image readers must return the pixels in [height, width, channel] order with
# the values the file holds. read_image_any() used to scramble every pixel when
# reading through magick (as.integer() on a magick bitmap is already H x W x C).
# A non-square, asymmetric test image catches transposition, which a square one
# can hide.

.test_png <- function() {
  set.seed(5)
  img <- array(runif(30 * 50 * 3), dim = c(30, 50, 3))
  img[1:5, 1:10, 1] <- 1; img[1:5, 1:10, 2:3] <- 0   # red block top-left
  p <- tempfile(fileext = ".png"); png::writePNG(img, p)
  list(path = p, img = png::readPNG(p))
}

test_that("read_image_any returns the PNG's pixels in H x W x C order", {
  skip_if_not_installed("magick")
  t <- .test_png()
  a <- suppressMessages(read_image_any(t$path, verbose = FALSE))
  expect_identical(dim(a)[1:2], c(30L, 50L))
  expect_equal(a[, , 1:3], t$img, tolerance = 1e-6)
})

test_that("magick_bitmap_array matches png::readPNG for RGB and RGBA", {
  skip_if_not_installed("magick")
  t <- .test_png()
  m <- magick::image_read(t$path)
  for (ch in c("rgb", "rgba")) {
    a <- magick_bitmap_array(magick::image_data(m, channels = ch))
    expect_identical(dim(a), c(30L, 50L, nchar(ch)))
    expect_equal(a[, , 1:3] / 255, t$img, tolerance = 1e-6)
  }
})

test_that("the preview downscaler keeps orientation", {
  skip_if_not_installed("magick")
  t <- .test_png()
  img <- array(as.integer(round(t$img * 255)), dim = dim(t$img))
  img <- img[rep(1:30, each = 4), rep(1:50, each = 4), ]   # 120 x 200
  prev <- make_preview_image(img, max_dim = 100L)
  expect_identical(dim(prev)[1:2], c(60L, 100L))
  expect_gt(mean(prev[1:8, 1:16, 1]), 200)                # red block stays top-left
  expect_lt(mean(prev[1:8, 1:16, 3]), 60)
})

test_that("the viewer reads large run images at display size, keeping the orientation", {
  skip_if_not_installed("magick")
  env <- new.env(parent = globalenv())
  withr::with_dir(file.path(REPO_ROOT, "shiny_app"), suppressWarnings(suppressMessages(
    source("global.R", local = env, encoding = "UTF-8"))))
  t <- .test_png()                                     # 30 x 50, red block top-left
  big <- suppressMessages(env$load_image_raster(t$path, max_dim = 25))
  expect_identical(dim(big)[1:2], c(15L, 25L))         # aspect kept
  expect_equal(attr(big, "resized_from_px"), c(50, 30))
  expect_gt(big[1, 1, 1], 0.8); expect_lt(big[1, 1, 2], 0.2)   # still red top-left
  # Placement code sees the original size and full coverage
  cov <- .raster_coverage(big)
  expect_equal(c(cov$orig_w, cov$orig_h, cov$fx, cov$fy), c(50, 30, 1, 1))
  # Small enough: read as before
  expect_identical(dim(suppressMessages(env$load_image_raster(t$path, max_dim = 100)))[1:2],
                   c(30L, 50L))
})
