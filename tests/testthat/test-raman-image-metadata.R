# WITec image placement file (R/measure_raman_placement.R), downsizing large
# images for the run (canonicalize_instrument_image, R/utils.R), and the
# viewer honouring a stage-frame centre (raman_image_extent_from_config,
# shiny_app/global.R).

local(sys.source(file.path(REPO_ROOT, "R", "measure_raman_placement.R"), envir = globalenv()))

.witec_txt <- function(dir, name, lines = c("Width: 13000", "Height: 13000",
                                            "X: -6011.3134765625", "Y: 12824.1048088074",
                                            "Z: 18.0499992370605", "Rotation: 0")) {
  p <- file.path(dir, name); writeLines(lines, p); p
}

test_that("placement file: top-left corner in the stage frame -> centre", {
  d <- withr::local_tempdir()
  w <- read_witec_image_placement(.witec_txt(d, "img.txt"))
  expect_equal(w$width_um, 13000)
  expect_equal(w$x_left_um, -6011.3134765625)
  expect_equal(w$center_x_um, -6011.3134765625 + 6500)
  expect_equal(w$center_y_um, 12824.1048088074 - 6500)
  expect_equal(w$rotation_deg, 0)
})

test_that("placement file: rotated or incomplete files are ignored", {
  d <- withr::local_tempdir()
  expect_null(suppressMessages(read_witec_image_placement(
    .witec_txt(d, "rot.txt", c("Width: 10", "Height: 10", "X: 0", "Y: 0", "Rotation: 90")))))
  expect_null(suppressMessages(read_witec_image_placement(
    .witec_txt(d, "bad.txt", c("Width: 10", "Height: 10", "X: 0")))))
  expect_message(read_witec_image_placement(.witec_txt(d, "notes.txt", "hello")),
                 "not a WITec image placement file")
})

test_that("placement file is found next to the image, also for a _small copy", {
  d <- withr::local_tempdir()
  txt <- .witec_txt(d, "Sample_B_IAEA_260928.txt")
  file.create(file.path(d, c("Sample_B_IAEA_260928.bmp", "Sample_B_IAEA_260928_small.bmp")))
  expect_equal(find_witec_image_placement(file.path(d, "Sample_B_IAEA_260928.bmp")), txt)
  expect_equal(find_witec_image_placement(file.path(d, "Sample_B_IAEA_260928_small.bmp")), txt)
  # Another name: found by its content (the only placement file in the folder)
  expect_equal(suppressMessages(find_witec_image_placement(file.path(d, "Other.bmp"))), txt)
  # NA switches discovery off; an explicit path is used as given
  expect_null(find_witec_image_placement(file.path(d, "Sample_B_IAEA_260928.bmp"),
                                         list(raman_image_metadata = NA)))
  other <- .witec_txt(d, "elsewhere.txt")
  expect_equal(find_witec_image_placement(file.path(d, "x.bmp"),
                                          list(raman_image_metadata = other)), other)
})

test_that("placement file by content: other .txt files are ignored, the best of several wins", {
  d <- withr::local_tempdir()
  .witec_txt(d, "notes.txt", "measured on Tuesday")                       # not a placement
  .witec_txt(d, "export settings.txt", c("Width: 13000", "Height: 13000", "X: 0"))  # incomplete
  file.create(file.path(d, "Raman image 1.bmp"))
  expect_null(suppressMessages(find_witec_image_placement(file.path(d, "Raman image 1.bmp"))))
  good <- .witec_txt(d, "Sample_A_IAEA.txt")                          # top-left (-6011, 12824)
  expect_equal(suppressMessages(find_witec_image_placement(file.path(d, "Raman image 1.bmp"))), good)
  # Two placement files: the one that puts the particles on the image
  off <- .witec_txt(d, "old export.txt",
                    c("Width: 13000", "Height: 13000", "X: 50000", "Y: 50000", "Rotation: 0"))
  parts <- data.frame(x_um = c(0, 100, -2000), y_um = c(6000, 5000, 3000))
  expect_equal(suppressMessages(find_witec_image_placement(
    file.path(d, "Raman image 1.bmp"), particles = parts)), good)
})

test_that("viewer: a stage-frame centre is used as is, never negated", {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  src <- readLines(gpath, warn = FALSE)
  i <- grep("^raman_image_extent_from_config <- function", src)[1]
  j <- i + which(src[(i + 1):length(src)] == "}")[1]
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(src[i:j], collapse = "\n")), envir = env)
  # Particles near y = 0: both readings of the centre contain all of them.
  x <- c(-10, 0, 10); y <- c(-10, 0, 10)
  cfg <- list(raman_image_width_um = 100, raman_image_height_um = 100,
              raman_image_center_x_um = 0, raman_image_center_y_um = 5)
  expect_true(env$raman_image_extent_from_config(cfg, x, y)$y_negated)   # tie -> negated
  cfg$raman_image_y_frame <- "stage"
  ext <- env$raman_image_extent_from_config(cfg, x, y)
  expect_false(ext$y_negated)
  expect_equal(c(ext$ymin, ext$ymax), c(-45, 55))
})

test_that("large images are downsized for the run and not copied", {
  skip_if_not_installed("magick")
  d <- withr::local_tempdir()
  src <- file.path(d, "big.bmp")
  img <- magick::image_blank(1200, 800, "black")
  magick::image_write(img, src, format = "bmp")
  out <- file.path(d, "inputs")
  info <- suppressMessages(canonicalize_instrument_image(
    src, out, "raman", max_preview_px = 2000L, max_canonical_px = 300,
    copy_original_max_mb = 0.5))
  expect_true(isTRUE(info$canonical_downsized))
  expect_equal(c(info$canonical_width, info$canonical_height), c(300, 200))  # aspect kept
  expect_equal(c(info$orig_width, info$orig_height), c(1200, 800))
  expect_null(info$orig_dest)                         # 2.7 MB original not copied
  expect_false(any(grepl("_original", list.files(out))))
  expect_false(is.na(info$md5))                       # still identified
  # A small image is left alone
  small <- suppressMessages(canonicalize_instrument_image(
    src, file.path(d, "inputs2"), "raman", max_canonical_px = 5000))
  expect_null(small$canonical_downsized)
  expect_equal(small$canonical_width, 1200)
  expect_false(is.null(small$orig_dest))
})

test_that("a downsized image over the size limit is written as JPEG", {
  skip_if_not_installed("magick")
  d <- withr::local_tempdir()
  src <- file.path(d, "noise.png")
  set.seed(1)
  a <- array(stats::runif(1000 * 1000 * 3), c(1000, 1000, 3))
  magick::image_write(magick::image_read(a), src, format = "png")
  info <- suppressMessages(canonicalize_instrument_image(
    src, file.path(d, "inputs"), "raman", max_canonical_px = 900, max_canonical_mb = 0.5))
  expect_match(info$canonical_path, "\\.jpg$")
  expect_lte(file.info(info$canonical_path)$size / 1024^2, 0.5 + 1e-9)
  expect_true(info$canonical_jpeg_quality <= 95)
})
