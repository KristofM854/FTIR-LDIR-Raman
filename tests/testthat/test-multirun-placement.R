# Multi-Run image-placement helpers: each must reproduce the corresponding
# single-instrument tab's extent from metadata. Pure functions -> unit-testable.

.load_placement_helpers <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  env <- new.env(parent = globalenv())
  txt <- paste(readLines(gpath, warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  env
}

test_that("FTIR/Bruker placement spans the raw particle extent", {
  env <- .load_placement_helpers()
  ext <- env$place_image_particle_extent(c(10, 50, 30, NA), c(-5, 20, 8, 100))
  # NA y pairs with NA-free x? both filtered independently by is.finite.
  expect_equal(ext$xmin, 10); expect_equal(ext$xmax, 50)
  expect_equal(ext$ymin, -5); expect_equal(ext$ymax, 100)
  expect_null(env$place_image_particle_extent(numeric(0), numeric(0)))
})

test_that("LDIR placement mirrors the scan-circle formula", {
  env <- .load_placement_helpers()
  meta <- data.frame(ldir_cx_px = 100, ldir_cy_px = 100,
                     ldir_scale_um_per_px = 2, ldir_image_width_px = 200,
                     ldir_image_height_px = 200, stringsAsFactors = FALSE)
  ext <- env$place_image_ldir_meta(meta)
  # xmin=-cx*s, xmax=(w-cx)*s, ymin=(cy-h)*s, ymax=cy*s
  expect_equal(ext$xmin, -200); expect_equal(ext$xmax, 200)
  expect_equal(ext$ymin, -200); expect_equal(ext$ymax, 200)
  # Missing calibration -> NULL (falls back in the viewer).
  expect_null(env$place_image_ldir_meta(data.frame(instrument = "ldir")))
})

test_that("Raman placement uses WITec extent with Y auto-detection", {
  env <- .load_placement_helpers()
  meta <- data.frame(raman_image_width_um = 1000, raman_image_height_um = 1000,
                     raman_image_center_x_um = 500, raman_image_center_y_um = 500,
                     stringsAsFactors = FALSE)
  # Particles cluster at positive Y -> the un-negated interpretation wins.
  x <- c(200, 400, 600, 800); y <- c(200, 400, 600, 800)
  ext <- env$place_image_raman_meta(meta, x, y)
  expect_equal(ext$xmin, 0);    expect_equal(ext$xmax, 1000)
  expect_equal(ext$ymin, 0);    expect_equal(ext$ymax, 1000)
  # No metadata -> NULL.
  expect_null(env$place_image_raman_meta(data.frame(instrument = "raman"), x, y))
})

test_that("Raman µm-per-pixel tier scales the image by the raster dims", {
  env <- .load_placement_helpers()
  meta <- data.frame(raman_um_per_px = 10, stringsAsFactors = FALSE)
  raw  <- array(0, dim = c(100, 200, 3))          # nrow=100, ncol=200
  ext  <- env$place_image_raman_umpx(meta, c(0, 1000), c(0, 500), raw)
  # centre = particle mean (500,250); half = px * upp / 2
  expect_equal(ext$xmin, -500); expect_equal(ext$xmax, 1500)   # 200*10/2 = 1000
  expect_equal(ext$ymin, -250); expect_equal(ext$ymax, 750)    # 100*10/2 = 500
  # No scale and no TIFF -> NULL (caller falls back to particle-extent fit).
  expect_null(env$place_image_raman_umpx(data.frame(a = 1), c(0, 1), c(0, 1), raw))
})

test_that("Raman dispatch prefers WITec (P1) then µm/px (P2)", {
  env <- .load_placement_helpers()
  raw <- array(0, dim = c(100, 200, 3))
  witec <- data.frame(raman_image_width_um = 1000, raman_image_height_um = 1000,
                      raman_image_center_x_um = 500, raman_image_center_y_um = 500)
  p1 <- env$place_image_multirun("raman", witec, c(200, 800), c(200, 800), raw)
  expect_equal(p1$xmax, 1000)                      # WITec extent, not µm/px
  umpx <- data.frame(raman_um_per_px = 10)
  p2 <- env$place_image_multirun("raman", umpx, c(0, 1000), c(0, 500), raw)
  expect_equal(p2$xmax, 1500)                      # falls through to µm/px
})

test_that("dispatch routes each instrument to the right placement", {
  env <- .load_placement_helpers()
  x <- c(0, 100); y <- c(0, 100)
  expect_equal(env$place_image_multirun("ftir_perkin", NULL, x, y)$xmax, 100)
  expect_equal(env$place_image_multirun("ftir_bruker", NULL, x, y)$xmax, 100)
  expect_null(env$place_image_multirun("nonsense", NULL, x, y))
})
