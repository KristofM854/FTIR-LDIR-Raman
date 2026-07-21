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

test_that("dispatch routes each instrument to the right placement", {
  env <- .load_placement_helpers()
  x <- c(0, 100); y <- c(0, 100)
  expect_equal(env$place_image_multirun("ftir_perkin", NULL, x, y)$xmax, 100)
  expect_equal(env$place_image_multirun("ftir_bruker", NULL, x, y)$xmax, 100)
  expect_null(env$place_image_multirun("nonsense", NULL, x, y))
})
