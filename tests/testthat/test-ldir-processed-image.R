# =============================================================================
# test-ldir-processed-image.R — LDIR analyzed (processed) particle-overlay image
# =============================================================================
# Covers the processed-image coordinate path added to 01c_ingest_ldir.R:
#   * find_ldir_processed_image()           — locate the "_analyzed" overlay
#   * extract_ldir_processed_image_coords() — segment coloured blobs -> centroids
# and documents the fallback contract (NULL -> optical extraction is used).
#
# helper-setup.R sources utils.R + 00_config.R (for read_image_any / make_config /
# log_message); the module under test is sourced here into globalenv so its
# functions resolve those helpers lexically.

local({
  rdir <- file.path(REPO_ROOT, "R")
  sys.source(file.path(rdir, "01c_ingest_ldir.R"), envir = globalenv())
})

# R < 4.4 has no base `%||%`; the sourced module relies on it. Provide it
# globally only when missing (a no-op on the R-release CI runner).
if (!exists("%||%", mode = "function")) {
  assign("%||%", function(a, b) if (is.null(a)) b else a, envir = globalenv())
}

# --- Fixtures ----------------------------------------------------------------

# Write a synthetic RGBA-free PNG: black background with solid coloured blobs.
# `blobs` is a list of list(rows=, cols=, channel=) where channel in 1:3 (RGB).
.write_synthetic_overlay <- function(path, h = 100L, w = 100L, blobs = list()) {
  arr <- array(0, dim = c(h, w, 3L))
  for (b in blobs) arr[b$rows, b$cols, b$channel] <- 1
  png::writePNG(arr, path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# find_ldir_processed_image()
# ---------------------------------------------------------------------------

test_that("find_ldir_processed_image locates the suffixed overlay next to the optical image", {
  d <- withr::local_tempdir()
  optical <- file.path(d, "sampleA.png")
  proc    <- file.path(d, "sampleA_analyzed.png")
  file.create(optical, proc)

  cfg <- list(ldir_processed_image_suffix = "_analyzed")
  found <- find_ldir_processed_image(optical, cfg)

  expect_false(is.null(found))
  expect_equal(normalizePath(found), normalizePath(proc))
})

test_that("find_ldir_processed_image matches any image extension for the overlay", {
  d <- withr::local_tempdir()
  optical <- file.path(d, "run1.png")
  proc    <- file.path(d, "run1_analyzed.tif")   # different extension than optical
  file.create(optical, proc)

  found <- find_ldir_processed_image(optical, list(ldir_processed_image_suffix = "_analyzed"))
  expect_equal(normalizePath(found), normalizePath(proc))
})

test_that("find_ldir_processed_image returns NULL when no overlay exists (fallback path)", {
  d <- withr::local_tempdir()
  optical <- file.path(d, "lonely.png")
  file.create(optical)   # only the optical image is present

  expect_null(find_ldir_processed_image(optical, list(ldir_processed_image_suffix = "_analyzed")))
})

test_that("find_ldir_processed_image returns NULL when the suffix is NULL (disabled)", {
  d <- withr::local_tempdir()
  optical <- file.path(d, "sampleB.png")
  proc    <- file.path(d, "sampleB_analyzed.png")
  file.create(optical, proc)

  expect_null(find_ldir_processed_image(optical, list(ldir_processed_image_suffix = NULL)))
})

test_that("find_ldir_processed_image honours a custom suffix from make_config", {
  d <- withr::local_tempdir()
  optical <- file.path(d, "sampleC.png")
  proc    <- file.path(d, "sampleC_analyzed.png")
  file.create(optical, proc)

  # make_config() default suffix is "_analyzed"
  cfg <- make_config()
  expect_equal(cfg$ldir_processed_image_suffix, "_analyzed")
  expect_equal(normalizePath(find_ldir_processed_image(optical, cfg)),
               normalizePath(proc))
})

# ---------------------------------------------------------------------------
# extract_ldir_processed_image_coords()
# ---------------------------------------------------------------------------

test_that("extract_ldir_processed_image_coords segments coloured blobs into centroids", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "overlay.png")
  # Green blob: rows/cols 20:30 (centre ~25,25); Blue blob: rows 60:70, cols 70:80.
  .write_synthetic_overlay(path, 100L, 100L, list(
    list(rows = 20:30, cols = 20:30, channel = 2L),   # green
    list(rows = 60:70, cols = 70:80, channel = 3L)    # blue
  ))

  # scale = 1 µm/px keeps the geometry easy to reason about; coordinates are
  # circle-centred (origin = image centre 50,50) and y increases upward.
  cfg <- list(ldir_processed_image_min_brightness = 30L,
              ldir_min_blob_area_px               = 5L,
              ldir_image_scale_um_per_px          = 1.0)

  res <- extract_ldir_processed_image_coords(path, config = cfg)
  df  <- res$particles

  # Returns particles + a circle_info calibration object.
  expect_type(res, "list")
  expect_true(all(c("particles", "circle_info") %in% names(res)))
  expect_true(isTRUE(res$circle_info$detected))
  expect_equal(res$circle_info$scale_um_per_px, 1)

  # Required column set (superset OK) + one row per blob
  req <- c("x_um", "y_um", "area_um2", "feret_max_um", "major_um",
           "minor_um", "aspect_ratio", "coord_source")
  expect_true(all(req %in% names(df)))
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2L)
  expect_true(all(df$coord_source == "processed_image"))

  # Circle-centred coords (image centre = 50,50; y up), scale = 1.
  # Green centroid (col~25, row~25) -> x = 25-50 = -25, y = 50-25 = +25.
  # Blue  centroid (col~75, row~65) -> x = 75-50 = +25, y = 50-65 = -15.
  green_row <- which.min((df$x_um - (-25))^2 + (df$y_um - 25)^2)
  blue_row  <- which.min((df$x_um - 25)^2   + (df$y_um - (-15))^2)
  expect_true(green_row != blue_row)
  expect_equal(df$x_um[green_row], -25, tolerance = 2)
  expect_equal(df$y_um[green_row],  25, tolerance = 2)
  expect_equal(df$x_um[blue_row],   25, tolerance = 2)
  expect_equal(df$y_um[blue_row],  -15, tolerance = 2)

  # Pixel centroids preserved for diagnostics/overlay.
  expect_equal(df$centroid_px_x[green_row], 25, tolerance = 2)
  expect_equal(df$centroid_px_y[green_row], 25, tolerance = 2)

  # Area ~ pixel count (11 x 11 = 121 px) with unit scale.
  expect_equal(df$area_um2[green_row], 121, tolerance = 5)
  expect_equal(df$area_um2[blue_row],  121, tolerance = 5)

  # Sizes are finite and positive.
  expect_true(all(df$feret_max_um > 0))
  expect_true(all(df$minor_um > 0))
  expect_true(all(is.finite(df$aspect_ratio)))
})

test_that("extract_ldir_processed_image_coords applies the µm/px scale factor", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "scaled.png")
  .write_synthetic_overlay(path, 100L, 100L, list(
    list(rows = 40:50, cols = 40:50, channel = 2L)   # single green blob
  ))

  df1 <- extract_ldir_processed_image_coords(
    path, config = list(ldir_image_scale_um_per_px = 1.0))$particles
  df2 <- extract_ldir_processed_image_coords(
    path, config = list(ldir_image_scale_um_per_px = 2.0))$particles

  expect_equal(nrow(df1), 1L)
  expect_equal(nrow(df2), 1L)
  # Linear measures scale by 2, area by 4 (coords are centred, so they also
  # scale linearly about the origin).
  expect_equal(df2$x_um,         df1$x_um * 2,         tolerance = 1e-6)
  expect_equal(df2$y_um,         df1$y_um * 2,         tolerance = 1e-6)
  expect_equal(df2$feret_max_um, df1$feret_max_um * 2, tolerance = 1e-6)
  expect_equal(df2$area_um2,     df1$area_um2 * 4,     tolerance = 1e-6)
  # Aspect ratio is scale-invariant.
  expect_equal(df2$aspect_ratio, df1$aspect_ratio, tolerance = 1e-6)
})

test_that("touching same-colour blobs stay separate via dominant-channel masks", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "touch.png")
  # A green and a blue square sharing an edge (cols 50/51) must remain 2 blobs
  # because colour quantization processes each channel independently.
  .write_synthetic_overlay(path, 100L, 100L, list(
    list(rows = 40:60, cols = 30:50, channel = 2L),   # green
    list(rows = 40:60, cols = 51:70, channel = 3L)    # blue (adjacent)
  ))

  df <- extract_ldir_processed_image_coords(path, config = list())$particles
  expect_equal(nrow(df), 2L)
})

test_that("extract_ldir_processed_image_coords returns an empty frame for a blank image", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "blank.png")
  .write_synthetic_overlay(path, 40L, 40L, list())   # all black

  res <- extract_ldir_processed_image_coords(path, config = list())
  df  <- res$particles
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  expect_true(all(c("x_um", "y_um", "coord_source") %in% names(df)))
  # circle_info is still returned so callers never see NA/NULL placement fields.
  expect_true(is.numeric(res$circle_info$cx_px))
  expect_true(isTRUE(res$circle_info$scale_um_per_px > 0))
})
