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

test_that("brightness segmentation keeps a white-cored blob whole (no drop/duplicate)", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "whitecore.png")
  # One particle: blue body with a WHITE core. Under dominant-channel masks the
  # white core (R-dominant) and blue body (B-dominant) split into two co-located
  # fragments; brightness mode must keep it as a single blob.
  arr <- array(0, dim = c(80L, 80L, 3L))
  arr[30:50, 30:50, 3] <- 0.85          # blue body
  arr[30:50, 30:50, 2] <- 0.55
  arr[37:43, 37:43, 1:3] <- 1           # white core
  png::writePNG(arr, path)

  # Brightness (default): one blob.
  res_b <- extract_ldir_processed_image_coords(path, config = list())
  expect_equal(nrow(res_b$particles), 1L)
  # Mean colour recorded for future material use.
  expect_true(all(c("mean_r", "mean_g", "mean_b") %in% names(res_b$particles)))
  expect_gt(res_b$particles$mean_b[1], res_b$particles$mean_r[1])   # bluish

  # Legacy color-channel mode splits the same particle into >1 fragment.
  res_c <- extract_ldir_processed_image_coords(
    path, config = list(ldir_processed_segmentation = "color_channel"))
  expect_gt(nrow(res_c$particles), 1L)
})

test_that("color_channel mode keeps touching different-colour blobs separate", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "touch.png")
  # A green and a blue square sharing an edge (cols 50/51).
  .write_synthetic_overlay(path, 100L, 100L, list(
    list(rows = 40:60, cols = 30:50, channel = 2L),   # green
    list(rows = 40:60, cols = 51:70, channel = 3L)    # blue (adjacent)
  ))

  # Legacy per-channel mode separates them by colour...
  df_c <- extract_ldir_processed_image_coords(
    path, config = list(ldir_processed_segmentation = "color_channel"))$particles
  expect_equal(nrow(df_c), 2L)

  # ...while brightness mode (default) merges the touching pair into one blob.
  df_b <- extract_ldir_processed_image_coords(path, config = list())$particles
  expect_equal(nrow(df_b), 1L)
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

# ---------------------------------------------------------------------------
# Shape descriptors on processed blobs
# ---------------------------------------------------------------------------

test_that("extract_ldir_processed_image_coords reports invariant shape descriptors", {
  skip_if_not_installed("png")

  d    <- withr::local_tempdir()
  path <- file.path(d, "shapes.png")
  # A near-square blob (low eccentricity) and a thin bar (high eccentricity).
  .write_synthetic_overlay(path, 100L, 120L, list(
    list(rows = 20:40, cols = 20:40, channel = 2L),   # square
    list(rows = 60:62, cols = 20:100, channel = 2L)   # thin horizontal bar
  ))

  df <- extract_ldir_processed_image_coords(path, config = list())$particles
  expect_equal(nrow(df), 2L)
  expect_true(all(c("eccentricity", "circularity", "solidity") %in% names(df)))

  # Descriptors on sane scales.
  expect_true(all(df$eccentricity >= 0 & df$eccentricity <= 1))
  expect_true(all(df$circularity  >  0 & df$circularity  <= 1))
  expect_true(all(df$solidity     >  0 & df$solidity     <= 1))

  # Identify the bar by its larger feret; it must be far more eccentric and
  # less circular than the square.
  bar    <- which.max(df$feret_max_um)
  square <- setdiff(seq_len(2), bar)
  expect_gt(df$eccentricity[bar], df$eccentricity[square])
  expect_gt(df$eccentricity[bar], 0.9)
  expect_lt(df$circularity[bar],  df$circularity[square])
  # Both blobs are convex -> high solidity.
  expect_gt(min(df$solidity), 0.9)
})

# ---------------------------------------------------------------------------
# relabel_ldir_low_hqi()
# ---------------------------------------------------------------------------

.ldir_ingest_like <- function() {
  data.frame(
    particle_id = paste0("A", 1:5),
    material    = c("Polyethylene terephtalate (PET)", "Polyamide (PA)",
                    "Tire", "Polyethylene terephthalate", "Polyether sulphone (PES)"),
    quality     = c(0.97, 0.84, 0.68, 0.90, 0.647),
    feret_max_um = c(50, 40, 30, 20, 10),
    stringsAsFactors = FALSE
  )
}

test_that("relabel_ldir_low_hqi relabels sub-threshold particles as unknown, keeps all rows", {
  df  <- .ldir_ingest_like()
  out <- relabel_ldir_low_hqi(df, list(ldir_hqi_unknown_threshold = 0.85))

  expect_equal(nrow(out), 5L)                       # nothing dropped
  # quality < 0.85 -> "unknown"; A2 (0.84), A3 (0.68), A5 (0.647).
  expect_equal(out$material,
               c("Polyethylene terephtalate (PET)", "unknown", "unknown",
                 "Polyethylene terephthalate", "unknown"))
  # Raw identification preserved.
  expect_equal(out$identification_raw, df$material)
})

test_that("relabel_ldir_low_hqi is a no-op when the threshold is NULL or 0", {
  df <- .ldir_ingest_like()
  expect_equal(relabel_ldir_low_hqi(df, list(ldir_hqi_unknown_threshold = NULL))$material,
               df$material)
  expect_equal(relabel_ldir_low_hqi(df, list(ldir_hqi_unknown_threshold = 0))$material,
               df$material)
})

# ---------------------------------------------------------------------------
# join_ldir_coords() shape-fingerprint disambiguation
# ---------------------------------------------------------------------------

test_that("shape fingerprint resolves a swap that size/rank alone cannot", {
  # Two particles with IDENTICAL size (so area/feret/rank are uninformative)
  # but distinct shapes. The image blobs are listed in REVERSED order, so only
  # the shape fingerprint can recover the correct assignment.
  excel <- data.frame(
    particle_id  = c("P1", "P2"),
    x_um = NA_real_, y_um = NA_real_,
    area_um2     = c(1000, 1000),
    feret_max_um = c(50, 50),
    major_um     = c(50, 50), minor_um = c(25, 25),
    aspect_ratio = c(2, 2),
    eccentricity = c(0.10, 0.90),   # P1 round-ish, P2 elongated
    circularity  = c(0.90, 0.50),
    solidity     = c(0.98, 0.80),
    stringsAsFactors = FALSE
  )
  image <- data.frame(
    particle_id  = c("B1", "B2"),
    x_um = c(100, 200), y_um = c(0, 0),   # B1 has P2's shape, B2 has P1's shape
    area_um2     = c(1000, 1000),
    feret_max_um = c(50, 50),
    major_um     = c(50, 50), minor_um = c(25, 25),
    eccentricity = c(0.90, 0.10),
    circularity  = c(0.50, 0.90),
    solidity     = c(0.80, 0.98),
    coord_source = "processed_image",
    stringsAsFactors = FALSE
  )

  cfg <- make_config(); cfg$ldir_join_weight_shape <- 1   # shape opt-in (default 0)
  joined <- join_ldir_coords(excel, image, config = cfg)

  # P1 (round) must take B2 (x=200); P2 (elongated) must take B1 (x=100).
  expect_equal(joined$x_um[joined$particle_id == "P1"], 200)
  expect_equal(joined$x_um[joined$particle_id == "P2"], 100)
})

test_that("area-based rank prevents the elongated-particle swap (A1/A2 regression)", {
  # A1: large AREA, moderate feret. A2: smaller area but LARGER feret (elongated,
  # like the real A2 with aspect 0.34). Ranking the image by feret while the
  # Excel is area-sorted put A1's row on A2's blob. With rank_metric="area"
  # (default) both sides order by area, so identity is preserved.
  excel <- data.frame(
    particle_id  = c("A1", "A2"),
    x_um = NA_real_, y_um = NA_real_,
    area_um2     = c(44100, 34775),   # A1 larger area
    feret_max_um = c(349, 428),       # ...but A2 larger feret
    major_um     = c(349, 428), minor_um = c(217, 147),
    eccentricity = c(0.52, 0.95),
    stringsAsFactors = FALSE
  )
  image <- data.frame(
    particle_id  = c("B_A2", "B_A1"),               # deliberately shuffled
    x_um = c(200, 100), y_um = c(0, 0),             # B_A1 (x=100) is A1's blob
    area_um2     = c(34775, 44100),
    feret_max_um = c(428, 349),
    major_um     = c(428, 349), minor_um = c(147, 217),
    eccentricity = c(0.95, 0.52),
    coord_source = "processed_image",
    stringsAsFactors = FALSE
  )

  joined <- join_ldir_coords(excel, image, config = make_config())
  # A1 must take its own blob (x=100), not A2's (x=200).
  expect_equal(joined$x_um[joined$particle_id == "A1"], 100)
  expect_equal(joined$x_um[joined$particle_id == "A2"], 200)
})

test_that("size gate stops shape from overriding a clear size difference", {
  # P1 (big) and P2 (small) differ ~26% in area, so size alone matches them
  # correctly. Their shapes are crossed so the shape term *wants* to swap them.
  # w_rank/w_ar are zeroed to isolate the size-gate vs shape interaction.
  excel <- data.frame(
    particle_id  = c("P1", "P2"),
    x_um = NA_real_, y_um = NA_real_,
    area_um2     = c(4400, 3400),
    feret_max_um = c(55, 48),
    major_um     = c(55, 48), minor_um = c(27, 24),
    eccentricity = c(0.10, 0.90),
    stringsAsFactors = FALSE
  )
  image <- data.frame(
    particle_id  = c("B1", "B2"),
    x_um = c(100, 200), y_um = c(0, 0),
    area_um2     = c(4400, 3400),     # B1 big, B2 small (matches P1, P2 by size)
    feret_max_um = c(55, 48),
    major_um     = c(55, 48), minor_um = c(27, 24),
    eccentricity = c(0.90, 0.10),     # ...but shapes are crossed
    coord_source = "processed_image",
    stringsAsFactors = FALSE
  )

  cfg <- make_config()
  cfg$ldir_join_weight_shape <- 1   # shape opt-in (default 0)
  cfg$ldir_join_weight_rank <- 0
  cfg$ldir_join_weight_ar   <- 0

  # Gate ON (default): size wins -> P1 keeps the big blob B1 (x=100).
  gated <- join_ldir_coords(excel, image, config = cfg)
  expect_equal(gated$x_um[gated$particle_id == "P1"], 100)
  expect_equal(gated$x_um[gated$particle_id == "P2"], 200)

  # Gate OFF: the crossed shape overrides size and forces the wrong swap,
  # demonstrating the gate is what prevents it.
  cfg_nogate <- cfg
  cfg_nogate$ldir_join_shape_size_gate <- 0
  ungated <- join_ldir_coords(excel, image, config = cfg_nogate)
  expect_equal(ungated$x_um[ungated$particle_id == "P1"], 200)
})

test_that("shape term uses eccentricity by default, not circularity/solidity", {
  # Identical sizes (gate = 1) so only shape decides. Eccentricity and
  # circularity point at OPPOSITE assignments; the default must follow
  # eccentricity, and an explicit circularity override must follow circularity.
  excel <- data.frame(
    particle_id  = c("P1", "P2"),
    x_um = NA_real_, y_um = NA_real_,
    area_um2     = c(1000, 1000), feret_max_um = c(50, 50),
    major_um     = c(50, 50), minor_um = c(25, 25),
    eccentricity = c(0.10, 0.90),
    circularity  = c(0.90, 0.10),
    stringsAsFactors = FALSE
  )
  image <- data.frame(
    particle_id  = c("B1", "B2"),
    x_um = c(100, 200), y_um = c(0, 0),
    area_um2     = c(1000, 1000), feret_max_um = c(50, 50),
    major_um     = c(50, 50), minor_um = c(25, 25),
    eccentricity = c(0.10, 0.90),   # by ecc: B1 ~ P1, B2 ~ P2 -> P1:x=100
    circularity  = c(0.10, 0.90),   # by circ: B2 ~ P1        -> P1:x=200
    coord_source = "processed_image",
    stringsAsFactors = FALSE
  )
  cfg <- make_config(); cfg$ldir_join_weight_shape <- 1
  cfg$ldir_join_weight_rank <- 0; cfg$ldir_join_weight_ar <- 0

  # Shape descriptor default (eccentricity): P1 -> B1 (x=100).
  d1 <- join_ldir_coords(excel, image, config = cfg)
  expect_equal(d1$x_um[d1$particle_id == "P1"], 100)

  # Force circularity: opposite assignment, P1 -> B2 (x=200).
  cfg_circ <- cfg; cfg_circ$ldir_join_shape_descriptors <- c("circularity")
  d2 <- join_ldir_coords(excel, image, config = cfg_circ)
  expect_equal(d2$x_um[d2$particle_id == "P1"], 200)
})

test_that("join_ldir_coords is unaffected when image blobs lack shape descriptors", {
  # Optical-path image_df (no eccentricity/circularity/solidity): the shape
  # term must stay inert and matching proceeds on size/rank as before.
  excel <- data.frame(
    particle_id  = c("P1", "P2"),
    x_um = NA_real_, y_um = NA_real_,
    area_um2     = c(4000, 1000),
    feret_max_um = c(80, 40),
    major_um     = c(80, 40), minor_um = c(40, 20),
    aspect_ratio = c(2, 2),
    eccentricity = c(0.1, 0.9), circularity = c(0.9, 0.5), solidity = c(0.98, 0.8),
    stringsAsFactors = FALSE
  )
  image <- data.frame(
    particle_id  = c("B1", "B2"),
    x_um = c(100, 200), y_um = c(0, 0),
    area_um2     = c(4000, 1000),
    feret_max_um = c(80, 40),
    major_um     = c(80, 40), minor_um = c(40, 20),
    coord_source = "circle_calibrated",
    stringsAsFactors = FALSE
  )
  joined <- join_ldir_coords(excel, image, config = make_config())
  # Larger P1 -> larger blob B1 (x=100); smaller P2 -> B2 (x=200).
  expect_equal(joined$x_um[joined$particle_id == "P1"], 100)
  expect_equal(joined$x_um[joined$particle_id == "P2"], 200)
})
