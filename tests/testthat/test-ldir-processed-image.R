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

# ---------------------------------------------------------------------------
# apply_ldir_coord_swaps()
# ---------------------------------------------------------------------------

.joined_like <- function() {
  data.frame(
    particle_id      = c("A19", "A20", "A28", "A29"),
    area_um2         = c(10062, 10025, 1375, 1287),   # Excel-intrinsic (must NOT move)
    material         = c("PET", "PET", "unknown", "unknown"),
    x_um             = c(100, 200, 300, 400),
    y_um             = c(11, 22, 33, 44),
    image_area_um2   = c(10000, 9900, 1300, 1250),
    image_feret_um   = c(130, 117, 45, 44),
    coord_match_cost = c(0.1, 0.2, 0.3, 0.4),
    coord_source     = "processed_image",
    stringsAsFactors = FALSE
  )
}

test_that("apply_ldir_coord_swaps exchanges only the image-assigned fields", {
  df  <- .joined_like()
  out <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("A19", "A29"))))

  a19 <- out[out$particle_id == "A19", ]
  a29 <- out[out$particle_id == "A29", ]
  # Coordinate + image fields swapped...
  expect_equal(a19$x_um, 400); expect_equal(a19$y_um, 44)
  expect_equal(a29$x_um, 100); expect_equal(a29$y_um, 11)
  expect_equal(a19$image_area_um2, 1250)
  expect_equal(a29$image_area_um2, 10000)
  expect_equal(a19$coord_match_cost, 0.4)
  # ...but Excel-intrinsic fields stay put.
  expect_equal(a19$area_um2, 10062)
  expect_equal(a29$area_um2, 1287)
  expect_equal(a19$material, "PET")
})

test_that("apply_ldir_coord_swaps rotates a 3-cycle (directed, order-independent)", {
  df  <- .joined_like()   # x_um: A19=100, A20=200, A28=300, A29=400
  # Each row c(id_clarity, id_R): id_clarity receives id_R's coordinate.
  cyc <- list(c("A19", "A20"), c("A20", "A28"), c("A28", "A19"))
  out <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = cyc))
  expect_equal(out$x_um[out$particle_id == "A19"], 200)   # <- A20
  expect_equal(out$x_um[out$particle_id == "A20"], 300)   # <- A28
  expect_equal(out$x_um[out$particle_id == "A28"], 100)   # <- A19
  expect_equal(out$x_um[out$particle_id == "A29"], 400)   # untouched

  # Order-independent: shuffling the rows gives the same result.
  out2 <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = rev(cyc)))
  expect_equal(out2$x_um, out$x_um)
})

test_that("a single row auto-completes into a swap", {
  df  <- .joined_like()
  out <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("A19", "A29"))))
  expect_equal(out$x_um[out$particle_id == "A19"], 400)   # <- A29
  expect_equal(out$x_um[out$particle_id == "A29"], 100)   # <- A19 (chain closed)
})

test_that("apply_ldir_coord_swaps aborts an ambiguous repeated-column instruction", {
  df <- .joined_like()
  # A19 named as id_clarity in two rows -> ambiguous -> abort, df unchanged.
  out <- apply_ldir_coord_swaps(
    df, list(ldir_coord_swaps = list(c("A19", "A20"), c("A19", "A28"))))
  expect_identical(out, df)
})

test_that("apply_ldir_coord_swaps is a no-op / warns on empty or bad input", {
  df <- .joined_like()
  expect_identical(apply_ldir_coord_swaps(df, list(ldir_coord_swaps = NULL)), df)
  expect_identical(apply_ldir_coord_swaps(df, list()), df)
  # Unknown id -> skipped, data unchanged.
  out <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("A19", "NOPE"))))
  expect_identical(out, df)
  # Malformed entry (not a pair) -> skipped.
  out2 <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("A19", "A20", "A28"))))
  expect_identical(out2, df)
})

# ---------------------------------------------------------------------------
# Coord-swaps CSV sidecar (auto-discovery + loading)
# ---------------------------------------------------------------------------

test_that("load_ldir_coord_swaps reads the sidecar next to the LDIR Excel", {
  d       <- withr::local_tempdir()
  xlsx    <- file.path(d, "Sample.xlsx")            # need not exist for lookup
  sidecar <- file.path(d, "Sample_coord_swaps.csv")
  writeLines(c("id_a,id_b,note",
               "A19,A29,near-tie",
               "A20,A28,"), sidecar)

  cfg   <- make_config()
  pairs <- load_ldir_coord_swaps(xlsx, cfg)
  expect_type(pairs, "list")
  expect_equal(length(pairs), 2L)
  expect_equal(pairs[[1]], c("A19", "A29"))   # `note` column ignored
  expect_equal(pairs[[2]], c("A20", "A28"))

  # And it round-trips through apply_ldir_coord_swaps.
  df  <- .joined_like()
  cfg$ldir_coord_swaps <- pairs
  out <- apply_ldir_coord_swaps(df, cfg)
  expect_equal(out$x_um[out$particle_id == "A19"], 400)
})

test_that("load_ldir_coord_swaps returns NULL when absent or disabled", {
  d    <- withr::local_tempdir()
  xlsx <- file.path(d, "NoSidecar.xlsx")
  expect_null(load_ldir_coord_swaps(xlsx, make_config()))     # no file

  # Suffix disabled.
  sidecar <- file.path(d, "NoSidecar_coord_swaps.csv")
  writeLines(c("id_a,id_b", "A1,A2"), sidecar)
  cfg <- make_config(); cfg$ldir_coord_swaps_suffix <- NULL
  expect_null(load_ldir_coord_swaps(xlsx, cfg))
})

test_that("load_ldir_coord_swaps accepts id_clarity/id_R headers", {
  d       <- withr::local_tempdir()
  xlsx    <- file.path(d, "Named.xlsx")
  sidecar <- file.path(d, "Named_coord_swaps.csv")
  writeLines(c("id_clarity,id_R,note", "A4,A3,", "A18,A17,keep"), sidecar)

  pairs <- suppressWarnings(load_ldir_coord_swaps(xlsx, make_config()))
  expect_equal(length(pairs), 2L)
  expect_equal(pairs[[1]], c("A4", "A3"))
  expect_equal(pairs[[2]], c("A18", "A17"))
})

test_that("load_ldir_coord_swaps falls back to the first two columns without id_a/id_b headers", {
  d       <- withr::local_tempdir()
  xlsx    <- file.path(d, "Fallback.xlsx")
  sidecar <- file.path(d, "Fallback_coord_swaps.csv")
  writeLines(c("from,to", "A5,A7", " , ", "A8,A9"), sidecar)   # blank row dropped

  pairs <- load_ldir_coord_swaps(xlsx, make_config())
  expect_equal(length(pairs), 2L)
  expect_equal(pairs[[1]], c("A5", "A7"))
  expect_equal(pairs[[2]], c("A8", "A9"))
})

test_that("load_ldir_coord_swaps auto-detects ';' and tab separators", {
  # A European-locale Excel writes ";"-delimited CSV. Read with sep="," that
  # collapses to a single column, and the whole sidecar was silently dropped —
  # the pipeline then ran with NO corrections at all.
  for (sep in c(";", "\t")) {
    d       <- withr::local_tempdir()
    xlsx    <- file.path(d, "Sep.xlsx")
    sidecar <- file.path(d, "Sep_coord_swaps.csv")
    writeLines(c(paste("id_clarity", "id_R", "note", sep = sep),
                 paste("A4", "A3", "", sep = sep),
                 paste("A18", "A17", "", sep = sep)), sidecar)

    pairs <- load_ldir_coord_swaps(xlsx, make_config())
    expect_equal(length(pairs), 2L, info = sep)
    expect_equal(pairs[[1]], c("A4", "A3"), info = sep)
    expect_equal(pairs[[2]], c("A18", "A17"), info = sep)
  }
})

test_that("load_ldir_coord_swaps strips a UTF-8 BOM from the header", {
  d       <- withr::local_tempdir()
  xlsx    <- file.path(d, "Bom.xlsx")
  sidecar <- file.path(d, "Bom_coord_swaps.csv")
  con <- file(sidecar, open = "wb")
  writeBin(charToRaw("﻿id_clarity,id_R\nA4,A3\n"), con)
  close(con)

  pairs <- load_ldir_coord_swaps(xlsx, make_config())
  expect_equal(length(pairs), 1L)
  expect_equal(pairs[[1]], c("A4", "A3"))
})

test_that("the full 10-row correction table applies as one global permutation", {
  # Regression for the reported 37/38/39 case: the table must relabel
  # 37 -> 39, 39 -> 38, 38 -> 37 (not a 37<->38 swap leaving 39 alone).
  ids <- paste0("A", 1:40)
  df  <- data.frame(particle_id = ids, x_um = seq_along(ids) * 100,
                    y_um = seq_along(ids) * 100, stringsAsFactors = FALSE)
  swaps <- list(c("A4","A3"),  c("A39","A37"), c("A38","A39"), c("A37","A38"),
                c("A19","A29"), c("A28","A20"), c("A18","A17"),
                c("A10","A9"), c("A9","A8"),   c("A6","A7"))

  out   <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = swaps))
  # owner[i] = the id whose ORIGINAL coordinate now sits on row i
  owner <- paste0("A", out$x_um / 100)
  relabel <- stats::setNames(out$particle_id, owner)   # R_old -> R_new

  expect_equal(unname(relabel[c("A37", "A39", "A38")]), c("A39", "A38", "A37"))
  # The two independent near-tie pairs stay independent 2-swaps.
  expect_equal(unname(relabel[c("A29", "A19")]), c("A19", "A29"))
  expect_equal(unname(relabel[c("A20", "A28")]), c("A28", "A20"))
  # The 8/9/10 open chain closes into a 3-cycle.
  expect_equal(unname(relabel[c("A9", "A8", "A10")]), c("A10", "A9", "A8"))
  # Everything not named in the table is untouched.
  named <- paste0("A", c(3,4,6,7,8,9,10,17,18,19,20,28,29,37,38,39))
  expect_true(all(owner[!ids %in% named] == ids[!ids %in% named]))
  # Order-independent.
  out2 <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = rev(swaps)))
  expect_equal(out2$x_um, out$x_um)
})

test_that("apply_ldir_coord_swaps tolerates case / separator / bare-number ids", {
  df  <- .joined_like()
  out <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("a19", " A_29 "))))
  expect_equal(out$x_um[out$particle_id == "A19"], 400)
  expect_equal(out$x_um[out$particle_id == "A29"], 100)

  out2 <- apply_ldir_coord_swaps(df, list(ldir_coord_swaps = list(c("19", "29"))))
  expect_equal(out2$x_um[out2$particle_id == "A19"], 400)
})
