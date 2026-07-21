# Multi-Run Raman placement (docs/multirun_image_placement_plan.md, plan v3):
# raman_image_width/height/center_x/center_y_um is a PER-DATASET WITec value
# baked into each pipeline run's manifest.json at run time. tools/
# reproducibility.R must resolve it from the matching run's manifest (matched
# by the raw Raman file's MD5), not from whatever happens to currently sit in
# R/00_config.R — hand-copying that value is what went stale and produced a
# ~2.5x wrong image scale in the Multi-Run viewer.

test_that("resolve_raman_calibration_from_manifests matches by raw-file MD5 and picks the latest run", {
  skip_if_not_installed("jsonlite")
  output_root <- withr::local_tempdir()
  raw_dir     <- withr::local_tempdir()

  raw_file <- file.path(raw_dir, "scan.csv")
  writeLines("particle data", raw_file)

  cfg <- make_config()
  cfg$raman_image_width_um    <- 5133.1
  cfg$raman_image_height_um   <- 5068.0
  cfg$raman_image_center_x_um <- 2296.3
  cfg$raman_image_center_y_um <- 2355.7
  write_manifest(file.path(output_root, "2026-07-09_1"), run_id = "2026-07-09_1",
                 config = cfg, input_paths = list(raman = raw_file))

  cfg$raman_image_width_um    <- 5237.9
  cfg$raman_image_height_um   <- 5171.5
  cfg$raman_image_center_x_um <- 2314.0
  cfg$raman_image_center_y_um <- 2354.0
  write_manifest(file.path(output_root, "2026-07-15_1"), run_id = "2026-07-15_1",
                 config = cfg, input_paths = list(raman = raw_file))

  # Force a deterministic run order independent of wall-clock timing between
  # the two write_manifest() calls above.
  bump_timestamp <- function(run_id, ts) {
    p <- file.path(output_root, run_id, "00_manifest", "manifest.json")
    m <- jsonlite::fromJSON(p, simplifyVector = TRUE)
    m$timestamp <- ts
    writeLines(jsonlite::toJSON(m, pretty = TRUE, auto_unbox = TRUE,
                                null = "null", na = "null"), p)
  }
  bump_timestamp("2026-07-09_1", "2026-07-09T11:20:44")
  bump_timestamp("2026-07-15_1", "2026-07-15T10:32:25")

  res <- resolve_raman_calibration_from_manifests(raw_file, output_root)
  expect_equal(res$run_id, "2026-07-15_1")
  expect_equal(res$width_um, 5237.9)
  expect_equal(res$height_um, 5171.5)
  expect_equal(res$center_x_um, 2314.0)
  expect_equal(res$center_y_um, 2354.0)
})

test_that("resolve_raman_calibration_from_manifests ignores runs whose raman input is a different file", {
  skip_if_not_installed("jsonlite")
  output_root <- withr::local_tempdir()
  raw_dir     <- withr::local_tempdir()

  raw_file   <- file.path(raw_dir, "scan_a.csv")
  other_file <- file.path(raw_dir, "scan_b.csv")
  writeLines("scan A content", raw_file)
  writeLines("scan B content, different", other_file)

  cfg <- make_config()
  cfg$raman_image_width_um    <- 9999
  cfg$raman_image_height_um   <- 9999
  cfg$raman_image_center_x_um <- 0
  cfg$raman_image_center_y_um <- 0
  write_manifest(file.path(output_root, "run_b"), run_id = "run_b",
                 config = cfg, input_paths = list(raman = other_file))

  expect_null(resolve_raman_calibration_from_manifests(raw_file, output_root))
})

test_that("resolve_raman_calibration_from_manifests skips runs with an incomplete calibration snapshot", {
  skip_if_not_installed("jsonlite")
  output_root <- withr::local_tempdir()
  raw_dir     <- withr::local_tempdir()
  raw_file <- file.path(raw_dir, "scan.csv")
  writeLines("particle data", raw_file)

  cfg <- make_config()
  cfg$raman_image_width_um  <- 5000
  cfg$raman_image_height_um <- NULL   # incomplete -> dropped from config_snapshot
  write_manifest(file.path(output_root, "run_incomplete"), run_id = "run_incomplete",
                 config = cfg, input_paths = list(raman = raw_file))

  expect_null(resolve_raman_calibration_from_manifests(raw_file, output_root))
})

test_that("resolve_raman_calibration_from_manifests returns NULL for a nonexistent file or output dir", {
  expect_null(resolve_raman_calibration_from_manifests("/no/such/file.csv", withr::local_tempdir()))
  expect_null(resolve_raman_calibration_from_manifests(tempfile(), "/no/such/output/dir"))
})
