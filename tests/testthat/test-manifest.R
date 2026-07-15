# Provenance: the run manifest must record config + environment and round-trip
# through JSON (the reproducibility record a Zenodo deposit relies on).

test_that("write_manifest records provenance that round-trips through JSON", {
  skip_if_not_installed("jsonlite")
  run_dir <- withr::local_tempdir()

  cfg  <- make_config()
  path <- write_manifest(run_dir, run_id = "unit_test_run", config = cfg,
                         stage = "started")
  expect_true(file.exists(path))

  m <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  expect_equal(m$run_id, "unit_test_run")
  expect_equal(m$stage, "started")
  expect_true(nzchar(m$r_version))
  expect_true(nzchar(m$platform))

  # Whitelisted config parameters are captured in the snapshot.
  expect_equal(m$config_snapshot$ldir_scan_diameter_um, 13000)
  expect_equal(m$config_snapshot$raman_image_width_um, cfg$raman_image_width_um)
})
