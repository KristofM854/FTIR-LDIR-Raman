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

test_that("the manifest is written outside a git checkout", {
  # An archived copy (e.g. the Zenodo zip) has no .git: the commit is NA and
  # the manifest must still be written.
  d <- tempfile(); dir.create(d)
  old <- setwd(d); on.exit(setwd(old), add = TRUE)
  expect_true(is.na(current_git_commit()))
  p <- suppressMessages(write_manifest(d, "t", make_config(), list(), stage = "started"))
  m <- jsonlite::fromJSON(p)
  expect_null(m$git_commit)
  expect_identical(m$stage, "started")
})

test_that("the run folder records every parameter and the R session", {
  d <- tempfile(); dir.create(d)
  cfg <- make_config(); cfg$match_dist_threshold_um <- 77
  files <- write_run_environment(d, cfg)
  full <- jsonlite::fromJSON(files[1])
  expect_equal(full$match_dist_threshold_um, 77)
  expect_true(all(c("align_seed", "icp_max_iterations", "match_lambda_area") %in% names(full)))
  expect_true(any(grepl("^R version", readLines(files[2]))))
})

test_that("the config snapshot names exist in the config", {
  # Snapshot keys that are not config keys are silently dropped -- the
  # manifest then claims to record a parameter it never does.
  src <- paste(readLines(file.path(REPO_ROOT, "R", "utils.R")), collapse = "\n")
  block <- regmatches(src, regexpr("cfg_keys <- c\\([^)]*\\)", src))
  keys <- regmatches(block, gregexpr('"[a-z_0-9]+"', block))[[1]]
  keys <- gsub('"', "", keys)
  expect_true(length(keys) > 10)
  expect_identical(setdiff(keys, names(make_config())), character(0))
})
