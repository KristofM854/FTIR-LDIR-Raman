# Regression test on the bundled example dataset (examples/dataset_A): one
# filter, all four instruments. Runs the complete pipeline in a subprocess and
# compares a fingerprint of the run (row counts, every matched pair and its
# distance, transforms, agreement tier counts) with the stored reference
# (tests/testthat/fixtures/dataset_A_fingerprint.csv).
#
# A change that moves any of these numbers fails here. When the change is
# intended, regenerate the reference and review its diff:
#   Rscript examples/run_example_A.R /tmp/exA
#   Rscript tools/regression_fingerprint.R /tmp/exA/<run> tests/testthat/fixtures/dataset_A_fingerprint.csv

test_that("the example dataset reproduces the reference results", {
  skip_if(nzchar(Sys.getenv("SKIP_E2E")), "SKIP_E2E set")
  for (p in c("readxl", "magick", "reticulate", "clue", "RANN", "png", "jsonlite"))
    skip_if_not_installed(p)
  py_ok <- tryCatch({
    sys.source(file.path(REPO_ROOT, "R", "utils_python.R"), envir = environment())
    suppressMessages(setup_python_detector()); TRUE
  }, error = function(e) FALSE)
  skip_if_not(py_ok, "Python with numpy/scipy/Pillow not available")

  out <- withr::local_tempdir()
  rscript <- file.path(R.home("bin"), "Rscript")
  log <- file.path(out, "run.log")
  status <- withr::with_dir(REPO_ROOT, system2(rscript,
    c(file.path("examples", "run_example_A.R"), shQuote(out)),
    stdout = log, stderr = log))
  expect_identical(status, 0L, info = paste(tail(readLines(log), 20), collapse = "\n"))
  run <- list.dirs(out, recursive = FALSE)
  expect_length(run, 1L)

  sys.source(file.path(REPO_ROOT, "tools", "regression_fingerprint.R"), envir = environment())
  got <- regression_fingerprint(run)
  ref <- utils::read.csv(file.path(REPO_ROOT, "tests", "testthat", "fixtures",
                                   "dataset_A_fingerprint.csv"), stringsAsFactors = FALSE)

  key <- function(d) paste(d$kind, d$key)
  expect_setequal(key(got), key(ref))
  m <- merge(ref, got, by = c("kind", "key"), suffixes = c(".ref", ".got"))
  tol <- c(n_rows = 0, agreement = 0, pair = 0.01, transform = 1e-3)
  bad <- m[abs(m$value.ref - m$value.got) > tol[m$kind], ]
  expect(nrow(bad) == 0,
         paste0(nrow(bad), " value(s) differ from the reference:\n",
                paste(utils::capture.output(print(utils::head(bad, 20))), collapse = "\n")))
})
