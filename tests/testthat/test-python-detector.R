# Runs the Python detector's own tests (inst/python/test_particle_detector.py:
# synthetic LDIR-like images with known particle positions) with the same
# interpreter the pipeline uses, so `Rscript run_tests.R` and CI cover them.

test_that("Python particle detector tests pass", {
  skip_if_not_installed("reticulate")
  py_ok <- tryCatch({
    sys.source(file.path(REPO_ROOT, "R", "utils_python.R"), envir = environment())
    suppressMessages(setup_python_detector()); TRUE
  }, error = function(e) FALSE)
  skip_if_not(py_ok, "Python with numpy/scipy/Pillow not available")
  py <- reticulate::py_config()$python
  out <- suppressWarnings(system2(py, c("-I", "-m", "unittest", "discover", "-s",
                                        shQuote(file.path(REPO_ROOT, "inst", "python")),
                                        "-p", shQuote("test_*.py")),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  expect(status == 0L, paste(c("Python tests failed:", tail(out, 30)), collapse = "\n"))
  expect_true(any(grepl("^OK", out)))
})
