#!/usr/bin/env Rscript
# Run the project's test suite:
#   Rscript run_tests.R
#
# Requires: testthat, RANN, clue (the packages the tested modules touch).
# The suite sources the R/ modules directly — no installation step needed.

if (!requireNamespace("testthat", quietly = TRUE))
  stop("testthat is required: install.packages('testthat')")

library(testthat)
res <- test_dir(
  file.path("tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE
)
