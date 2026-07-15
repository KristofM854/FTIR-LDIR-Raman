# B4: try_or() must surface failures in the log, not swallow them silently.

test_that("try_or returns the value on success without logging", {
  expect_silent(res <- try_or(sqrt(16), default = -1, what = "sqrt"))
  expect_equal(res, 4)
})

test_that("try_or logs a warning and returns the default on error", {
  expect_message(
    res <- try_or(stop("boom"), default = -1, what = "deliberate step"),
    "deliberate step failed: boom"
  )
  expect_equal(res, -1)
})
