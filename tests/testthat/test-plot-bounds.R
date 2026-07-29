# Viewport guards (shiny_app/app.R): .brush_bounds() and sanitize_bounds().
#
# Regression for "non-finite location and/or size for viewport" — grid aborts
# the ENTIRE plot (the image appears to vanish) when coord_fixed() is handed a
# limit that is NA/Inf, or a zero-width range under expand = FALSE. The
# reported trigger was a brush drag so small it was effectively a click.

.load_bounds_helpers <- function() {
  apath <- file.path(REPO_ROOT, "shiny_app", "app.R")
  skip_if_not(file.exists(apath), "shiny_app/app.R not found")
  src <- readLines(apath, warn = FALSE)
  i <- grep("^  \\.MIN_SPAN <- ", src)[1]
  j <- grep("^  observeEvent\\(input\\$ftir_brush,", src)[1]
  skip_if(is.na(i) || is.na(j), "viewport guards not found in app.R")
  env <- new.env(parent = globalenv())
  eval(parse(text = gsub("^  ", "", paste(src[i:(j - 1)], collapse = "\n"))),
       envir = env)
  env
}

.finite_box <- function(b) {
  all(is.finite(c(b$x, b$y))) && diff(b$x) > 0 && diff(b$y) > 0
}

test_that("a click-sized brush never becomes a zoom", {
  env <- .load_bounds_helpers()
  expect_null(env$.brush_bounds(list(xmin = 5, xmax = 5, ymin = 2, ymax = 9)))
  expect_null(env$.brush_bounds(list(xmin = 1, xmax = 9, ymin = 2, ymax = 2)))
  expect_null(env$.brush_bounds(list(xmin = 5, xmax = 5, ymin = 2, ymax = 2)))
  expect_null(env$.brush_bounds(NULL))
  expect_null(env$.brush_bounds(list(xmin = NA, xmax = 9, ymin = 2, ymax = 9)))
})

test_that("a genuine brush is passed through unchanged", {
  env <- .load_bounds_helpers()
  expect_identical(env$.brush_bounds(list(xmin = 1, xmax = 9, ymin = 2, ymax = 8)),
                   list(x = c(1, 9), y = c(2, 8)))
})

test_that("sanitize_bounds replaces non-finite limits with the fallback", {
  env <- .load_bounds_helpers()
  # All-NA coordinates make min()/max(na.rm = TRUE) return +Inf / -Inf.
  expect_true(.finite_box(env$sanitize_bounds(list(x = c(Inf, -Inf), y = c(0, 10)))))
  # Image placement metadata that could not be resolved yields NA extents.
  expect_true(.finite_box(env$sanitize_bounds(list(x = c(NA, NA), y = c(0, 10)))))
  expect_true(.finite_box(env$sanitize_bounds(NULL)))
  expect_true(.finite_box(env$sanitize_bounds(list(x = c(0, 10)))))          # y missing
  expect_true(.finite_box(env$sanitize_bounds(list(x = c(1, 2, 3), y = c(0, 1)))))
})

test_that("sanitize_bounds widens a degenerate axis about its centre", {
  env <- .load_bounds_helpers()
  b <- env$sanitize_bounds(list(x = c(500, 500), y = c(0, 10)))
  expect_true(.finite_box(b))
  expect_equal(mean(b$x), 500)     # stays where the user was looking
  expect_equal(b$y, c(0, 10))      # the healthy axis is untouched
})

test_that("sanitize_bounds normalises reversed limits and passes valid ones", {
  env <- .load_bounds_helpers()
  expect_identical(env$sanitize_bounds(list(x = c(9, 1), y = c(8, 2))),
                   list(x = c(1, 9), y = c(2, 8)))
  good <- list(x = c(0, 10000), y = c(0, 10000))
  expect_identical(env$sanitize_bounds(good), good)
})

test_that("sanitize_bounds honours a per-viewer fallback", {
  env <- .load_bounds_helpers()
  ldir_fb <- list(x = c(-7000, 7000), y = c(-7000, 7000))
  expect_identical(env$sanitize_bounds(NULL, ldir_fb), ldir_fb)
})

test_that("the guarded limits actually render where the raw ones abort", {
  # The end-to-end check: this is the exact failure from the bug report.
  skip_if_not_installed("ggplot2")
  env <- .load_bounds_helpers()
  df <- data.frame(x = c(1, 5, 9), y = c(2, 5, 8))
  draw <- function(b) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point() +
         ggplot2::coord_fixed(xlim = b$x, ylim = b$y, expand = FALSE)
    f <- tempfile(fileext = ".png")
    grDevices::png(f, 400, 300)
    on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
    tryCatch({ print(p); TRUE }, error = function(e) FALSE)
  }
  zero_width <- list(x = c(5, 5), y = c(2, 8))
  expect_false(suppressWarnings(draw(zero_width)))              # reproduces it
  expect_true(suppressWarnings(draw(env$sanitize_bounds(zero_width))))
})
