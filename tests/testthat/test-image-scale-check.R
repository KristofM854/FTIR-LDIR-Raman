# image_scale_from_particle_area() (shiny_app/global.R).
#
# The Raman micrograph is placed from per-dataset WITec values held in
# 00_config.R. Carrying a previous scan's values into a new run draws the image
# at the wrong size, and raman_image_extent_from_config() does not catch it:
# that guard only rejects an extent which fails to CONTAIN the particles, so a
# stale extent that is merely too large still passes.
#
# The analysed particle areas give an independent handle on µm/px — in a
# dark-field image the bright pixels should cover the reported physical area.
# It is coarse (the count moves with the brightness threshold), so it is only
# ever used to flag a gross mismatch.

.load_scale_helper <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  src <- readLines(gpath, warn = FALSE)
  i <- grep("^image_scale_from_particle_area <- function", src)[1]
  skip_if(is.na(i), "helper not found in global.R")
  j <- i + which(src[(i + 1):length(src)] == "}")[1]
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(src[i:j], collapse = "\n")), envir = env)
  env
}

# Dark-field frame: n_bright lit pixels on a dark background.
.fake_image <- function(side = 200, n_bright = 400, channels = 3L) {
  m <- matrix(0.02, side, side)
  m[seq_len(n_bright)] <- 0.9
  if (is.na(channels)) return(m)
  array(rep(m, channels), dim = c(side, side, channels))
}

test_that("the estimate recovers a known scale", {
  env <- .load_scale_helper()
  # 400 bright px representing 3600 µm² total => 3 µm/px.
  img <- .fake_image(n_bright = 400)
  est <- env$image_scale_from_particle_area(img, rep(3600 / 4, 4))
  expect_equal(est, 3, tolerance = 1e-9)
})

test_that("it works on a single-channel image too", {
  env <- .load_scale_helper()
  img <- .fake_image(n_bright = 400, channels = NA)
  expect_equal(env$image_scale_from_particle_area(img, 3600), 3, tolerance = 1e-9)
})

test_that("it takes the brightest channel, not the first", {
  env <- .load_scale_helper()
  side <- 200
  dark <- matrix(0.02, side, side)
  lit  <- dark; lit[seq_len(400)] <- 0.9
  img  <- array(c(dark, dark, lit), dim = c(side, side, 3))   # blue-only signal
  expect_equal(env$image_scale_from_particle_area(img, 3600), 3, tolerance = 1e-9)
})

test_that("it declines when the image cannot support an estimate", {
  env <- .load_scale_helper()
  expect_null(env$image_scale_from_particle_area(NULL, 1000))
  # Nothing bright.
  expect_null(env$image_scale_from_particle_area(.fake_image(n_bright = 0), 1000))
  # Not dark-field: more than 20% lit, so the count means nothing.
  expect_null(env$image_scale_from_particle_area(
    .fake_image(side = 100, n_bright = 4000), 1000))
  # No usable areas.
  expect_null(env$image_scale_from_particle_area(.fake_image(), numeric(0)))
  expect_null(env$image_scale_from_particle_area(.fake_image(), c(NA, NA)))
  expect_null(env$image_scale_from_particle_area(.fake_image(), -5))
})

test_that("a stale per-dataset extent is flagged, ordinary slop is not", {
  # Mirrors the ratio test the viewer applies (warn outside 0.5x .. 2x).
  env <- .load_scale_helper()
  img <- .fake_image(side = 200, n_bright = 400)             # => 3 µm/px truth
  est <- env$image_scale_from_particle_area(img, 3600)
  warns <- function(width_um) {
    r <- (width_um / 200) / est
    r > 2 || r < 0.5
  }
  expect_true(warns(600 * 2.5))    # stale extent 2.5x too large — the real bug
  expect_true(warns(600 / 2.5))    # and the too-small direction
  expect_false(warns(600))         # correct
  expect_false(warns(600 * 1.3))   # threshold-driven slop stays quiet
  expect_false(warns(600 * 0.8))
})
