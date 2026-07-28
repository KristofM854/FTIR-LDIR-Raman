# View-orientation helpers (shiny_app/global.R): the display-only rotation +
# mirror that brings an instrument's native tab into the Raman orientation.
# Pure geometry -> unit-testable without booting Shiny.

.load_view_helpers <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  env <- new.env(parent = globalenv())
  txt <- paste(readLines(gpath, warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  env
}

.dihedral <- function() {
  out <- list()
  for (fl in c(FALSE, TRUE))
    for (d in c(0L, 90L, -90L, 180L))
      out[[length(out) + 1L]] <- list(deg = d, flip = fl)
  out
}

test_that("each rotation is undone by its inverse", {
  env <- .load_view_helpers()
  x <- c(1, 3, -2, 0); y <- c(2, -1, 5, 7)
  for (d in c(0L, 90L, -90L, 180L)) {
    b <- env$view_transform_xy(x, y, d)
    f <- env$view_transform_xy(b$x, b$y, -d)
    expect_equal(f$x, x, info = paste("deg", d))
    expect_equal(f$y, y, info = paste("deg", d))
  }
})

test_that("the mirror flips Y only and is its own inverse", {
  env <- .load_view_helpers()
  x <- c(1, 3, -2); y <- c(2, -1, 5)
  m <- env$view_transform_xy(x, y, 0L, TRUE)
  expect_equal(m$x, x)
  expect_equal(m$y, -y)
  m2 <- env$view_transform_xy(m$x, m$y, 0L, TRUE)
  expect_equal(m2$x, x); expect_equal(m2$y, y)
})

test_that("an X mirror is a Y mirror plus 180 deg — so one axis suffices", {
  # This is why the UI offers a single "Mirror (flip Y)" checkbox: the four
  # rotations times that one mirror already span all eight orientations.
  env <- .load_view_helpers()
  x <- c(1, 3, -2); y <- c(2, -1, 5)
  a <- env$view_transform_xy(x, y, 180L, TRUE)
  expect_equal(a$x, -x)
  expect_equal(a$y,  y)
})

test_that("the transformed extent still bounds the transformed corners", {
  env <- .load_view_helpers()
  ext <- list(xmin = 0, xmax = 10, ymin = -4, ymax = 6)
  g <- expand.grid(x = c(ext$xmin, ext$xmax), y = c(ext$ymin, ext$ymax))
  for (cd in .dihedral()) {
    p <- env$view_transform_xy(g$x, g$y, cd$deg, cd$flip)
    e <- env$view_transform_extent(ext, cd$deg, cd$flip)
    lbl <- paste("deg", cd$deg, "flip", cd$flip)
    expect_equal(range(p$x), c(e$xmin, e$xmax), info = lbl)
    expect_equal(range(p$y), c(e$ymin, e$ymax), info = lbl)
  }
})

test_that("raster transforms are involutions and preserve storage mode", {
  env <- .load_view_helpers()
  r <- array(seq_len(2 * 3 * 3), dim = c(2, 3, 3))   # integer RGB-ish raster
  expect_identical(env$view_transform_raster(env$view_transform_raster(r, 0L, TRUE),
                                             0L, TRUE), r)
  expect_identical(env$view_transform_raster(env$view_transform_raster(r, 180L),
                                             180L), r)
  expect_identical(Reduce(function(a, .) env$rotate_raster_view(a, 90L), 1:4, r), r)
  # 90 deg swaps the pixel dimensions; a mirror leaves them alone.
  expect_identical(dim(env$view_transform_raster(r, 90L))[1:2], dim(r)[2:1])
  expect_identical(dim(env$view_transform_raster(r, 0L, TRUE)), dim(r))
  # 2D (greyscale) rasters take the same path.
  expect_identical(env$flip_raster_view(env$flip_raster_view(r[, , 1], TRUE), TRUE),
                   r[, , 1])
})

test_that("particles stay on their background pixel through a view transform", {
  # The invariant the feature rests on: rotating/mirroring the scene must move
  # the points, the extent and the raster together, so a particle drawn on a
  # blob is still on that blob afterwards.
  env <- .load_view_helpers()
  nr <- 7; nc <- 11
  img <- list(raster = matrix(seq_len(nr * nc), nrow = nr, ncol = nc),
              xmin = -300, xmax = 800, ymin = -100, ymax = 600)

  pixel_at <- function(im, x, y) {           # raster row 1 is the top (ymax)
    nrw <- nrow(im$raster); ncl <- ncol(im$raster)
    cc <- min(max(floor((x - im$xmin) / (im$xmax - im$xmin) * ncl) + 1, 1), ncl)
    rr <- min(max(floor((im$ymax - y) / (im$ymax - im$ymin) * nrw) + 1, 1), nrw)
    im$raster[rr, cc]
  }

  set.seed(42)
  px <- runif(200, img$xmin, img$xmax)
  py <- runif(200, img$ymin, img$ymax)
  before <- mapply(function(a, b) pixel_at(img, a, b), px, py)

  for (cd in .dihedral()) {
    e <- env$view_transform_extent(img, cd$deg, cd$flip)
    im2 <- list(raster = env$view_transform_raster(img$raster, cd$deg, cd$flip),
                xmin = e$xmin, xmax = e$xmax, ymin = e$ymin, ymax = e$ymax)
    p <- env$view_transform_xy(px, py, cd$deg, cd$flip)
    after <- mapply(function(a, b) pixel_at(im2, a, b), p$x, p$y)
    expect_identical(after, before,
                     info = paste("deg", cd$deg, "flip", cd$flip))
  }
})

test_that("auto_view_dihedral recovers every orientation, translation aside", {
  env <- .load_view_helpers()
  set.seed(1)
  n <- 40; bx <- runif(n, 0, 1000); by <- runif(n, 0, 700)
  for (cd in .dihedral()) {
    p <- env$view_transform_xy(bx, by, cd$deg, cd$flip)
    got <- env$auto_view_dihedral(bx, by, p$x + 250, p$y - 90)
    lbl <- paste("deg", cd$deg, "flip", cd$flip)
    expect_equal(got$deg, cd$deg, info = lbl)
    expect_equal(isTRUE(got$flip), cd$flip, info = lbl)
  }
})

test_that("auto_view_dihedral tells a 180 deg rotation from a Y flip", {
  # The whole reason the mirror exists: on an asymmetric cloud these are
  # different scenes, and no rotation can turn one into the other.
  env <- .load_view_helpers()
  set.seed(7)
  n <- 40; bx <- runif(n, 0, 1000); by <- runif(n, 0, 700)
  rot <- env$view_transform_xy(bx, by, 180L, FALSE)
  mir <- env$view_transform_xy(bx, by, 0L,   TRUE)
  g_rot <- env$auto_view_dihedral(bx, by, rot$x, rot$y)
  g_mir <- env$auto_view_dihedral(bx, by, mir$x, mir$y)
  expect_equal(g_rot$deg, 180L); expect_false(isTRUE(g_rot$flip))
  expect_equal(g_mir$deg, 0L);   expect_true(isTRUE(g_mir$flip))
})

test_that("auto detection falls back to identity when the clouds disagree", {
  env <- .load_view_helpers()
  set.seed(3)
  a <- data.frame(x = runif(30, 0, 1000), y = runif(30, 0, 1000))
  b <- data.frame(x = runif(30, 0, 1000), y = runif(30, 0, 1000))
  got <- env$auto_view_dihedral(a$x, a$y, b$x, b$y)
  expect_equal(got$deg, 0L)
  expect_false(isTRUE(got$flip))
  # Too few points to decide -> identity, not a guess.
  tiny <- env$auto_view_dihedral(c(1, 2), c(1, 2), c(1, 2), c(1, 2))
  expect_equal(tiny$deg, 0L)
  expect_false(isTRUE(tiny$flip))
})

test_that("the LDIR wrapper is unchanged: rotations only, bare integer", {
  env <- .load_view_helpers()
  set.seed(11)
  n <- 30; bx <- runif(n, 0, 900); by <- runif(n, 0, 600)
  for (d in c(0L, 90L, -90L, 180L)) {
    p <- env$rotate_xy_view(bx, by, d)
    got <- env$ldir_auto_view_rotation(bx, by, p$x + 40, p$y + 15)
    expect_identical(got, d, info = paste("deg", d))
  }
  # A mirrored reference is not a rotation — must not be reported as one.
  m <- env$view_transform_xy(bx, by, 0L, TRUE)
  expect_identical(env$ldir_auto_view_rotation(bx, by, m$x, m$y), 0L)
})

test_that("auto detection stays fast on a large particle cloud", {
  # Scoring is O(n_src * n_ref) per candidate over 8 candidates; the helper
  # thins each cloud so a big FTIR run cannot stall the tab.
  env <- .load_view_helpers()
  set.seed(5)
  n <- 4000; bx <- runif(n, 0, 5000); by <- runif(n, 0, 5000)
  p <- env$view_transform_xy(bx, by, 180L, FALSE)
  el <- system.time(got <- env$auto_view_dihedral(bx, by, p$x + 30, p$y - 10))[["elapsed"]]
  expect_equal(got$deg, 180L)
  expect_lt(el, 20)
})
