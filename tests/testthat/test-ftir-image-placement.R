# FTIR image placement (R/ftir_image_placement.R): the image must sit under its
# particles to within half a pixel, independent of which particles are shown,
# and the right way up. Scenes come from tools/ftir_synthetic.R, whose geometry
# is known exactly, so placement is checked against ground truth.

sys.source(file.path(REPO_ROOT, "tools", "ftir_synthetic.R"), envir = environment())

# Where a placement draws the image pixel that truly holds each particle,
# minus where the particle is. In um; compare against half a pixel.
.placement_err <- function(pl, sc) {
  tr <- sc$truth; p <- sc$particles
  W <- ncol(sc$image); H <- nrow(sc$image)
  u <- (p$x_um - tr$xmin) / (tr$xmax - tr$xmin) * W      # true edge px
  v <- (tr$ymax - p$y_um) / (tr$ymax - tr$ymin) * H
  # pl$raster may be row-flipped; annotation_raster draws its row 1 at ymax.
  flipped <- !identical(pl$raster, sc$image)
  if (flipped) v <- H - v
  px <- pl$xmin + u * (pl$xmax - pl$xmin) / W
  py <- pl$ymax - v * (pl$ymax - pl$ymin) / H
  sqrt((px - p$x_um)^2 + (py - p$y_um)^2)
}
.edge_err <- function(pl, tr)
  abs(unlist(pl[c("xmin", "xmax", "ymin", "ymax")]) - unlist(tr[c("xmin", "xmax", "ymin", "ymax")]))

.hard_scene <- function(...) do.call(make_synthetic_ftir,
                                     utils::modifyList(synthetic_ftir_scenes()$hard, list(...)))

test_that("blob detector finds every rendered particle, colour bar excluded", {
  sc <- .hard_scene()
  b <- detect_image_blobs(sc$image)
  expect_equal(nrow(b), nrow(sc$particles))
  expect_identical(attr(b, "polarity"), "bright")
})

test_that("non-square scan with margins + colour bar: positions recovered within half a pixel", {
  sc <- .hard_scene()
  p  <- sc$particles
  pl <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  expect_identical(pl$method, "registration")
  half_px <- 0.5 * min(sc$truth$px_x, sc$truth$px_y)
  expect_lt(max(.placement_err(pl, sc)), half_px)
  expect_true(all(.edge_err(pl, sc$truth) < half_px))
  # The image is genuinely non-square and offset from the origin.
  expect_gt((pl$xmax - pl$xmin) / (pl$ymax - pl$ymin), 1.4)
  expect_lt(pl$xmin, 0)
})

test_that("the legacy particle-bounding-box placement fails the same check", {
  # Guards the test itself: stretching the image onto min/max of the
  # particles (the old FTIR tab) misplaces this scene by many pixels.
  sc <- make_synthetic_ftir()
  p  <- sc$particles
  legacy <- list(raster = sc$image, xmin = min(p$x_um), xmax = max(p$x_um),
                 ymin = min(p$y_um), ymax = max(p$y_um))
  err <- .placement_err(legacy, sc)
  expect_gt(max(err), 50 * sc$truth$px_x)
  fixed <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  expect_lt(max(.placement_err(fixed, sc)), 0.5 * sc$truth$px_x)
})

test_that("image extent does not change when the particle filter changes", {
  sc <- .hard_scene()
  p  <- sc$particles
  full <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  half_px <- 0.5 * min(sc$truth$px_x, sc$truth$px_y)
  subsets <- list(top20 = order(-p$feret_max_um)[1:20],
                  big   = which(p$feret_max_um >= 100),
                  left  = which(p$x_um < 2500))
  for (nm in names(subsets)) {
    i  <- subsets[[nm]]
    pl <- place_ftir_image(sc$image, p$x_um[i], p$y_um[i], p$feret_max_um[i])
    expect_identical(pl$method, "registration", label = nm)
    expect_true(all(abs(unlist(pl[c("xmin", "xmax", "ymin", "ymax")]) -
                        unlist(full[c("xmin", "xmax", "ymin", "ymax")])) < half_px),
                label = nm)
  }
  # A stored registration is reused verbatim, whatever particles are passed.
  reused <- place_ftir_image(sc$image, p$x_um[1:3], p$y_um[1:3],
                             registration = full$registration)
  expect_equal(unlist(reused[c("xmin", "xmax", "ymin", "ymax")]),
               unlist(full[c("xmin", "xmax", "ymin", "ymax")]))
})

test_that("stored placement round-trips through JSON", {
  sc <- make_synthetic_ftir(n = 40, seed = 11)
  p  <- sc$particles
  pl <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  f  <- tempfile(fileext = ".json")
  write_ftir_image_placement(pl, f)
  reg <- read_ftir_image_placement(f)
  again <- place_ftir_image(sc$image, numeric(0), numeric(0), registration = reg)
  expect_equal(unlist(again[c("xmin", "xmax", "ymin", "ymax")]),
               unlist(pl[c("xmin", "xmax", "ymin", "ymax")]), tolerance = 1e-9)
})

test_that("Y orientation: an asymmetric pattern lands the right way up", {
  # One big particle near the top-left, small ones elsewhere: upside down or
  # mirrored, the big blob would land on a small particle.
  pts <- rbind(c(600, 4300), c(4200, 600), c(4300, 4200), c(2500, 1500),
               c(1200, 900), c(3500, 3000), c(800, 2600), c(2900, 4400))
  sc <- make_synthetic_ftir(n = nrow(pts), pins = pts, seed = 2)
  # Pinned particles render at Feret 60; paint the first one at 400 um.
  tr <- sc$truth
  c0 <- (600 - tr$xmin) / tr$px_x; r0 <- (tr$ymax - 4300) / tr$px_y; rr <- 200 / tr$px_x
  for (r in seq_len(nrow(sc$image))) {
    cs <- which(((seq_len(ncol(sc$image)) - 0.5 - c0)^2 + (r - 0.5 - r0)^2) <= rr^2)
    if (length(cs)) { sc$image[r, cs, 1] <- 1; sc$image[r, cs, 2] <- 0.85; sc$image[r, cs, 3] <- 0.1 }
  }
  sc$particles$feret_max_um[1] <- 400
  p <- sc$particles

  pl <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  expect_identical(pl$method, "registration")
  expect_lt(max(.placement_err(pl, sc)), 0.5 * tr$px_x)
  # Raster row 1 is drawn at ymax: the big blob is in the top rows and the
  # big particle has the largest y.
  g <- (pl$raster[, , 1] + pl$raster[, , 2] + pl$raster[, , 3]) / 3
  top_rows <- seq_len(round(0.3 * nrow(g)))
  big_r <- which.max(rowSums(g > 0.6 & col(g) < 0.3 * ncol(g)))
  expect_true(big_r %in% top_rows)

  # The same image exported upside down (row 1 = smallest Y) is detected and
  # flipped back, rather than drawn mirrored.
  up <- sc; up$image <- sc$image[nrow(sc$image):1, , ]
  pl_up <- place_ftir_image(up$image, p$x_um, p$y_um, p$feret_max_um)
  expect_identical(pl_up$method, "registration")
  expect_gt(pl_up$registration$by, 0)
  expect_equal(pl_up$raster, sc$image)
  expect_lt(max(abs(unlist(pl_up[c("xmin", "xmax", "ymin", "ymax")]) -
                    unlist(tr[c("xmin", "xmax", "ymin", "ymax")]))), 0.5 * tr$px_x)
})

test_that("config extent (P1) wins; unregistrable image falls back, flagged", {
  sc <- make_synthetic_ftir(n = 30, seed = 3)
  p  <- sc$particles
  cfg <- list(ftir_image_width_um = 6000, ftir_image_height_um = 5500)
  pl <- place_ftir_image(sc$image, p$x_um, p$y_um, cfg = cfg)
  expect_identical(pl$method, "config_extent")
  expect_equal(unlist(pl[c("xmin", "xmax", "ymin", "ymax")]),
               c(xmin = 0, xmax = 6000, ymin = 0, ymax = 5500))
  # Bruker fields are read with their own prefix only.
  expect_identical(place_ftir_image(sc$image, p$x_um, p$y_um, cfg = cfg,
                                    prefix = "ftir_bruker_image")$method, "registration")

  blank <- array(0.2, dim = c(300, 300, 3))
  fb <- place_ftir_image(blank, p$x_um, p$y_um)
  expect_identical(fb$method, "particle_extent_approx")
  # Too few particles to register: also the flagged fallback.
  few <- make_synthetic_ftir(n = 4, pins = NULL, seed = 4)
  expect_identical(place_ftir_image(few$image, few$particles$x_um,
                                    few$particles$y_um)$method,
                   "particle_extent_approx")
})

test_that("a downsampled raster (viewer upload path) keeps the same placement", {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath))
  skip_if_not_installed("ggplot2")
  suppressPackageStartupMessages(library(ggplot2))
  env <- new.env(parent = globalenv())
  txt <- paste(readLines(gpath, warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)

  sc <- .hard_scene()
  p  <- sc$particles
  full <- place_ftir_image(sc$image, p$x_um, p$y_um, p$feret_max_um)
  small <- env$downsample_raster(sc$image, max_dim = 700)   # k = 3, drops remainder
  pl <- place_ftir_image(small, p$x_um, p$y_um, p$feret_max_um)
  k <- ceiling(ncol(sc$image) / 700)
  # Accurate to a fraction of one DOWNSAMPLED pixel, and the extent matches
  # the original image minus the dropped remainder columns/rows.
  expect_lt(max(abs(c(pl$xmin - full$xmin, pl$ymax - full$ymax))), 0.5 * k * sc$truth$px_x)
  expect_equal(pl$xmax - pl$xmin, ncol(small) * k * sc$truth$px_x, tolerance = 0.5 * k * sc$truth$px_x)
  # The stored full-resolution registration also applies to the small raster.
  re <- place_ftir_image(small, numeric(0), numeric(0), registration = full$registration)
  expect_lt(abs(re$xmin - full$xmin), 1e-6)
  expect_equal(re$xmax - re$xmin, ncol(small) * k * full$registration$bx, tolerance = 1e-6)

  # Report path: a report-filtered df with the full place_df gives the
  # viewer's placement, and true-size circles build without error.
  png_path <- tempfile(fileext = ".png"); png::writePNG(sc$image, png_path)
  d_all <- data.frame(x_orig = p$x_um, y_orig = p$y_um, x = p$x_um, y = p$y_um,
                      feret_max = p$feret_max_um, match_status = "unmatched",
                      particle_id = p$particle_id)
  d_flt <- d_all[d_all$feret_max >= 150, ]
  img <- env$report_instrument_image("ftir_perkin", d_flt, list(), png_path,
                                     place_df = d_all)
  expect_equal(unlist(img[c("xmin", "xmax", "ymin", "ymax")]),
               unlist(full[c("xmin", "xmax", "ymin", "ymax")]), tolerance = 1e-6)
  g <- env$make_scatter(d_flt, img, list(x = c(img$xmin, img$xmax), y = c(img$ymin, img$ymax)),
                        "t", true_size = TRUE)
  expect_s3_class(g, "ggplot")
  expect_silent(ggplot2::ggplot_build(g))
})
