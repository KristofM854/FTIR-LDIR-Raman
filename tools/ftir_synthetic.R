# =============================================================================
# ftir_synthetic.R -- synthetic FTIR export (particle CSV + absorption PNG)
# =============================================================================
# Builds a scene whose geometry is known exactly, so placement can be checked
# against ground truth rather than against itself. Used by
# tools/diag_ftir_overlay.R and tests/testthat/test-ftir-image-placement.R.
#
# Geometry: the data region of the PNG covers the scan rectangle `scan`
# (um) with `data_px` pixels; `margins` (px) of white border are added around
# it, plus an optional colour bar in the right margin. Particle coordinates
# are written as the CSV's "Coord. [um]" column ("[x;y]"), Y pointing up.
# =============================================================================

make_synthetic_ftir <- function(n = 150,
                                scan = c(xmin = 0, xmax = 5000, ymin = 0, ymax = 5000),
                                data_px = c(w = 960, h = 960),
                                margins = c(left = 0, right = 0, top = 0, bottom = 0),
                                colorbar = FALSE,
                                particle_box = c(100, 4950, 150, 4580),
                                pins = rbind(c(4950, 150), c(4870, 4580), c(100, 2600)),
                                feret_range = c(25, 300),
                                seed = 1L) {
  set.seed(seed)
  px_x <- (scan[["xmax"]] - scan[["xmin"]]) / data_px[["w"]]
  px_y <- (scan[["ymax"]] - scan[["ymin"]]) / data_px[["h"]]
  W <- margins[["left"]] + data_px[["w"]] + margins[["right"]]
  H <- margins[["top"]]  + data_px[["h"]] + margins[["bottom"]]

  # Particles: pinned extremes + random, with a minimum separation so blobs
  # do not merge. Feret log-uniform over feret_range.
  feret <- exp(runif(n, log(feret_range[1]), log(feret_range[2])))
  pts <- matrix(NA_real_, n, 2)
  npin <- if (is.null(pins)) 0L else nrow(pins)
  if (npin > 0) { pts[seq_len(npin), ] <- pins; feret[seq_len(npin)] <- 60 }
  i <- npin
  tries <- 0L
  while (i < n && tries < 200000L) {
    tries <- tries + 1L
    p <- c(runif(1, particle_box[1], particle_box[2]),
           runif(1, particle_box[3], particle_box[4]))
    k <- i + 1L
    if (i > 0) {
      d <- sqrt((pts[seq_len(i), 1] - p[1])^2 + (pts[seq_len(i), 2] - p[2])^2)
      if (any(d < (feret[seq_len(i)] + feret[k]) / 2 + 3 * max(px_x, px_y))) next
    }
    pts[k, ] <- p; i <- k
  }
  pts <- pts[seq_len(i), , drop = FALSE]; feret <- feret[seq_len(i)]

  # Pixel (edge) coordinates of each particle centre in the full PNG.
  col_c <- margins[["left"]] + (pts[, 1] - scan[["xmin"]]) / px_x
  row_c <- margins[["top"]]  + (scan[["ymax"]] - pts[, 2]) / px_y

  # Render: white margins, dark-blue map background with mild noise, particles
  # as bright discs (yellow/red, like a false-colour absorbance map).
  img <- array(1, dim = c(H, W, 3))
  dr <- margins[["top"]]  + seq_len(data_px[["h"]])
  dc <- margins[["left"]] + seq_len(data_px[["w"]])
  bg <- c(0.10, 0.12, 0.45)
  for (ch in 1:3)
    img[dr, dc, ch] <- pmin(1, pmax(0, bg[ch] + rnorm(length(dr) * length(dc), 0, 0.01)))
  fg <- c(1.00, 0.85, 0.10)
  for (j in seq_len(nrow(pts))) {
    rx <- feret[j] / 2 / px_x; ry <- feret[j] / 2 / px_y
    cs <- max(1, floor(col_c[j] - rx)):min(W, ceiling(col_c[j] + rx))
    rs <- max(1, floor(row_c[j] - ry)):min(H, ceiling(row_c[j] + ry))
    inside <- outer(((rs - 0.5) - row_c[j]) / ry, ((cs - 0.5) - col_c[j]) / rx,
                    function(a, b) a^2 + b^2 <= 1)
    for (ch in 1:3) {
      m <- img[rs, cs, ch]; m[inside] <- fg[ch]; img[rs, cs, ch] <- m
    }
  }
  if (isTRUE(colorbar) && margins[["right"]] >= 30) {
    cb_c <- margins[["left"]] + data_px[["w"]] + 10 + seq_len(15)
    cb_r <- margins[["top"]] + seq_len(data_px[["h"]])
    ramp <- seq(1, 0, length.out = length(cb_r))
    img[cb_r, cb_c, 1] <- ramp; img[cb_r, cb_c, 2] <- 0.5 * ramp
    img[cb_r, cb_c, 3] <- 1 - ramp
    img[cb_r, cb_c[c(1, 15)], ] <- 0      # frame lines
  }

  # Whole-PNG extent in um (edge convention).
  truth <- list(
    xmin = scan[["xmin"]] - margins[["left"]] * px_x,
    xmax = scan[["xmax"]] + margins[["right"]] * px_x,
    ymin = scan[["ymin"]] - margins[["bottom"]] * px_y,
    ymax = scan[["ymax"]] + margins[["top"]] * px_y,
    px_x = px_x, px_y = px_y, W = W, H = H)

  particles <- data.frame(
    particle_id = sprintf("P%03d", seq_len(nrow(pts))),
    x_um = pts[, 1], y_um = pts[, 2], feret_max_um = feret,
    stringsAsFactors = FALSE)
  list(image = img, particles = particles, truth = truth)
}

#' Write the scene in the instrument's export format: a comma-separated CSV
#' with "Coord. [um]" = "[x;y]" (what ingest_ftir() parses) and a PNG.
write_synthetic_ftir <- function(scene, dir, stem = "synthetic") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- scene$particles
  csv <- data.frame(
    "Identifier"      = p$particle_id,
    "Group"           = "PE",
    "Coord. [um]"     = sprintf("[%.1f;%.1f]", p$x_um, p$y_um),
    "Area on map [um2]" = pi * (p$feret_max_um / 2)^2,
    "Major dim [um]"  = p$feret_max_um,
    "Minor dim [um]"  = p$feret_max_um,
    "Feret min [um]"  = p$feret_max_um,
    "Max AAU score"   = 0.9,
    check.names = FALSE)
  csv_path <- file.path(dir, paste0(stem, ".csv"))
  png_path <- file.path(dir, paste0("Average Abs.( ", stem, " ).png"))
  utils::write.csv(csv, csv_path, row.names = FALSE)
  png::writePNG(scene$image, png_path)
  list(csv = csv_path, png = png_path)
}

#' The two scenes the diagnostic and tests use.
#' "reported": square 5000 um scan from (0,0), no margins, particles spanning
#'             X 100-4950 / Y 150-4580 -- the geometry read off the screenshot.
#' "hard":     non-square 6000 x 4000 um scan offset from the origin, white
#'             margins on every side and a colour bar.
synthetic_ftir_scenes <- function() {
  list(
    reported = list(),
    hard = list(n = 120,
                scan = c(xmin = -500, xmax = 5500, ymin = 250, ymax = 4250),
                data_px = c(w = 1500, h = 1000),
                margins = c(left = 60, right = 140, top = 40, bottom = 70),
                colorbar = TRUE,
                particle_box = c(-200, 5100, 500, 4000),
                pins = NULL, seed = 7L))
}
