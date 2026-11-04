# A short FTIR-type particle list (the Lumos lists only library-matched
# particles) cannot fix its pose against a dense Raman field; the instrument's
# own image shows every particle. align_and_match(image_cloud) fits the pose on
# the image's spots when the particles do not determine it. Also: the Bruker
# image extent from the OPUS values (opus_ftir_bruker_image) and the pose
# limits (no mirror, scale near 1).

local({
  for (m in c("02_prefilter.R", "03_normalize.R", "03b_landmark_align.R",
              "05_transform.R", "08b_material_map.R", "08_agreement.R",
              "align_and_match.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
})

.mk <- function(id, x, y, size, mat, q) data.frame(
  particle_id = id, x_um = x, y_um = y, feret_max_um = size, feret_min_um = size * 0.6,
  major_um = size, minor_um = size * 0.6, area_um2 = size^2 * 0.5, material = mat,
  quality = q, stringsAsFactors = FALSE)

# Dense Raman field; a 12-particle Bruker list rotated 177 deg and shifted; a
# striped photo in the Bruker frame (7.9 um/px) showing 80% of all particles
.spot_scene <- function() {
  set.seed(11)
  n <- 900; A <- 12000
  tx <- stats::runif(n, 0, A); ty <- stats::runif(n, 0, A)
  sz <- pmin(400, 20 + stats::rexp(n, 1 / 40))
  raman <- .mk(paste0("R", 1:n), tx, ty, sz,
               sample(c("Polyethylene", "Cellulose", "Polyamide", "PET"), n, TRUE), 90)
  th <- 177 * pi / 180
  bx <- cos(th) * tx - sin(th) * ty + 13000; by <- sin(th) * tx + cos(th) * ty + 12500
  sel <- sample(which(sz < 60), 12)
  bruker <- .mk(paste0("B", 1:12), bx[sel] + stats::rnorm(12, 0, 5),
                by[sel] + stats::rnorm(12, 0, 5), sz[sel], "Protein", 0.8)
  px <- 7.9; W <- 2048; H <- 1949; x0 <- -1500; y0 <- -1500
  img <- matrix(0.8, H, W) + 0.08 * sin(seq_len(H) / 25)
  for (i in which(stats::runif(n) < 0.8)) {
    cx <- (bx[i] - x0) / px; cy <- H - (by[i] - y0) / px; r <- max(1.5, sz[i] / 2 / px)
    rr <- max(1, floor(cy - r)):min(H, ceiling(cy + r))
    cc <- max(1, floor(cx - r)):min(W, ceiling(cx + r))
    if (!length(rr) || !length(cc) || min(rr) > H || min(cc) > W) next
    d <- outer((rr - cy)^2, (cc - cx)^2, "+") <= r^2
    img[rr, cc][d] <- 0.3
  }
  list(raman = raman, bruker = bruker, raw = array(img, c(H, W, 3)),
       pl = list(xmin = x0, xmax = x0 + W * px, ymin = y0, ymax = y0 + H * px,
                 flip_rows = FALSE, approximate = FALSE))
}

test_that("a pose the particle list cannot fix is fitted on the image's spots", {
  sc <- .spot_scene()
  cfg <- make_config(); cfg$output_dir <- withr::local_tempdir()
  ref <- suppressMessages(prepare_reference(sc$raman, cfg, "raman"))
  spots <- image_spot_cloud(sc$raw, sc$pl, sc$bruker$x_um, sc$bruker$y_um, max_spots = 600)
  expect_gt(nrow(spots), 200)
  without <- suppressMessages(align_and_match(sc$bruker, ref, cfg, "FTIR (Bruker)", "ftir_bruker"))
  expect_false(without$alignment_determined)
  with <- suppressMessages(align_and_match(sc$bruker, ref, cfg, "FTIR (Bruker)", "ftir_bruker",
                                           image_cloud = spots))
  expect_true(with$alignment_determined)
  expect_match(with$method, "image_spots")
  expect_equal(with$icp$params$rotation_deg, -177, tolerance = 0.5 / 177)
  expect_equal(with$icp$params$scale, 1, tolerance = 0.01)
  expect_equal(with$match$match_stats$n_matched, 12)
  # An approximate placement gives no spots: its frame is not the table's
  expect_null(image_spot_cloud(sc$raw, utils::modifyList(sc$pl, list(approximate = TRUE)),
                               sc$bruker$x_um, sc$bruker$y_um))
})

test_that("the Bruker pose search is limited to no mirror and a scale near 1", {
  cfg <- apply_pose_limits(make_config(), "ftir_bruker")
  expect_false(cfg$ransac_allow_mirror)
  expect_equal(c(cfg$align_scale_min, cfg$align_scale_max), c(0.97, 1.03))
  expect_equal(c(cfg$icp_scale_min, cfg$icp_scale_max), c(0.97, 1.03))
  expect_false(apply_pose_limits(make_config(), "ftir")$ransac_allow_mirror)    # PerkinElmer too
  expect_identical(apply_pose_limits(make_config(), "ldir")$ransac_allow_mirror,
                   make_config()$ransac_allow_mirror)
  # Without mirrored poses the search is half as large
  set.seed(3); rx <- stats::runif(200, 0, 1e4); ry <- stats::runif(200, 0, 1e4)
  p2 <- pose_significance(12, 30, rx, ry, 50)
  p1 <- pose_significance(12, 30, rx, ry, 50, mirror = FALSE)
  expect_lt(p2$p, 1)
  expect_equal(p1$p, p2$p / 2)
})

test_that("the OPUS values place the Bruker image, and stale values are refused", {
  cfg <- list(opus_ftir_bruker_image = list(
    pixel_size_x_um = 7.907902, pixel_size_y_um = 7.908996, pos_x_px = 6110.180843,
    pos_y_px = 604.290875, size_x_px = 2048, size_y_px = 1949,
    map_origin_x_um = 49564, map_origin_y_um = 6050.38))
  x <- c(264, 11753, 595); y <- c(3798, 7013, 73)
  e <- ftir_extent_from_opus(cfg, "ftir_bruker_image", 2048, 1949, x, y)
  expect_equal(unlist(e), c(xmin = -1245.3, xmax = 14950.1, ymin = -1271.0, ymax = 14143.6),
               tolerance = 1e-4)
  expect_message(expect_null(ftir_extent_from_opus(cfg, "ftir_bruker_image", 3000, 2000, x, y)),
                 "another scan")
  expect_message(expect_null(ftir_extent_from_opus(cfg, "ftir_bruker_image", 2048, 1949,
                                                   x + 30000, y)), "inside")
  # place_ftir_image uses it before any registration, also for a downsized raster
  raw <- array(0.5, c(974, 1024, 3)); attr(raw, "resized_from_px") <- c(2048, 1949)
  pl <- place_ftir_image(raw, x, y, cfg = cfg, prefix = "ftir_bruker_image")
  expect_identical(pl$method, "opus_metadata")
  expect_equal(pl$xmin, e$xmin)
})
