# Anchor selection and landmark selection must never hand the aligner an empty
# or degenerate particle set.
#
# Background: a Polycarbonate-dominated sample produced zero Raman alignment
# anchors (align_raman_materials is PET/PP). nrow == 0 made the centroid NaN,
# every x_norm/y_norm became NaN, and the failure only surfaced ~40 lines later
# as "NA/NaN/Inf in foreign function call (arg 1)" from RANN::nn2 inside the
# RANSAC coarse search. These tests pin the three guards that stop that.

local({
  rdir <- file.path(REPO_ROOT, "R")
  for (m in c("03_normalize.R", "03b_landmark_align.R")) {
    sys.source(file.path(rdir, m), envir = globalenv())
  }
})

.cfg <- make_config()

.particles <- function(size, material = "Polycarbonate", n = length(size)) {
  data.frame(
    particle_id  = paste0("P_", seq_len(n)),
    x_um         = seq_len(n) * 100,
    y_um         = seq_len(n) * 50,
    feret_max_um = size,
    major_um     = size,
    minor_um     = size / 2,
    material     = rep_len(material, n),
    stringsAsFactors = FALSE
  )
}

# ---- select_material_anchors ------------------------------------------------

test_that("the material filter is honoured when it leaves enough anchors", {
  df  <- .particles(rep(50, 10), material = c(rep("PET", 6), rep("Polycarbonate", 4)))
  out <- select_material_anchors(df, c("PET", "Polypro"), min_count = 4, "FTIR")
  expect_equal(nrow(out), 6L)
  expect_true(all(out$material == "PET"))
})

test_that("a filter that would empty the anchor set is dropped, not applied", {
  # The PC sample: nothing matches PET/PP at all.
  df  <- .particles(rep(50, 9), material = "Polycarbonate")
  out <- select_material_anchors(df, c("Polyethylene terephtalate", "Polypropylene"),
                                 min_count = 4, "Raman")
  expect_equal(nrow(out), 9L)
  expect_false(any(is.na(out$x_um)))
})

test_that("a filter leaving fewer than min_count anchors is dropped", {
  df  <- .particles(rep(50, 10), material = c(rep("PET", 3), rep("Polycarbonate", 7)))
  out <- select_material_anchors(df, "PET", min_count = 4, "FTIR")
  expect_equal(nrow(out), 10L)   # 3 < 4, so keep everything
})

test_that("no configured materials means no filtering", {
  df <- .particles(rep(50, 5))
  expect_equal(nrow(select_material_anchors(df, NULL, 4, "FTIR")), 5L)
  expect_equal(nrow(select_material_anchors(df, character(0), 4, "FTIR")), 5L)
})

test_that("NA materials never slip into the anchor set", {
  df <- .particles(rep(50, 8), material = c(rep("PET", 5), rep(NA_character_, 3)))
  out <- select_material_anchors(df, "PET", min_count = 4, "FTIR")
  expect_equal(nrow(out), 5L)
  expect_false(any(is.na(out$material)))
})

# ---- normalize_coordinates guards -------------------------------------------

test_that("an empty cloud is rejected with a clear message, not a NaN centroid", {
  ok <- .particles(rep(50, 5))
  expect_error(normalize_coordinates(ok, ok[0, ], normalize_scale = FALSE),
               "centroid source is empty")
  expect_error(normalize_coordinates(ok[0, ], ok, normalize_scale = FALSE),
               "centroid source is empty")
})

test_that("an all-NA coordinate column is rejected", {
  ok <- .particles(rep(50, 5))
  bad <- ok
  bad$x_um <- NA_real_
  expect_error(normalize_coordinates(ok, bad, normalize_scale = FALSE),
               "entirely NA")
})

test_that("a healthy pair still normalizes to zero-mean coordinates", {
  a <- .particles(rep(50, 5))
  b <- .particles(rep(50, 7))
  res <- normalize_coordinates(a, b, normalize_scale = FALSE)
  expect_equal(mean(res$ftir$x_norm), 0)
  expect_equal(mean(res$raman$y_norm), 0)
  expect_false(any(is.na(res$ftir$x_norm)))
})

# ---- apply_normalization ----------------------------------------------------

test_that("apply_normalization reproduces normalize_coordinates on the same set", {
  a   <- .particles(rep(50, 6))
  res <- normalize_coordinates(a, a, normalize_scale = FALSE)
  out <- apply_normalization(a, res$ftir_centroid, res$ftir_scale)
  expect_equal(out$x_norm, res$ftir$x_norm)
  expect_equal(out$y_norm, res$ftir$y_norm)
})

test_that("apply_normalization puts a subset in the SAME frame as the parent", {
  # This is the property that makes a full-cloud centroid safe to use for the
  # material anchor subset: the subset's coordinates must not shift.
  a      <- .particles(rep(50, 10))
  res    <- normalize_coordinates(a, a, normalize_scale = FALSE)
  subset <- a[c(2, 5, 9), ]
  out    <- apply_normalization(subset, res$ftir_centroid, res$ftir_scale)
  expect_equal(out$x_norm, res$ftir$x_norm[c(2, 5, 9)])
})

test_that("apply_normalization survives a degenerate scale", {
  a <- .particles(rep(50, 4))
  expect_false(any(is.na(apply_normalization(a, c(0, 0), 0)$x_norm)))
  expect_false(any(is.na(apply_normalization(a, c(0, 0), NA)$x_norm)))
})

# ---- Tier 2 aligner interchangeability --------------------------------------

test_that("global_register_align is a drop-in for ransac_align in Tier 2", {
  # Tier 2 picks whichever aligner pairs more particles, so both must return the
  # same shape: main.R reads $transform, $params$scale/$rotation_deg/$reflected
  # and $n_inliers off whichever one wins.
  cfg <- make_config()
  cfg$ransac_coarse_step_deg <- 10L
  cfg$ransac_n_iterations    <- 100L

  src <- make_cloud(18, seed = 3)
  names(src) <- c("x_norm", "y_norm")
  src$feret_max_um <- seq(40, 200, length.out = 18)
  M  <- make_similarity(s = 1, deg = 18, tx = 120, ty = -75)
  tr <- apply_transform_points(src$x_norm, src$y_norm, M)
  ref <- data.frame(x_norm = tr$x_transformed, y_norm = tr$y_transformed,
                    feret_max_um = src$feret_max_um)

  r <- ransac_align(src, ref, cfg)
  g <- global_register_align(src, ref, cfg, allow_mirror = cfg$ransac_allow_mirror)

  for (res in list(r, g)) {
    expect_true(is.matrix(res$transform))
    expect_equal(dim(res$transform), c(3L, 3L))
    expect_false(any(is.na(res$transform)))
    expect_true(is.numeric(res$params$scale))
    expect_true(is.numeric(res$params$rotation_deg))
    expect_true(is.logical(res$params$reflected))
    expect_true(is.numeric(res$n_inliers))
  }
})

# ---- end-to-end: the original crash must not be reachable -------------------

test_that("a PC-only sample aligns instead of crashing in RANN::nn2", {
  # Reproduces the reported failure end to end: FTIR carries the anchor material
  # (PET), Raman is pure Polycarbonate, so the PET/PP filter empties the Raman
  # anchor set. Previously: NaN centroid -> NaN x_norm -> RANN::nn2 rejects the
  # reference matrix with "NA/NaN/Inf in foreign function call (arg 1)".
  cfg <- make_config()
  cfg$ransac_coarse_step_deg <- 10L
  cfg$ransac_n_iterations    <- 100L

  ftir <- .particles(seq(260, 120, length.out = 12),
                     material = c(rep("PET", 6), rep("Polycarbonate", 6)))
  ftir$x_um <- c(120, 800, 1500, 2300, 3100, 3600, 4200, 700, 1900, 2800, 3900, 4600)
  ftir$y_um <- c(400, 1200, 900, 2100, 1600, 3000, 2400, 3300, 3800, 800, 4100, 2900)

  # Raman sees the same particles, shifted and rotated, all identified as PC
  M  <- make_similarity(s = 1, deg = 12, tx = 250, ty = -180)
  tr <- apply_transform_points(ftir$x_um, ftir$y_um, M)
  raman <- ftir
  raman$material <- "Polycarbonate"
  raman$x_um <- tr$x_transformed
  raman$y_um <- tr$y_transformed

  ftir_anchors  <- select_material_anchors(ftir,  cfg$align_ftir_materials,  4, "FTIR")
  raman_anchors <- select_material_anchors(raman, cfg$align_raman_materials, 4, "Raman")
  expect_equal(nrow(raman_anchors), nrow(raman))  # filter dropped, not applied

  norm <- normalize_coordinates(ftir, raman, normalize_scale = FALSE)
  expect_false(any(is.na(norm$raman_centroid)))

  ftir_n  <- apply_normalization(ftir_anchors,  norm$ftir_centroid,  norm$ftir_scale)
  raman_n <- apply_normalization(raman_anchors, norm$raman_centroid, norm$raman_scale)
  expect_false(any(is.na(cbind(raman_n$x_norm, raman_n$y_norm))))

  res <- ransac_align(ftir_n, raman_n, cfg)
  expect_false(any(is.na(res$transform)))
  expect_gt(res$n_inliers, 0)
})
