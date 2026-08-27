# Task 3: shared alignment primitives (R/align_helpers.R), extracted from
# global_register_align(). The characterization test pins the exact pre-refactor
# output so the extraction is provably behaviour-preserving.

test_that("align_pose_xy applies rotation, scale and mirror", {
  p <- align_pose_xy(90, 1, FALSE, 1, 0)          # 90 deg: (1,0) -> (0,1)
  expect_equal(p$x, 0, tolerance = 1e-12)
  expect_equal(p$y, 1, tolerance = 1e-12)

  p2 <- align_pose_xy(0, 2, FALSE, 3, -4)          # scale 2, no rotation
  expect_equal(p2$x, 6); expect_equal(p2$y, -8)

  pm <- align_pose_xy(0, 1, TRUE, 3, 4)            # mirror flips the y cross-term
  expect_equal(pm$x, 3); expect_equal(pm$y, -4)
})

test_that("align_one_to_one counts mutual matches, each point used once", {
  px <- c(0, 100, 9999); py <- c(0, 0, 0)
  qx <- c(0, 100);       qy <- c(0, 0)
  expect_equal(align_one_to_one(px, py, qx, qy, tol = 10), 2L)
  expect_equal(align_one_to_one(px, py, qx + 500, qy, tol = 10), 0L)
})

test_that("align_rspan is the robust 5-95 percentile span, floored above zero", {
  expect_equal(align_rspan(c(0, 100)),
               unname(diff(stats::quantile(c(0, 100), c(0.05, 0.95)))))
  expect_gt(align_rspan(rep(5, 10)), 0)
})

test_that("global_register_align reproduces the pre-refactor baseline", {
  set.seed(5)
  src <- data.frame(x_norm = runif(20, 0, 1000), y_norm = runif(20, 0, 1000))
  M   <- make_similarity(s = 1, deg = 30, tx = 150, ty = 50)
  tr  <- apply_transform_points(src$x_norm, src$y_norm, M)
  ref <- data.frame(x_norm = tr$x_transformed, y_norm = tr$y_transformed)
  cfg <- make_config(); cfg$ransac_coarse_step_deg <- 10L

  g <- global_register_align(src, ref, cfg)
  expect_equal(g$n_inliers, 20)

  # Re-pinned when the scale sweep was narrowed to the plausible band (the old
  # seq(0.2, 1.3, 0.05) offered the search collapsed poses, which is how a real
  # run ended up at scale 0.2336). The grid points moved slightly, so this
  # baseline moved with them.
  #
  # The new values are NEARER ground truth on every parameter -- the synthetic
  # case is scale 1.0, 30 deg, tx 150, ty 50:
  #   scale 0.9908 -> 0.9930 ; tx 152.05 -> 151.57 ; ty 55.66 -> 54.32
  # so this is an accuracy improvement, not a drift to be tolerated.
  baseline <- matrix(c(0.85997694,   0.49650792, 0,
                      -0.49650792,   0.85997694, 0,
                       151.56799790, 54.31971715, 1), nrow = 3)
  expect_equal(unname(g$transform), baseline, tolerance = 1e-6)
})
