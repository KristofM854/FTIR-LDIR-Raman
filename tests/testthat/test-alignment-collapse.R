# =============================================================================
# test-alignment-collapse.R -- regression for the LDIR<->Raman scale collapse
# =============================================================================
# Real geometry from run 2026-08-27_11, which produced a transform of
# scale 0.2336 / rotation -18.03 deg where the truth is scale ~1.05 /
# rotation ~+90 deg. The LDIR cloud was shrunk 4.3x into a blob and the
# overlay was unusable.
#
# The failure is not in the pose search -- that finds the correct answer. It is
# that nearest-neighbour RMS REWARDS collapse: packing the source into a dense
# region of the target lowers RMS as the fit gets more wrong. Measured here:
#   correct pose   scale 1.06  RMS 301 um  31 one-to-one inliers
#   collapsed pose scale 0.23  RMS  67 um  18 one-to-one inliers
# So RMS-based selection picks the wrong transform and inlier count does not.

ldir_x <- c(166.2, -307.0, 818.8, -913.1, -762.0, 1649.8, 27.1, -1881.1, -1564.1, -1500.2, -415.2, 527.4, -858.9, -1016.3, 773.9, -1687.8, 93.7, -1404.6, 1228.1, 266.7, 1852.2, -1835.6, 1664.3, 1936.5, 738.8, 2095.8, 1354.1, 1348.2, 1004.5, 1954.2, 1135.1, 1107.3, -2171.4, -2091.4, 1344.1, 1057.0, 2095.4, 1880.3, 1232.9, 837.9, 814.8, 953.4, 894.4, 1274.0, 865.3, 1301.0, 31.5, 2181.1, 1593.7, -1620.9, 1483.9, -578.2, 850.3, -1923.6, 2213.0, 904.8, 1104.4, 748.2, -1566.1, 2306.0, 252.3, -114.9, -13.6)
ldir_y <- c(1460.7, 1627.7, -2002.3, -404.6, 2222.0, 1963.9, 735.8, 1226.8, -804.7, -757.0, -1376.5, -872.9, -1587.8, 1484.1, -1193.5, -1214.9, -1186.1, 128.6, -757.0, 1884.4, -2062.2, 1805.8, -1250.0, -2269.1, 2222.9, 1644.3, -1600.2, 1571.4, -51.1, -1915.0, 2040.9, -504.9, -2307.8, -1658.2, -425.3, -485.3, 1371.4, 1339.0, -1012.6, -318.9, -455.1, -392.0, -1010.7, -1421.6, -397.6, 1616.8, 197.9, -83.6, -112.5, 242.1, -1106.0, -1880.8, -1259.8, -1359.0, -315.0, -490.1, -1836.0, -2137.1, -442.5, 214.0, -1412.8, 1630.0, -2021.3)
raman_x <- c(-2547.1, -1020.3, -3980.9, -1540.6, -484.1, -3067.1, -1763.5, -474.3, -1491.8, -4676.2, -4389.5, -3854.8, -105.8, 113.0, -595.5, -1936.7, -1219.2, -1847.8, -148.0, -261.2, -654.7, -4643.9, -1125.4, -4620.8, -1551.2, -1955.2, -1216.6, -2008.0, -1825.2, -4620.8, -4629.7, -782.1, -7.4, -4659.7, -1874.3, -1809.9, -2291.8, -1872.9, -3990.8, -969.4, -1803.9, -1772.7, -1779.3, 208.6, -3702.3)
raman_y <- c(4739.2, 3115.0, 3689.6, 3323.5, 65.1, 2296.2, 3463.5, 4501.8, 3598.4, 1384.9, 4021.1, 2445.7, 4289.9, 4382.3, 3744.2, 3212.9, 3607.3, 3728.6, 3180.0, 4394.4, 3627.5, 1395.0, 3879.0, 1430.0, 3181.7, 3176.2, 3249.3, 3163.2, 3319.5, 1435.8, 1441.2, 3659.3, 3104.1, 1420.7, 3215.8, 3415.9, 3286.8, 3302.7, 1907.3, 4075.4, 3148.5, 3246.1, 3393.4, -4.7, 2199.5)

.centre <- function(v) v - mean(v)
lx <- .centre(ldir_x); ly <- .centre(ldir_y)
rx <- .centre(raman_x); ry <- .centre(raman_y)
TOL <- 200

# Score a pose the way the pipeline does: align_score_pose() recovers the
# translation by VOTING over pairwise offsets, then counts one-to-one inliers.
# This matters -- with a naive centroid-matched translation the collapse can
# out-score the truth, because centroid matching is a poor offset for the
# correct pose. The criterion is only collapse-proof when each candidate is
# given its own best translation, which is exactly what the sweep does.
.score <- function(deg, s, mir = FALSE)
  align_score_pose(deg, s, mir, lx, ly, rx, ry, TOL)$n

test_that("one-to-one inlier count prefers the correct pose over the collapse", {
  n_correct   <- .score(90, 1.05)
  n_collapsed <- .score(-18, 0.234)
  expect_gt(n_correct, n_collapsed)
})

test_that("RMS does NOT separate them -- it prefers the collapse", {
  # Pinned so nobody reintroduces an RMS-based acceptance gate.
  rms_at <- function(deg, s) {
    r <- align_score_pose(deg, s, FALSE, lx, ly, rx, ry, TOL)
    p <- align_pose_xy(deg, s, FALSE, lx, ly)
    px <- p$x + r$tx; py <- p$y + r$ty
    d <- sqrt(outer(rx, px, "-")^2 + outer(ry, py, "-")^2)
    sqrt(mean(apply(d, 2, min)^2))
  }
  expect_lt(rms_at(-18, 0.234), rms_at(90, 1.05))
})

test_that("collapsing cannot beat the true pose once each gets its best offset", {
  n_base <- .score(90, 1.05)
  for (s in c(0.5, 0.25, 0.1, 0.05))
    expect_lte(.score(90, s), n_base)
})

test_that("global_register_align recovers the correct pose, not the collapse", {
  cfg <- list(ransac_inlier_dist_um = TOL, ransac_coarse_step_deg = 10)
  src <- data.frame(x_norm = lx, y_norm = ly)
  ref <- data.frame(x_norm = rx, y_norm = ry)
  res <- global_register_align(src, ref, cfg, allow_mirror = FALSE)
  expect_gt(res$params$scale, 0.75)
  expect_lt(res$params$scale, 1.4)
  expect_lt(abs(abs(res$params$rotation_deg) - 90), 15)
})

test_that("icp_refine clamps a collapsing scale back into the plausible band", {
  # Feed ICP a deliberately wrong starting pose; unconstrained it walks the
  # scale down (measured: 1.05 -> 0.58 from deg=0). The guard must hold it.
  cfg <- list(icp_max_iterations = 30, icp_convergence_thresh = 0.01,
              icp_max_pair_dist_um = 500, ransac_allow_mirror = FALSE,
              icp_trim_pct = 0.2)
  M0 <- matrix(c(1.05, 0, 0, 0, 1.05, 0, 0, 0, 1), nrow = 3)   # deg 0, s 1.05
  out <- icp_refine(data.frame(x_norm = lx, y_norm = ly),
                    data.frame(x_norm = rx, y_norm = ry), M0, cfg)
  expect_gte(out$params$scale, 0.8 - 1e-9)
  expect_lte(out$params$scale, 1.25 + 1e-9)
})
