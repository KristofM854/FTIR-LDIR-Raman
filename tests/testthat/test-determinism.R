# B3: alignment must be reproducible run-to-run, and must not perturb the
# caller's global RNG (the two properties the local-seed + restore pattern buys).

# Determinism is independent of how exhaustive the search is, so shrink the grid
# and iteration count to keep these tests fast (production defaults are 1 deg /
# 2000 iterations, which is needlessly slow for a reproducibility assertion).
fast_config <- function() {
  cfg <- make_config()
  cfg$ransac_coarse_step_deg <- 10L
  cfg$ransac_n_iterations    <- 100L
  cfg
}

test_that("ransac_align is deterministic across runs", {
  ftir <- make_cloud(15, seed = 11)
  names(ftir) <- c("x_norm", "y_norm")
  M  <- make_similarity(s = 1, deg = 20, tx = 200, ty = -100)
  tr <- apply_transform_points(ftir$x_norm, ftir$y_norm, M)
  raman <- data.frame(x_norm = tr$x_transformed, y_norm = tr$y_transformed)
  cfg <- fast_config()

  r1 <- ransac_align(ftir, raman, cfg)
  r2 <- ransac_align(ftir, raman, cfg)
  expect_identical(r1$transform, r2$transform)
  expect_identical(r1$n_inliers, r2$n_inliers)
})

test_that("global_register_align is deterministic and RNG-safe", {
  src <- make_cloud(20, seed = 5)
  names(src) <- c("x_norm", "y_norm")
  M  <- make_similarity(s = 1, deg = 30, tx = 150, ty = 50)
  tr <- apply_transform_points(src$x_norm, src$y_norm, M)
  ref <- data.frame(x_norm = tr$x_transformed, y_norm = tr$y_transformed)
  cfg <- fast_config()

  g1 <- global_register_align(src, ref, cfg)
  g2 <- global_register_align(src, ref, cfg)
  expect_identical(g1$transform, g2$transform)

  # The aligner seeds a *local* RNG stream and restores .Random.seed on exit,
  # so an identical global draw must bracket a call to it unchanged.
  set.seed(123); before <- runif(3)
  set.seed(123); invisible(global_register_align(src, ref, cfg)); after <- runif(3)
  expect_identical(before, after)
})

test_that("align_seed controls the sampler (different seed is still stable)", {
  src <- make_cloud(40, seed = 9)
  names(src) <- c("x_norm", "y_norm")
  M  <- make_similarity(s = 1, deg = 22, tx = 80, ty = -60)
  tr <- apply_transform_points(src$x_norm, src$y_norm, M)
  ref <- data.frame(x_norm = tr$x_transformed, y_norm = tr$y_transformed)

  cfg_a <- fast_config(); cfg_a$align_seed <- 1L
  cfg_b <- fast_config(); cfg_b$align_seed <- 999L
  # Each seed is internally reproducible.
  expect_identical(global_register_align(src, ref, cfg_a)$transform,
                   global_register_align(src, ref, cfg_a)$transform)
  expect_identical(global_register_align(src, ref, cfg_b)$transform,
                   global_register_align(src, ref, cfg_b)$transform)
})
