# Transform math: the coordinate primitives every alignment path relies on.

test_that("similarity transform round-trips exactly", {
  src <- make_cloud(12, seed = 42)
  M   <- make_similarity(s = 1.4, deg = 25, tx = 300, ty = -150)
  dst <- apply_transform_points(src$x, src$y, M)

  est <- estimate_similarity_transform(src$x, src$y,
                                       dst$x_transformed, dst$y_transformed)
  expect_equal(est$scale,        1.4,  tolerance = 1e-8)
  expect_equal(est$rotation_deg, 25,   tolerance = 1e-8)
  expect_equal(est$tx,           300,  tolerance = 1e-6)
  expect_equal(est$ty,          -150,  tolerance = 1e-6)
  expect_lt(est$residual_rms, 1e-6)
  expect_false(est$reflected)
})

test_that("identity transform leaves points unchanged", {
  src <- make_cloud(5, seed = 7)
  out <- apply_transform_points(src$x, src$y, diag(3))
  expect_equal(out$x_transformed, src$x)
  expect_equal(out$y_transformed, src$y)
})

test_that("estimate_similarity_transform detects a reflection", {
  src   <- make_cloud(10, seed = 3)
  M_ref <- make_similarity(s = 1, deg = 40, tx = 10, ty = 20, reflect = TRUE)
  dst   <- apply_transform_points(src$x, src$y, M_ref)

  est <- estimate_similarity_transform(src$x, src$y,
                                       dst$x_transformed, dst$y_transformed,
                                       allow_reflection = TRUE)
  expect_true(est$reflected)
  expect_lt(est$residual_rms, 1e-6)
})
