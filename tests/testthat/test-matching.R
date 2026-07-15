# Particle matching: optimal 1-to-1 assignment and the dedicated LDIR gate.

test_that("Hungarian matching pairs coincident particles 1-to-1", {
  shared <- data.frame(x = c(0, 100, 200, 300), y = c(0, 50, 100, 150))
  src <- data.frame(
    particle_id  = c("S1", "S2", "S3", "S4", "S5"),
    x_aligned    = c(shared$x, 5000), y_aligned = c(shared$y, 5000),
    feret_max_um = c(10, 20, 30, 40, 50), area_um2 = c(100, 400, 900, 1600, 2500),
    major_um     = c(12, 22, 32, 42, 52), minor_um = c(8, 18, 28, 38, 48),
    stringsAsFactors = FALSE
  )
  ref <- data.frame(
    particle_id  = c("R1", "R2", "R3", "R4", "R5"),
    x_norm       = c(shared$x, -5000), y_norm = c(shared$y, -5000),
    feret_max_um = c(10, 20, 30, 40, 50), area_um2 = c(100, 400, 900, 1600, 2500),
    major_um     = c(12, 22, 32, 42, 52), minor_um = c(8, 18, 28, 38, 48),
    stringsAsFactors = FALSE
  )

  res <- match_particles(src, ref, make_config(),
                         src_label = "ftir", ref_label = "raman")

  expect_equal(res$match_stats$n_matched, 4L)
  # 1-to-1: no reference particle assigned twice.
  expect_equal(anyDuplicated(res$matched$raman_idx), 0L)
  # The four coincident particles pair up; the far-flung pair does not.
  expect_setequal(res$matched$ftir_particle_id, c("S1", "S2", "S3", "S4"))
})

test_that("LDIR pairing uses the looser LDIR distance gate", {
  # One pair separated by 200 um: outside the 100 um FTIR gate, inside the
  # 250 um LDIR gate (config$match_dist_threshold_ldir_um).
  src <- data.frame(particle_id = "L1", x_aligned = 0, y_aligned = 0,
                    feret_max_um = 30, area_um2 = 900, major_um = 30, minor_um = 25,
                    stringsAsFactors = FALSE)
  ref <- data.frame(particle_id = "R1", x_norm = 200, y_norm = 0,
                    feret_max_um = 30, area_um2 = 900, major_um = 30, minor_um = 25,
                    stringsAsFactors = FALSE)
  cfg <- make_config()
  cfg$ldir_force_complete_match <- FALSE   # honour the distance gate

  ftir_res <- match_particles(src, ref, cfg, src_label = "ftir", ref_label = "raman")
  expect_equal(ftir_res$match_stats$n_matched, 0L)

  ldir_res <- match_particles(src, ref, cfg, src_label = "ldir", ref_label = "raman")
  expect_equal(ldir_res$match_stats$n_matched, 1L)
})
