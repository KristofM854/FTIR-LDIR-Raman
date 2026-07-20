# Intra-instrument reproducibility engine: registration, spatial linking,
# consensus particles, and concordance/accuracy metrics.

local({
  rdir <- file.path(REPO_ROOT, "R")
  for (m in c("08b_material_map.R", "reproducibility.R"))
    sys.source(file.path(rdir, m), envir = globalenv())
})

# Rigidly move a cloud (rotation + translation + small jitter).
.move <- function(df, deg, tx, ty, jitter = 1.5, seed = 1) {
  set.seed(seed)
  th <- deg * pi / 180
  x <- df$x_um * cos(th) - df$y_um * sin(th) + tx + rnorm(nrow(df), 0, jitter)
  y <- df$x_um * sin(th) + df$y_um * cos(th) + ty + rnorm(nrow(df), 0, jitter)
  df$x_um <- x; df$y_um <- y; df
}

.base_filter <- function(n = 12, seed = 7) {
  set.seed(seed)
  data.frame(
    particle_id  = paste0("p", seq_len(n)),
    x_um = runif(n, 0, 6000), y_um = runif(n, 0, 6000),
    area_um2 = 1000, major_um = 50, minor_um = 40,
    feret_max_um = runif(n, 30, 200),
    material = "Polyethylene terephthalate", quality = 0.9,
    stringsAsFactors = FALSE)
}

test_that("rigid registration recovers a re-seated cloud", {
  a <- .base_filter()
  b <- .move(a, deg = 3, tx = 150, ty = -90)
  reg <- repro_rigid_register(b, a, gate = 400)
  # Aligned b should sit on top of a (within the injected jitter).
  resid <- sqrt((reg$x_aligned - a$x_um)^2 + (reg$y_aligned - a$y_um)^2)
  expect_lt(median(resid), 5)
})

test_that("consensus links particles and surfaces detection + material flips", {
  run1 <- .base_filter()
  run2 <- .move(run1, deg = 3,  tx = 120, ty = -80, seed = 2)
  # Run 3: re-seated the other way, particle p5 not detected, p3 misclassified.
  run3 <- .move(run1[-5, ], deg = -2, tx = -70, ty = 140, seed = 3)
  run3$material[run3$particle_id == "p3"] <- "Polypropylene"

  ref_fam <- classify_family_vec("Polyethylene terephthalate")
  res <- run_reproducibility(list(run1, run2, run3),
                             gate = 50, reference_family = ref_fam)
  s <- res$summary

  # 12 physical particles; p5 seen in 2 of 3 runs, the rest in all 3.
  expect_equal(s$n_consensus, 12L)
  expect_equal(s$detected_in_all, 11L)
  expect_equal(unname(s$detection_breakdown["2"]), 1L)
  expect_equal(unname(s$detection_breakdown["3"]), 11L)

  # One material flip (p3) among the 12 multi-run particles -> 11/12 concordant.
  expect_equal(round(s$material_concordance, 4), round(11 / 12, 4))
  # Accuracy dips below 1 because of the single wrong call, but stays high.
  expect_true(s$accuracy_vs_reference < 1 && s$accuracy_vs_reference > 0.9)

  # Registration is tight, so positional jitter is small.
  expect_lt(s$median_pos_jitter_um, 10)

  # The consensus row for p3 is flagged discordant.
  p3 <- res$consensus[!res$consensus$material_concordant, ]
  expect_equal(nrow(p3), 1L)
})

test_that("run_reproducibility needs at least two runs", {
  expect_error(run_reproducibility(list(.base_filter())))
})

test_that("long table has one row per detected (particle x run) with repro flags", {
  run1 <- .base_filter()
  run2 <- .move(run1, deg = 3,  tx = 120, ty = -80, seed = 2)
  run3 <- .move(run1[-5, ], deg = -2, tx = -70, ty = 140, seed = 3)
  res  <- run_reproducibility(list(run1, run2, run3), gate = 50)
  long <- repro_long_table(res)

  # 12 in all-3 (11) contribute 3 rows, the 2-run particle contributes 2:
  # 11*3 + 1*2 = 35 detections total.
  expect_equal(nrow(long), 35L)
  expect_setequal(unique(long$run), 1:3)
  expect_true(all(c("consensus_id", "x_aligned", "material_family",
                    "n_runs_detected", "material_concordant") %in% names(long)))
  # Every consensus particle's row count equals its n_runs_detected.
  by_c <- tapply(long$run, long$consensus_id, length)
  expect_true(all(by_c == long$n_runs_detected[match(names(by_c), long$consensus_id)]))
})
