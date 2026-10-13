# Tests for tools/validation.R on synthetic data with known answers: the null
# model finds few chance matches for a true correspondence and many for a
# dense field; the LDIR join check detects a scrambled join.

local({
  for (m in c("08b_material_map.R", "08_agreement.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
  sys.source(file.path(REPO_ROOT, "tools", "validation.R"), envir = globalenv())
})

.field <- function(n, span, seed) {
  set.seed(seed)
  data.frame(particle_id = paste0("P", seq_len(n)),
             x = stats::runif(n, 0, span), y = stats::runif(n, 0, span),
             feret_max_um = round(stats::runif(n, 25, 300)), stringsAsFactors = FALSE)
}

.pair_from_field <- function(f, jitter = 2) {
  ref <- data.frame(particle_id = paste0("R", seq_len(nrow(f))), x_norm = f$x, y_norm = f$y,
                    feret_max_um = f$feret_max_um, feret_min_um = f$feret_max_um / 2,
                    area_um2 = f$feret_max_um^2 / 2, material = "PET", stringsAsFactors = FALSE)
  src <- data.frame(particle_id = paste0("S", seq_len(nrow(f))),
                    x_aligned = f$x + stats::rnorm(nrow(f), 0, jitter),
                    y_aligned = f$y + stats::rnorm(nrow(f), 0, jitter),
                    feret_max_um = f$feret_max_um, feret_min_um = f$feret_max_um / 2,
                    area_um2 = f$feret_max_um^2 / 2, material = "PET", stringsAsFactors = FALSE)
  list(src = src, ref = ref, src_label = "ftir", ref_label = "raman")
}

test_that("null model: sparse field -> chance far below observed; dense field -> near Poisson", {
  cfg <- make_config()
  sparse <- .pair_from_field(.field(60, 5000, 1))
  nm <- suppressMessages(null_model_matches(sparse, cfg, n_iter = 40))
  expect_equal(nm$observed, 60)
  expect_lt(nm$chance_rotation_mean, 10)
  expect_lt(nm$chance_shift_mean, 10)
  # Poisson: 60 * (1 - exp(-60 / 5000^2 * pi * 100^2)) ~= 4.4
  expect_equal(nm$chance_poisson, 60 * (1 - exp(-60 / 5000^2 * pi * 100^2)), tolerance = 0.05)

  dense <- .pair_from_field(.field(400, 2000, 2))
  nd <- suppressMessages(null_model_matches(dense, cfg, n_iter = 20))
  # 400 particles in 2 x 2 mm with a 100 um gate: most would pair by chance
  expect_gt(nd$chance_shift_mean, 200)
  expect_gt(nd$chance_share_pct, 50)
})

test_that("null model leaves the caller's RNG state unchanged", {
  set.seed(7); before <- .Random.seed
  suppressMessages(null_model_matches(.pair_from_field(.field(30, 4000, 3)), make_config(), n_iter = 5))
  expect_identical(.Random.seed, before)
})

test_that("gate sensitivity: baseline gate keeps every baseline pair", {
  cfg <- make_config()
  p <- .pair_from_field(.field(50, 5000, 4))
  p$matched <- suppressMessages(match_particles(p$src, p$ref, cfg))$matched
  g <- suppressMessages(sensitivity_gate(p, cfg, gates = c(5, 100)))
  expect_equal(g$baseline_pairs_kept_pct[g$gate_um == 100], 100)
  expect_lt(g$n_matched[g$gate_um == 5], g$n_matched[g$gate_um == 100])
})

test_that("LDIR join check: correct join is significant, scrambled join is not", {
  set.seed(5)
  n <- 60
  fam <- sample(c("PET", "PP", "PE", "PS"), n, replace = TRUE)
  size <- round(stats::runif(n, 25, 300))
  mk <- function(ldir_idx) list(ldir = list(ldir_raman_match = list(matched = data.frame(
    ldir_particle_id = paste0("L", ldir_idx), raman_particle_id = paste0("R", seq_len(n)),
    match_distance = 5, ldir_material = fam[ldir_idx], raman_material = fam,
    ldir_feret_max_um = size[ldir_idx] * stats::runif(n, 0.9, 1.1), raman_feret_max_um = size,
    stringsAsFactors = FALSE))))
  good <- ldir_join_check(mk(seq_len(n)), n_iter = 200)$summary
  expect_true(all(good$p_value < 0.01))
  expect_true(good$informative[good$statistic == "family_agreement"])
  bad <- ldir_join_check(mk(sample.int(n)), n_iter = 200)$summary
  expect_true(all(bad$p_value > 0.01))
  # A single-polymer sample: the family statistic is flagged uninformative
  fam[] <- "PET"
  one <- ldir_join_check(mk(seq_len(n)), n_iter = 50)$summary
  expect_false(one$informative[one$statistic == "family_agreement"])
})
