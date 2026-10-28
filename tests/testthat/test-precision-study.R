# Precision study after ISO 5725-2 (R/precision_study.R): counts per filter,
# n replicates on each of p instruments.

for (f in c("utils.R", "particle_contract.R", "00_config.R", "01_ingest.R", "01c_ingest_ldir.R",
            "01_readers.R", "08b_material_map.R", "08_agreement.R", "align_helpers.R",
            "ftir_segments.R",
            "precision_study.R"))
  sys.source(file.path(REPO_ROOT, "R", f), envir = globalenv())

test_that("variance components equal the one-way ANOVA", {
  set.seed(3)
  y <- c(rpois(3, 50), rpois(3, 60), rpois(3, 45), rpois(3, 52))
  lab <- rep(c("ftir", "ftir_bruker", "raman", "ldir"), each = 3)
  r <- iso5725_stratum(y, lab)
  a <- stats::anova(stats::lm(y ~ factor(lab)))
  msw <- a$`Mean Sq`[2]; msb <- a$`Mean Sq`[1]
  expect_equal(r$stats$s_r, sqrt(msw))
  expect_equal(r$stats$s_L, sqrt(max(0, (msb - msw) / 3)))
  expect_equal(r$stats$s_R^2, r$stats$s_r^2 + r$stats$s_L^2)
  expect_equal(r$stats$r_limit, 2.8 * r$stats$s_r)
  expect_lt(r$stats$s_r_low, r$stats$s_r); expect_gt(r$stats$s_r_high, r$stats$s_r)
  # s_L is 0 when the instruments differ less than their replicates
  r0 <- iso5725_stratum(c(10, 20, 30, 20, 21, 19), rep(c("a", "b"), each = 3))
  expect_equal(r0$stats$s_L, 0); expect_equal(r0$stats$s_R, r0$stats$s_r)
})

test_that("Mandel indicator levels match their null distributions (p = 4, n = 3)", {
  cr <- mandel_critical(4, 3)
  expect_equal(c(cr$h_5, cr$h_1), c(1.42, 1.49), tolerance = 0.01)   # ISO 5725-2 Table 6
  set.seed(1)
  k <- replicate(40000, { s2 <- apply(matrix(stats::rnorm(12), 3), 2, stats::var)
                          sqrt(s2[1] / mean(s2)) })
  expect_equal(unname(stats::quantile(k, c(0.95, 0.99))), c(cr$k_5, cr$k_1), tolerance = 0.02)
})

test_that("counts per stratum, zeros filled, quality and size cutoffs applied", {
  cfg <- make_config()
  mk <- function(feret, mat, q) data.frame(particle_id = seq_along(feret), feret_max_um = feret,
                                           material = mat, quality = q, stringsAsFactors = FALSE)
  parts <- list(
    "ftir|1"  = mk(c(10, 30, 60, 150, 300), c("PET", "PET", "Polypropylene", "Cellulose", "PET"),
                   c(0.9, 0.9, 0.9, 0.5, 0.95)),
    "ftir|2"  = mk(c(30, 60, 150), c("PET", "PET", "Protein"), c(0.9, 0.9, 0.9)),
    "raman|1" = mk(c(25, 70), c("PET", "Polypropylene"), c(90, 60)),
    "raman|2" = mk(c(25, 70, 80), c("PET", "PET", "Polypropylene"), c(90, 80, 75)))
  des <- read_precision_design(data.frame(instrument = c("FTIR (PerkinElmer)", "ftir", "Raman", "raman"),
                                          replicate = c(1, 2, 1, 2), file = NA))
  cnt <- precision_counts(des, cfg, particles = parts)
  g <- function(inst, rep, grp, sc = "all sizes")
    cnt$count[cnt$instrument == inst & cnt$replicate == rep & cnt$group == grp & cnt$size_class == sc]
  expect_equal(g("ftir", "1", "All particles"), 4)          # the 10 um particle is below 20 um
  expect_equal(g("ftir", "1", "PET"), 2)
  expect_equal(g("ftir", "1", "Plastics"), 3)               # cellulose fails AAU 0.7
  expect_equal(g("ftir", "1", "PET", ">200"), 1)
  expect_equal(g("raman", "1", "PP"), 0)                    # HQI 60 < 70
  expect_equal(g("raman", "2", "Protein"), 0)               # zero filled
  st <- precision_iso5725(cnt)
  expect_setequal(unique(st$instruments$instrument), c("ftir", "raman"))
  expect_true(all(c("s_r", "s_L", "s_R", "cv_R_pct") %in% names(st$iso)))
})

test_that("the design names instruments flexibly and rejects duplicates", {
  d <- read_precision_design(data.frame(instrument = c("PerkinElmer", "Bruker", "LDIR", "raman"),
                                        replicate = 1, file = "x"))
  expect_equal(d$instrument, c("ftir", "ftir_bruker", "ldir", "raman"))
  expect_error(read_precision_design(data.frame(instrument = c("ldir", "ldir"), replicate = 1,
                                                file = "x")), "twice")
})
