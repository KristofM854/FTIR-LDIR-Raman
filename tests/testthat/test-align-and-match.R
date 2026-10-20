# Tests for R/align_and_match.R: the Raman reference frame, align_and_match()
# on a synthetic pair with a known transform, cross-pair matching and the
# triplet quality gate.

local({
  for (m in c("02_prefilter.R", "03_normalize.R", "03b_landmark_align.R",
              "05_transform.R", "08b_material_map.R", "08_agreement.R",
              "align_and_match.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
})

# A Raman field and an FTIR view of the same particles: rotated 30 deg,
# shifted, with 1 um jitter, missing a few particles and with a few extra.
.synthetic_pair <- function() {
  set.seed(42)
  n <- 80
  raman <- data.frame(
    particle_id = paste0("R", seq_len(n)),
    x_um = runif(n, 0, 4000), y_um = runif(n, 0, 4000),
    feret_max_um = round(runif(n, 25, 300)), stringsAsFactors = FALSE)
  raman$feret_min_um <- raman$feret_max_um / 2
  raman$major_um <- raman$feret_max_um; raman$minor_um <- raman$feret_min_um
  raman$area_um2 <- raman$feret_max_um * raman$feret_min_um
  raman$material <- sample(c("PET", "PP", "PE"), n, replace = TRUE)
  raman$quality <- 90
  raman$source_file <- "raman.csv"

  keep <- 6:n                                   # FTIR misses the first five
  th <- 30 * pi / 180
  x <- raman$x_um[keep]; y <- raman$y_um[keep]
  ftir <- raman[keep, ]
  ftir$particle_id <- paste0("F", keep)
  ftir$x_um <- cos(th) * x - sin(th) * y + 1500 + rnorm(length(keep), 0, 1)
  ftir$y_um <- sin(th) * x + cos(th) * y - 700  + rnorm(length(keep), 0, 1)
  ftir$quality <- 0.9
  extra <- ftir[1:3, ]
  extra$particle_id <- paste0("X", 1:3)
  extra$x_um <- c(-3000, -3100, -3200); extra$y_um <- c(9000, 9100, 9200)
  list(raman = raman, ftir = rbind(ftir, extra))
}

.cfg <- function() {
  cfg <- make_config()
  cfg$align_ftir_materials <- NULL
  cfg$align_raman_materials <- NULL
  cfg
}

test_that("prepare_raman_reference: frame centred on the full cleaned cloud", {
  p <- .synthetic_pair()
  ref <- suppressMessages(prepare_raman_reference(p$raman, .cfg()))
  expect_equal(ref$centroid, c(mean(p$raman$x_um), mean(p$raman$y_um)))
  expect_equal(ref$scale, 1)
  expect_equal(mean(ref$clean$x_norm), 0)
  expect_equal(nrow(ref$for_match), nrow(p$raman))
  # Every set is in the same frame
  i <- match(ref$for_transform$particle_id, ref$clean$particle_id)
  expect_equal(ref$for_transform$x_norm, ref$clean$x_norm[i])
})

test_that("cloud_normalization: scale only when asked, errors on empty input", {
  d <- data.frame(x_um = c(0, 2), y_um = c(0, 0))
  expect_equal(cloud_normalization(d)$scale, 1)
  expect_equal(cloud_normalization(d, TRUE)$scale, 1)   # RMS distance 1
  expect_equal(cloud_normalization(data.frame(x_um = c(0, 4), y_um = 0), TRUE)$scale, 2)
  expect_error(cloud_normalization(d[0, ]), "is empty")
  expect_error(cloud_normalization(data.frame(x_um = NA_real_, y_um = 1)), "entirely NA")
})

test_that("align_and_match recovers a known rotation + shift and matches the shared particles", {
  p <- .synthetic_pair()
  cfg <- .cfg()
  ref <- suppressMessages(prepare_raman_reference(p$raman, cfg))
  res <- suppressMessages(align_and_match(p$ftir, ref, cfg, "FTIR (test)"))
  expect_false(res$skipped)
  # FTIR -> Raman is the inverse: rotation -30 deg, scale 1
  expect_equal(res$icp$params$scale, 1, tolerance = 0.01)
  expect_equal(res$icp$params$rotation_deg, -30, tolerance = 0.5)
  m <- res$match$matched
  # Every matched pair is a true correspondence, and nearly all are found
  expect_true(all(sub("^F", "R", m$ftir_particle_id) == m$raman_particle_id))
  expect_gte(nrow(m), 70)
  expect_false(any(grepl("^X", m$ftir_particle_id)))
  expect_lt(stats::median(m$match_distance), 5)
  # All pairs are PET/PP/PE with identical labels: all Exact
  expect_equal(res$agreement$tiered_rates$n_exact, nrow(m))
})

test_that("align_and_match skips a table with too few particles", {
  p <- .synthetic_pair()
  ref <- suppressMessages(prepare_raman_reference(p$raman, .cfg()))
  res <- suppressMessages(align_and_match(p$ftir[1:2, ], ref, .cfg(), "FTIR (test)"))
  expect_true(res$skipped)
  expect_equal(nrow(res$prefiltered), 2)
  expect_null(res$match)
})

test_that("match_aligned_pair matches two aligned tables, NULL when one is missing", {
  a <- data.frame(particle_id = c("a1", "a2"), x_aligned = c(0, 500), y_aligned = c(0, 0),
                  feret_max_um = 50, feret_min_um = 20, area_um2 = 800, material = "PET")
  b <- a; b$particle_id <- c("b1", "b2"); b$x_aligned <- b$x_aligned + 3
  m <- suppressMessages(match_aligned_pair(a, b, make_config(), "ftir_bruker", "ftir_perkin"))
  expect_equal(m$match_stats$n_matched, 2)
  # LDIR on the reference side also gets the LDIR gate (75 um): 65 um apart
  b$x_aligned <- a$x_aligned + 65
  expect_equal(suppressMessages(match_aligned_pair(a, b, make_config(), "ftir_bruker",
                                                   "ftir_perkin"))$match_stats$n_matched, 0)
  expect_equal(suppressMessages(match_aligned_pair(a, b, make_config(), "ftir_bruker",
                                                   "ldir"))$match_stats$n_matched, 2)
  expect_null(match_aligned_pair(a, NULL, make_config(), "ftir_bruker", "ldir"))
  expect_null(match_aligned_pair(a, a[0, ], make_config(), "ftir_bruker", "ldir"))
})

test_that("build_triplets: low-quality and unidentified sides are not scored", {
  fm <- list(matched = data.frame(
    raman_particle_id = c("r1", "r2", "r3", "r4"), ftir_particle_id = c("f1", "f2", "f3", "f4"),
    ftir_material = c("PET", "PET", "PET", "PP"), ftir_quality = c(0.9, 0.5, 0.9, 0.9),
    match_distance = 5))
  lm <- list(matched = data.frame(
    raman_particle_id = c("r1", "r2", "r3", "r4"), ldir_particle_id = c("l1", "l2", "l3", "l4"),
    ldir_material = c("PET", "PET", "unknown", "PP"), ldir_quality = c(0.95, 0.95, 0.95, 0.5),
    match_distance = 8))
  raman <- data.frame(particle_id = c("r1", "r2", "r3", "r4"),
                      material = c("PET", "PET", "PET", "PE"), quality = c(90, 90, 50, 90))
  cfg <- make_config()
  t <- suppressMessages(build_triplets(fm, lm, raman, cfg))
  t <- t[order(t$raman_particle_id), ]
  expect_equal(t$n_instruments_scored, c(3, 2, 1, 2))
  expect_equal(t$material_consensus,
               c("Full agreement",                            # all good
                 "LDIR+Raman agree",                          # FTIR below 0.70
                 "Not scored (quality / no identification)",  # Raman HQI 50, LDIR unknown
                 "No agreement"))                             # LDIR below 0.85; PP vs PE
  expect_equal(t$n_instrument_agreement, c(3, 1, 0, 0))
  expect_equal(nrow(build_triplets(fm, list(matched = lm$matched[0, ]), raman, cfg)), 0)
})

test_that("load.R lists every pipeline module, and each exists", {
  e <- new.env()
  sys.source(file.path(REPO_ROOT, "R", "load.R"), envir = e)
  expect_true(all(file.exists(file.path(REPO_ROOT, "R", e$pipeline_modules))))
  expect_true(all(c("run_pipeline.R", "align_and_match.R", "ldir_stage.R",
                    "particle_contract.R", "01_readers.R") %in% e$pipeline_modules))
  expect_false(anyDuplicated(e$pipeline_modules) > 0)
})

test_that("landmarks are trusted only with enough one-to-one landmark pairs", {
  p <- .synthetic_pair()
  cfg <- .cfg()
  ref <- suppressMessages(prepare_raman_reference(p$raman, cfg))
  src <- apply_normalization(prefilter_ftir(p$ftir, 0, 0),
                             c(mean(p$ftir$x_um), mean(p$ftir$y_um)))
  lm <- suppressMessages(landmark_align(src, ref$for_transform, cfg, src_label = "FTIR"))
  expect_true(lm$confident)
  expect_gte(lm$n_one_to_one, cfg$landmark_min_pairs)
  cfg$landmark_min_pairs <- lm$n_one_to_one + 1
  lm2 <- suppressMessages(landmark_align(src, ref$for_transform, cfg, src_label = "FTIR"))
  expect_false(lm2$confident)
})

test_that("the Tier 2 race is scored on one-to-one pairs", {
  a <- data.frame(x_norm = c(0, 100, 200), y_norm = 0)
  b <- data.frame(x_norm = c(0, 100), y_norm = 0)
  I <- diag(3)
  cfg <- list(ransac_inlier_dist_um = 10)
  expect_equal(tier2_pairs(I, a, b, cfg), 2)         # the third point has no partner
  collapse <- I; collapse[1, 1] <- collapse[2, 2] <- 0   # every point onto (0, 0)
  expect_equal(tier2_pairs(collapse, a, b, cfg), 1)  # many-to-one would count 3
  expect_equal(tier2_pairs(NULL, a, b, cfg), -1)
})

test_that("triplets flag LDIR and FTIR particles that are not each other's direct match", {
  fm <- list(matched = data.frame(raman_particle_id = c("r1", "r2"), ftir_particle_id = c("f1", "f2"),
                                  ftir_material = "PET", ftir_quality = 0.9, match_distance = 5))
  lm <- list(matched = data.frame(raman_particle_id = c("r1", "r2"), ldir_particle_id = c("l1", "l2"),
                                  ldir_material = "PET", ldir_quality = 0.95, match_distance = 5))
  raman <- data.frame(particle_id = c("r1", "r2"), material = "PET", quality = 90)
  direct <- list(matched = data.frame(ldir_particle_id = c("l1", "l2"), ftir_particle_id = c("f1", "f9")))
  t <- suppressMessages(build_triplets(fm, lm, raman, make_config(), ldir_ftir_match = direct))
  t <- t[order(t$raman_particle_id), ]
  expect_equal(t$direct_ldir_ftir_match, c(TRUE, FALSE))
  expect_false("direct_ldir_ftir_match" %in%
               names(suppressMessages(build_triplets(fm, lm, raman, make_config()))))
})

test_that("composite matches never share a reference fragment", {
  # Two adjacent source particles whose search circles both cover fragments R1-R3
  src <- data.frame(particle_id = c("S_big", "S_small"), x_aligned = c(0, 30), y_aligned = 0,
                    major_um = c(200, 120), area_um2 = c(3000, 1500))
  ref <- data.frame(particle_id = paste0("R", 1:5), x_norm = c(10, 20, 30, 60, 70), y_norm = 0,
                    area_um2 = 600)
  cm <- suppressMessages(find_composite_matches(src, ref, make_config()))
  ids <- unlist(strsplit(cm$ref_ids, ";"))
  expect_false(anyDuplicated(ids) > 0)
  expect_equal(cm$src_particle_id[1], "S_big")       # source order kept
})
