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

test_that("LDIR pairing uses its own distance gate", {
  # One pair separated by 35 um: inside the 50 um FTIR gate, outside the
  # 20 um LDIR gate (config$match_dist_threshold_ldir_um, ~3 sigma of the
  # LDIR positional error).
  src <- data.frame(particle_id = "L1", x_aligned = 0, y_aligned = 0,
                    feret_max_um = 30, area_um2 = 900, major_um = 30, minor_um = 25,
                    stringsAsFactors = FALSE)
  ref <- data.frame(particle_id = "R1", x_norm = 35, y_norm = 0,
                    feret_max_um = 30, area_um2 = 900, major_um = 30, minor_um = 25,
                    stringsAsFactors = FALSE)
  cfg <- make_config()
  cfg$ldir_force_complete_match <- FALSE   # honour the distance gate

  ftir_res <- match_particles(src, ref, cfg, src_label = "ftir", ref_label = "raman")
  expect_equal(ftir_res$match_stats$n_matched, 1L)

  ldir_res <- match_particles(src, ref, cfg, src_label = "ldir", ref_label = "raman")
  expect_equal(ldir_res$match_stats$n_matched, 0L)
  ref$x_norm <- 15
  ldir_res <- match_particles(src, ref, cfg, src_label = "ldir", ref_label = "raman")
  expect_equal(ldir_res$match_stats$n_matched, 1L)
})

# Regression for the partner-robbing bug: ldir_force_complete_match = TRUE drops
# the spatial gate and minimises the GLOBAL SUM of pairing costs, so a
# partner-less particle can steal a good particle's correct match. The gate
# (force_complete = FALSE) must leave partner-less particles unmatched instead.
.robbing_fixture <- function() {
  mk_src <- function(id, x, y) data.frame(
    particle_id = id, x_aligned = x, y_aligned = y,
    area_um2 = 1000, feret_max_um = 50, major_um = 50, minor_um = 50,
    stringsAsFactors = FALSE)
  mk_ref <- function(id, x, y) data.frame(
    particle_id = id, x_norm = x, y_norm = y,
    area_um2 = 1000, feret_max_um = 50, major_um = 50, minor_um = 50,
    stringsAsFactors = FALSE)
  list(
    src = rbind(mk_src("S_good", 0, 0), mk_src("S_orphan", 0, 300)),
    # R_perfect coincides with S_good; the only other refs are far away.
    ref = rbind(mk_ref("R_perfect", 0, 0), mk_ref("R_far", 0, 100000))
  )
}

test_that("distance gate leaves a partner-less LDIR particle unmatched", {
  fx  <- .robbing_fixture()
  cfg <- make_config()
  cfg$match_dist_threshold_ldir_um <- 250
  cfg$ldir_force_complete_match    <- FALSE

  res  <- match_particles(fx$src, fx$ref, cfg, src_label = "ldir", ref_label = "raman")
  pair <- setNames(res$matched$raman_particle_id, res$matched$ldir_particle_id)

  # S_good keeps its coincident partner; S_orphan (nearest ref 300 um away,
  # outside the gate) is reported unmatched rather than force-assigned.
  expect_equal(unname(pair["S_good"]), "R_perfect")
  expect_false("S_orphan" %in% names(pair))
  expect_equal(res$match_stats$n_matched, 1L)
})

test_that("force_complete would instead force the orphan onto a spurious partner", {
  fx  <- .robbing_fixture()
  cfg <- make_config()
  cfg$match_dist_threshold_ldir_um <- 250
  cfg$ldir_force_complete_match    <- TRUE   # the disproven setting

  res  <- match_particles(fx$src, fx$ref, cfg, src_label = "ldir", ref_label = "raman")
  # Both particles get forced into pairs — the orphan lands on the far ref, a
  # spurious ~100 mm match. This is exactly the behaviour the default disables.
  expect_equal(res$match_stats$n_matched, 2L)
  orphan <- res$matched[res$matched$ldir_particle_id == "S_orphan", ]
  expect_gt(orphan$match_distance, 1000)
})

test_that("matching per group of nearby particles equals one dense assignment", {
  # The matcher solves each connected group of candidate pairs on its own; the
  # optimum must be the one of a single dense problem (BIG = no pair).
  set.seed(11)
  for (it in 1:15) {
    ns <- sample(5:120, 1); nr <- sample(5:150, 1); W <- sample(c(400, 1500), 1)
    mk <- function(n) data.frame(particle_id = paste0("p", seq_len(n)),
                                 x = runif(n, 0, W), y = runif(n, 0, W),
                                 feret_max_um = runif(n, 10, 200), area_um2 = runif(n, 50, 5000),
                                 major_um = runif(n, 10, 200), minor_um = runif(n, 5, 50))
    a <- mk(ns); b <- mk(nr)
    a$x_aligned <- a$x; a$y_aligned <- a$y; b$x_norm <- b$x; b$y_norm <- b$y
    cfg <- make_config()
    m <- suppressMessages(match_particles(a, b, cfg))$matched
    # Dense reference: the same cost as the matcher, BIG outside the gate
    n <- max(ns, nr); C <- matrix(1e12, n, n)
    asp <- function(d, i) d$major_um[i] / max(d$minor_um[i], 1)
    for (i in seq_len(ns)) for (j in seq_len(nr)) {
      d <- sqrt((a$x[i] - b$x[j])^2 + (a$y[i] - b$y[j])^2)
      thr <- max(cfg$match_dist_threshold_um,
                 cfg_value(cfg, "match_adaptive_dist_factor") * a$major_um[i])
      if (d > thr || sum(sqrt((b$x - a$x[i])^2 + (b$y - a$y[i])^2) < d) >= 20) next
      fr <- a$feret_max_um[i] / b$feret_max_um[j]
      if (max(fr, 1 / fr) > cfg$match_max_size_ratio) next
      sa <- asp(a, i); sb <- asp(b, j)
      C[i, j] <- d + 0.3 * abs(log(a$area_um2[i] / b$area_um2[j])) +
        0.3 * abs(log(a$feret_max_um[i] / b$feret_max_um[j])) +
        0.1 * abs(sa - sb) / max(sa, sb, 1)
    }
    s <- as.integer(clue::solve_LSAP(C))
    ok <- seq_len(ns)[s[seq_len(ns)] <= nr & C[cbind(seq_len(ns), s[seq_len(ns)])] < 1e12]
    expect_equal(nrow(m), length(ok), info = paste("problem", it))
    expect_equal(sum(m$match_score), sum(C[cbind(ok, s[ok])]), tolerance = 1e-9,
                 info = paste("problem", it))
  }
})

test_that("pairs whose sizes differ by more than match_max_size_ratio are not made", {
  # A 148 um particle next to a 6 um one (the PC reference filter's A5) is no
  # pair, however close; within the ratio the pair is made.
  src <- data.frame(particle_id = "S1", x_aligned = 0, y_aligned = 0,
                    feret_max_um = 148, area_um2 = 9000, major_um = 148, minor_um = 60)
  ref <- data.frame(particle_id = "R1", x_norm = 5, y_norm = 0,
                    feret_max_um = 6, area_um2 = 20, major_um = 6, minor_um = 4)
  cfg <- make_config()
  for (meth in c("hungarian", "greedy")) {
    cfg$match_method <- meth
    expect_equal(suppressMessages(match_particles(src, ref, cfg))$match_stats$n_matched, 0L,
                 info = meth)
    ref$feret_max_um <- 40
    expect_equal(suppressMessages(match_particles(src, ref, cfg))$match_stats$n_matched, 1L,
                 info = meth)
    ref$feret_max_um <- 6
  }
  cfg$match_max_size_ratio <- Inf
  expect_equal(suppressMessages(match_particles(src, ref, cfg))$match_stats$n_matched, 1L)
})

test_that("the positional error fit recovers sigma among chance neighbours", {
  set.seed(5)
  W <- 13000; n_ref <- 2400; n_true <- 100; s <- 7
  rx <- runif(n_ref, 0, W); ry <- runif(n_ref, 0, W)
  # 100 true partners offset by N(0, 7 um), 300 sources with no partner
  sx <- c(rx[1:n_true] + rnorm(n_true, 0, s), runif(300, 0, W))
  sy <- c(ry[1:n_true] + rnorm(n_true, 0, s), runif(300, 0, W))
  pe <- .position_error(sx, sy, rx, ry, n_ref / W^2, 75)
  expect_gt(pe$sigma, 5.5); expect_lt(pe$sigma, 8.5)
  # No true partners: no peak, no sigma
  pe0 <- .position_error(runif(400, 0, W), runif(400, 0, W), rx, ry, n_ref / W^2, 75)
  expect_true(is.na(pe0$sigma) || pe0$sigma > 15)
})

test_that("pair_confidence: near pairs are true, far ones in a dense field are not", {
  fit <- list(sigma = 7, w = 0.6, rho = 1.5e-5, window = 75)
  p <- pair_confidence(c(0, 5, 15, 30, 60), fit)
  expect_true(all(diff(p) < 0))                  # falls with distance
  expect_gt(p[2], 0.95); expect_lt(p[5], 0.01)
  # A large particle's centroid is less certain: same distance, higher p
  expect_gt(pair_confidence(30, fit, major = 400, adaptive_factor = 0.15),
            pair_confidence(30, fit, major = 40, adaptive_factor = 0.15))
  expect_true(all(is.na(pair_confidence(c(1, 2), list(sigma = NA)))))
})

test_that("the LDIR gate is set per run from the fitted positional error", {
  set.seed(8)
  W <- 6000; n_ref <- 1500
  ref <- data.frame(particle_id = paste0("R", 1:n_ref), x_norm = runif(n_ref, 0, W),
                    y_norm = runif(n_ref, 0, W), feret_max_um = 40, area_um2 = 1000,
                    major_um = 40, minor_um = 30)
  i <- 1:200                                    # true partners with sigma = 4 um
  src <- data.frame(particle_id = paste0("L", i), x_aligned = ref$x_norm[i] + rnorm(200, 0, 4),
                    y_aligned = ref$y_norm[i] + rnorm(200, 0, 4), feret_max_um = 40,
                    area_um2 = 1000, major_um = 40, minor_um = 30)
  cfg <- make_config()
  msgs <- character(0)
  res <- withCallingHandlers(match_particles(src, ref, cfg, src_label = "ldir", ref_label = "raman"),
                             message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                                     invokeRestart("muffleMessage") })
  expect_gt(res$match_stats$chance_gate_um, 10); expect_lt(res$match_stats$chance_gate_um, 15)
  expect_true(any(grepl("x sigma", msgs)))
  cfg$match_gate_from_sigma <- FALSE
  res2 <- suppressMessages(match_particles(src, ref, cfg, src_label = "ldir", ref_label = "raman"))
  expect_equal(res2$match_stats$chance_gate_um, 20)
})
