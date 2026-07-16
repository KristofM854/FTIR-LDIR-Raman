# LDIR<->Raman acceptance-gate classification (viewer).
#
# Regression guard for the "40/40 LDIR matched" miscount: with
# ldir_force_complete_match = TRUE the pipeline pairs EVERY LDIR particle with a
# Raman particle regardless of distance, so the raw matched CSV always looks
# fully matched. The viewer must reclassify forced over-gate pairs (match_distance
# above the acceptance gate) as unmatched so the summary, overlay and tables agree.
#
# The functions under test live in shiny_app/global.R, which normally boots the
# Shiny app (library()/relative source()). We load only the pure helpers by
# neutralising those side effects and stubbing the one external dependency.

.load_global_helpers <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  env <- new.env(parent = globalenv())
  env$classify_family_vec <- function(x) x
  # Deterministic manifest so the gate resolves to a known value.
  env$load_run_manifest <- function(run_dir)
    list(config_snapshot = list(match_dist_threshold_ldir_um = 250))
  txt <- paste(readLines(gpath, warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  env
}

.synthetic_ldir_match <- function(distances) {
  n <- length(distances)
  data.frame(
    ldir_particle_id  = paste0("L", seq_len(n)),
    raman_particle_id = paste0("R", seq_len(n)),
    ldir_x_aligned = seq_len(n) * 100, ldir_y_aligned = seq_len(n) * 100,
    ldir_x_um = seq_len(n), ldir_y_um = seq_len(n),
    ldir_area_um2 = 1, ldir_major_um = 1, ldir_minor_um = 1, ldir_feret_max_um = 1,
    ldir_material = "PET", ldir_quality = 1,
    match_id = seq_len(n), match_score = 1,
    match_distance = distances,
    stringsAsFactors = FALSE
  )
}

test_that("annotate_ldir_gate flags only within-gate pairs as genuine", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 120, 200, 240, 10, 400, 600, 900))
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")

  expect_equal(data$ldir_match_gate_um, 250)
  expect_equal(data$ldir_raman_matched$within_gate,
               c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("build_instrument_dfs reports over-gate forced pairs as unmatched LDIR", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 120, 200, 240, 10, 400, 600, 900))
  data <- env$annotate_ldir_gate(
    list(ldir_raman_matched = m, unmatched_ldir = NULL,
         matched = NULL, unmatched_raman = NULL), "dummy")

  res <- env$build_instrument_dfs(data)
  # 5 within gate -> matched; 3 over gate -> unmatched; all 8 still present.
  expect_equal(sum(res$ldir$match_status == "matched"), 5L)
  expect_equal(sum(res$ldir$match_status == "unmatched"), 3L)
  expect_equal(nrow(res$ldir), 8L)
  # Over-gate rows must not claim a Raman partner.
  expect_false(any(res$ldir$matched_to_raman[res$ldir$match_status == "unmatched"]))
})

test_that("genuine-pair helpers gate the Raman->LDIR flag", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 400))  # R1 genuine, R2 over-gate
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")

  # Raw unmatched_raman CSV schema (x_norm/x_um/feret_max_um), as loaded.
  raman <- data.frame(particle_id = c("R1", "R2", "R3"),
                      x_norm = 0, y_norm = 0, x_um = 0, y_um = 0,
                      area_um2 = 1, major_um = 1, minor_um = 1, feret_max_um = 1,
                      material = "PET", quality = 1, stringsAsFactors = FALSE)
  data$matched <- NULL
  data$unmatched_raman <- raman
  res <- env$build_instrument_dfs(data)
  flag <- res$raman$matched_to_ldir[match(c("R1", "R2", "R3"), res$raman$particle_id)]
  expect_equal(flag, c(TRUE, FALSE, FALSE))
})

test_that("legacy match frames without match_distance stay fully matched", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 400))
  m$match_distance <- NULL  # simulate a pre-gate export
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")
  expect_true(all(data$ldir_raman_matched$within_gate))
})

test_that("regate_ldir retunes within_gate to an explicit gate", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 120, 200, 240, 10, 400, 600, 900))
  d <- list(ldir_raman_matched = m)

  expect_equal(sum(env$regate_ldir(d, 100)$ldir_raman_matched$within_gate), 2L)
  expect_equal(sum(env$regate_ldir(d, 250)$ldir_raman_matched$within_gate), 5L)
  expect_equal(sum(env$regate_ldir(d, 1000)$ldir_raman_matched$within_gate), 8L)

  # Stored gate follows the override.
  expect_equal(env$regate_ldir(d, 500)$ldir_match_gate_um, 500)

  # Invalid gate is a no-op (no within_gate written, data untouched).
  out <- env$regate_ldir(d, NA)
  expect_false("within_gate" %in% names(out$ldir_raman_matched))
})

test_that("ldir_acceptance_gate falls back to 250 when unrecorded", {
  env <- .load_global_helpers()
  # Override the manifest loader to return no gate keys.
  env$load_run_manifest <- function(run_dir) list(config_snapshot = list())
  expect_equal(env$ldir_acceptance_gate("dummy"), 250)
})
