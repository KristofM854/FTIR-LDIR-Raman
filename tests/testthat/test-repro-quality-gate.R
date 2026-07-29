# Multi-Run identification-quality gate (shiny_app/app.R, repro_filtered_pts()).
#
# The gate is an INCLUSION criterion, not a display toggle like run visibility:
# dropping a run's call for a particle means that particle was not identified
# in that run at the chosen threshold. So the per-consensus flags
# (n_runs_detected, material_concordant) are recomputed from the surviving
# rows. Skipping that recompute would leave the rings describing the
# unfiltered analysis while the plot showed the filtered one — a particle whose
# only good call was dropped would still be ringed as reproduced everywhere.
#
# The logic is inside server(), so these tests pin the transformation itself.

.q_gate <- function(pts, qr) {
  if (!is.null(qr) && length(qr) == 2 && "quality" %in% names(pts)) {
    keep <- is.na(pts$quality) | (pts$quality >= qr[1] & pts$quality <= qr[2])
    if (!all(keep)) {
      pts <- pts[keep, , drop = FALSE]
      if (nrow(pts) == 0) return(pts)
      cid <- as.character(pts$consensus_id)
      det <- tapply(pts$run, pts$consensus_id, function(r) length(unique(r)))
      pts$n_runs_detected <- as.integer(det[cid])
      con <- tapply(pts$material_family, pts$consensus_id,
                    function(f) length(unique(f[!is.na(f)])) <= 1L)
      pts$material_concordant <- as.logical(con[cid])
    }
  }
  pts
}

# Mirrors the viewer's status rule.
.q_status <- function(p, n_runs) {
  ifelse(!as.logical(p$material_concordant), "discordant",
         ifelse(p$n_runs_detected < n_runs, "missing", "consensus"))
}

# Three runs. c1 reproduces everywhere but run 3's call is low quality; c2 is
# discordant and the odd one out is the low-quality call; c3 carries no quality
# at all (an ingester that never reported it).
.q_fixture <- function() {
  data.frame(
    consensus_id    = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    run             = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
    quality         = c(.95, .92, .30, .90, .91, .25, NA, NA, NA),
    material_family = c("PP", "PP", "PP", "PP", "PP", "PET", "PE", "PE", "PE"),
    n_runs_detected = 3L,
    material_concordant = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
                            TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE)
}

test_that("a wide-open slider is a no-op", {
  pts <- .q_fixture()
  expect_identical(.q_gate(pts, c(0, 1)), pts)
  expect_identical(.q_gate(pts, NULL), pts)
})

test_that("particles below the threshold stop counting as detected", {
  g <- .q_gate(.q_fixture(), c(0.5, 1))
  c1 <- g[g$consensus_id == 1, ]
  expect_equal(nrow(g), 7L)                       # two low-quality rows gone
  expect_equal(unique(c1$n_runs_detected), 2L)    # was 3
  expect_true(all(.q_status(c1, 3) == "missing")) # ring follows the filter
})

test_that("a disagreement carried only by a filtered-out call clears", {
  g  <- .q_gate(.q_fixture(), c(0.5, 1))
  c2 <- g[g$consensus_id == 2, ]
  expect_true(all(c2$material_concordant))          # PET call was the bad one
  expect_true(all(.q_status(c2, 3) == "missing"))   # not "discordant" any more
})

test_that("rows without a quality value are kept and left alone", {
  # An instrument or run that never reported quality must not vanish just
  # because the slider moved.
  g  <- .q_gate(.q_fixture(), c(0.5, 1))
  c3 <- g[g$consensus_id == 3, ]
  expect_equal(nrow(c3), 3L)
  expect_equal(unique(c3$n_runs_detected), 3L)
  expect_true(all(.q_status(c3, 3) == "consensus"))
})

test_that("both ends of the slider gate", {
  pts <- .q_fixture()
  expect_equal(nrow(.q_gate(pts, c(2, 3))), 3L)      # only the NA group left
  expect_equal(nrow(.q_gate(pts, c(0, 0.5))), 5L)    # two low + three NA
})

test_that("the gate can empty the selection, and says so upstream", {
  pts <- .q_fixture()
  pts <- pts[!is.na(pts$quality), ]                  # drop the NA-quality group
  expect_equal(nrow(.q_gate(pts, c(2, 3))), 0L)
})

test_that("the reproducibility points table carries quality through", {
  # The viewer can only offer the slider if the writer emits the column.
  rpath <- file.path(REPO_ROOT, "R", "reproducibility.R")
  skip_if_not(file.exists(rpath))
  src <- paste(readLines(rpath, warn = FALSE), collapse = "\n")
  expect_match(src, "quality\\s*=\\s*if \\(\"quality\" %in% names\\(d\\)\\)")
})
