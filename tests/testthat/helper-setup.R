# =============================================================================
# helper-setup.R — sourced automatically by testthat before every test file.
# =============================================================================
# The project is run by sourcing main.R (not installed as a package), so the
# tests source the individual R/ modules under test directly. None of these
# modules has source-time side effects beyond defining functions (the reticulate
# block in 00_config.R is guarded by requireNamespace and is a no-op when
# reticulate is absent).

.find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(d, "main.R")) && dir.exists(file.path(d, "R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate repo root (expected main.R + R/) upward from ", start)
}

REPO_ROOT <- .find_repo_root()

local({
  rdir <- file.path(REPO_ROOT, "R")
  # Order matters only for %||% / helpers used at source time; these four are
  # self-contained given base R + RANN + clue.
  for (m in c("utils.R", "00_config.R", "04_ransac.R", "07_match.R")) {
    sys.source(file.path(rdir, m), envir = globalenv())
  }
})

# --- Synthetic fixtures ------------------------------------------------------

#' 3x3 homogeneous similarity transform (scale s, rotation deg, translation).
make_similarity <- function(s, deg, tx, ty, reflect = FALSE) {
  th <- deg * pi / 180
  build_transform_matrix(s * cos(th), s * sin(th), tx, ty, reflect)
}

#' Deterministic pseudo-random 2D point cloud, columns x/y.
make_cloud <- function(n, seed = 1L, span = 1000) {
  set.seed(seed)
  data.frame(x = runif(n, 0, span), y = runif(n, 0, span))
}
