# End-to-end: run_pipeline() on a synthetic filter measured by all four
# instruments, every one read through the generic column-mapping reader with a
# KNOWN transform from the Raman frame (rotation, scale, shift, units, y
# orientation). Ground truth for every number: the recovered transforms, which
# particle pairs, the agreement tiers and exclusions, the triplets. Then the
# same run again: the output folder must be byte-identical (apart from
# timestamps and the run's own path).

local({
  sys.source(file.path(REPO_ROOT, "R", "load.R"), envir = globalenv())
  load_pipeline(root = REPO_ROOT)
})

.truth <- function() {
  set.seed(2026)
  n <- 120
  t <- data.frame(k = seq_len(n),
                  x = stats::runif(n, 0, 5000), y = stats::runif(n, 0, 5000),
                  size = round(stats::runif(n, 30, 300)))
  t$material <- c("Polyethylene", "Polypropylene", "PET", "Polystyrene")[(t$k - 1) %% 4 + 1]
  t
}

# Designed differences between FTIR (PerkinElmer) and Raman, all on particles
# both instruments see (k <= 100):
#   D_family  FTIR "HDPE" where Raman says Polyethylene      -> Family tier
#   D_dis     FTIR "Polystyrene" where Raman Polypropylene   -> Disagree
#   Q_low     Raman HQI 50 (below raman_hqi_threshold 70)    -> excluded
#   U         Raman "unknown"                                -> excluded
.sets <- list(D_family = c(1, 5, 9, 13, 17, 21, 25, 29),
              D_dis    = c(2, 6, 10, 14, 18),
              Q_low    = c(3, 7, 11),
              U        = c(4, 8))

.transform_pts <- function(x, y, deg, s, tx, ty) {
  th <- deg * pi / 180
  list(x = s * (cos(th) * x - sin(th) * y) + tx, y = s * (sin(th) * x + cos(th) * y) + ty)
}

.write_inputs <- function(dir) {
  t <- .truth(); S <- .sets
  set.seed(7)
  jit <- function(n) stats::rnorm(n, 0, 1.5)
  # Raman (reference frame): sees everything
  raman <- data.frame(ID = paste0("R", t$k), Xum = t$x, Yum = t$y, Feret = t$size,
                      Mat = t$material, HQI = 92)
  raman$HQI[t$k %in% S$Q_low] <- 50
  raman$Mat[t$k %in% S$U] <- "unknown"
  utils::write.csv(raman, file.path(dir, "raman.csv"), row.names = FALSE)

  # FTIR (PerkinElmer): k <= 100 plus 4 particles outside the Raman field;
  # rotated 178 deg, scale 0.997, shifted; written in MILLIMETRES
  f <- t[t$k <= 100, ]
  p <- .transform_pts(f$x, f$y, 178, 0.997, 3500, 2500)
  ftir <- data.frame(ID = paste0("F", f$k), X_mm = (p$x + jit(nrow(f))) / 1000,
                     Y_mm = (p$y + jit(nrow(f))) / 1000, L_mm = f$size / 1000,
                     Polymer = f$material, AAU = 0.9)
  ftir$Polymer[f$k %in% S$D_family] <- "HDPE"
  ftir$Polymer[f$k %in% S$D_dis] <- "Polystyrene"
  ftir <- rbind(ftir, data.frame(ID = paste0("FX", 1:4), X_mm = 20 + 0:3, Y_mm = 20,
                                 L_mm = 0.1, Polymer = "PET", AAU = 0.9))
  utils::write.csv(ftir, file.path(dir, "ftir.csv"), row.names = FALSE)

  # FTIR (Bruker): k in 11..110, rotated 3 deg, scale 1.002, image-style y (down)
  b <- t[t$k >= 11 & t$k <= 110, ]
  p <- .transform_pts(b$x, b$y, 3, 1.002, -400, 150)
  bruker <- data.frame(ID = paste0("B", b$k), x = p$x + jit(nrow(b)), y = -(p$y + jit(nrow(b))),
                       size = b$size, material = b$material, score = 0.9)
  utils::write.csv(bruker, file.path(dir, "bruker.csv"), row.names = FALSE)

  # LDIR, coordinates in the table: k in 21..120, rotated -90 deg, scale 1.03
  l <- t[t$k >= 21, ]
  p <- .transform_pts(l$x, l$y, -90, 1.03, 6000, 1000)
  ldir <- data.frame(ID = paste0("L", l$k), x = p$x + jit(nrow(l)), y = p$y + jit(nrow(l)),
                     size = l$size * 1.03, Identification = l$material, Quality = 0.95)
  utils::write.csv(ldir, file.path(dir, "ldir.csv"), row.names = FALSE)

  # Bruker's mapping as a JSON file (the other three as R lists)
  jsonlite::write_json(list(
    columns = list(particle_id = "ID", x_um = "x", y_um = "y", feret_max_um = "size",
                   material = "material", quality = "score"),
    y_axis = "down", quality_metric = "AAU"), file.path(dir, "bruker_map.json"),
    auto_unbox = TRUE)
}

.config <- function(dir, out) {
  cfg <- make_config(ftir_path = file.path(dir, "ftir.csv"),
                     raman_path = file.path(dir, "raman.csv"), output_dir = out)
  cfg$ftir_bruker_path <- file.path(dir, "bruker.csv")
  cfg$ldir_path <- file.path(dir, "ldir.csv")
  cfg$raman_column_map <- list(columns = list(particle_id = "ID", x_um = "Xum", y_um = "Yum",
    feret_max_um = "Feret", material = "Mat", quality = "HQI"), quality_metric = "HQI")
  cfg$ftir_column_map <- list(columns = list(particle_id = "ID", x_um = "X_mm", y_um = "Y_mm",
    feret_max_um = "L_mm", material = "Polymer", quality = "AAU"),
    coordinate_unit_um = 1000, size_unit_um = 1000, quality_metric = "AAU")
  cfg$ftir_bruker_column_map <- file.path(dir, "bruker_map.json")
  cfg$ldir_column_map <- list(columns = list(particle_id = "ID", x_um = "x", y_um = "y",
    feret_max_um = "size", material = "Identification", quality = "Quality"))
  cfg
}

# The diagnostics PDF cannot encode "\u2013" in a non-UTF-8 locale (CI, POSIX):
# a cosmetic warning per plot title, muffled here and nothing else.
.run <- function(dir, out) withCallingHandlers(
  suppressMessages(run_pipeline(.config(dir, out), report = FALSE)),
  warning = function(w) if (grepl("conversion failure", conditionMessage(w)))
    invokeRestart("muffleWarning"))

.pair_ok <- function(m, a, b) {
  ka <- sub("^[A-Z]+", "", m[[paste0(a, "_particle_id")]])
  kb <- sub("^[A-Z]+", "", m[[paste0(b, "_particle_id")]])
  ka == kb
}

test_that("synthetic four-instrument run recovers the known transforms, pairs and tiers", {
  dir <- withr::local_tempdir(); out <- withr::local_tempdir()
  .write_inputs(dir)
  r <- .run(dir, out)
  S <- .sets

  # Transforms: each instrument -> Raman is the inverse of what was applied
  pe <- r$ftir$icp$params
  expect_equal(pe$scale, 1 / 0.997, tolerance = 2e-3)
  expect_equal(((pe$rotation_deg + 180) %% 360) - 180, -178, tolerance = 0.2)
  br <- r$ftir_bruker$icp$params
  expect_equal(br$scale, 1 / 1.002, tolerance = 2e-3)
  expect_equal(br$rotation_deg, -3, tolerance = 0.2)
  expect_false(isTRUE(br$reflected))           # y-down was undone by the reader
  ld <- extract_transform_params(r$ldir$ldir_icp$transform)
  expect_equal(ld$scale, 1 / 1.03, tolerance = 5e-3)

  # Pairs: every shared particle found, every pair a true correspondence.
  # The Raman particles below HQI 70 (Q_low: k 3, 7, 11) are not matched
  # (match_quality_filter): FTIR k 1..100 -> 97 pairs, Bruker k 11..110 -> 99.
  m <- r$ftir$match$matched
  expect_equal(nrow(m), 100 - length(S$Q_low))   # FTIR extras unmatched
  expect_true(all(.pair_ok(m, "ftir", "raman")))
  mb <- r$ftir_bruker$match$matched
  expect_equal(nrow(mb), 100 - sum(S$Q_low >= 11))
  expect_true(all(.pair_ok(mb, "ftir", "raman")))
  ml <- r$ldir$ldir_ref_match$matched
  expect_equal(nrow(ml), 100)
  expect_true(all(.pair_ok(ml, "ldir", "raman")))
  bp <- r$cross$ftir_bruker_ftir_perkin$matched
  expect_equal(nrow(bp), 90)                   # k 11..100
  expect_true(all(.pair_ok(bp, "ftir_bruker", "ftir_perkin")))

  # Agreement FTIR (PerkinElmer)-Raman, exactly as designed
  tr <- r$ftir$agreement$tiered_rates
  expect_equal(tr$n_pairs_matched, 100 - length(S$Q_low))
  expect_equal(tr$n_excluded_quality_b, 0)                # not matched in the first place
  expect_false(any(m$raman_particle_id %in% paste0("R", S$Q_low)))
  expect_equal(tr$n_excluded_unidentified, length(S$U))
  expect_equal(tr$n_family, length(S$D_family))
  expect_equal(tr$n_disagree, length(S$D_dis))
  expect_equal(tr$n_exact, 100 - length(unlist(S)))
  expect_equal(tr$n_unresolved, 0)

  # Triplets: Raman particles matched by both FTIR (PerkinElmer) and LDIR
  tp <- r$ldir$triplets
  expect_equal(nrow(tp), 80)                   # k 21..100
  low <- tp$raman_particle_id %in% paste0("R", S$Q_low)
  expect_true(all(tp$n_instruments_scored[!low] == 3))

  # Contract, files, join check
  for (k in c("ftir", "ftir_bruker", "raman", "ldir"))
    expect_silent(validate_particles(r$raw[[k]], what = k))
  run <- r$config$output_dir
  for (f in c("05_matches/matched_ftir_perkin_raman.csv", "05_matches/matched_ftir_bruker_raman.csv",
              "05_matches/matched_ldir_raman.csv", "05_matches/triplets_3way.csv",
              "04_alignment/transform_params_ftir_perkin_raman.txt",
              "07_diagnostics/ldir_join_check.csv", "00_manifest/config_full.json"))
    expect_true(file.exists(file.path(run, f)), info = f)
  jc <- utils::read.csv(file.path(run, "07_diagnostics/ldir_join_check.csv"))
  expect_lt(jc$p_value[jc$statistic == "size_spearman"], 0.01)
})

test_that("the same run twice gives byte-identical outputs", {
  dir <- withr::local_tempdir(); o1 <- withr::local_tempdir(); o2 <- withr::local_tempdir()
  .write_inputs(dir)
  r1 <- .run(dir, o1)$config$output_dir
  r2 <- .run(dir, o2)$config$output_dir
  files <- sort(list.files(r1, recursive = TRUE))
  expect_identical(files, sort(list.files(r2, recursive = TRUE)))
  # Lines that legitimately differ: timestamps, the run's own path, git state
  volatile <- "Generated:|timestamp|updated_at|created|output_dir|debug_dir|run_dir|run_id|Sys.time|git_|R_session|Date|date:"
  differ <- character(0)
  for (f in files) {
    a <- file.path(r1, f); b <- file.path(r2, f)
    if (grepl("\\.(pdf|png|jpg)$", f)) {
      if (!identical(tools::md5sum(a)[[1]], tools::md5sum(b)[[1]])) differ <- c(differ, f)
      next
    }
    la <- readLines(a, warn = FALSE); lb <- readLines(b, warn = FALSE)
    la <- gsub(r1, "<RUN>", la, fixed = TRUE); lb <- gsub(r2, "<RUN>", lb, fixed = TRUE)
    keep <- function(l) l[!grepl(volatile, l)]
    if (!identical(keep(la), keep(lb))) differ <- c(differ, f)
  }
  # PDFs embed their creation time; everything else must match
  expect_identical(setdiff(differ, grep("\\.pdf$", files, value = TRUE)), character(0))
})

test_that("every combination of two to four instruments runs and pairs correctly", {
  dir <- withr::local_tempdir()
  .write_inputs(dir)
  path_key <- c(raman = "raman_path", ftir = "ftir_path",
                ftir_bruker = "ftir_bruker_path", ldir = "ldir_path")
  # Shared particles of each pair (k ranges in .write_inputs)
  k <- list(raman = 1:120, ftir = 1:100, ftir_bruker = 11:110, ldir = 21:120)
  roles <- names(path_key)
  combos <- unlist(lapply(2:4, function(n) combn(roles, n, simplify = FALSE)),
                   recursive = FALSE)
  expect_length(combos, 11)
  # The viewer's loaders (the report reads the run the same way)
  viewer <- new.env(parent = globalenv())
  withr::with_dir(REPO_ROOT, suppressWarnings(
    source(file.path("shiny_app", "global.R"), local = viewer, encoding = "UTF-8")))
  for (cmb in combos) {
    what <- paste(cmb, collapse = "+")
    out <- withr::local_tempdir()
    cfg <- .config(dir, out)
    for (r in setdiff(roles, cmb)) cfg[path_key[[r]]] <- list(NULL)
    r <- withCallingHandlers(
      suppressMessages(run_pipeline(cfg, report = FALSE)),
      warning = function(w) if (grepl("conversion failure", conditionMessage(w)))
        invokeRestart("muffleWarning"))
    ref <- REFERENCE_ORDER[REFERENCE_ORDER %in% cmb][1]
    expect_identical(r$reference, ref, info = what)
    # Every other instrument is paired with the reference, every pair true
    expect_setequal(names(r$pairs),
                    vapply(setdiff(cmb, ref), function(s) pair_labels(s, ref)$stem, ""))
    for (p in r$pairs) {
      m <- p$match$matched
      src <- names(ROLE_FILE)[ROLE_FILE == p$labels$src_file]
      # Every true pair found. Particles of one instrument outside the other's
      # field can make a chance pair (LDIR, 150 um gate, without Raman: about
      # one expected at this density)
      # Raman particles below the HQI cutoff are not matched
      n_true <- length(setdiff(intersect(k[[src]], k[[ref]]),
                               if (ref == "raman") .sets$Q_low))
      ok <- .pair_ok(m, p$labels$src_prefix, p$labels$ref_prefix)
      expect_equal(sum(ok), n_true, info = paste(what, p$labels$stem))
      expect_lte(sum(!ok), 2, label = paste(what, p$labels$stem, "chance pairs"))
      expect_true(file.exists(file.path(r$config$output_dir, "05_matches",
                                        paste0("matched_", p$labels$stem, ".csv"))),
                  info = paste(what, p$labels$stem))
    }
    # Every instrument in the shared frame for the viewer
    for (s in cmb)
      expect_true(file.exists(file.path(r$config$output_dir, "04_alignment",
                                        paste0("aligned_", ROLE_FILE[[s]], ".csv"))),
                  info = paste(what, s))
    expect_true(file.exists(file.path(r$config$output_dir, "07_diagnostics",
                                      "match_statistics.txt")), info = what)
    # The viewer finds the run and shows exactly these instruments, all
    # particles of each
    expect_true(viewer$is_run_dir(r$config$output_dir), info = what)
    dfs <- viewer$build_instrument_dfs(
      viewer$load_run_data(list(dir = r$config$output_dir, format = "subdir")))
    for (s in roles) {
      n <- if (is.null(dfs[[s]])) 0L else nrow(dfs[[s]])
      n_expected <- if (s %in% cmb) length(k[[s]]) + if (s == "ftir") 4L else 0L else 0L
      expect_equal(n, n_expected, info = paste(what, "viewer", s))
    }
  }
})

test_that("fewer than two instruments is an error before anything is written", {
  out <- withr::local_tempdir()
  cfg <- make_config(raman_path = "raman.csv", output_dir = out)
  expect_error(run_pipeline(cfg, report = FALSE), "At least two instruments")
  expect_length(list.files(out), 0)
  expect_error(choose_reference("ldir"), "At least two")
  expect_error(choose_reference(character(0)), "none")
})
