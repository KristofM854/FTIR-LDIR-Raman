# PDF report content layer (shiny_app/global.R).
#
# The two table builders are the SAME code path the Summary tab's HTML tables
# use, so these tests pin the numbers the report and the tab both show. The page
# builders are checked for the properties that matter when a page is printed to
# a device: they return something drawable, and a bad page cannot take the whole
# report down with it.

.load_report_helpers <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  skip_if_not_installed("ggplot2")
  env <- new.env(parent = globalenv())
  txt <- paste(readLines(gpath, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  # The table builders classify materials; that lives in the pipeline module.
  sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = env)
  env
}

.pdf_particles <- function(materials, n = length(materials),
                           feret = seq(20, 200, length.out = n),
                           status = "matched") {
  data.frame(material     = rep_len(materials, n),
             feret_max    = rep_len(feret, n),
             match_status = rep_len(status, n),
             stringsAsFactors = FALSE)
}

# ---- report_plastics_table --------------------------------------------------

test_that("the plastics table counts families per device and totals them", {
  env <- .load_report_helpers()
  devices <- list(
    "FTIR (PerkinElmer)" = .pdf_particles(c("PET", "PET", "PP")),
    "Raman"              = .pdf_particles(c("Polyethylene terephtalate", "Polypropylene")))
  tbl <- env$report_plastics_table(devices)

  expect_true(all(c("Family", "Category", "FTIR (PerkinElmer)", "Raman") %in% names(tbl)))
  expect_identical(tbl$Family[nrow(tbl)], "Total")

  pet <- tbl[tbl$Family == "PET", ]
  expect_equal(pet$`FTIR (PerkinElmer)`, 2L)
  # Cross-instrument nomenclature is canonicalised, so PET lands in one row.
  expect_equal(pet$Raman, 1L)
  expect_equal(tbl[tbl$Family == "PP", ]$`FTIR (PerkinElmer)`, 1L)

  # Total is the per-device sum over the families shown.
  tot <- tbl[nrow(tbl), ]
  expect_equal(tot$`FTIR (PerkinElmer)`, 3L)
  expect_equal(tot$Raman, 2L)
})

test_that("the plastics table drops empty devices and excludes Unknown", {
  env <- .load_report_helpers()
  tbl <- env$report_plastics_table(list(
    "FTIR (PerkinElmer)" = .pdf_particles(c("PET", "PET")),
    "Raman"              = NULL,
    "LDIR"               = .pdf_particles(character(0), n = 0)))
  expect_false("Raman" %in% names(tbl))
  expect_false("LDIR"  %in% names(tbl))

  # summarise_plastics() excludes Unknown, so the Total must not count it.
  tbl2 <- env$report_plastics_table(list(
    A = .pdf_particles(c("PET", "definitely not a polymer", "PP"))))
  expect_false("Unknown" %in% tbl2$Family)
  expect_equal(tbl2$A[tbl2$Family == "Total"], 2L)
})

test_that("the plastics table orders synthetic families before semi-synthetic", {
  env <- .load_report_helpers()
  tbl <- env$report_plastics_table(list(A = .pdf_particles(c("Cellulose", "PET"))))
  body <- tbl[tbl$Family != "Total", ]
  expect_equal(body$Category, c("Synthetic", "Semi-synthetic"))
})

test_that("no classifiable data yields a zero-row table, not an error", {
  env <- .load_report_helpers()
  expect_equal(nrow(env$report_plastics_table(list())), 0L)
  expect_equal(nrow(env$report_plastics_table(list(A = NULL, B = NULL))), 0L)
})

# ---- report_size_stats_table -----------------------------------------------

test_that("size statistics match a hand computation and carry units", {
  env <- .load_report_helpers()
  dfs <- list(ftir = .pdf_particles(rep("PET", 4), feret = c(10, 20, 30, 40),
                                    status = c("matched", "matched", "unmatched", "matched")))
  tbl <- env$report_size_stats_table(dfs)
  expect_equal(nrow(tbl), 1L)
  expect_equal(tbl$Instrument, "FTIR (PerkinElmer)")
  expect_equal(tbl$Total, 4L)
  # No passes_cutoff column: every particle counts as above the cutoff
  expect_equal(tbl[["Above quality cutoff"]], 4L)
  expect_equal(tbl[["Below quality cutoff"]], 0L)
  dfs$ftir$passes_cutoff <- c(TRUE, TRUE, FALSE, TRUE)
  tbl <- env$report_size_stats_table(dfs)
  expect_equal(tbl[["Above quality cutoff"]], 3L)
  expect_equal(tbl[["Below quality cutoff"]], 1L)
  expect_match(tbl$Mean,   "^25\\b")
  expect_match(tbl$Median, "^25\\b")
  expect_match(tbl$Range,  "^10")
  # Units survive as real micro signs, not mojibake.
  expect_true(all(grepl("µm$", c(tbl$Mean, tbl$Median, tbl$Range))))
  expect_true(grepl("–", tbl$Range))
})

test_that("size statistics skip instruments with no data", {
  env <- .load_report_helpers()
  tbl <- env$report_size_stats_table(list(
    ftir  = .pdf_particles(rep("PET", 3)),
    raman = NULL,
    ldir  = .pdf_particles(character(0), n = 0)))
  expect_equal(tbl$Instrument, "FTIR (PerkinElmer)")
  expect_equal(nrow(env$report_size_stats_table(list())), 0L)
  # FTIR (Bruker) has its own row
  tbl <- env$report_size_stats_table(list(ftir_bruker = .pdf_particles(rep("PET", 2))))
  expect_equal(tbl$Instrument, "FTIR (Bruker)")
})

test_that("matching quality table: one row per pair from the run's files", {
  env <- .load_report_helpers()
  run <- withr::local_tempdir()
  dir.create(file.path(run, "07_diagnostics")); dir.create(file.path(run, "05_matches"))
  dir.create(file.path(run, "06_agreement"))
  utils::write.csv(data.frame(pair = c("ldir_raman", "ftir_bruker_ftir_perkin"), n_matched = c(4, 2),
                              gate_um = c(75, 50), p_chance_partner = 0.05,
                              chance_pairs_expected = c(1, 0.1), chance_share = c(0.25, 0.05),
                              alignment_determined = c(TRUE, FALSE),
                              position_sigma_um = c(7.3, NA)),
                   file.path(run, "07_diagnostics", "chance_pairs.csv"), row.names = FALSE)
  utils::write.csv(data.frame(match_distance = c(2, 4, 6, 8)),
                   file.path(run, "05_matches", "matched_ldir_raman.csv"), row.names = FALSE)
  utils::write.csv(data.frame(tier = c("Exact", "Family", "Disagree", "Unresolved")),
                   file.path(run, "06_agreement", "agreement_pairwise_ldir_raman.csv"),
                   row.names = FALSE)
  t <- env$report_matching_table(run)
  expect_equal(nrow(t), 2)
  expect_equal(t$Pair, c("LDIR \u2013 Raman", "FTIR (Bruker) \u2013 FTIR (PE)"))
  expect_equal(t[["Dist. med / p90"]][1], "5.0 / 7.4")
  expect_equal(t[["By chance"]][1], "1.0 (25%)")
  expect_equal(t$Scored[1], "4")
  expect_equal(t[["Family agr. (95% CI)"]][1], "67% (21-94%)")   # 2 of 3 resolved
  expect_equal(t$Scored[2], "\u2013")                  # no agreement table
  expect_equal(t[["Sigma (um)"]], c("7.3", "\u2013"))
  expect_equal(t[["Family agr. (95% CI)"]][2], "not determined")
  expect_equal(t[["Family agr., p_true>=0.9"]][1], "-")  # no p_true column
  expect_equal(nrow(env$report_matching_table(withr::local_tempdir())), 0)
})

test_that("detection table: one row per direction, a cell per size class", {
  env <- .load_report_helpers()
  run <- withr::local_tempdir()
  dir.create(file.path(run, "07_diagnostics"))
  d <- expand.grid(size_class = c("20-50", "all sizes"), group = c("all", "plastics"),
                   stringsAsFactors = FALSE)
  d <- rbind(cbind(found_by = "ldir", of = "raman", d), cbind(found_by = "raman", of = "ldir", d))
  d$n_base <- 10; d$n_found <- c(5, 8, 5, 8, 10, 10, 10, 10)
  d$rate_pct <- 100 * d$n_found / d$n_base
  utils::write.csv(d, file.path(run, "07_diagnostics", "detection_rates.csv"), row.names = FALSE)
  t <- env$report_detection_table(run, "all")
  expect_equal(t$Direction, c("LDIR finds Raman", "Raman finds LDIR"))
  expect_equal(t[["20-50 um"]][1], "50% (24-76%) n=10")
  expect_equal(t[["all sizes"]][2], "100% (72-100%) n=10")
  expect_equal(length(env$report_detection_pages(run)), 2)
  expect_equal(nrow(env$report_detection_table(withr::local_tempdir())), 0)
})

# ---- page builders ----------------------------------------------------------

test_that("text pages keep monospace alignment and wrap only long lines", {
  env <- .load_report_helpers()
  # strwrap() collapses runs of spaces; the padded label columns must survive.
  aligned <- "Run                 : 2026-03-09_3"
  expect_identical(env$.report_wrap(aligned), aligned)
  long <- paste(rep("word", 60), collapse = " ")
  expect_gt(length(env$.report_wrap(long, width = 40)), 1L)
  expect_identical(env$.report_wrap(""), "")
  expect_length(env$.report_wrap(character(0)), 0L)
  expect_s3_class(env$report_text_page("T", aligned, "sub"), "ggplot")
})

test_that("table pages are drawable and degrade to text when empty", {
  env <- .load_report_helpers()
  skip_if_not_installed("gridExtra")
  df <- data.frame(A = 1:3, B = c("x", "y", "z"))
  pg <- env$report_table_page(df, "Title", "caption")
  expect_true(inherits(pg, "gtable") || inherits(pg, "grob") || inherits(pg, "ggplot"))
  # Empty input must still produce a page rather than NULL or an error.
  expect_s3_class(env$report_table_page(df[0, ], "Title"), "ggplot")
  expect_s3_class(env$report_table_page(NULL, "Title"), "ggplot")
})

test_that("figure and grid pages pass plots through and tolerate NULL", {
  env <- .load_report_helpers()
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(ggplot2::aes(x, y))
  expect_s3_class(env$report_figure_page(p, "Title", "caption"), "ggplot")
  expect_null(env$report_figure_page(NULL, "Title"))
  skip_if_not_installed("gridExtra")
  expect_false(is.null(env$report_grid_page(list(p, NULL, p), "Grid", "cap")))
  expect_null(env$report_grid_page(list(NULL, NULL), "Grid"))
  expect_null(env$report_grid_page(list(), "Grid"))
})

# ---- write_report_pdf -------------------------------------------------------

test_that("a multi-page PDF is written with one page per element", {
  env <- .load_report_helpers()
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(ggplot2::aes(x, y))
  f <- tempfile(fileext = ".pdf"); on.exit(unlink(f), add = TRUE)
  n <- env$write_report_pdf(list(env$report_text_page("A", "body"), p, p), f)
  expect_equal(n, 3L)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 1000)
  # A PDF, not an empty or truncated file.
  expect_identical(rawToChar(readBin(f, "raw", 4L)), "%PDF")
})

test_that("NULL pages are skipped and an empty report still produces a PDF", {
  env <- .load_report_helpers()
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(ggplot2::aes(x, y))
  f <- tempfile(fileext = ".pdf"); on.exit(unlink(f), add = TRUE)
  expect_equal(env$write_report_pdf(list(NULL, p, NULL), f), 1L)
  g <- tempfile(fileext = ".pdf"); on.exit(unlink(g), add = TRUE)
  # No pages -> a one-page "nothing to report" PDF rather than a corrupt file.
  expect_equal(env$write_report_pdf(list(), g), 1L)
  expect_identical(rawToChar(readBin(g, "raw", 4L)), "%PDF")
})

test_that("one failing page does not cost the reader the rest of the report", {
  env <- .load_report_helpers()
  good <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(ggplot2::aes(x, y))
  bad  <- structure(list(), class = "ggplot")   # errors when printed
  f <- tempfile(fileext = ".pdf"); on.exit(unlink(f), add = TRUE)
  n <- suppressMessages(env$write_report_pdf(list(good, bad, good), f))
  expect_equal(n, 3L)                            # placeholder page, not a gap
  expect_gt(file.size(f), 1000)
  expect_identical(rawToChar(readBin(f, "raw", 4L)), "%PDF")
})

test_that("every material family gets a colour, so no stack is drawn invisible", {
  env <- .load_report_helpers()
  fams <- c(env$synthetic_families, env$semi_synthetic_families, env$natural_families,
            env$inorganic_families, "Unknown", "Something new")
  expect_false(anyNA(env$family_colour(fams)))
  cts <- list(LDIR = table(c("Natural", "Chitin", "Chitin", "Protein", "Mineral", "PE")))
  p <- env$build_material_barplot_gg(cts)
  b <- ggplot2::ggplot_build(p)
  expect_false(anyNA(b$data[[1]]$fill))
})

test_that("the title page lists input file names in full, within the page width", {
  env <- .load_report_helpers()
  long <- paste0(strrep("Comparstic2026 BORA2025 240314 Faanui SW1L R4 A ", 2), "260727_lumos.bmp")
  m <- list(inputs = list(
    ftir = list(basename = "Sample A 260805_stitched.csv", format = "CSV", size_bytes = 1400),
    ftir_bruker_image = list(basename = long, format = "BMP", size_bytes = 1.5e7, width = 2048, height = 1949)),
    config_snapshot = list(ftir_segment_tables = paste0("/x/Sample A 260805M", 1:6, ".csv")))
  l <- env$report_input_file_lines(m)
  expect_true(all(nchar(l) <= 92))
  expect_false(any(grepl("...", l, fixed = TRUE)))
  expect_identical(gsub("\\s+", " ", paste(trimws(l[grepl("Comparstic|lumos", l)]), collapse = " ")),
                   gsub("\\s+", " ", long))
  expect_true(any(grepl("stitched from 6 segments: Sample A 260805M1-M6.csv", l, fixed = TRUE)))
})

test_that("a wide table is fitted to the page", {
  env <- .load_report_helpers()
  skip_if_not_installed("gridExtra")
  cell <- "82% (77-86%) n=280, corr. 64%"
  df <- data.frame(Direction = "Raman finds LDIR", a = cell, b = cell, c = cell, d = cell, e = cell,
                   f = cell, g = cell)
  g <- env$report_table_page(df, "Detection")
  tg <- g$grobs[[1]]$grobs[[1]]
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_lte(grid::convertWidth(sum(tg$widths), "in", valueOnly = TRUE), 10.2)
})


# ---- Summary: alignment check and agreement (from a run folder) -------------

.fake_agreement_run <- function() {
  skip_if_not_installed("withr")
  suppressPackageStartupMessages(library(ggplot2))
  d <- withr::local_tempdir(.local_envir = parent.frame())
  for (s in c("04_alignment", "05_matches", "07_diagnostics")) dir.create(file.path(d, s))
  w <- function(x, f) utils::write.csv(x, file.path(d, f), row.names = FALSE)
  w(data.frame(instrument = c("ldir", "ftir_perkin"), reference = "Raman",
               method = c("landmark_ransac", "ransac"), scale = 1, rotation_deg = c(-178.4, 2),
               reflected = FALSE, icp_rms_um = c(14.1, 30), pose_p = c(4e-22, 0.2),
               determined = c(TRUE, FALSE), n_matched = c(2, 1)),
    "04_alignment/alignment_summary.csv")
  w(data.frame(pair = c("ldir_raman", "ftir_perkin_raman"), alignment_determined = c(TRUE, FALSE),
               n_matched = c(2, 1), gate_um = 20, p_chance_partner = 0.01,
               chance_pairs_expected = 0.1, chance_share = 0.05, position_sigma_um = c(7.1, NA)),
    "07_diagnostics/chance_pairs.csv")
  w(data.frame(found_by = c("raman", "ldir"), of = c("ldir", "raman"), group = "all",
               size_class = "all sizes", n_base = c(4, 3), n_found = c(2, 2),
               rate_pct = c(50, 66.7), ci_low_pct = c(15, 20), ci_high_pct = c(85, 94),
               n_found_true = 2, rate_true_pct = 50, gate_um = 20, n_detected = c(3, 2),
               detected_pct = c(75, 66.7), det_ci_low_pct = 30, det_ci_high_pct = 95,
               det_chance_expected = 0.1, detected_corrected_pct = 75),
    "07_diagnostics/detection_rates.csv")
  w(data.frame(of = "ldir", found_by = "raman", particle_id = c("A1", "A2", "A3", "A4"),
               passes = c(TRUE, TRUE, TRUE, FALSE),
               status = c("matched", "matched", "partner_below_cutoff", "this_below_cutoff")),
    "05_matches/partner_status.csv")
  # LDIR A1-A3 pass (quality >= 0.85), A4 not; Raman R1, R2 pass, R3 below 70
  w(data.frame(particle_id = paste0("A", 1:4), quality = c(0.9, 0.95, 0.9, 0.5),
               feret_max_um = 100, material = c("PET", "PE", "PP", "PET")),
    "04_alignment/aligned_ldir.csv")
  w(data.frame(particle_id = paste0("R", 1:3), quality = c(90, 80, 50), feret_max_um = 100,
               material = c("Polyethylene terephthalate", "Cellulose", "PP")),
    "04_alignment/aligned_raman.csv")
  w(data.frame(ldir_particle_id = c("A1", "A2"), ldir_material = c("PET", "PE"),
               raman_particle_id = c("R1", "R2"),
               raman_material = c("Polyethylene terephthalate", "Cellulose"),
               within_gate = TRUE),
    "05_matches/matched_ldir_raman.csv")
  # Pairs through an undetermined alignment are left out
  w(data.frame(ftir_particle_id = "F1", ftir_material = "PET", raman_particle_id = "R3",
               raman_material = "PP"), "05_matches/matched_ftir_perkin_raman.csv")
  d
}

test_that("summary tables read the alignment, pairs and partner status of a run", {
  env <- .load_report_helpers()
  run <- .fake_agreement_run()
  al <- env$summary_alignment_table(run)
  expect_equal(al$Instrument, c("LDIR", "FTIR (PerkinElmer)"))
  expect_equal(al$Verdict, c("fixed by the data", "NOT DETERMINED"))
  expect_equal(al[[8]][1], "7.1")                           # positional error sigma
  pm <- env$summary_pair_matrix(run)
  expect_equal(unname(pm$identified_pct["LDIR", "Raman"]), 50)
  expect_equal(unname(pm$detected_pct["LDIR", "Raman"]), 75)
  expect_true(is.na(pm$identified_pct["LDIR", "LDIR"]))
  # Only the determined pair is used
  expect_named(env$run_pair_matches(run), "ldir_raman")
  tb <- env$summary_partner_status_counts(run, passing_only = TRUE)
  expect_equal(sum(tb$n), 3L)
  expect_equal(tb$n[tb$status == "matched"], 2L)
  expect_s3_class(env$plot_partner_status(tb), "ggplot")
})

test_that("particle groups count each combination of instruments once", {
  env <- .load_report_helpers()
  run <- .fake_agreement_run()
  g <- env$summary_particle_groups(run)
  # Pairs A1-R1 and A2-R2; A3 alone (LDIR); A4 and R3 below their cutoff
  expect_equal(g$n[g$combo == "Raman + LDIR"], 2L)
  expect_equal(g$n[g$combo == "LDIR"], 1L)
  expect_false("Raman" %in% g$combo)
  # A2 (PE) - R2 (Cellulose) has a synthetic member; A3 (PP) is synthetic too
  gs <- env$summary_particle_groups(run, synthetic_only = TRUE)
  expect_equal(sum(gs$n), 3L)
  expect_s3_class(env$plot_particle_groups(g), "ggplot")
})

test_that("material concordance cross-tabulates the families of a pair", {
  env <- .load_report_helpers()
  run <- .fake_agreement_run()
  expect_equal(unname(env$summary_concordance_choices(run)), "ldir_raman")
  tb <- env$summary_concordance(run, "ldir_raman")
  expect_equal(sum(tb), 2L)
  expect_equal(names(dimnames(tb)), c("LDIR", "Raman"))
  expect_equal(tb["PET", "PET"], 1L)
  pg <- env$report_agreement_pages(run)
  expect_true(length(pg$main) >= 3)
  expect_length(pg$appendix, 1L)
})

test_that("the supporting-files table lists every side file of the run's instruments", {
  env <- .load_report_helpers()
  d <- withr::local_tempdir()
  dir.create(file.path(d, "00_manifest"))
  utils::write.csv(data.frame(kind = "ldir_sam", file = "run_(1).sam", path = "x/run_(1).sam",
                              used_for = "exact join"),
                   file.path(d, "00_manifest", "supporting_files.csv"), row.names = FALSE)
  sf <- env$supporting_files_table(d, c("raman", "ldir"))
  expect_setequal(sf$kind, c("raman_placement", "ldir_analyzed_image", "ldir_sam", "ldir_coord_swaps"))
  expect_equal(sf$file[sf$kind == "ldir_sam"], "run_(1).sam")
  expect_true(is.na(sf$file[sf$kind == "raman_placement"]))
})
