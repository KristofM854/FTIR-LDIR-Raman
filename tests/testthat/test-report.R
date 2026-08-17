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
  expect_equal(tbl$Instrument, "FTIR")
  expect_equal(tbl$Total, 4L)
  expect_equal(tbl$Matched, 3L)
  expect_equal(tbl$Unmatched, 1L)
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
  expect_equal(tbl$Instrument, "FTIR")
  expect_equal(nrow(env$report_size_stats_table(list())), 0L)
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
