# =============================================================================
# test-input-classify.R -- recognising input files by their content
# =============================================================================
# main.R's one-dialog input (classify_input_file() / assign_input_roles() in
# R/00b_file_input.R). Two folders:
#   * dataset A copied under neutral names (t1.csv, i1.png, ...), so only the
#     content can decide;
#   * a header-only copy of the Faanui field sample: the real file names and
#     the columns, sheets and image headers of its exports (no particle data).

.load_classifier <- function() {
  skip_if_not_installed("readxl"); skip_if_not_installed("png")
  for (m in c("measure_raman_placement.R", "ftir_segments.R", "00b_file_input.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
}

# Answers questions from a named list (question pattern -> answer) and
# records every question asked
.asker <- function(answers = list()) {
  log <- character(0)
  f <- function(q, choices) {
    log <<- c(log, q)
    for (p in names(answers)) if (grepl(p, q, fixed = TRUE)) return(answers[[p]])
    "not used"
  }
  list(ask = f, asked = function() log)
}

.dataset_a_neutral <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  a <- file.path(REPO_ROOT, "examples", "dataset_A")
  map <- c(FTIR_perkin.csv = "t1.csv", FTIR_bruker.csv = "t2.csv", Raman_data.csv = "t3.csv",
           LDIR_data.xlsx = "t4.xlsx", LDIR_coord_swaps.csv = "t5.csv",
           FTIR_image.png = "i1.png", LDIR_image.png = "i2.png",
           LDIR_image_analyzed.png = "i3.png", Raman_image_small.png = "i4.png",
           README.md = "notes.md")
  for (f in names(map)) file.copy(file.path(a, f), file.path(d, map[[f]]))
  d
}

# A BMP header only (54 bytes): what classify_input_file() reads
.write_bmp_header <- function(path, w, h, bits, ppm) {
  i32 <- function(v) writeBin(as.integer(v), raw(), size = 4, endian = "little")
  i16 <- function(v) writeBin(as.integer(v), raw(), size = 2, endian = "little")
  writeBin(c(charToRaw("BM"), i32(54), i32(0), i32(54), i32(40), i32(w), i32(h), i16(1),
             i16(bits), i32(0), i32(0), i32(ppm), i32(ppm), i32(0), i32(0)), path)
}

.faanui_headers <- function() {
  skip_if_not_installed("writexl")
  d <- withr::local_tempdir(.local_envir = parent.frame())
  f <- function(x) file.path(d, x)
  ftir_cols <- c("Identifier", "Coord. [pixels]", "Coord. [um]", "Max AAU score",
                 "Spectrum yielding score", "Group", "Group ID", "No. of pixels [-]",
                 "Area on map [um2]", "Major dim [um]", "Minor dim [um]",
                 "Feret min [um]", "Volume [um3]", "Mass [ng]")
  raman_cols <- c("Particle Name", "Material", "HQI", "RamanSignal", "FluorescenceSignal",
                  "Area [um2]", "Aspect Ratio", "Feret Max [um]", "Feret Min [um]",
                  "Visual Center Point X [um]", "Visual Center Point Y [um]")
  ldir_cols <- c("#", "Id", "Width (um)", "Height (um)", "Diameter (um)",
                 "Aspect Ratio", "Area (um2)", "Identification", "Quality", "Is Valid")
  writeLines(paste(ftir_cols, collapse = ";"),
             f("Comparstic2026 BORA2025 240314 Faanui SW1L R4 A 260727_lumos.csv"), useBytes = TRUE)
  empty <- function(cols) as.data.frame(stats::setNames(rep(list(character(0)), length(cols)), cols),
                                        check.names = FALSE)
  writexl::write_xlsx(list(Data = empty(ftir_cols)),
                      f("Comparstic2026 BORA2025 240314 Faanui SW1L R4 A 260805_perkin.xlsx"))
  writexl::write_xlsx(list(Info = data.frame(a = "Instrument Name", b = "Agilent 8700 LDIR"),
                           Particles = empty(ldir_cols)),
                      f("Comparstic2026 BORA2025 240314 FAANUI SW1L R4 A Au 260721.xlsx"))
  for (lib in c("IAEA", "STJapan"))
    writeLines(paste(raman_cols, collapse = ";"),
               f(paste0("Comparstic2026 Bora2025 240314 Faanui SW1L R4 A ", lib, " 260924.csv")), useBytes = TRUE)
  writeLines(c("Width: 13000", "Height: 13000", "X: -6500", "Y: 6500", "Z: 0", "Rotation: 0"),
             f("Comparstic2026 Bora2025 240314 Faanui SW1L R4 A IAEA 260924.txt"))
  .write_bmp_header(f("Comparstic2026 BORA2025 240314 Faanui SW1L R4 260805_perkin.bmp"), 2000, 2000, 24, 2000)
  .write_bmp_header(f("Comparstic2026 BORA2025 240314 Faanui SW1L R4 A 260727_lumos.bmp"), 2048, 1949, 32, 3780)
  .write_bmp_header(f("Comparstic2026 Bora2025 240314 Faanui SW1L R4 A IAEA 260928.bmp"), 22512, 22512, 24, 0)
  # LDIR image and its analyzed image (black with a few coloured particles)
  raw_img <- array(runif(60 * 60 * 3, 0, 0.2), c(60, 60, 3))
  ana <- array(0, c(60, 60, 3)); ana[10:15, 10:15, ] <- rep(c(166, 206, 227) / 255, each = 36)
  ana[30:33, 40:44, ] <- rep(c(251, 154, 153) / 255, each = 20)
  png::writePNG(raw_img, f("Comparstic2026 BORA2025 240314 FAANUI SW1L R4 A Au 260721.png"))
  png::writePNG(ana, f("Comparstic2026 BORA2025 240314 FAANUI SW1L R4 A Au 260721_analyzed.png"))
  # .sam: the Clarity name at byte 112, as in the field file
  writeBin(c(as.raw(c(2, 0, 0, 0, 1, 0, 0, 0, 7, 0, 0, 0, 0x1b, 0, 0, 0)), raw(95),
             charToRaw("Agilent.Clarity.Models"), raw(64)),
           f("Comparstic2026 BORA2025 240314 FAANUI SW1L R4 Au 260721 _(1).sam"))
  d
}

test_that("each kind of file is recognised from its content", {
  .load_classifier()
  d <- .dataset_a_neutral()
  k <- function(x) classify_input_file(file.path(d, x))
  expect_equal(k("t1.csv")$kind, "ftir_table"); expect_true(is.na(k("t1.csv")$role))
  expect_equal(k("t3.csv")$kind, "raman_table")
  expect_equal(k("t4.xlsx")$kind, "ldir_table")
  expect_equal(k("t5.csv")$kind, "ldir_coord_swaps")
  expect_equal(k("i3.png")$kind, "ldir_analyzed_image")
  expect_equal(k("i2.png")$kind, "image")                  # the raw LDIR scan is not a palette image
  expect_equal(k("notes.md")$kind, "unknown")
  f <- .faanui_headers()
  kf <- function(p) classify_input_file(list.files(f, p, full.names = TRUE))
  expect_equal(kf("lumos\\.bmp$")$role, "ftir_bruker")
  expect_equal(kf("IAEA 260928\\.bmp$")$role, "raman")
  expect_equal(kf("\\.sam$")$kind, "ldir_sam")
  expect_equal(kf("\\.txt$")$kind, "raman_placement")
  expect_equal(kf("perkin\\.xlsx$")$kind, "ftir_table")
})

test_that("dataset A under neutral names: content decides, the rest is asked", {
  .load_classifier()
  d <- .dataset_a_neutral()
  a <- .asker(list("FTIR table t1.csv" = "FTIR (PerkinElmer)",
                   "image is i1.png" = "FTIR (PerkinElmer)", "image is i4.png" = "Raman"))
  r <- assign_input_roles(list.files(d, full.names = TRUE), a$ask)
  p <- lapply(r$paths, basename)
  expect_equal(p$ftir_file, "t1.csv"); expect_equal(p$ftir_bruker_file, "t2.csv")
  expect_equal(p$raman_file, "t3.csv"); expect_equal(p$ldir_file, "t4.xlsx")
  expect_equal(p$ldir_image, "i2.png"); expect_equal(p$ftir_image, "i1.png")
  expect_equal(p$raman_image, "i4.png")
  # The analyzed image is paired with its raw image by size, so its path is passed on
  expect_equal(basename(r$config$ldir_processed_image_path), "i3.png")
  expect_equal(basename(r$config$ldir_coord_swaps_path), "t5.csv")
  # Asked once for the FTIR tables (the second follows) and once per plain PNG
  expect_length(a$asked(), 3)
  expect_true(any(grepl("not used .*notes.md", r$lines)))
})

test_that("the Faanui folder: one question, which Raman table", {
  .load_classifier()
  f <- .faanui_headers()
  a <- .asker(list("Several Raman tables" =
                     "Comparstic2026 Bora2025 240314 Faanui SW1L R4 A IAEA 260924.csv"))
  r <- assign_input_roles(list.files(f, full.names = TRUE), a$ask)
  p <- lapply(r$paths, basename)
  expect_equal(a$asked(), "Several Raman tables: which one is this run's?")
  expect_match(p$ftir_bruker_file, "lumos\\.csv$")       # paired with the 32-bit Lumos image
  expect_match(p$ftir_bruker_image, "lumos\\.bmp$")
  expect_match(p$ftir_file, "perkin\\.xlsx$")
  expect_match(p$ftir_image, "perkin\\.bmp$")
  expect_match(p$raman_file, "IAEA 260924\\.csv$")
  expect_match(p$raman_image, "260928\\.bmp$")           # a different date than its .txt
  expect_match(p$ldir_file, "Au 260721\\.xlsx$")
  expect_match(p$ldir_image, "Au 260721\\.png$")
  expect_match(basename(r$config$ldir_sam_path), "\\.sam$")
  expect_match(basename(r$config$raman_image_metadata), "\\.txt$")
  # The raw LDIR image has the analyzed image's name: found by the pipeline itself
  expect_null(r$config$ldir_processed_image_path)
  expect_true(any(grepl("^not used .*STJapan", r$lines)))
})

test_that("two candidates for a side file are asked; a side file without its table is not used", {
  .load_classifier()
  f <- .faanui_headers()
  file.copy(list.files(f, "\\.sam$", full.names = TRUE), file.path(f, "second.sam"))
  a <- .asker(list("Several Raman" = "x", "Several LDIR project files" = "second.sam"))
  r <- assign_input_roles(list.files(f, full.names = TRUE), a$ask)
  expect_equal(basename(r$config$ldir_sam_path), "second.sam")
  # Without an LDIR table, the .sam and LDIR images are listed as not used
  g <- .faanui_headers()
  unlink(list.files(g, "Au 260721\\.xlsx$", full.names = TRUE))
  r <- assign_input_roles(list.files(g, full.names = TRUE), .asker(list("Several Raman" = "x"))$ask)
  expect_null(r$config$ldir_sam_path); expect_null(r$paths$ldir_image)
  expect_true(any(grepl("no LDIR table", r$lines)))
})

test_that("a segmented PerkinElmer scan is one FTIR (PerkinElmer) input", {
  .load_classifier()
  d <- withr::local_tempdir()
  a <- file.path(REPO_ROOT, "examples", "dataset_A")
  for (k in 1:3) {
    file.copy(file.path(a, "FTIR_perkin.csv"), file.path(d, sprintf("scan_M%d.csv", k)))
    .write_bmp_header(file.path(d, sprintf("scan_M%d.bmp", k)), 500, 500, 24, 0)
  }
  file.copy(file.path(a, "FTIR_bruker.csv"), file.path(d, "other.csv"))
  file.copy(file.path(a, "Raman_data.csv"), file.path(d, "r.csv"))
  a2 <- .asker()
  r <- assign_input_roles(list.files(d, full.names = TRUE), a2$ask)
  expect_equal(basename(r$paths$ftir_file), "scan_M1.csv")
  expect_equal(basename(r$paths$ftir_image), "scan_M1.bmp")
  # The other FTIR table is then the Bruker one, without a question
  expect_equal(basename(r$paths$ftir_bruker_file), "other.csv")
  expect_length(a2$asked(), 0)
})
