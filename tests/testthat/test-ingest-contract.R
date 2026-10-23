# Reader input contract (R/01_ingest.R): an export that does not meet it must
# fail with a message naming the problem -- never be read into wrong numbers.

sys.source(file.path(REPO_ROOT, "R", "01_ingest.R"), envir = globalenv())

.ftir_header <- paste0("Identifier,Coord. [pixels],Coord. [\u00b5m],Max AAU score,",
                       "Group,Area on map [\u00b5m\u00b2],Major dim [\u00b5m],",
                       "Minor dim [\u00b5m],Feret min [\u00b5m]")
.ftir_rows <- c("MP_0,[15;264],[94;1650],0.944,Polycarbonate,15898.4,184.1,110,125.6",
                "MP_1,[60;678],[375;4238],0.821,Polycarbonate,5585.9,124.2,57.3,83.7")

.write_csv <- function(lines) {
  p <- tempfile(fileext = ".csv")
  writeLines(enc2utf8(lines), p, useBytes = TRUE)
  p
}

test_that("a well-formed FTIR export reads in micrometres", {
  d <- suppressMessages(ingest_ftir(.write_csv(c(.ftir_header, .ftir_rows))))
  expect_equal(d$x_um, c(94, 375))
  expect_equal(d$y_um, c(1650, 4238))
  expect_equal(d$quality, c(0.944, 0.821))
})

test_that("the pixel coordinate column is never used as micrometres", {
  hdr <- sub("Coord. \\[\u00b5m\\]", "Position", .ftir_header)
  expect_error(suppressMessages(ingest_ftir(.write_csv(c(hdr, .ftir_rows)))),
               "no micrometre coordinate column")
})

test_that("an export with no particle rows fails clearly", {
  expect_error(suppressMessages(ingest_ftir(.write_csv(.ftir_header))),
               "contains no particle rows")
})

test_that("repeated FTIR identifiers (MP_-1 for every non-plastic) are renumbered, the export kept", {
  # The Spotlight export labels every non-plastic MP_-1
  rows <- c(.ftir_rows, sub("^\"?MP_[0-9]+\"?", "MP_-1", .ftir_rows[2]),
            sub("^\"?MP_[0-9]+\"?", "MP_-1", .ftir_rows[2]))
  d <- suppressMessages(ingest_ftir(.write_csv(c(.ftir_header, rows))))
  expect_equal(d$particle_id, paste0("P_", seq_along(rows)))
  expect_equal(d$identifier_export[length(rows)], "MP_-1")
  expect_false(anyDuplicated(d$particle_id) > 0)
  # Unique identifiers are left alone
  d1 <- suppressMessages(ingest_ftir(.write_csv(c(.ftir_header, .ftir_rows))))
  expect_null(d1$identifier_export)
  expect_false(any(startsWith(d1$particle_id, "P_")))
})

test_that("a stitched table keeps its P_ numbers and its exported identifiers", {
  hdr <- sub("^Identifier,", "Identifier,Identifier (export),", .ftir_header)
  d <- suppressMessages(ingest_ftir(.write_csv(c(hdr, paste0("P_", 1:2, ",MP_-1,", sub("^[^,]*,", "", .ftir_rows))))))
  expect_equal(d$particle_id, c("P_1", "P_2"))
  expect_equal(d$identifier_export, c("MP_-1", "MP_-1"))
  expect_equal(d$x_um, c(94, 375))
})

test_that("a semicolon export with decimal commas reads the same numbers", {
  eu <- c(gsub(",", ";", .ftir_header),
          '"MP_0";"[15;264]";"[94,5;1650,25]";"0,944";"Polycarbonate";"15898,4";"184,1";"110";"125,6"')
  d <- suppressMessages(ingest_ftir(.write_csv(eu)))
  expect_equal(c(d$x_um, d$y_um, d$quality, d$area_um2), c(94.5, 1650.25, 0.944, 15898.4))
})

test_that("as_number reads decimal commas and leaves points alone", {
  expect_equal(as_number(c("1,5", "2.5", "-0,25", "3", NA, "abc")),
               c(1.5, 2.5, -0.25, 3, NA, NA))
})

test_that("a Raman export with duplicate IDs fails", {
  r <- c("Particle Name,Visual Center Point X [um],Visual Center Point Y [um],Material,HQI",
         "R1,10,20,PET,80", "R1,30,40,PE,90")
  expect_error(suppressMessages(ingest_raman(.write_csv(r))), "duplicated particle ID.*R1")
})
