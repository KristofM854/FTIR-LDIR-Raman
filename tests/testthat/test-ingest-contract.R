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

test_that("a Lumos export with \"[x,y]\" coordinates is read, with a warning", {
  # Seen in Bruker Lumos exports: ";"-separated, decimal comma, and the
  # coordinates written with a comma between x and y
  hdr <- gsub(",", ";", .ftir_header)
  rows <- c("MP_-1;[47,677];[264,3798];0,752;Protein;440,6;34,2;16,4;22,4",
            "MP_0;[106,779];[595,4370];0,681;PA;1007,1;61,4;20,9;32,0")
  msgs <- character(0)
  d <- withCallingHandlers(ingest_ftir(.write_csv(c(hdr, rows))),
                           message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                                   invokeRestart("muffleMessage") })
  expect_equal(d$x_um, c(264, 595))
  expect_equal(d$y_um, c(3798, 4370))
  expect_equal(d$quality, c(0.752, 0.681))
  expect_true(any(grepl("[x,y]", msgs, fixed = TRUE)))
  # A decimal comma inside "[x;y]" stays a decimal comma
  expect_equal(parse_ftir_coordinates(c("[94,5;1650,2]", "[264,3798]"))$x_um, c(94.5, 264))
})


# AUDIT L-07: "Area" must not resolve to "Convex Area" when that column comes first
test_that("find_column prefers a column that starts with the pattern", {
  df <- data.frame(a = 1, b = 2, c = 3, check.names = FALSE)
  names(df) <- c("Convex Area [um2]", "Area [um2]", "Feret Max [um]")
  expect_equal(find_column(df, "Area"), "Area [um2]")
  expect_equal(find_column(df, "Feret Max"), "Feret Max [um]")
  expect_equal(find_column(df, "Convex"), "Convex Area [um2]")
  names(df)[2] <- "Total surface"
  expect_equal(find_column(df, "Area"), "Convex Area [um2]")   # substring as the fallback
})
