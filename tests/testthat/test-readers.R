# Tests for the reader entry point (R/01_readers.R): the generic
# column-mapping reader, the per-instrument material tables, read_instrument()
# for every role, and the rule that vendor column names stay in the readers.

local({
  for (m in c("01_ingest.R", "01c_ingest_ldir.R", "08b_material_map.R", "01_readers.R"))
    sys.source(file.path(REPO_ROOT, "R", m), envir = globalenv())
})

.write_csv <- function(df, sep = ",") {
  p <- tempfile(fileext = ".csv")
  utils::write.table(df, p, sep = sep, row.names = FALSE, quote = TRUE)
  p
}

# A lab's own export: millimetres, image-style y (down), decimal comma, its
# own column names.
.lab_table <- function() data.frame(
  `Part no` = c("A1", "A2", "A3"),
  `X [mm]`  = c("1,5", "2,0", "0,25"),
  `Y [mm]`  = c("0,5", "1,0", "3,0"),
  `Len`     = c("0,10", "0,05", "0,2"),
  `Wid`     = c("0,05", "0,02", "0,1"),
  `Area`    = c("0,004", "0,001", "0,016"),
  `Polymer` = c("PE", "Polypropylene", "my lab name"),
  `Score`   = c("80", "65", "91"),
  check.names = FALSE, stringsAsFactors = FALSE)

.lab_map <- list(
  columns = list(particle_id = "Part no", x_um = "X [mm]", y_um = "Y [mm]",
                 feret_max_um = "Len", feret_min_um = "Wid", area_um2 = "Area",
                 material = "Polymer", quality = "Score"),
  coordinate_unit_um = 1000, size_unit_um = 1000, y_axis = "down",
  quality_metric = "Lab score")

test_that("generic reader applies units, y orientation and decimal commas", {
  p <- .write_csv(.lab_table(), sep = ";")
  df <- read_mapped_particles(p, .lab_map)
  expect_equal(df$particle_id, c("A1", "A2", "A3"))
  expect_equal(df$x_um, c(1500, 2000, 250))
  expect_equal(df$y_um, c(-500, -1000, -3000))           # y down -> negated
  expect_equal(df$feret_max_um, c(100, 50, 200))
  expect_equal(df$major_um, df$feret_max_um)             # filled from feret_max
  expect_equal(df$minor_um, df$feret_min_um)
  expect_equal(df$area_um2, c(4000, 1000, 16000))        # mm^2 -> um^2
  expect_equal(df$quality, c(80, 65, 91))
  expect_equal(unique(df$source_file), basename(p))
})

test_that("generic reader reads the mapping from JSON and Excel files", {
  skip_if_not_installed("writexl")
  x <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Particles = .lab_table()), x)
  j <- tempfile(fileext = ".json")
  jsonlite::write_json(c(.lab_map, sheet = "Particles"), j, auto_unbox = TRUE)
  df <- read_mapped_particles(x, j)
  expect_equal(df$x_um, c(1500, 2000, 250))
})

test_that("generic reader: clear errors for a bad mapping", {
  p <- .write_csv(.lab_table())
  m <- .lab_map; m$columns$x_um <- "X (mm)"
  expect_error(read_mapped_particles(p, m), "column\\(s\\) not in .*: X \\(mm\\)")
  m <- .lab_map; m$columns$colour <- "Polymer"
  expect_error(load_column_map(m), "unknown field\\(s\\) colour")
  m <- .lab_map; m$columns$material <- NULL
  expect_error(load_column_map(m), "'material' is required")
  m <- .lab_map; m$y_axis <- "left"
  expect_error(read_mapped_particles(p, m), "y_axis must be")
  expect_error(load_column_map(list(sheet = 1)), "needs a 'columns' entry")
  dup <- .lab_table(); dup$`Part no` <- "A1"
  expect_error(read_mapped_particles(.write_csv(dup), .lab_map), "duplicated particle ID")
})

test_that("generic reader numbers rows when there is no ID column", {
  m <- .lab_map; m$columns$particle_id <- NULL
  df <- read_mapped_particles(.write_csv(.lab_table()), m)
  expect_equal(df$particle_id, c("1", "2", "3"))
})

test_that("material tables translate labels, keep the raw label, and reject duplicates", {
  tab <- data.frame(label = c("my lab name", "Polypropylene"),
                    material = c("Polystyrene", "PP"))
  expect_equal(translate_materials(c("My Lab Name ", "PE", NA), tab),
               c("Polystyrene", "PE", NA))
  expect_error(load_material_table("raman", list(raman_material_table =
    data.frame(label = c("a", "A"), material = c("x", "y")))), "lists a label twice")
  expect_error(load_material_table("raman", list(raman_material_table =
    data.frame(name = "a", material = "x"))), "'label' and 'material'")
  # The shipped tables exist for every role and are empty
  for (role in names(instrument_roles)) {
    t <- load_material_table(role, list())
    expect_named(t, c("label", "material"))
    expect_equal(nrow(t), 0)
  }
  # A user CSV with comment lines
  f <- tempfile(fileext = ".csv")
  writeLines(c("# my table", "label,material", "Hostalen GM 7040,Polyethylene"), f)
  expect_equal(load_material_table("ldir", list(ldir_material_table = f))$material,
               "Polyethylene")
})

test_that("read_instrument with a column map returns a validated canonical table", {
  p <- .write_csv(.lab_table())
  tab <- data.frame(label = "my lab name", material = "Polystyrene")
  cfg <- list(raman_path = p, raman_column_map = .lab_map, raman_material_table = tab,
              sample_id = "filter 7")
  df <- read_instrument("raman", cfg)
  expect_true(all(canonical_particle_columns %in% names(df)))
  expect_equal(unique(df$instrument), "raman")
  expect_equal(unique(df$quality_metric), "Lab score")
  expect_equal(unique(df$sample_id), "filter 7")
  expect_equal(df$identification_raw, c("PE", "Polypropylene", "my lab name"))
  expect_equal(df$material, c("PE", "Polypropylene", "Polystyrene"))
  expect_equal(df$material_family, c("PE", "PP", "PS"))
})

test_that("read_instrument: LDIR via a column map keeps coordinates and relabels low quality", {
  t <- .lab_table(); t$Score <- c("0,9", "0,5", "0,95")
  p <- .write_csv(t)
  m <- .lab_map; m$quality_metric <- NULL
  df <- read_instrument("ldir", list(ldir_path = p, ldir_column_map = m,
                                     ldir_hqi_unknown_threshold = 0.85))
  expect_equal(df$x_um, c(1500, 2000, 250))
  expect_equal(df$material, c("PE", "unknown", "my lab name"))
  expect_equal(df$identification_raw[2], "Polypropylene")
  expect_equal(unique(df$quality_metric), "LDIR quality")   # role default
})

test_that("read_instrument: built-in readers on the example dataset", {
  D <- file.path(REPO_ROOT, "examples", "dataset_A")
  skip_if_not(dir.exists(D))
  cfg <- make_config()
  cfg$ftir_path <- file.path(D, "FTIR_perkin.csv")
  cfg$ftir_bruker_path <- file.path(D, "FTIR_bruker.csv")
  cfg$raman_path <- file.path(D, "Raman_data.csv")
  cfg$ldir_path <- file.path(D, "LDIR_data.xlsx")
  n <- c(ftir = 33, ftir_bruker = 42, raman = 43, ldir = 40)
  for (role in names(n)) {
    df <- read_instrument(role, cfg)
    expect_equal(nrow(df), n[[role]], info = role)
    expect_equal(unique(df$instrument), role)
    expect_true(all(c("material_polymer", "material_family") %in% names(df)))
  }
  expect_true(all(is.na(read_instrument("ldir", cfg)$x_um)))
  expect_error(read_instrument("ftir", list()), "no input file configured")
  expect_error(read_instrument("sem", cfg), "Unknown instrument role")
})

test_that("vendor column names appear only in the instrument readers", {
  # Column headers of the built-in export formats. Downstream code must read
  # the canonical columns, never these.
  vendor <- c("Max AAU score", "Area on map", "Major dim", "Minor dim",
              "Feret min", "Feret Min", "hqi_value", "Is Valid", "Identification",
              "Aspect Ratio", "Eccentricity", "Circularity", "Solidity", "Perimeter")
  readers <- c("R/01_ingest.R", "R/01c_ingest_ldir.R")
  files <- c(list.files(file.path(REPO_ROOT, "R"), "\\.R$", full.names = TRUE),
             file.path(REPO_ROOT, "main.R"),
             list.files(file.path(REPO_ROOT, "shiny_app"), "\\.R$", full.names = TRUE))
  files <- files[!sub(paste0(REPO_ROOT, "/"), "", files, fixed = TRUE) %in% readers]
  hits <- character(0)
  for (f in files) {
    code <- sub("#.*$", "", readLines(f, warn = FALSE))     # ignore comments
    for (v in vendor) {
      i <- grep(paste0("[\"'`]", v, "[\"'`]"), code)
      if (length(i)) hits <- c(hits, sprintf("%s:%d %s", basename(f), i, v))
    }
    # The fuzzy column lookups are reader tools
    i <- grep("\\b(find_column|safe_numeric|safe_col_numeric)\\(", code)
    if (length(i)) hits <- c(hits, sprintf("%s:%d column lookup", basename(f), i))
  }
  expect_identical(hits, character(0))
})

test_that("an Excel export is read from its only sheet, whatever its name", {
  skip_if_not_installed("writexl")
  d <- data.frame(Identifier = c("MP_0", "MP_1"), `Coord. [um]` = c("[100;200]", "[300;400]"),
                  `Max AAU score` = c(0.9, 0.8), Group = c("PET", "PP"),
                  `Area on map [um2]` = c(400, 900), `Major dim [um]` = c(30, 40),
                  `Minor dim [um]` = c(20, 25), check.names = FALSE)
  f <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Data = d), f)                 # not "Long_Table"
  x <- suppressMessages(read_data_file(f, sheet = "Long_Table"))
  expect_equal(nrow(x), 2)
  writexl::write_xlsx(list(Data = d, Other = d), f)      # ambiguous: say which sheets exist
  expect_error(suppressMessages(read_data_file(f, sheet = "Long_Table")), "'Data', 'Other'")
})
