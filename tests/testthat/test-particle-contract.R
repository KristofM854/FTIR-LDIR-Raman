# Tests for the internal data contract (R/particle_contract.R): the canonical
# particle table, the canonical image object and their validators.

.contract_df <- function(n = 3) {
  as_canonical_particles(data.frame(
    particle_id = paste0("P", seq_len(n)),
    x_um = seq_len(n) * 10, y_um = seq_len(n) * 5,
    feret_max_um = 20, feret_min_um = 10, major_um = 20, minor_um = 10,
    area_um2 = 150, material = "PE", quality = 0.9, source_file = "a.csv",
    stringsAsFactors = FALSE), "ftir", "AAU")
}

test_that("as_canonical_particles adds the contract columns and keeps values", {
  raw <- data.frame(particle_id = c("a", "b"), x_um = 1:2, y_um = 3:4,
                    material = "PP", extra = 7, stringsAsFactors = FALSE)
  df <- as_canonical_particles(raw, "raman", "HQI", sample_id = "S1")
  expect_true(all(canonical_particle_columns %in% names(df)))
  expect_equal(names(df)[seq_along(raw)], names(raw))   # reader columns stay first
  expect_equal(df$x_um, 1:2)
  expect_equal(df$extra, c(7, 7))
  expect_equal(unique(df$instrument), "raman")
  expect_equal(unique(df$quality_metric), "HQI")
  expect_equal(unique(df$sample_id), "S1")
  expect_true(all(is.na(df$area_um2)))
  # Columns the reader already set are not overwritten
  raw$instrument <- "custom"
  expect_equal(unique(as_canonical_particles(raw, "raman")$instrument), "custom")
})

test_that("validate_particles accepts a valid table", {
  df <- .contract_df()
  expect_identical(validate_particles(df), df)
})

test_that("validate_particles lists every problem at once", {
  df <- .contract_df()
  df$particle_id[2] <- "P1"
  df$feret_max_um[1] <- -1
  df$quality <- as.character(df$quality)
  df$area_um2 <- NULL
  err <- tryCatch(validate_particles(df, what = "Test table"), error = conditionMessage)
  expect_match(err, "Test table does not meet the particle contract")
  expect_match(err, "missing column\\(s\\): area_um2")
  expect_match(err, "duplicated particle_id: P1")
  expect_match(err, "'quality' is not numeric")
  expect_match(err, "negative values in 'feret_max_um'")
})

test_that("validate_particles: coordinates required unless require_coords = FALSE", {
  df <- .contract_df()
  df$x_um <- NA_real_; df$y_um <- NA_real_
  expect_error(validate_particles(df), "no particle has coordinates")
  expect_silent(validate_particles(df, require_coords = FALSE))
  df$x_um[1] <- Inf
  expect_error(validate_particles(df, require_coords = FALSE), "infinite coordinate")
  df <- .contract_df(); df$particle_id[1] <- ""
  expect_error(validate_particles(df), "empty particle_id")
  expect_error(validate_particles(.contract_df()[0, ]), "has no rows")
})

test_that("validate_particles reports partly missing coordinates without failing", {
  df <- .contract_df(); df$x_um[1] <- NA
  expect_message(validate_particles(df), "1 of 3 particles have no coordinates")
})

test_that("image object: pixel <-> um mapping, y up", {
  img <- new_image_object("ldir", "x.png", 100, 80, pixel_size_um = 2,
                          origin_px = c(50, 40), provenance = "scan_circle")
  expect_equal(unname(img$extent_um), c(-100, 100, -80, 80))
  um <- image_pixels_to_um(c(50, 60, 50), c(40, 40, 30), img)
  expect_equal(um$x_um, c(0, 20, 0))
  expect_equal(um$y_um, c(0, 0, 20))          # a row higher up = larger y
  px <- image_um_to_pixels(um$x_um, um$y_um, img)
  expect_equal(px$col_px, c(50, 60, 50))
  expect_equal(px$row_px, c(40, 40, 30))
})

test_that("image object from an extent reproduces the extent, both y orientations", {
  ext <- c(xmin = -10, xmax = 190, ymin = 5, ymax = 105)
  up <- image_object_from_extent("ftir", "f.png", 400, 200, ext, "config")
  expect_equal(unname(up$extent_um), unname(ext))
  expect_equal(up$pixel_size_um, c(0.5, 0.5))
  # Pixel (0, 0) is the top-left corner: (xmin, ymax)
  expect_equal(unlist(image_pixels_to_um(0, 0, up)), c(x_um = -10, y_um = 105))
  down <- image_object_from_extent("ftir", "f.png", 400, 200, ext, "registration",
                                   y_axis = "down")
  expect_equal(unname(down$extent_um), unname(ext))
  expect_equal(unlist(image_pixels_to_um(0, 0, down)), c(x_um = -10, y_um = 5))
  expect_equal(unlist(image_pixels_to_um(400, 200, down)), c(x_um = 190, y_um = 105))
})

test_that("validate_image rejects broken geometry and unknown provenance", {
  img <- new_image_object("raman", "r.png", 10, 10, 1, c(0, 10), "config")
  bad <- img
  bad$pixel_size_um <- c(1, -1)
  bad$provenance <- "guess"
  bad$y_axis <- "sideways"
  err <- tryCatch(validate_image(bad), error = conditionMessage)
  expect_match(err, "pixel_size_um must be two positive numbers")
  expect_match(err, "unknown provenance 'guess'")
  expect_match(err, "y_axis must be 'up' or 'down'")
  expect_error(new_image_object("raman", "r.png", 0, 10, 1, c(0, 0), "config"),
               "width_px must be > 0")
})

test_that("write_image_objects writes one record per image", {
  d <- tempfile(); dir.create(d)
  img <- new_image_object("ldir", "l.png", 10, 10, 2, c(5, 5), "excel_area_scale")
  p <- write_image_objects(d, list(ldir = img, ftir = NULL))
  rec <- jsonlite::fromJSON(p)
  expect_named(rec, "ldir")
  expect_equal(rec$ldir$provenance, "excel_area_scale")
  expect_equal(rec$ldir$extent_um, c(-10, 10, -10, 10))
})
