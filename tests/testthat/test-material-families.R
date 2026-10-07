# Material family classification (R/08b_material_map.R).
#
# The family rules were rebuilt against the full S.T. Japan Raman library
# listings (L600xx). Those listings are licensed and git-ignored, so nothing
# here reads them: the names below are the cases that went wrong before (PEG
# counted as PE, "POLY(ETHYLENE TEREPHTHALATE)" unrecognised, nylon matched by
# "STEARAMIDE", ...) plus the instrument names the pipeline already relied on,
# and the library index is exercised with a small synthetic one.

sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = globalenv())

fam <- function(x, index = NULL) classify_family_vec(x, index = index)

test_that("instrument names used so far keep their family", {
  expect_identical(
    fam(c("Polyethylene terephtalate (PET)", "PET", "Polypro",
          "Polypropylene (PP) Lab bench", "Polypropylene", "Polycarbonate (PC)",
          "Nylon 12", "HDPE", "LDPE", "Polystyrene", "PVC", "Cellulose acetate",
          "Polyamide (naturally occurring)", "Polyurethane", "PTFE")),
    c("PET", "PET", "PP", "PP", "PP", "PC", "PA", "PE", "PE", "PS", "PVC",
      "Cellulose", "Natural", "PU", "PTFE"))
})

test_that("upper/lower case and bracket spelling do not change the family", {
  expect_identical(
    fam(c("polyethylene", "POLYETHYLENE", "PolyEthylene", "Poly(ethylene)",
          "POLY(ETHYLENE TEREPHTHALATE)", "poly(ethylene terephthalate)",
          "POLY(VINYL CHLORIDE)", "Poly(methyl methacrylate) 120000")),
    c("PE", "PE", "PE", "PE", "PET", "PET", "PVC", "PMMA"))
})

test_that("polymer-sounding non-plastics are not counted as plastics", {
  expect_identical(
    fam(c("POLYETHYLENE GLYCOL 4000", "PEG | MACROGOL 6000 | POLYETHYLENE GLYCOL 6000",
          "BRIJ (R) 78 | POLYETHYLENE GLYCOL OCTADECYL ETHER",
          "POLY(PROPYLENE GLYCOL) 500", "STEARAMIDE | OCTADECANAMIDE",
          "DIOCTYL PHTHALATE | BIS(2-ETHYLHEXYL) PHTHALATE", "IRGANOX 1010",
          "ANTIPYRINE | PHENAZONE | PHENYLONE")),
    c("Additive", "Additive", "Additive", "Additive", "Additive", "Additive",
      "Additive", "Unknown"))
  # polyesters / polyethers that merely start with "poly(ethylene|propylene ..."
  expect_identical(
    fam(c("POLYETHYLENEIMINE", "POLY(PROPYLENE ADIPATE)", "POLY(ETHYLENE SUCCINATE)")),
    rep("Other polymer", 3))
  # acrylic monomers are reagents, not acrylic plastic
  expect_identical(fam(c("BUTYL ACRYLATE", "ACRYLIC ACID")), c("Unknown", "Unknown"))
  # natural-product names that contain a pattern by accident
  expect_identical(fam(c("L-TYROSINE", "SPECTINOMYCIN DIHYDROCHLORIDE")),
                   c("Unknown", "Unknown"))
})

test_that("copolymers, blends and masterbatches go to the right family", {
  expect_identical(
    fam(c("Ethylene/Vinyl Acetate Copolymer 80:20",
          "COLOR MASTERBATCH POLY[ETHYLENE-CO-(VINYL ACETATE)] + 45% YELLOW PIGMENT #1",
          "ELVAX 265, COPOLYMER EVA TYPE",
          "COLOR MASTERBATCH POLYETHYLENE + 30% RED PIGMENT #1",
          "POLYPROPYLENE WITH 20% TALC",
          "STYRENE/ACRYLONITRILE COPOLYMER", "POLY(ACRYLONITRILE:BUTADIENE:STYRENE) #1",
          "POLY(STYRENE:BUTADIENE)", "POLY(STYRENE-ETHYLENE-BUTYLENE)",
          "KELTAN 512, COPOLYMER EPDM TYPE", "POLYSTYRENE HIGH IMPACT",
          "NYLON 6,T | POLYTRIMETHYL HEXAMETHYLENE TEREPHTHALAMIDE",
          "Poly(1,4-butylene Terephthalate)", "POLYACETAL",
          "EPOXY RESIN, BISPHENOL A + EPICHLOROHYDRIN", "POWDER COATING H50 | EPOXYPOLYESTER")),
    c("EVA", "EVA", "EVA", "PE", "PP", "ABS", "ABS", "Rubber", "Rubber", "Rubber",
      "PS", "PA", "PET", "POM", "Epoxy", "Epoxy"))
})

test_that("minerals, fillers and pigments get their own coarse families", {
  expect_identical(
    fam(c("Carbonate", "Sand", "CALCIUM CARBONATE", "TALC", "KAOLIN", "BARIUM SULFATE",
          "TITANIUM(IV) OXIDE ANATASE", "PIGMENT WHITE 6 BASED | TITANIUM DIOXIDE",
          "IRON OXIDE", "PIGMENT BLUE 15")),
    c(rep("Mineral", 6), rep("Pigment", 4)))
})

test_that("every family name classifies back to itself", {
  all_fams <- c(names(default_polymer_families), "Inorganic", "Other compound")
  expect_identical(fam(all_fams), all_fams)
  expect_identical(fam(toupper(all_fams)), all_fams)
})

test_that("categories cover every family and keep plastics apart", {
  cats <- classify_category_vec(c(names(default_polymer_families), "Inorganic",
                                  "Other compound", "Unknown"))
  expect_true(all(cats %in% material_category_levels))
  expect_false(any(cats[-length(cats)] == "Unknown"))
  expect_identical(classify_category_vec(c("EVA", "Other polymer", "Additive",
                                           "Pigment", "Mineral", "Inorganic",
                                           "Other compound")),
                   c("Synthetic", "Synthetic", "Additive/Pigment", "Additive/Pigment",
                     "Inorganic", "Inorganic", "Other"))
})

# --- Spectral-library index --------------------------------------------------

.toy_index <- function() {
  .prepare_library_index(data.frame(
    library  = c("L60002", "L60002", "L60019", "L60020", "L60020"),
    entry_id = 1:5,
    entry    = c("HOSTALEN GM 7040 | POLYETHYLENE",
                 "LUTENSOL XYZ | POLYETHYLENE GLYCOL ETHER",
                 "SODIUM CHLORIDE | HALITE",
                 "ACETAMINOPHEN | PARACETAMOL",
                 "KOHJIN VCI | CYANOMETHANE"),
    stringsAsFactors = FALSE))
}

test_that("trade names resolve through their library synonyms", {
  idx <- .toy_index()
  expect_identical(fam("Hostalen GM 7040"), "Unknown")          # no index
  expect_identical(fam(c("Hostalen GM 7040", "hostalen gm-7040", "HOSTALEN GM 7040 | POLYETHYLENE"),
                       index = idx), c("PE", "PE", "PE"))
  expect_identical(fam("Lutensol XYZ", index = idx), "Additive")
})

test_that("recognised non-plastic entries are labelled, unknown names stay Unknown", {
  idx <- .toy_index()
  expect_identical(fam(c("Halite", "Paracetamol", "something else"), index = idx),
                   c("Inorganic", "Other compound", "Unknown"))
  # A complete entry line is matched as a whole, not via a synonym it shares
  # with an unrelated entry.
  expect_identical(fam("KOHJIN VCI | CYANOMETHANE", index = idx), "Other compound")
})

test_that("the index can be switched off and is off for the test suite", {
  expect_null(spectral_library_dir())
  expect_null(get_spectral_library_index())
  withr_dir <- tempfile("no_such_dir")
  old <- options(ftir.spectral_library_dir = withr_dir)
  on.exit(options(old), add = TRUE)
  expect_null(spectral_library_dir())
})

test_that("a library listing PDF is parsed, wrapped entries re-joined", {
  skip_if_not_installed("pdftools")
  dir <- tempfile("spectral_libraries")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  pdf_path <- file.path(dir, "L69999 Test Raman Spectra Database.pdf")

  lines <- c(
    "S.T.Japan-Europe GmbH          Test Database          Raman Spectra",
    "L69999",
    "ACRYLIC ACID",
    paste0("BRIJ (R) 78 | POLYETHYLENE GLYCOL OCTADECYL ETHER | POLYOXYETHYLENE (20) ",
           "MONO-"),
    "1-OCTADECYL ETHER",   # wrapped tail: out of alphabetical order
    "HOSTALEN GM 7040 | POLYETHYLENE",
    "TALC",
    "contact@stjapan.de          www.stjapan.de          1 of 1",
    "All information within this document is the property of S.T.Japan-Europe GmbH"
  )
  grDevices::pdf(pdf_path, width = 11, height = 8.5, family = "Courier")
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, length(lines) + 1))
  graphics::text(0.01, rev(seq_along(lines)), lines, adj = 0, cex = 0.8)
  grDevices::dev.off()

  parsed <- parse_spectral_library_pdf(pdf_path)
  expect_identical(unique(parsed$library), "L69999")
  expect_identical(parsed$entry, c(
    "ACRYLIC ACID",
    paste0("BRIJ (R) 78 | POLYETHYLENE GLYCOL OCTADECYL ETHER | POLYOXYETHYLENE (20) ",
           "MONO-1-OCTADECYL ETHER"),
    "HOSTALEN GM 7040 | POLYETHYLENE",
    "TALC"))

  idx <- load_spectral_library_index(dir)
  expect_true(file.exists(file.path(dir, "library_index.csv")))
  expect_identical(fam(c("Hostalen GM 7040", "Brij (R) 78"), index = idx),
                   c("PE", "Additive"))
})
