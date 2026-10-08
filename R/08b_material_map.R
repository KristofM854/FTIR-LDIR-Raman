# =============================================================================
# 08b_material_map.R — Cross-instrument material equivalence mapping
# =============================================================================
#
# FTIR, Raman, and LDIR identify materials using different spectral libraries
# and nomenclature. This module provides:
#   1. A configurable regex mapping from raw names → canonical types
#   2. A material family classification (PET → PET, Nylon 6 → PA, talc →
#      Mineral, PEG → Additive, ...)
#   3. A tiered agreement scorer: Exact / Family / Filler/pigment / Disagree
#   4. A category classifier: Synthetic / Semi-synthetic / Natural/Organic /
#      Additive/Pigment / Inorganic / Other / Unknown
#   5. An optional lookup in the local spectral-library index (the S.T. Japan
#      Raman library listings shipped with the WITec instrument), so trade
#      names and synonyms resolve to the same family as their generic name
#      ("HOSTALEN GM 7040" → its synonym "POLYETHYLENE" → PE).
#
# The mapping tables are stored in config (material_map_ftir, material_map_raman,
# material_map_ldir). When NULL, the pipeline uses built-in defaults below.
#
# Spectral-library index
# ----------------------
# The library listings are licensed and must not be committed: they live in the
# git-ignored `spectral_libraries/` folder at the repository root (one PDF per
# library, e.g. "L60035 Microplastics and Related Compounds.pdf"). On first use
# the PDFs are parsed (needs the `pdftools` package) into
# `spectral_libraries/library_index.csv`, which is reused until a PDF changes.
# Without the folder everything still works: classification then relies on the
# name rules alone, exactly as before the index existed.
# =============================================================================


# ---- Default material family definitions ----
# Each family maps a canonical family name to a character vector of patterns
# (case-insensitive, Perl regex). Patterns are tested against a flattened form
# of the name (upper case, brackets removed, "POLY(" joined to the monomer, so
# "Poly(ethylene terephthalate)" is seen as "POLYETHYLENE TEREPHTHALATE").
# Families are tried IN ORDER and the first match wins — so specific families
# sit above the generic ones they would otherwise fall into (PET above PE, the
# plastics above "Additive", everything above the "Other polymer" catch-all).

# Names that start like a commodity polymer but are something else: PEG and
# other polyether surfactants, polyethyleneimine, and the "POLY(ETHYLENE
# ADIPATE)"-type polyesters. Used as a negative lookahead after PE / PP.
.poly_not_polyolefin <- paste0(
  "(?!.*(GLYCOL|OXIDE\\b|IMINE|ETHER|ETHOXYL|TEREPH|NAPHTH|CHLORINAT|CHLOROSUL|",
  "ADIPATE|AZELATE|FUMARATE|MALEATE|PHTHALATE|SUCCINATE|SEBACATE|GLUTARATE|",
  "TETRAHYDRO|TERTAHYDRO|HEXACHLOR|ENDOMETHYLENE|DIMETHACRYLATE|CARBONATE|SULFIDE))"
)

default_polymer_families <- list(
  # --- Synthetic polymers -----------------------------------------------------
  PET = c(
    "^PET$", "^PETG$", "^PBT$",
    "\\b(POLY)?ETHYLENE\\s*TEREPH(?!.*AMIDE)",  # also catches the "terephtalate" typo
    "^POLY.*TEREPHT?HALATE",               # PBT and other terephthalate polyesters
    "(?<!EPOXY)POLYESTER"
  ),
  PP = c(
    "^PP$", "^BOPP$", "^CPP$", "^POLYPRO$",
    paste0("POLYPROPYLENE?+(?!\\s*-?\\s*(GLYCOL|OXIDE\\b|ADIPATE|FUMARATE|MALEATE|",
           "ISOPHTHALATE|PHTHALATE|SUCCINATE|SEBACATE|HEXACHLOR|TETRAHYDRO|",
           "TERTAHYDRO|CARBONATE))")
  ),
  EVA = c(
    "^EVA$",
    "ETHYLENE\\W*(CO\\W*)?VINYL\\s*ACETATE",
    "\\bEVA\\b(?!L)",                      # "COPOLYMER EVA TYPE", not EVAL (EVOH)
    "\\(EVA\\)"
  ),
  PE = c(
    "^PE$", "^(HD|LD|LLD|MD|VLD|UHMW)PE$", "^PE[- ]?(HD|LD|LLD|MD|UHMW)$",
    paste0("POLYETHYLENE(?!\\s*-\\s*CO\\b|\\s*:)", .poly_not_polyolefin),
    "(HIGH|LOW|MEDIUM)[ .-]?DENSITY[ -]?POLYETHYLENE"
  ),
  ABS = c(
    "^ABS$", "^SAN$", "^ASA$", "^AES$",
    "ACRYLONITRILE.*BUTADIENE.*STYRENE",
    "STYRENE\\W*(CO\\W*)?ACRYLONITRILE",
    "ACRYLONITRILE\\W*(CO\\W*)?STYRENE",
    "ACRYLONITRILE\\W*(CO\\W*)?ETHYLENE\\W*(CO\\W*)?PROPYLENE\\W*(CO\\W*)?STYRENE",
    "COPOLYMER (ABS|SAN) TYPE"
  ),
  PS = c(
    "^PS$", "^EPS$", "^XPS$", "^HIPS$",
    "POLYSTYRENE(?!.*(SULFON|\\bCO\\b|:|BLOCK|-B-|BUTADIENE|ISOPRENE|BUTYLENE))",
    "^STYRENE POLYMER$", "STYROFOAM"
  ),
  PVC = c(
    "^U?PVC$", "^VINYON$",
    "POLYVINYL\\s*CHLORIDE",
    "VINYL\\s*CHLORIDE\\s*/\\s*VINYL\\s*ACETATE",
    "VINYL\\s*ACETATE\\s*:\\s*VINYL\\s*CHLORIDE"
  ),
  PA = c(
    "^PA$", "^PA\\s?\\d", "\\bNYLON",
    "^POLYAMIDE(?!.*(NATURAL|OCCURRING|IMIDE))",  # PA but NOT "naturally occurring"
    "POLY\\s*CAPROLACTAM", "POLYLAURYLLACTAM", "POLYUNDECANOAMIDE",
    "POLYHEXAMETHYLENE\\s*(ADIPAMIDE|SEBACAMIDE|DODECANEDIAMIDE|NONANEDIAMIDE|AZELAMIDE)",
    "TEREPHTHALAMIDE", "PHENYLENE\\s*ISOPHT?HALAMIDE", "\\bARAMIDS?\\b"
  ),
  PC = c(
    "^PC$",
    "POLYCARBONATE",
    "BISPHENOL\\s*-?\\s*A\\s*CARBONATE"
  ),
  PMMA = c(
    "^PMMA$", "^ACRYLIC$", "PLEXIGLAS", "PERSPEX",
    "POLY\\s*METHYL\\s*METHACRYLATE(?!.*(:|\\bCO\\b|BUTYL|STYRENE))",
    "^METHYL METHACRYLATE POLYMER$"
  ),
  PU = c(
    "^PUR?$",
    "POLYURETHAN",
    "POLYISOCYANATE", "DIISOCYANATE BASED", "POLYPHENYL ISOCYANATE"
  ),
  PTFE = c(
    "^PTFE$",
    "POLYTETRAFLUOROETHYLENE(?!\\s*-\\s*CO\\b)",
    "TEFLON"
  ),
  POM = c(
    "^POM$", "POLYACETAL", "POLYOXYMETHYLENE", "ACETAL\\s*(BASED\\s*)?(CO)?POLYMER",
    "DELRIN"
  ),
  Epoxy = c(
    "^EPOXY", "EPOXY\\s*RESIN", "EPOXYPOLYESTER", "ARALDITE",
    "BISPHENOL\\s*-?\\s*A\\W*(CO\\W*)?EPICHLOR"
  ),
  Rubber = c(
    "RUBBER", "^SBR$", "^NBR$", "^EPDM$", "^EPM$", "^SBS$", "^SEBS$", "^IIR$",
    "\\bTIRES?\\b", "\\bTYRES?\\b",
    "COPOLYMER (EPDM|SBR) TYPE", "POLYMER PB TYPE",
    "POLYBUTADIENE", "POLYISOPRENE", "POLYCHLOROPRENE", "NEOPRENE", "POLYISOBUTYLENE",
    "ISOBUTYLENE\\W*(CO\\W*)?ISOPRENE",
    "STYRENE\\W*(CO\\W*)?(BUTADIENE|ISOPRENE|ETHYLENE\\W*BUTYLENE)",
    "BUTADIENE\\W*(CO\\W*)?(STYRENE|ACRYLONITRILE)",
    "ACRYLONITRILE\\W*(CO\\W*)?BUTADIENE",
    "ETHYLENE\\W*(CO\\W*)?PROPYLENE(?!.*GLYCOL)",
    "\\bLATEX\\b", "ELASTOMER"
  ),
  # --- Semi-synthetic -----------------------------------------------------------
  Cellulose = c(
    "Cellulos",     # matches Cellulose, Cellulosic, Cellulose Acetate, CAB
    "^CAB$", "CARMELLOSE",
    "Rayon", "Viscose", "Lyocell"
  ),
  Acrylate = c(
    # Acrylic polymers only — the monomers (butyl acrylate, acrylic acid, ...)
    # are reagents, not particles of plastic.
    "^ACRYL(ATE|AMIDE)S?$", "^POLY.*ACRYL", "ACRYL.*(CO|TER)POLYMER",
    "(CO|TER)POLYMER.*ACRYL", "ACRYLATE BASED", "IONOMER",
    "EUDRAGIT", "CARBOPOL", "ACRITAMER"
  ),
  # --- Natural / organic ------------------------------------------------------------
  Protein = c(
    "Protein", "Keratin", "Silk", "^Wool$",
    "COLLAGEN", "\\bGELATINE?\\b", "CASEIN", "ALBUMIN", "FIBROIN"
  ),
  Chitin = c(
    "Chitin", "Chitosan"
  ),
  Natural = c(
    "naturally\\s*occurring",   # covers "Polyamide (naturally occurring)" from LDIR
    "Cotton",
    "^Natural$", "^Organic$", "^Biofilm$",
    "^(?!GLUE).*\\bWOOD\\b", "LIGNIN", "STARCH", "POLYSACCHARID", "\\bPECTINE?\\b", "POLYGALACTURON",
    "CARRAGEENAN", "\\bAGAR\\b", "ALGINATE", "CAROTEN", "CHLOROPHYLL",
    "BEESWAX", "CARNAUBA", "SHELLAC", "COLOPHONY", "\\bROSIN\\b", "\\bAMBER\\b",
    "\\bCOPAL\\b", "\\bDAMMAR\\b", "BALSAM", "PINE RESIN", "GUAJAC", "\\bGUM\\b"
  ),
  # --- Additives and colourants (found in/on plastics, not a polymer) -------------
  Additive = c(
    # polyether surfactants / PEG — the main false "PE" hit in the libraries
    "POLYETHYLENE\\s*-?\\s*GLYCOL", "^PEG\\b", "MACROGOL", "POLYOXYETHYLENE",
    "POLYSORBATE", "\\bTWEEN\\b", "\\bTRITON\\b", "\\bBRIJ\\b", "\\bIGEPAL\\b",
    "POLOXAMER", "POLYPROPYLENE\\s*-?\\s*GLYCOL", "POLYETHYLENE\\s*OXIDE\\b",
    # plasticisers
    "^(?!POLY).*PHTHALATE", "^(?!POLY).*\\b(DI|BIS).*ADIPATE", "TRIMELLITATE",
    "(TRIBUTYL|TRIETHYL|ACETYL TRIBUTYL) CITRATE", "(TRICRESYL|TRIPHENYL) PHOSPHATE",
    "CHLORINATED PARAFFIN", "PLASTICI[SZ]ER",
    # antioxidants, UV and heat stabilisers
    "IRGANOX", "IRGAFOS", "TINUVIN", "CHIMASSORB", "^BHT$",
    "BUTYLATED HYDROXYTOLUENE", "ANTIOXIDANT", "UV ABSORBER",
    # slip agents, lubricants, flame retardants
    "ERUCAMIDE", "OLEAMIDE", "STEARAMIDE", "^(?!POLY).*STEARATE",
    "DECABROMO", "TETRABROMOBISPHENOL", "HEXABROMOCYCLODODECANE", "FLAME RETARDANT"
  ),
  Pigment = c(
    "PIGMENT", "\\bDYES?\\b",
    "^(ACID|BASIC|DIRECT|DISPERSE|REACTIVE|SOLVENT|VAT|MORDANT|FOOD) (RED|BLUE|YELLOW|GREEN|BLACK|BROWN|VIOLET|ORANGE|WHITE)\\b",
    "PHTHALOCYANIN", "QUINACRIDON", "ULTRAMARINE", "PRUSSIAN BLUE", "INDIGO",
    "ALIZARIN", "RHODAMINE",
    "TITANIUM\\s*DIOXIDE", "TITANIUM\\s*\\(?IV\\)?\\s*OXIDE", "^TIO2$", "ANATASE", "RUTILE",
    "CARBON BLACK", "LAMP\\s*BLACK", "IVORY BLACK",
    "IRON\\s*\\(?III\\)?\\s*OXIDE", "IRON OXIDE", "FERRIC OXIDE", "HA?EMATITE",
    "MAGNETITE", "GOETHITE", "OCHRE",
    "ZINC OXIDE", "ZINC WHITE", "LITHOPONE",
    "CHROME (YELLOW|GREEN|OXIDE)", "CADMIUM (RED|YELLOW|SULFIDE|SELENIDE)"
  ),
  # --- Minerals and fillers ---------------------------------------------------------
  Mineral = c(
    "^Carbonate$", "^Sand$",
    "CALCIUM CARBONATE", "MAGNESIUM CARBONATE", "BASIC.*CARBONATE", "CARBONATE,? BASIC",
    "CALCITE", "ARAGONITE", "DOLOMITE", "MAGNESITE", "\\bCHALK\\b", "LIMESTONE",
    "QUARTZ", "\\bSILICA\\b", "SILICON\\s*(\\(?IV\\)?\\s*)?(DI)?OXIDE", "CRISTOBALITE",
    "^GLASS", "KAOLIN", "\\bTALC\\b", "\\bMICA\\b", "MUSCOVITE", "BIOTITE",
    "(?<!HYPO)\\bCHLORITE\\b", "FELDSPAR", "ALBITE", "ORTHOCLASE", "MONTMORILLONITE",
    "BENTONITE", "ILLITE", "\\bCLAY\\b", "ZEOLITE", "(?<!ORTHO)SILICATE",
    "GYPSUM", "CALCIUM SULFATE", "BARIUM SULFATE", "BARYTE", "BARITE", "ANHYDRITE",
    "APATITE", "CALCIUM PHOSPHATE", "DIATOM"
  ),
  # --- Any other synthetic polymer (catch-all, keep LAST) ---------------------------
  "Other polymer" = c(
    paste0("^POLY(?!(SORBATE|OXYETHYLENE|GLYC|SACCHAR|GALACTURON|MYXIN|PHOSPH|",
           "SULFIDE|CHLORO(TER|BI)PHENYL|AMINE|ALKYLBENZENE|POX))"),
    "COPOLYMER", "TERPOLYMER", "\\bRESIN\\b",
    "^(PVDF|PVDC|PVA|PVAC|PVB|PVOH|EVOH|PEEK|PPS|PSU|PES|PEI|PLA|PCL|PBS|PHB|PHA|PBAT|PAN|PVP|PPO|PPE|PEO|PDMS|PI)$",
    "SILICONE", "SILOXANE", "DIMETHICONE"
  )
)


# ---- Category classification ----

synthetic_families    <- c("PET", "PP", "PE", "PS", "PVC", "PA", "PC",
                           "PMMA", "PU", "PTFE", "ABS", "Rubber",
                           "EVA", "POM", "Epoxy", "Other polymer")
semi_synthetic_families <- c("Cellulose", "Acrylate")
natural_families      <- c("Protein", "Chitin", "Natural")
additive_families     <- c("Additive", "Pigment")
inorganic_families    <- c("Mineral", "Inorganic")
other_compound_families <- c("Other compound")

# Reporting order of the categories (tables, summaries).
material_category_levels <- c("Synthetic", "Semi-synthetic", "Natural/Organic",
                              "Additive/Pigment", "Inorganic", "Other", "Unknown")

# Families assigned only from the spectral-library index (no name rule can
# produce them): a recognised entry of the inorganics library that no rule
# above claims is "Inorganic"; any other recognised entry is "Other compound"
# (a known reference substance — a drug, solvent, reagent — rather than an
# unidentified name, which stays "Unknown").
.inorganic_library_ids <- c("L60019")


# ---- Name normalisation ----

#' Normalise a material name for matching: upper case, typographic dashes
#' replaced by "-", runs of whitespace collapsed.
#'
#' @param x Character vector
#' @return Character vector
.ascii_dashes <- function(x) gsub("[\u2010-\u2015\u2212]", "-", x, perl = TRUE)

normalize_material_key <- function(x) {
  x <- enc2utf8(as.character(x))
  x[is.na(x)] <- ""
  # U+2010..U+2015 and U+2212: hyphen, non-breaking hyphen, figure dash, en/em
  # dash, horizontal bar, minus — the PDFs and some exports use them for "-".
  x <- .ascii_dashes(x)
  x <- toupper(x)
  trimws(gsub("\\s+", " ", x))
}

# Form the family patterns are tested against: brackets dropped and "POLY ("
# joined to its monomer, so "Poly(vinyl chloride)" reads "POLYVINYL CHLORIDE".
.flatten_for_rules <- function(x) {
  x <- normalize_material_key(x)
  x <- gsub("[][(){}]", " ", x)
  x <- trimws(gsub("\\s+", " ", x))
  gsub("\\bPOLY\\s+(?=[A-Z])", "POLY", x, perl = TRUE)
}

# Lookup key for the library index: only letters and digits survive, so case,
# spacing, brackets and punctuation differences never prevent a match.
.library_lookup_key <- function(x) gsub("[^A-Z0-9]", "", normalize_material_key(x))


# ---- Spectral-library index ----

.material_index_env <- new.env(parent = emptyenv())

#' Locate the spectral-library folder.
#'
#' Uses option `ftir.spectral_library_dir` when set (set it to NA to disable the
#' index, as the test suite does); otherwise looks for `spectral_libraries/` in
#' the working directory and up to two levels above it (main.R runs from the
#' repo root, the Shiny app from shiny_app/, tests from tests/testthat/).
#'
#' @return Directory path, or NULL when none is found / the index is disabled.
spectral_library_dir <- function() {
  opt <- getOption("ftir.spectral_library_dir")
  if (!is.null(opt)) {
    if (length(opt) != 1 || is.na(opt) || !nzchar(opt) || !dir.exists(opt)) return(NULL)
    return(normalizePath(opt, winslash = "/"))
  }
  cands <- c("spectral_libraries", file.path("..", "spectral_libraries"),
             file.path("..", "..", "spectral_libraries"))
  hit <- cands[dir.exists(cands)]
  if (length(hit) == 0) return(NULL)
  normalizePath(hit[[1]], winslash = "/")
}

#' Parse one S.T. Japan library listing (PDF) into entries.
#'
#' Each line of the listing is one library entry; synonyms of the same entry are
#' separated by " | ". Page headers/footers and the closing
#' legal notice are dropped. Entries too long for
#' the page are wrapped onto the next line; a line is re-joined to the one above
#' when that one reaches the right margin and the line breaks the listing's
#' alphabetical order.
#'
#' @param pdf_path Path to the PDF
#' @return data.frame(library, entry) — one row per entry, synonyms unsplit
parse_spectral_library_pdf <- function(pdf_path) {
  if (!requireNamespace("pdftools", quietly = TRUE))
    stop("Reading the spectral-library PDFs needs the 'pdftools' package: ",
         "install.packages('pdftools')")
  lib_id <- regmatches(basename(pdf_path), regexpr("L\\d{5}", basename(pdf_path)))
  if (length(lib_id) == 0) lib_id <- tools::file_path_sans_ext(basename(pdf_path))

  pages <- pdftools::pdf_text(pdf_path)
  lines <- unlist(strsplit(pages, "\n", fixed = TRUE), use.names = FALSE)
  lines <- trimws(.ascii_dashes(enc2utf8(lines)))
  junk <- !nzchar(lines) |
    grepl("S\\.T\\.Japan|stjapan\\.de|^L\\d{5}$|\\d+ of \\d+$", lines)
  lines <- lines[!junk]
  # The listing ends with the vendor's copyright notice and disclaimer.
  legal <- grep("^All information within this document|copyright", lines,
                ignore.case = TRUE)
  if (length(legal) > 0) lines <- lines[seq_len(legal[1] - 1L)]
  if (length(lines) == 0) return(data.frame(library = character(0), entry = character(0)))

  wrap_at <- 0.8 * max(nchar(lines))
  # Locale-independent "a sorts before b", on letters and digits only (the
  # listing's own ordering of brackets and punctuation is not byte order).
  key <- function(x) utf8ToInt(gsub("[^A-Z0-9]", "", toupper(x)))
  before <- function(a, b) {
    ia <- key(a); ib <- key(b)
    n <- min(length(ia), length(ib))
    d <- which(ia[seq_len(n)] != ib[seq_len(n)])[1]
    if (is.na(d)) length(ia) < length(ib) else ia[d] < ib[d]
  }
  out <- character(0)
  for (i in seq_along(lines)) {
    ln <- lines[i]
    k <- length(out)
    # A wrapped tail follows a line that reached the right margin and does not
    # fit the alphabetical order between the entry above and the line below.
    # (Some listings clip long names instead of wrapping; there the next line
    # is a regular entry and stays in order.)
    if (k > 0 && nchar(lines[i - 1L]) >= wrap_at &&
        (before(ln, out[k]) || (i < length(lines) && before(lines[i + 1L], ln)))) {
      sep <- if (grepl("-$", lines[i - 1L])) "" else " "
      out[k] <- paste0(out[k], sep, ln)
    } else {
      out <- c(out, ln)
    }
  }
  data.frame(library = lib_id, entry = out, stringsAsFactors = FALSE)
}

#' Build (or refresh) the spectral-library index CSV from the PDFs in `dir`.
#'
#' @param dir Folder holding the library PDFs
#' @param out_csv Where to write the index
#' @return The index data.frame(library, entry_id, entry), invisibly
build_spectral_library_index <- function(dir = spectral_library_dir(),
                                         out_csv = file.path(dir, "library_index.csv")) {
  if (is.null(dir) || !dir.exists(dir)) stop("Spectral-library folder not found")
  pdfs <- list.files(dir, pattern = "\\.pdf$", full.names = TRUE, ignore.case = TRUE)
  if (length(pdfs) == 0) stop("No library PDFs in ", dir)
  idx <- do.call(rbind, lapply(sort(pdfs), parse_spectral_library_pdf))
  idx$entry_id <- seq_len(nrow(idx))
  idx <- idx[, c("library", "entry_id", "entry")]
  utils::write.csv(idx, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(idx)
}

# Turn the entry table into the in-memory lookup: key → entry ids, and each
# entry's synonyms.
.prepare_library_index <- function(idx) {
  syn <- strsplit(idx$entry, "\\s*\\|\\s*")
  syn <- lapply(syn, function(s) unique(trimws(s[nzchar(trimws(s))])))
  rows <- rep(seq_len(nrow(idx)), lengths(syn))
  keys <- .library_lookup_key(unlist(syn, use.names = FALSE))
  # The full line is also a key, for exports that print the whole entry.
  keys <- c(keys, .library_lookup_key(idx$entry))
  rows <- c(rows, seq_len(nrow(idx)))
  ok <- nzchar(keys)
  list(
    library  = idx$library,
    synonyms = syn,
    by_key   = split(rows[ok], keys[ok])
  )
}

#' Load the spectral-library index, building or refreshing the CSV as needed.
#'
#' @param dir Folder holding the PDFs (and the cached library_index.csv)
#' @return Prepared index (list), or NULL when unavailable
load_spectral_library_index <- function(dir = spectral_library_dir()) {
  if (is.null(dir)) return(NULL)
  csv  <- file.path(dir, "library_index.csv")
  pdfs <- list.files(dir, pattern = "\\.pdf$", full.names = TRUE, ignore.case = TRUE)
  stale <- length(pdfs) > 0 &&
    (!file.exists(csv) || any(file.mtime(pdfs) > file.mtime(csv)))
  idx <- NULL
  if (stale && requireNamespace("pdftools", quietly = TRUE)) {
    idx <- tryCatch(build_spectral_library_index(dir, csv), error = function(e) {
      message("Spectral-library index: could not parse the PDFs (",
              conditionMessage(e), ")")
      NULL
    })
  } else if (stale && !file.exists(csv)) {
    message("Spectral-library index: PDFs found in ", dir, " but the 'pdftools' ",
            "package is not installed, so they are ignored. ",
            "install.packages('pdftools') to use them.")
  }
  if (is.null(idx) && file.exists(csv)) {
    idx <- utils::read.csv(csv, stringsAsFactors = FALSE, encoding = "UTF-8")
  }
  if (is.null(idx) || nrow(idx) == 0) return(NULL)
  message("Spectral-library index: ", nrow(idx), " entries from ",
          length(unique(idx$library)), " libraries (", dir, ")")
  .prepare_library_index(idx)
}

#' The spectral-library index in use (loaded once per session, on first use).
#'
#' @return Prepared index, or NULL when no library folder is available
get_spectral_library_index <- function() {
  dir <- spectral_library_dir()
  dir_key <- if (is.null(dir)) "" else dir
  if (!identical(.material_index_env$dir, dir_key)) {
    .material_index_env$index <- load_spectral_library_index(dir)
    .material_index_env$dir   <- dir_key
  }
  .material_index_env$index
}

#' Find a material name in the spectral-library index.
#'
#' @param name Character scalar
#' @param index Prepared index (see get_spectral_library_index())
#' @return List(synonyms, libraries) of every matching entry, or NULL
lookup_library_entry <- function(name, index = get_spectral_library_index()) {
  if (is.null(index)) return(NULL)
  # The whole name first: an export that prints the complete entry line must
  # not be mixed up with another entry that merely shares one synonym.
  key  <- .library_lookup_key(name)
  rows <- if (nzchar(key)) index$by_key[[key]]
  if (is.null(rows)) {
    parts <- strsplit(name, "\\s*\\|\\s*")[[1]]
    rows <- unique(unlist(index$by_key[intersect(.library_lookup_key(parts),
                                                  names(index$by_key))],
                          use.names = FALSE))
  }
  if (length(rows) == 0) return(NULL)
  list(synonyms  = unique(unlist(index$synonyms[rows], use.names = FALSE)),
       libraries = unique(index$library[rows]))
}


# ---- Family classification ----

# One combined regex per family (patterns OR-ed), compiled lazily per family list.
.family_regex <- function(families) {
  vapply(families, function(p) paste0("(?:", p, ")", collapse = "|"), character(1))
}

.first_family <- function(strings, fam_regex) {
  if (length(strings) == 0) return(NA_character_)
  for (fam in names(fam_regex)) {
    if (any(grepl(fam_regex[[fam]], strings, ignore.case = TRUE, perl = TRUE)))
      return(fam)
  }
  NA_character_
}

#' Classify a raw material name into a material family
#'
#' 1. A name that already is a family name ("PE", "Other polymer") is kept.
#' 2. The family rules are tried on the name (each "|"-separated part).
#' 3. If none matches and the spectral-library index is available, the rules are
#'    tried on every synonym of the matching library entry — that is how trade
#'    names ("LUPOLEN 6021 D") resolve to their polymer.
#' 4. A recognised library entry that no rule claims becomes "Inorganic" (from
#'    the inorganics library) or "Other compound".
#'
#' @param name Character scalar — raw material name
#' @param families Named list of character vectors (family → regex patterns)
#' @param index Spectral-library index (NULL = name rules only)
#' @return Character scalar — family name, or "Unknown"
classify_family <- function(name, families = default_polymer_families,
                            index = get_spectral_library_index()) {
  if (is.na(name) || !nzchar(trimws(name))) return("Unknown")
  fam_regex <- .family_regex(families)
  .classify_one(name, families, fam_regex, index)
}

.classify_one <- function(name, families, fam_regex, index) {
  key <- normalize_material_key(name)
  all_fams <- c(names(families), "Inorganic", "Other compound")
  same <- all_fams[toupper(all_fams) == key]
  if (length(same) > 0) return(same[[1]])

  parts <- trimws(strsplit(name, "|", fixed = TRUE)[[1]])
  fam <- .first_family(.flatten_for_rules(parts[nzchar(parts)]), fam_regex)
  if (!is.na(fam)) return(fam)

  hit <- lookup_library_entry(name, index)
  if (is.null(hit)) return("Unknown")
  fam <- .first_family(.flatten_for_rules(hit$synonyms), fam_regex)
  if (!is.na(fam)) return(fam)
  if (any(hit$libraries %in% .inorganic_library_ids)) "Inorganic" else "Other compound"
}


#' Vectorized family classification
#'
#' Each distinct name is classified once.
#'
#' @param names Character vector of raw material names
#' @param families Named list of character vectors
#' @param index Spectral-library index (NULL = name rules only)
#' @return Character vector of family names
classify_family_vec <- function(names, families = default_polymer_families,
                                index = get_spectral_library_index()) {
  names <- as.character(names)
  if (length(names) == 0) return(character(0))
  fam_regex <- .family_regex(families)
  u <- unique(names)
  fam_u <- vapply(u, function(n) {
    if (is.na(n) || !nzchar(trimws(n))) "Unknown"
    else .classify_one(n, families, fam_regex, index)
  }, character(1), USE.NAMES = FALSE)
  fam_u[match(names, u)]
}


#' Classify a material family into a reporting category
#'
#' @param family Character scalar — material family name
#' @return One of `material_category_levels`
classify_category <- function(family) {
  if (family %in% synthetic_families)      return("Synthetic")
  if (family %in% semi_synthetic_families) return("Semi-synthetic")
  if (family %in% natural_families)        return("Natural/Organic")
  if (family %in% additive_families)       return("Additive/Pigment")
  if (family %in% inorganic_families)      return("Inorganic")
  if (family %in% other_compound_families) return("Other")
  "Unknown"
}


#' Vectorized category classification
#'
#' @param families Character vector of material family names
#' @return Character vector of categories
classify_category_vec <- function(families) {
  vapply(families, classify_category, character(1), USE.NAMES = FALSE)
}


# Families that are typically compounded INTO a plastic: fillers (Mineral:
# CaCO3, talc, kaolin, BaSO4, silica), pigments (TiO2, iron oxides, organic
# pigments) and additives. When one instrument reports such a family and the
# other a plastic, both can be right about the same particle — e.g. Raman sees
# the strongly scattering TiO2 while FTIR sees the polymer around it.
filler_pigment_families <- c("Mineral", "Pigment", "Additive")

#' Is a pair "plastic on one side, filler/pigment/additive on the other"?
#'
#' @param fam_a,fam_b Character vectors of material families (same length)
#' @return Logical vector
is_filler_pigment_pair <- function(fam_a, fam_b) {
  plastic <- c(synthetic_families, semi_synthetic_families)
  (fam_a %in% plastic & fam_b %in% filler_pigment_families) |
    (fam_b %in% plastic & fam_a %in% filler_pigment_families)
}

#' Vectorized tiered agreement scoring
#'
#' Tiers, best first:
#'   Exact          — same family and the same name
#'   Family         — same family (e.g. "HDPE" vs "Polyethylene")
#'   Filler/pigment — one side a plastic, the other a filler, pigment or
#'                    additive (see filler_pigment_families): compatible, not
#'                    a confirmation of the polymer
#'   Disagree       — anything else
#'
#' @param ftir_names Character vector of raw FTIR material names
#' @param raman_names Character vector of raw Raman material names
#' @param families Polymer family definitions
#' @return Character vector: "Exact", "Family", "Filler/pigment" or "Disagree"
score_tiered_agreement_vec <- function(ftir_names, raman_names,
                                       families = default_polymer_families) {
  ftir_fam  <- classify_family_vec(ftir_names, families)
  raman_fam <- classify_family_vec(raman_names, families)

  result <- rep("Disagree", length(ftir_names))
  result[is_filler_pigment_pair(ftir_fam, raman_fam)] <- "Filler/pigment"
  same_fam <- ftir_fam == raman_fam & ftir_fam != "Unknown"
  result[same_fam] <- "Family"

  # Exact: same family AND same canonical (uppercased trimmed) name
  ftir_canon  <- toupper(trimws(ftir_names))
  raman_canon <- toupper(trimws(raman_names))
  exact <- same_fam & (ftir_canon == raman_canon)
  result[exact] <- "Exact"

  result
}


#' Map raw material names to canonical types using a user-provided mapping
#'
#' Each entry in `mapping` is: canonical_name = "regex pattern".
#' The first matching pattern wins. Unmatched names pass through unchanged
#' (uppercased + trimmed).
#'
#' @param raw_names Character vector of raw material names
#' @param mapping Named list: names = canonical types, values = regex patterns
#' @return Character vector of canonical material names
map_materials <- function(raw_names, mapping) {
  if (is.null(mapping) || length(mapping) == 0) {
    return(NULL)  # signal: use default normalization
  }

  result <- rep(NA_character_, length(raw_names))

  for (canonical in names(mapping)) {
    pattern <- mapping[[canonical]]
    hits <- grepl(pattern, raw_names, ignore.case = TRUE)
    # Only assign if not already mapped (first match wins)
    result[hits & is.na(result)] <- toupper(canonical)
  }

  # Anything unmatched: fall back to uppercase trimmed raw name
  unmapped <- is.na(result)
  if (any(unmapped)) {
    result[unmapped] <- toupper(trimws(raw_names[unmapped]))
  }

  result
}


#' Apply material mapping to matched pairs for agreement scoring
#'
#' If material maps are configured, uses them. Otherwise falls back to
#' the abbreviation-based normalize_material() from 08_agreement.R.
#'
#' @param ftir_materials Character vector of raw FTIR material names
#' @param raman_materials Character vector of raw Raman material names
#' @param config Configuration list with material_map_ftir / material_map_raman
#' @return List with ftir_canonical and raman_canonical character vectors
map_paired_materials <- function(ftir_materials, raman_materials, config) {
  ftir_mapped  <- map_materials(ftir_materials, config$material_map_ftir)
  raman_mapped <- map_materials(raman_materials, config$material_map_raman)

  # Fall back to abbreviation normalization if no mapping was configured
  if (is.null(ftir_mapped))  ftir_mapped  <- normalize_material(ftir_materials)
  if (is.null(raman_mapped)) raman_mapped <- normalize_material(raman_materials)

  list(
    ftir_canonical  = ftir_mapped,
    raman_canonical = raman_mapped
  )
}
