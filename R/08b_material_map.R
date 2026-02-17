# =============================================================================
# 08b_material_map.R — Cross-instrument material equivalence mapping
# =============================================================================
#
# FTIR, Raman, and LDIR identify materials using different spectral libraries
# and nomenclature. This module provides:
#   1. A configurable regex mapping from raw names → canonical types
#   2. A polymer family classification (PET → Polyester, PA → Polyamide, ...)
#   3. A tiered agreement scorer: Exact / Family / Disagree
#   4. A category classifier: Synthetic / Semi-synthetic / Natural / Unknown
#
# The mapping tables are stored in config (material_map_ftir, material_map_raman,
# material_map_ldir). When NULL, the pipeline uses built-in defaults below.
# =============================================================================


# ---- Default polymer family definitions ----
# Each family maps a canonical family name to a character vector of patterns
# (case-insensitive regex). First match wins.

default_polymer_families <- list(
  PET = c(
    "^PET$",
    "Polyethylene\\s*tereph",
    "Polyester"
  ),
  PP = c(
    "^PP$",
    "Polypro"
  ),
  PE = c(
    "^PE$", "^HDPE$", "^LDPE$", "^LLDPE$",
    "Polyethylene(?!.*tereph)",  # PE but NOT PET
    "High.density polyethylene",
    "Low.density polyethylene"
  ),
  PS = c(
    "^PS$",
    "Polystyrene"
  ),
  PVC = c(
    "^PVC$",
    "Polyvinyl\\s*chloride"
  ),
  PA = c(
    "^PA$", "^PA\\d", "^Nylon",
    "^Polyamide(?!.*natural)"   # PA but NOT "naturally occurring"
  ),
  PC = c(
    "^PC$",
    "Polycarbonate"
  ),
  PMMA = c(
    "^PMMA$",
    "Polymethyl\\s*methacrylate",
    "^Acrylic$"
  ),
  PU = c(
    "^PU$",
    "Polyurethane"
  ),
  PTFE = c(
    "^PTFE$",
    "Polytetrafluoroethylene",
    "Teflon"
  ),
  ABS = c(
    "^ABS$",
    "Acrylonitrile.*butadiene.*styrene"
  ),
  Rubber = c(
    "Rubber", "^SBR$", "^NBR$", "^EPDM$",
    "Tire", "Tyre"
  ),
  Cellulose = c(
    "Cellulos",     # matches Cellulose, Cellulosic, Cellulose Acetate, CAB
    "^CAB$",
    "Rayon", "Viscose", "Lyocell"
  ),
  Acrylate = c(
    "Polyacrylamide", "Polyacrylate",
    "Acrylamide", "Acrylate"
  ),
  Protein = c(
    "Protein", "Keratin", "Silk", "^Wool$"
  ),
  Chitin = c(
    "Chitin", "Chitosan"
  ),
  Natural = c(
    "naturally\\s*occurring",
    "^Carbonate$", "^Sand$",
    "Cotton",
    "^Natural$", "^Organic$", "^Biofilm$"
  )
)


# ---- Category classification ----

synthetic_families    <- c("PET", "PP", "PE", "PS", "PVC", "PA", "PC",
                           "PMMA", "PU", "PTFE", "ABS", "Rubber")
semi_synthetic_families <- c("Cellulose", "Acrylate")
natural_families      <- c("Protein", "Chitin", "Natural")


#' Classify a raw material name into a polymer family
#'
#' @param name Character scalar — raw material name
#' @param families Named list of character vectors (family → regex patterns)
#' @return Character scalar — family name, or "Unknown"
classify_family <- function(name, families = default_polymer_families) {
  if (is.na(name) || !nzchar(trimws(name))) return("Unknown")
  name_clean <- trimws(name)
  for (fam in names(families)) {
    for (pat in families[[fam]]) {
      if (grepl(pat, name_clean, ignore.case = TRUE, perl = TRUE)) {
        return(fam)
      }
    }
  }
  "Unknown"
}


#' Vectorized family classification
#'
#' @param names Character vector of raw material names
#' @param families Named list of character vectors
#' @return Character vector of family names
classify_family_vec <- function(names, families = default_polymer_families) {
  vapply(names, classify_family, character(1), families = families,
         USE.NAMES = FALSE)
}


#' Classify a polymer family into a reporting category
#'
#' @param family Character scalar — polymer family name
#' @return One of "Synthetic", "Semi-synthetic", "Natural/Organic", "Unknown"
classify_category <- function(family) {
  if (family %in% synthetic_families)      return("Synthetic")
  if (family %in% semi_synthetic_families) return("Semi-synthetic")
  if (family %in% natural_families)        return("Natural/Organic")
  "Unknown"
}


#' Vectorized category classification
#'
#' @param families Character vector of polymer family names
#' @return Character vector of categories
classify_category_vec <- function(families) {
  vapply(families, classify_category, character(1), USE.NAMES = FALSE)
}


#' Score tiered agreement between two material names
#'
#' Returns one of:
#'   "Exact"   — same canonical family AND same canonical name
#'   "Family"  — same family but different canonical name
#'   "Disagree" — different families
#'
#' @param ftir_name Raw FTIR material name
#' @param raman_name Raw Raman material name
#' @param families Polymer family definitions
#' @return Character scalar: "Exact", "Family", or "Disagree"
score_tiered_agreement <- function(ftir_name, raman_name,
                                   families = default_polymer_families) {
  ftir_fam  <- classify_family(ftir_name, families)
  raman_fam <- classify_family(raman_name, families)

  if (ftir_fam == "Unknown" || raman_fam == "Unknown") return("Disagree")
  if (ftir_fam != raman_fam) return("Disagree")

  # Same family — check if canonical names match exactly
  ftir_canon  <- toupper(trimws(ftir_name))
  raman_canon <- toupper(trimws(raman_name))
  if (ftir_canon == raman_canon) return("Exact")

  # Same family, different raw names: still "Family"
  "Family"
}


#' Vectorized tiered agreement scoring
#'
#' @param ftir_names Character vector of raw FTIR material names
#' @param raman_names Character vector of raw Raman material names
#' @param families Polymer family definitions
#' @return Character vector: "Exact", "Family", or "Disagree"
score_tiered_agreement_vec <- function(ftir_names, raman_names,
                                       families = default_polymer_families) {
  ftir_fam  <- classify_family_vec(ftir_names, families)
  raman_fam <- classify_family_vec(raman_names, families)

  result <- rep("Disagree", length(ftir_names))
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
