# =============================================================================
# 08b_material_map.R — Cross-instrument material equivalence mapping
# =============================================================================
#
# FTIR and Raman identify materials using different spectral libraries and
# nomenclature. This module provides a configurable mapping from each
# instrument's raw material names to a shared canonical vocabulary, so that
# agreement scoring compares like with like.
#
# The actual mapping tables (config$material_map_ftir, config$material_map_raman)
# should be populated by the user. When NULL, the pipeline falls back to the
# existing abbreviation-based normalization in 08_agreement.R.
# =============================================================================

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
