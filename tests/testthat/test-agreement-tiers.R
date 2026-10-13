# Tiered material agreement (R/08b_material_map.R, R/08_agreement.R).
#
# Primary agreement = family hit (Exact + Family); secondary = exact polymer.
# "Filler/pigment" is a plastic on one instrument and a filler, pigment or
# additive on the other (Raman sees the TiO2 in a particle whose polymer FTIR
# identifies): neither a confirmation nor a disagreement. "Unresolved" is a
# pair the names cannot settle (unrecognised name, or two different polymers
# in the "Other polymer" catch-all): flagged for inspection, outside the rates.
# Low-quality and unidentified particles never enter the denominator.

sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = globalenv())
sys.source(file.path(REPO_ROOT, "R", "08_agreement.R"), envir = globalenv())

test_that("tiers: exact, family, filler/pigment, unresolved, disagree", {
  ftir  <- c("PET", "HDPE", "Polypropylene", "Polystyrene", "TALC",
             "PA", "PE", "Unknown", "PE", "Calcium carbonate")
  raman <- c("PET", "Polyethylene", "TITANIUM(IV) OXIDE ANATASE", "CALCIUM CARBONATE",
             "Polypropylene", "Irganox 1010", "Nylon 6", "TALC", "Sodium chloride",
             "Kaolin")
  expect_identical(score_tiered_agreement_vec(ftir, raman),
                   c("Exact", "Family", "Filler/pigment", "Filler/pigment",
                     "Filler/pigment", "Filler/pigment", "Disagree", "Unresolved",
                     "Unresolved", "Family"))
})

test_that("the same polymer under different spellings is an exact hit", {
  a <- c("PET", "Polypro", "Polycarbonate (PC)", "Polyvinyl alcohol", "Nylon 6")
  b <- c("Polyethylene terephtalate", "Polypropylene", "PC", "PVOH", "Polyamide 6")
  expect_identical(score_tiered_agreement_vec(a, b), rep("Exact", 5))
})

test_that("exact polymer and its upstream family are both reported", {
  x <- c("Polybutylene terephthalate", "Chlorinated polyethylene", "HDPE", "PET")
  expect_identical(canonical_polymer_vec(x), c("PBT", "CPE", "HDPE", "PET"))
  expect_identical(classify_family_vec(x), c("PET", "PE", "PE", "PET"))
  # PBT vs PET: family hit, not exact
  expect_identical(score_tiered_agreement_vec("PBT", "PET"), "Family")
  d <- annotate_materials(data.frame(material = c(x, "unknown", NA)))
  expect_identical(d$material_polymer, c("PBT", "CPE", "HDPE", "PET", "UNKNOWN", "UNKNOWN"))
  expect_identical(d$material_family, c("PET", "PE", "PE", "PET", "Unknown", "Unknown"))
})

test_that("two different polymers in the catch-all are unresolved, not a family hit", {
  expect_identical(score_tiered_agreement_vec("Polyvinyl alcohol", "Silicone"), "Unresolved")
  expect_identical(score_tiered_agreement_vec("Silicone", "PE"), "Disagree")
  expect_false(family_agrees("Other polymer", "Other polymer"))
  expect_false(family_agrees("Unknown", "Unknown"))
  expect_true(family_agrees("PE", "PE"))
})

.pairs <- function(ftir, raman, fq = 0.9, rq = 90, ...) {
  data.frame(match_id = seq_along(ftir), ftir_material = ftir, raman_material = raman,
             ftir_quality = fq, raman_quality = rq, match_distance = 5,
             ftir_feret_max_um = 120, ..., stringsAsFactors = FALSE)
}

test_that("rates are over resolved pairs; family is primary, exact secondary", {
  m <- .pairs(c("PET", "PP", "PP", "PE", "PS", "Silicone"),
              c("PET", "TITANIUM(IV) OXIDE ANATASE", "Polypropylene", "Nylon 6",
                "Polystyrene", "Polyvinyl alcohol"))
  ag <- suppressMessages(analyze_agreement(list(matched = m)))
  tr <- ag$tiered_rates
  expect_identical(c(tr$n_total, tr$n_unresolved, tr$n_resolved), c(6L, 1L, 5L))
  expect_identical(c(tr$n_exact, tr$n_family, tr$n_filler, tr$n_disagree), c(3L, 0L, 1L, 1L))
  expect_equal(tr$family_or_better_pct, 60)        # 3 of 5 resolved
  expect_equal(tr$exact_pct, 60)
  expect_equal(ag$agreement_rate, 0.6)
  expect_identical(ag$summary_df$tier,
                   c("Exact", "Filler/pigment", "Exact", "Disagree", "Exact", "Unresolved"))
  expect_true(all(c("polymer_a", "family_a", "agree_family", "agree_exact") %in% names(ag$summary_df)))
})

test_that("low-quality and unidentified particles never enter the denominator", {
  m <- .pairs(c("PET", "PET", "PET", "PET", "unknown"),
              c("PET", "PET", "PE",  "PET", "PET"),
              fq = c(0.9, 0.5, 0.9, 0.9, 0.9), rq = c(90, 90, 60, NA, 90))
  ag <- suppressMessages(analyze_agreement(list(matched = m), make_config()))
  tr <- ag$tiered_rates
  # 0.5 < FTIR 0.70; 60 < Raman 70; NA Raman HQI; "unknown" label
  expect_identical(tr$n_excluded_quality_a, 1L)
  expect_identical(tr$n_excluded_quality_b, 2L)
  expect_identical(tr$n_excluded_unidentified, 1L)
  expect_identical(tr$n_total, 1L)
  expect_equal(tr$family_or_better_pct, 100)
})

test_that("LDIR quality uses the Agilent unknown threshold; join warnings are counted", {
  m <- data.frame(match_id = 1:3, ldir_material = c("PET", "PET", "PET"),
                  raman_material = c("PET", "PET", "PET"),
                  ldir_quality = c(0.95, 0.80, 0.90), raman_quality = 90,
                  ldir_coord_join_warning = c(TRUE, FALSE, FALSE),
                  match_distance = 5, ldir_feret_max_um = 100)
  ag <- suppressMessages(analyze_agreement(list(matched = m), make_config(),
                                           instrument_a = "LDIR", instrument_b = "Raman"))
  expect_identical(ag$tiered_rates$n_excluded_quality_a, 1L)   # 0.80 < 0.85
  expect_identical(ag$tiered_rates$n_total, 2L)
  expect_identical(ag$tiered_rates$n_join_warning, 1L)        # kept, flagged
  expect_identical(ag$summary_df$coord_join_warning, c(TRUE, FALSE))
})

test_that("only plastic families pair with fillers", {
  expect_true(all(is_filler_pigment_pair(c("PE", "Cellulose", "Other polymer"),
                                         c("Pigment", "Mineral", "Additive"))))
  expect_false(any(is_filler_pigment_pair(
    c("Natural", "Unknown", "Mineral", "Inorganic", "Other compound"),
    c("Pigment", "Mineral", "Pigment", "PE", "PE"))))
})
