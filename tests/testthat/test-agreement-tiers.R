# Tiered material agreement (R/08b_material_map.R, R/08_agreement.R).
#
# "Filler/pigment" is the tier for a plastic on one instrument and a filler,
# pigment or additive on the other — typically Raman reporting the TiO2 or
# CaCO3 compounded into a particle whose polymer FTIR identifies. It must not
# be lumped in with "Disagree", nor counted as a confirmation of the polymer.

sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = globalenv())
sys.source(file.path(REPO_ROOT, "R", "08_agreement.R"), envir = globalenv())

test_that("plastic vs filler/pigment/additive gets its own tier, both ways round", {
  ftir  <- c("PET", "HDPE", "Polypropylene", "Polystyrene", "TALC",
             "PA", "PE", "Unknown", "PE", "Calcium carbonate")
  raman <- c("PET", "Polyethylene", "TITANIUM(IV) OXIDE ANATASE", "CALCIUM CARBONATE",
             "Polypropylene", "Irganox 1010", "Nylon 6", "TALC", "Sodium chloride",
             "Kaolin")
  expect_identical(score_tiered_agreement_vec(ftir, raman),
                   c("Exact", "Family", "Filler/pigment", "Filler/pigment",
                     "Filler/pigment", "Filler/pigment", "Disagree", "Disagree",
                     "Disagree", "Family"))
})

test_that("only plastic families pair with fillers", {
  expect_true(all(is_filler_pigment_pair(c("PE", "Cellulose", "Other polymer"),
                                         c("Pigment", "Mineral", "Additive"))))
  expect_false(any(is_filler_pigment_pair(
    c("Natural", "Unknown", "Mineral", "Inorganic", "Other compound"),
    c("Pigment", "Mineral", "Pigment", "PE", "PE"))))
})

test_that("analyze_agreement counts the filler/pigment tier separately", {
  matched <- data.frame(
    match_id       = 1:5,
    ftir_material  = c("PET", "PP", "PP", "PE", "PS"),
    raman_material = c("PET", "TITANIUM(IV) OXIDE ANATASE", "TALC", "Nylon 6", "Polystyrene"),
    raman_quality  = 900,
    match_distance = 5,
    ftir_feret_max_um = c(40, 60, 150, 250, 400),
    stringsAsFactors = FALSE)
  ag <- suppressMessages(analyze_agreement(list(matched = matched)))
  tr <- ag$tiered_rates
  expect_identical(c(tr$n_exact, tr$n_family, tr$n_filler, tr$n_disagree),
                   c(1L, 1L, 2L, 1L))
  expect_equal(tr$filler_pct, 40)
  expect_equal(tr$family_or_better_pct, 40)   # fillers do not confirm the polymer
  expect_identical(ag$summary_df$tier,
                   c("Exact", "Filler/pigment", "Filler/pigment", "Disagree", "Family"))
  pp <- ag$agreement_detail[ag$agreement_detail$material_a == "PP", ]
  expect_identical(c(pp$n_filler, pp$n_disagree), c(2L, 0L))
  expect_identical(sum(ag$size_analysis$n_disagree), 1L)
  expect_identical(sum(ag$size_analysis$n_filler), 2L)
})
