# Per-material show/hide checkbox groups (shiny_app/app.R).
#
# The logic lives inside server(), so it cannot be sourced directly. These
# tests pin the state machine it implements — the part that is easy to get
# wrong, because an empty checkbox group reads back as NULL and that NULL is
# ambiguous:
#
#   NULL + never populated  -> no filter yet   -> show everything
#   NULL + populated        -> all unticked    -> show nothing
#
# Collapsing those two would either hide every particle while a run loads, or
# silently ignore the user unticking the whole list.

# Mirrors mat_keep(): resolve a group to the families to keep (NULL = no filter).
.mat_keep <- function(sel, ready) {
  if (is.null(sel)) {
    if (!isTRUE(ready)) return(NULL)
    return(character(0))
  }
  sel
}

# Mirrors update_material_choices(): keep existing ticks, tick families new to
# this run, drop families it no longer has.
.pick_selected <- function(prev_choices, prev_sel, families) {
  if (is.null(prev_choices)) families
  else union(intersect(prev_sel, families), setdiff(families, prev_choices))
}

.apply <- function(fams, keep) if (is.null(keep)) fams else fams[fams %in% keep]

SAMPLE <- c("PE", "PE", "PP", "PET", "unknown")

test_that("an unpopulated group does not filter anything", {
  expect_null(.mat_keep(NULL, FALSE))
  expect_identical(.apply(SAMPLE, .mat_keep(NULL, FALSE)), SAMPLE)
})

test_that("a freshly populated group ticks every family present", {
  fams <- sort(unique(SAMPLE))
  sel  <- .pick_selected(NULL, NULL, fams)
  expect_identical(sel, fams)
  expect_identical(.apply(SAMPLE, .mat_keep(sel, TRUE)), SAMPLE)
})

test_that("unticking one family hides only that family", {
  fams <- sort(unique(SAMPLE))
  expect_identical(.apply(SAMPLE, .mat_keep(setdiff(fams, "PE"), TRUE)),
                   c("PP", "PET", "unknown"))
  # "unknown" is switchable like any polymer — the point of the feature.
  expect_identical(.apply(SAMPLE, .mat_keep(setdiff(fams, "unknown"), TRUE)),
                   c("PE", "PE", "PP", "PET"))
})

test_that("unticking everything shows nothing, not everything", {
  keep <- .mat_keep(NULL, TRUE)
  expect_identical(keep, character(0))
  expect_length(.apply(SAMPLE, keep), 0L)
})

test_that("switching runs keeps ticks and admits new families ticked", {
  fams <- sort(unique(SAMPLE))            # PE, PET, PP, unknown
  sel  <- setdiff(fams, "PE")            # user hid PE
  new  <- c("PE", "PP", "PET", "unknown", "PS")

  sel2 <- .pick_selected(fams, sel, new)
  expect_false("PE" %in% sel2)           # deliberate choice survives
  expect_true("PS" %in% sel2)            # new family is visible by default
  expect_true(all(c("PP", "PET", "unknown") %in% sel2))

  # A family the new run lacks simply drops out.
  sel3 <- .pick_selected(new, sel2, c("PP", "PET"))
  expect_identical(sort(sel3), c("PET", "PP"))
})

test_that("the two NULL states are distinguishable as cache keys", {
  # The overlay plot is cached; keying on the raw input would collide these
  # two states and serve a stale plot.
  expect_false(identical(.mat_keep(NULL, FALSE), .mat_keep(NULL, TRUE)))
})
