# Reference-material picker in tools/reproducibility.R.
#
# The reference used to be typed free-hand. Anything classify_family() did not
# recognise fell through to "Unknown" and the accuracy metric then scored
# nothing — a silent, spelling-dependent failure. It is now a numbered menu,
# like the instrument question, offering polymer FAMILY codes.
#
# The property that makes it typo-proof is pinned below: every code offered
# classifies back to itself, so the value handed to the accuracy metric is
# exactly the family the user pointed at.

.load_ref_picker <- function() {
  tpath <- file.path(REPO_ROOT, "tools", "reproducibility.R")
  skip_if_not(file.exists(tpath), "tools/reproducibility.R not found")
  src <- readLines(tpath, warn = FALSE)
  i <- grep("^\\.REFERENCE_MATERIALS <- list\\(", src)[1]
  j <- grep("^#' Prompt for instrument type", src)[1]
  skip_if(is.na(i) || is.na(j), "picker block not found")
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(src[i:(j - 1)], collapse = "\n")), envir = env)
  sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = env)
  env
}

# Drive the picker with a scripted queue of menu answers.
.run_picker <- function(env, keys) {
  q <- keys
  fake <- function(choices, title = NULL) { a <- q[1]; q <<- q[-1]; a }
  suppressMessages(env$select_reference_materials(.menu = fake))
}

.n_opts  <- function(env) length(env$.REFERENCE_MATERIALS)
.done_key <- function(env) .n_opts(env) + 1L

test_that("every offered polymer classifies back to itself", {
  # This is the whole point: the stored value goes straight through
  # classify_family_vec(), so a code that did not round-trip would silently
  # score against the wrong family — or none at all.
  env <- .load_ref_picker()
  codes <- vapply(env$.REFERENCE_MATERIALS, `[`, character(1), 1L)
  expect_identical(unname(env$classify_family_vec(codes)), codes)
  expect_false(any(env$classify_family_vec(codes) == "Unknown"))
})

test_that("the list covers the synthetic and semi-synthetic families", {
  # Those are exactly the families the accuracy metric can score, so anything
  # missing here would be unreachable from the menu.
  env <- .load_ref_picker()
  codes <- vapply(env$.REFERENCE_MATERIALS, `[`, character(1), 1L)
  expect_setequal(codes, c(env$synthetic_families, env$semi_synthetic_families))
  expect_false(anyDuplicated(codes) > 0)
  # Each entry is a code plus a human-readable label.
  expect_true(all(vapply(env$.REFERENCE_MATERIALS, length, integer(1)) == 2L))
  expect_true(all(nzchar(vapply(env$.REFERENCE_MATERIALS, `[`, character(1), 2L))))
})

test_that("the most common microplastic polymers come first", {
  env <- .load_ref_picker()
  codes <- vapply(env$.REFERENCE_MATERIALS, `[`, character(1), 1L)
  expect_identical(codes[1:6], c("PE", "PP", "PET", "PS", "PVC", "PA"))
})

test_that("a monotype filter is one number then DONE", {
  env <- .load_ref_picker()
  expect_identical(.run_picker(env, c(2L, .done_key(env))), "PP")
})

test_that("a mixed-polymer standard accumulates, in pick order", {
  env <- .load_ref_picker()
  expect_identical(.run_picker(env, c(1L, 2L, 3L, .done_key(env))),
                   c("PE", "PP", "PET"))
})

test_that("skipping yields NULL, so accuracy is simply not reported", {
  env <- .load_ref_picker()
  expect_null(.run_picker(env, 0L))                    # cancelled
  expect_null(.run_picker(env, .done_key(env)))        # DONE, nothing picked
  # Cancelling discards a part-built selection rather than half-applying it.
  expect_null(.run_picker(env, c(1L, 2L, 0L)))
})

test_that("picking an already-selected polymer unselects it", {
  env <- .load_ref_picker()
  d <- .done_key(env)
  expect_identical(.run_picker(env, c(2L, 2L, 1L, d)), "PE")   # mis-click undone
  expect_identical(.run_picker(env, c(5L, 5L, 5L, d)), "PVC")  # odd taps stay on
  expect_null(.run_picker(env, c(3L, 3L, d)))                  # back to none
})

test_that("the DONE entry sits after every polymer, whatever the list length", {
  env <- .load_ref_picker()
  seen <- NULL
  fake <- function(choices, title = NULL) { seen <<- choices; .done_key(env) }
  suppressMessages(env$select_reference_materials(.menu = fake))
  expect_length(seen, .n_opts(env) + 1L)
  expect_match(seen[length(seen)], "^DONE")
  expect_false(any(grepl("^DONE", seen[-length(seen)])))
})
