# =============================================================================
# test-file-picker-cancel.R -- escaping the required-file dialog
# =============================================================================
# A required slot (Raman) used to re-prompt forever when the operator pressed
# Cancel: the dialog reopened immediately with no way out, so the only escape
# was killing the R process. Cancelling twice in a row now aborts the run.
#
# The picker lives inside an `if (input_mode == "explicit")` block in main.R,
# so it is extracted by source text -- the same approach test-plot-bounds.R
# uses for the viewport guards.

.load_picker <- function(interactive_flag = TRUE) {
  mpath <- file.path(REPO_ROOT, "main.R")
  skip_if_not(file.exists(mpath), "main.R not found")
  src <- readLines(mpath, warn = FALSE)
  i <- grep("^  # Abort the run cleanly", src)[1]
  j <- grep("^  DATA_FILTER  <- ", src)[1]
  if (is.na(i) || is.na(j))
    stop("The file-picker block was not found in main.R: its marker lines changed; update this test.")
  env <- new.env(parent = globalenv())
  eval(parse(text = gsub("^  ", "", paste(src[i:(j - 1)], collapse = "\n"))),
       envir = env)
  env$.is_windows <- TRUE
  env$.last_dir   <- NULL
  env$interactive <- function() interactive_flag
  env
}

# Drive the picker with a queue of dialog results; "" means the user cancelled.
.run_picker <- function(responses, required, interactive_flag = TRUE) {
  env <- .load_picker(interactive_flag)
  q <- responses
  env$choose.files <- function(...) {
    v <- q[1]; q <<- q[-1]
    if (is.na(v)) "" else v
  }
  suppressMessages(tryCatch(
    env$.pick_file("test", required = required),
    pipeline_cancelled = function(c) "CANCELLED"))
}

test_that("an optional slot returns NULL on Cancel", {
  expect_null(.run_picker(c(""), required = FALSE))
})

test_that("a required slot still accepts a file", {
  f <- tempfile(fileext = ".csv"); file.create(f)
  expect_identical(.run_picker(c(f), required = TRUE), f)
})

test_that("one Cancel on a required slot re-prompts rather than aborting", {
  # The operator gets a second chance -- a mis-click must not kill the run.
  f <- tempfile(fileext = ".csv"); file.create(f)
  expect_identical(.run_picker(c("", f), required = TRUE), f)
})

test_that("two consecutive Cancels on a required slot abort the pipeline", {
  # The regression: this used to loop forever.
  expect_identical(.run_picker(c("", ""), required = TRUE), "CANCELLED")
})

test_that("a non-interactive session aborts on the first Cancel", {
  # Nobody is there to press Cancel a second time, so re-prompting would spin.
  expect_identical(.run_picker(c(""), required = TRUE,
                               interactive_flag = FALSE), "CANCELLED")
})

test_that("the cancel condition is classed, not a bare error", {
  env <- .load_picker()
  cond <- tryCatch(suppressMessages(env$.cancel_pipeline("test reason")),
                   condition = function(c) c)
  expect_s3_class(cond, "pipeline_cancelled")
  expect_match(conditionMessage(cond), "test reason")
})

# The whole explicit-input block: instrument list, slots, confirmation
.run_input_block <- function(instruments, files, confirm = "YES") {
  src <- readLines(file.path(REPO_ROOT, "main.R"), warn = FALSE)
  i <- grep("^  # Abort the run cleanly", src)[1]
  j <- grep("^\\} else if \\(input_mode == \"interactive\"\\)", src)[1]
  env <- new.env(parent = globalenv())
  env$.is_windows <- TRUE; env$.last_dir <- NULL
  env$interactive <- function() TRUE
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  env$select.list <- function(...) instruments
  q <- files; asked <- character(0)
  env$choose.files <- function(caption, ...) {
    asked <<- c(asked, caption); v <- q[1]; q <<- q[-1]; if (is.na(v)) "" else v
  }
  shown <- NULL
  env$winDialog <- function(type, message) { shown <<- message; confirm }
  txt <- gsub("utils::(select\\.list|winDialog)", "\\1", paste(src[i:(j - 1)], collapse = "\n"))
  out <- suppressMessages(tryCatch({ eval(parse(text = txt), envir = env); "RAN" },
                                   pipeline_cancelled = function(c) conditionMessage(c)))
  list(out = out, env = env, asked = asked, shown = shown)
}

test_that("only the chosen instruments' dialogs open, and the files are confirmed", {
  r <- .run_input_block(c("Raman", "LDIR"), c("r.csv", "", "l.csv", "l.png"))
  expect_identical(r$out, "RAN")
  expect_length(r$asked, 4)
  expect_true(all(grepl("^(Raman|LDIR)", r$asked)))
  expect_identical(r$env$raman_file, "r.csv"); expect_null(r$env$raman_image)
  expect_null(r$env$ftir_file)
  expect_match(r$shown, "l.png", fixed = TRUE)
})

test_that("Cancel on the instrument list or No at the confirmation ends the run at once", {
  expect_match(.run_input_block(character(0), character(0))$out, "no instruments selected")
  expect_match(.run_input_block("Raman", character(0))$out, "at least two")
  r <- .run_input_block(c("Raman", "LDIR"), c("wrong.csv", "", "l.csv", ""), confirm = "NO")
  expect_match(r$out, "not confirmed")
})
