# =============================================================================
# validate_example_A.R -- validation checks on the bundled example dataset
# =============================================================================
# Runs the pipeline on examples/dataset_A (without the report) and the checks
# in tools/validation.R: chance matches (null model), gate and alignment
# sensitivity, LDIR join check. From the repository root:
#
#   Rscript examples/validate_example_A.R            # -> output/validation_A/
#   Rscript examples/validate_example_A.R /tmp/val
#
# Takes a few minutes (the alignment sweep re-runs every alignment).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else file.path("output", "validation_A")
D <- file.path("examples", "dataset_A")
if (!dir.exists(D)) stop("Run from the repository root (examples/dataset_A not found).")

source("R/load.R"); load_pipeline()
source("tools/validation.R")

config <- make_config(ftir_path  = file.path(D, "FTIR_perkin.csv"),
                      raman_path = file.path(D, "Raman_data.csv"),
                      output_dir = file.path(out, "run"))
config$ftir_bruker_path <- file.path(D, "FTIR_bruker.csv")
config$ldir_path        <- file.path(D, "LDIR_data.xlsx")
config$ftir_image       <- file.path(D, "FTIR_image.png")
config$raman_image      <- file.path(D, "Raman_image_small.png")
config$ldir_image       <- file.path(D, "LDIR_image.png")

result <- run_pipeline(config, report = FALSE)
v <- validate_run(result, out_dir = out)

show <- function(title, d) { cat("\n==", title, "==\n"); print(d, row.names = FALSE) }
show("Chance matches (null model)", v$null_model)
show("Gate sensitivity", v$gate)
show("Alignment sensitivity", v$alignment)
if (!is.null(v$ldir_join)) show("LDIR join check", v$ldir_join$summary)
cat("\nCSV files in", out, "\n")
