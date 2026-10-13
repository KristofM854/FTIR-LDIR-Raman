# =============================================================================
# run_example_A.R -- run the full pipeline on the bundled example dataset
# =============================================================================
# One filter measured with all four instruments (see examples/dataset_A/README.md).
# From the repository root:
#
#   Rscript examples/run_example_A.R            # results in output/<date>_<n>/
#   Rscript examples/run_example_A.R /tmp/out   # results under /tmp/out/
#
# Needs the R packages in DESCRIPTION and Python with
# inst/python/requirements.txt (LDIR image analysis). Takes about 1-2 minutes.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
D <- file.path("examples", "dataset_A")
if (!dir.exists(D)) stop("Run from the repository root (examples/dataset_A not found).")

input_mode        <- "hardcoded"
output_dir        <- if (length(args) >= 1) args[1] else "output"
ftir_file         <- file.path(D, "FTIR_perkin.csv")
ftir_image        <- file.path(D, "FTIR_image.png")
ftir_bruker_file  <- file.path(D, "FTIR_bruker.csv")
raman_file        <- file.path(D, "Raman_data.csv")
raman_image       <- file.path(D, "Raman_image_small.png")
ldir_file         <- file.path(D, "LDIR_data.xlsx")
ldir_image        <- file.path(D, "LDIR_image.png")   # LDIR_image_analyzed.png and
                                                      # LDIR_coord_swaps.csv are found
                                                      # next to it automatically
source("main.R")
