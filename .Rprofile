# Activate the project's renv library when it has been initialised.
#
# First-time setup (once, on your machine — it needs the packages installed to
# capture their exact versions, which cannot be done in a headless CI image):
#   install.packages("renv")
#   renv::init(bare = TRUE)                 # generates renv/activate.R
#   renv::snapshot(type = "explicit")       # writes renv.lock from DESCRIPTION
#   git add renv.lock renv/activate.R .Rprofile && git commit
#
# Thereafter (and for anyone re-running the deposit):
#   renv::restore()                         # installs the exact locked versions
#
# Until renv/activate.R exists this is a harmless no-op, so the pipeline still
# runs against the system library.
if (file.exists("renv/activate.R")) source("renv/activate.R")
