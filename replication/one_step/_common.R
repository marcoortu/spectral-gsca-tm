## Common setup for all one_step scripts. Source this first.
## Run from the package root, e.g.
##   Rscript replication/one_step/00_noiseless_check.R
suppressMessages({
  # locate package root relative to this file's known location
  pkg_root <- normalizePath(file.path(getwd()))
  # if launched from elsewhere, allow override
  if (!file.exists(file.path(pkg_root, "DESCRIPTION")))
    pkg_root <- normalizePath("c:/Users/Utente/WORKSPACE/R/sgscatm")
  library(devtools)
  devtools::load_all(pkg_root, quiet = TRUE)
  source(file.path(pkg_root, "replication/simulation/sim_dgp.R"))
  source(file.path(pkg_root, "replication/simulation/sim_utils.R"))
  source(file.path(pkg_root, "replication/one_step/estimators.R"))
})
set.seed(1)
ONE_STEP_DIR <- file.path(pkg_root, "replication/one_step")
FIG_DIR      <- file.path(ONE_STEP_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

norm_F <- function(X) sqrt(mean(X^2))
