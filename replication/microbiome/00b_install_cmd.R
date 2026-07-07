#!/usr/bin/env Rscript
# Retry curatedMetagenomicData install alone, with full error surfacing.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
log <- function(...) cat(sprintf("[cmd] %s\n", sprintf(...)))
if (requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
  log("already installed"); quit(save = "no")
}
tryCatch(
  BiocManager::install("curatedMetagenomicData", update = FALSE, ask = FALSE),
  error = function(e) log("ERROR: %s", conditionMessage(e))
)
log("curatedMetagenomicData available: %s",
    requireNamespace("curatedMetagenomicData", quietly = TRUE))
