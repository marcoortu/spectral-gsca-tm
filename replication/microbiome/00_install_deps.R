#!/usr/bin/env Rscript
# ===================================================================
# Install dependencies for the microbiome verification run.
# Idempotent: skips anything already installed.
# ===================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)  # long downloads

log <- function(...) cat(sprintf("[install] %s\n", sprintf(...)))

ensure_cran <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      log("installing CRAN: %s", p)
      install.packages(p, quiet = TRUE)
    } else log("already have CRAN: %s", p)
  }
}

ensure_cran(c("BiocManager", "vegan", "Matrix", "ggplot2", "dplyr", "tidyr"))

bioc <- c("curatedMetagenomicData", "ALDEx2", "ANCOMBC",
          "TreeSummarizedExperiment", "SummarizedExperiment", "mia")
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    log("installing Bioc: %s", p)
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE, quiet = TRUE),
      error = function(e) log("FAILED Bioc %s: %s", p, conditionMessage(e))
    )
  } else log("already have Bioc: %s", p)
}

log("=== availability after install ===")
allpkgs <- c("BiocManager", "vegan", "curatedMetagenomicData", "ALDEx2",
             "ANCOMBC", "TreeSummarizedExperiment", "mia", "stm")
for (p in allpkgs) log("%-28s %s", p, requireNamespace(p, quietly = TRUE))
log("DONE")
