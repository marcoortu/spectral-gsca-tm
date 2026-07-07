#!/usr/bin/env Rscript
# Install SIAMCAT (ships the Zeller 2014 CRC metagenomic cohort as example data)
# as a robust real-data source, given the curatedMetagenomicData/rbiom conflict.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
log <- function(...) cat(sprintf("[siamcat] %s\n", sprintf(...)))
if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  log("installing SIAMCAT ...")
  tryCatch(BiocManager::install("SIAMCAT", update = FALSE, ask = FALSE),
           error = function(e) log("ERROR: %s", conditionMessage(e)))
} else log("already installed")
log("SIAMCAT available: %s", requireNamespace("SIAMCAT", quietly = TRUE))
if (requireNamespace("SIAMCAT", quietly = TRUE)) {
  suppressPackageStartupMessages(library(SIAMCAT))
  data("feat.crc.zeller", package = "SIAMCAT")
  data("meta.crc.zeller", package = "SIAMCAT")
  log("feat dim: %d x %d", nrow(feat.crc.zeller), ncol(feat.crc.zeller))
  log("meta dim: %d x %d", nrow(meta.crc.zeller), ncol(meta.crc.zeller))
  log("meta cols: %s", paste(colnames(meta.crc.zeller), collapse = ", "))
}
