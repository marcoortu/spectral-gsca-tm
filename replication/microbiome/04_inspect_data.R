#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SIAMCAT))
data("feat.crc.zeller", package = "SIAMCAT")
data("meta.crc.zeller", package = "SIAMCAT")
feat <- feat.crc.zeller; meta <- meta.crc.zeller
cat(sprintf("feat: %d x %d (features x samples)\n", nrow(feat), ncol(feat)))
cat("first 8 feature rownames:\n"); print(head(rownames(feat), 8))
cat("last 4 feature rownames:\n"); print(tail(rownames(feat), 4))
cat("\ncolsums (should be ~1 if rel abund):\n"); print(summary(colSums(feat)))
cat("\nmeta head:\n"); print(head(meta))
cat("\nGroup table:\n"); print(table(meta$Group))
cat("Gender table:\n"); print(table(meta$Gender))
cat(sprintf("Age missing: %d, BMI missing: %d, Gender missing: %d\n",
            sum(is.na(meta$Age)), sum(is.na(meta$BMI)), sum(is.na(meta$Gender))))
cat("\nAny known CRC taxa present?\n")
known <- c("Fusobacterium","Peptostreptococcus","Parvimonas","Gemella",
           "Porphyromonas","Solobacterium")
for (g in known) {
  hits <- grep(g, rownames(feat), value = TRUE, ignore.case = TRUE)
  cat(sprintf("  %-18s : %d matches; e.g. %s\n", g, length(hits),
              paste(head(hits,2), collapse=" | ")))
}
