#!/usr/bin/env Rscript
# ===================================================================
# Phase 2 prep — Zeller 2014 CRC metagenomic cohort (via SIAMCAT).
# Builds genus- and species-level compositions + covariate matrix.
#
# NOTE (logged deviation): the pre-registration named curatedMetagenomicData
# for a pooled multi-cohort set. That package fails to install on this system
# (Bioc 3.22 rbiom API break: 'unifrac' no longer exported, via the mia load
# chain). We substitute the Zeller 2014 French CRC cohort bundled in SIAMCAT —
# one of the same cohorts cMD pools and the canonical dataset establishing
# Fusobacterium/Peptostreptococcus CRC enrichment. Single-cohort, so no
# study/country adjustment is needed. This is disclosed, not silent.
# ===================================================================
suppressPackageStartupMessages(library(SIAMCAT))
data("feat.crc.zeller", package = "SIAMCAT")
data("meta.crc.zeller", package = "SIAMCAT")
feat <- as.matrix(feat.crc.zeller)     # 1754 species x 141 samples (rel abund)
meta <- meta.crc.zeller

# ---- parse genus from species label ----
species_to_genus <- function(x) {
  x <- sub("\\s*\\[[a-z]:[0-9]+\\]\\s*$", "", x)  # strip [h:NNN]/[u:NNN]
  x <- sub("^unnamed\\s+", "", x)                  # strip 'unnamed '
  vapply(strsplit(x, "\\s+"), `[`, character(1), 1) # first token = genus
}

# drop UNMAPPED / unclassified rows
keep <- !grepl("^UNMAPPED$|unclassified", rownames(feat), ignore.case = TRUE)
feat <- feat[keep, , drop = FALSE]
genus <- species_to_genus(rownames(feat))

# ---- samples x features ----
Xsp <- t(feat)                            # 141 x N_species (rel abund)
# aggregate species -> genus (sum relative abundances)
Xg  <- t(rowsum(t(Xsp), group = genus))   # 141 x N_genus
cat(sprintf("species: %d, genera: %d, samples: %d\n",
            ncol(Xsp), ncol(Xg), nrow(Xg)))

# ---- covariates: complete cases on Age, BMI, Gender, Group ----
meta$BMI <- suppressWarnings(as.numeric(as.character(meta$BMI)))
meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
cc <- complete.cases(meta[, c("Age","BMI","Gender","Group")])
cat(sprintf("complete-case samples: %d (dropped %d for missing BMI)\n",
            sum(cc), sum(!cc)))
meta <- meta[cc, ]; Xsp <- Xsp[cc, ]; Xg <- Xg[cc, ]

C <- cbind(
  study_condition = as.integer(meta$Group == "CRC"),   # 1 = CRC
  age             = as.numeric(scale(meta$Age)),
  bmi             = as.numeric(scale(meta$BMI)),
  gender          = as.integer(meta$Gender == "M")      # 1 = male
)
rownames(C) <- rownames(meta)

# ---- prevalence filter helper + renormalize ----
prev_filter <- function(X, min_prev = 0.10) {
  prev <- colMeans(X > 0)
  X <- X[, prev >= min_prev, drop = FALSE]
  X / rowSums(X)                          # renormalize to compositions
}

Xg10  <- prev_filter(Xg,  0.10)
Xg20  <- prev_filter(Xg,  0.20)
Xg30  <- prev_filter(Xg,  0.30)
Xsp10 <- prev_filter(Xsp, 0.10)
Xsp20 <- prev_filter(Xsp, 0.20)

cat(sprintf("genus  N at prev>=10/20/30%%: %d / %d / %d\n",
            ncol(Xg10), ncol(Xg20), ncol(Xg30)))
cat(sprintf("species N at prev>=10/20%%: %d / %d\n", ncol(Xsp10), ncol(Xsp20)))

# sparsity + library-size (from original rel-abund; report nonzero fraction)
sparsity_g  <- mean(Xg10  == 0)
sparsity_sp <- mean(Xsp10 == 0)
cat(sprintf("sparsity (frac zeros) genus10=%.3f  species10=%.3f\n",
            sparsity_g, sparsity_sp))

# known CRC-enriched genera present at genus10
known <- c("Fusobacterium","Peptostreptococcus","Parvimonas","Gemella",
           "Porphyromonas","Solobacterium")
cat("known CRC genera present at genus>=10%:\n")
print(intersect(known, colnames(Xg10)))

saveRDS(list(Xg10 = Xg10, Xg20 = Xg20, Xg30 = Xg30,
             Xsp10 = Xsp10, Xsp20 = Xsp20,
             Xg_full = Xg, Xsp_full = Xsp,
             C = C, meta = meta, genus_map = genus,
             known_crc = known),
        "output/phase2_data.rds")
cat("Saved output/phase2_data.rds\n")
