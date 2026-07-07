#!/usr/bin/env Rscript
# ===================================================================
# Phase 2 competitors + concordance. Gates: G6 (speed), G5 (concordance).
# Competitors on the SAME genus table + covariates: PERMANOVA, ALDEx2, STM,
# ANCOM-BC. All timed. Then map sgscatm disease direction to genus loadings
# and test overlap with known CRC taxa and with ALDEx2/PERMANOVA.
# ===================================================================
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
set.seed(2026)
D  <- readRDS("output/phase2_data.rds")
F2 <- readRDS("output/phase2_fit.rds")
Xg <- as.matrix(D$Xg10); C <- D$C; M <- nrow(Xg); N <- ncol(Xg)
covdf <- data.frame(study_condition = C[,"study_condition"],
                    age = C[,"age"], bmi = C[,"bmi"], gender = C[,"gender"])
timings <- list()

# ---- sgscatm timing (from fit script) ----
timings$sgscatm <- F2$t_fit + F2$t_se

# ---- PERMANOVA (vegan::adonis2 on Aitchison distance) ----
suppressPackageStartupMessages(library(vegan))
Xp  <- Xg + 0.5 / N; Xp <- Xp / rowSums(Xp)
clr <- log(Xp) - rowMeans(log(Xp))
Dait <- dist(clr, method = "euclidean")           # Aitchison distance
t0 <- proc.time()[3]
perm <- adonis2(Dait ~ study_condition + age + bmi + gender,
                data = covdf, by = "margin", permutations = 999)
timings$permanova <- proc.time()[3] - t0
cat("=== PERMANOVA (adonis2, margin) ===\n"); print(perm)
perm_p <- setNames(perm$`Pr(>F)`[seq_len(4)], rownames(perm)[seq_len(4)])

# ---- ALDEx2 (aldex.glm, 128 MC) ----
aldex_p <- NULL; aldex_sig <- character(0)
tryCatch({
  suppressPackageStartupMessages(library(ALDEx2))
  depth  <- 1e6
  counts <- round(sweep(Xg, 1, rowSums(Xg), "/") * depth)   # pseudo-counts
  reads  <- t(counts)                                        # features x samples
  storage.mode(reads) <- "integer"
  mm <- model.matrix(~ study_condition + age + bmi + gender, data = covdf)
  t0 <- proc.time()[3]
  clrx <- aldex.clr(reads, mm, mc.samples = 128, denom = "all", verbose = FALSE)
  glm  <- aldex.glm(clrx, mm, verbose = FALSE)
  timings$aldex2 <- proc.time()[3] - t0
  pcol <- grep("study_condition.*pval$", colnames(glm), value = TRUE)[1]
  qcol <- grep("study_condition.*(BH|holm)$", colnames(glm), value = TRUE)[1]
  aldex_p   <- setNames(glm[[pcol]], rownames(glm))
  aldex_sig <- rownames(glm)[glm[[qcol]] < 0.1]
  cat(sprintf("\nALDEx2: %d genera with BH<0.1 for study_condition\n", length(aldex_sig)))
  print(head(aldex_sig, 15))
}, error = function(e) { cat("ALDEx2 FAILED:", conditionMessage(e), "\n") })

# ---- STM (variational EM) ----
tryCatch({
  suppressPackageStartupMessages(library(stm))
  depth  <- 1e4
  counts <- round(sweep(Xg, 1, rowSums(Xg), "/") * depth)
  docs <- lapply(seq_len(M), function(i) {
    idx <- which(counts[i, ] > 0)
    rbind(as.integer(idx), as.integer(counts[i, idx]))
  })
  vocab <- colnames(Xg)
  t0 <- proc.time()[3]
  stm_fit <- stm(documents = docs, vocab = vocab, K = F2$Kstar,
                 prevalence = ~ study_condition + age + bmi + gender,
                 data = covdf, max.em.its = 75, init.type = "Spectral",
                 verbose = FALSE)
  ep <- estimateEffect(1:F2$Kstar ~ study_condition + age + bmi + gender,
                       stm_fit, metadata = covdf, uncertainty = "Global")
  timings$stm <- proc.time()[3] - t0
  cat(sprintf("\nSTM fit + estimateEffect done (%.2fs)\n", timings$stm))
}, error = function(e) { cat("STM FAILED:", conditionMessage(e), "\n"); timings$stm <<- NA })

# ---- ANCOM-BC2 ----
tryCatch({
  suppressPackageStartupMessages({ library(ANCOMBC); library(TreeSummarizedExperiment) })
  depth  <- 1e6
  counts <- round(sweep(Xg, 1, rowSums(Xg), "/") * depth)
  tse <- TreeSummarizedExperiment(
    assays = list(counts = t(counts)),
    colData = DataFrame(covdf))
  t0 <- proc.time()[3]
  anc <- ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "study_condition + age + bmi + gender",
                  p_adj_method = "BH", prv_cut = 0, verbose = FALSE)
  timings$ancombc <- proc.time()[3] - t0
  rr <- anc$res
  qc <- grep("q_study_condition", colnames(rr), value = TRUE)[1]
  ancombc_sig <- rr$taxon[rr[[qc]] < 0.1]
  cat(sprintf("\nANCOM-BC2: %d genera q<0.1 for study_condition (%.2fs)\n",
              length(ancombc_sig), timings$ancombc))
}, error = function(e) { cat("ANCOMBC FAILED:", conditionMessage(e), "\n"); timings$ancombc <<- NA })

# ================= G5 concordance =================
# disease direction in the fit's (varimax) score basis, project Psi -> genera
fit <- F2$fit
d_dir  <- fit$Bz["study_condition", ]                 # (K-1) vector
d_dir  <- d_dir / sqrt(sum(d_dir^2))
loadings <- as.numeric(t(fit$Psi) %*% d_dir)          # N genus loadings
names(loadings) <- colnames(Xg)
ord <- order(abs(loadings), decreasing = TRUE)
topn <- 25L
top_genera <- names(loadings)[ord[seq_len(topn)]]
cat(sprintf("\n=== sgscatm top-%d disease-direction genera ===\n", topn))
print(data.frame(genus = top_genera, loading = round(loadings[ord[seq_len(topn)]],3)),
      row.names = FALSE)

known <- D$known_crc
# hypergeometric: overlap of top-loading set with known CRC genera among all N
known_in_pool <- intersect(known, colnames(Xg))
hit <- length(intersect(top_genera, known_in_pool))
p_hyper <- phyper(hit - 1L, length(known_in_pool), N - length(known_in_pool),
                  topn, lower.tail = FALSE)
cat(sprintf("\nG5b: known CRC genera in pool = %d; in sgscatm top-%d = %d; hypergeom p = %.4g\n",
            length(known_in_pool), topn, hit, p_hyper))
cat("  overlap:", paste(intersect(top_genera, known_in_pool), collapse=", "), "\n")

# concordance sgscatm-top vs ALDEx2-significant
if (!is.null(aldex_p)) {
  aldex_top <- names(sort(aldex_p))[seq_len(min(topn, length(aldex_p)))]
  ov <- length(intersect(top_genera, aldex_top))
  cat(sprintf("\nG5: sgscatm-top vs ALDEx2-top-%d overlap = %d (Jaccard %.2f)\n",
              topn, ov, ov/(2*topn-ov)))
}

# G5a covariate-significance agreement
cat("\n=== G5a covariate significance ===\n")
comp <- data.frame(covariate = c("study_condition","age","bmi","gender"),
                   sgscatm_wald_p = F2$wald$p,
                   permanova_p = perm_p[c("study_condition","age","bmi","gender")])
print(comp, row.names = FALSE, digits = 3)

saveRDS(list(timings = timings, perm = perm, perm_p = perm_p,
             aldex_p = aldex_p, aldex_sig = aldex_sig,
             loadings = loadings, top_genera = top_genera,
             p_hyper = p_hyper, hit = hit, known_in_pool = known_in_pool,
             comp = comp),
        "output/phase2_competitors.rds")
cat("\n=== G6 timings (seconds) ===\n")
print(unlist(timings))
cat("\nSaved output/phase2_competitors.rds\n")
