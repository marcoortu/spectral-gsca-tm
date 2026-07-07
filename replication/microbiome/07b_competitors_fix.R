#!/usr/bin/env Rscript
# Fix ALDEx2 (two-group mode; glm mode hits an ALDEx2/R4 class bug) and STM
# (reindex vocab), and recompute G5 concordance with directional + SNR views.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
set.seed(2026)
D  <- readRDS("output/phase2_data.rds")
F2 <- readRDS("output/phase2_fit.rds")
CP <- readRDS("output/phase2_competitors.rds")
Xg <- as.matrix(D$Xg10); C <- D$C; M <- nrow(Xg); N <- ncol(Xg)
covdf <- as.data.frame(C)
timings <- CP$timings

# ---- ALDEx2 two-group (CRC vs control) ----
aldex_p <- NULL; aldex_res <- NULL
tryCatch({
  suppressPackageStartupMessages(library(ALDEx2))
  counts <- round(sweep(Xg,1,rowSums(Xg),"/")*1e6); reads <- t(counts)
  storage.mode(reads) <- "integer"
  conds <- ifelse(C[,"study_condition"]==1, "CRC","CTR")
  t0 <- proc.time()[3]
  ax <- aldex(reads, conds, mc.samples=128, test="t", effect=TRUE,
              denom="all", verbose=FALSE)
  timings$aldex2 <- proc.time()[3] - t0
  aldex_res <- ax
  aldex_p <- setNames(ax$wi.eBH, rownames(ax))     # BH-adjusted Wilcoxon
  aldex_sig <- rownames(ax)[ax$wi.eBH < 0.1]
  cat(sprintf("ALDEx2 two-group done (%.2fs); %d genera BH<0.1\n",
              timings$aldex2, length(aldex_sig)))
  print(head(rownames(ax)[order(ax$wi.eBH)], 15))
}, error=function(e) cat("ALDEx2 FAILED:", conditionMessage(e), "\n"))

# ---- STM (reindexed vocab) ----
tryCatch({
  suppressPackageStartupMessages(library(stm))
  counts <- round(sweep(Xg,1,rowSums(Xg),"/")*1e4)
  keep <- colSums(counts) > 0
  counts <- counts[, keep, drop=FALSE]; vocab <- colnames(Xg)[keep]
  docs <- lapply(seq_len(M), function(i){ idx <- which(counts[i,]>0)
    rbind(as.integer(idx), as.integer(counts[i,idx])) })
  ok <- vapply(docs, function(d) ncol(d)>0, logical(1))
  t0 <- proc.time()[3]
  stm_fit <- stm(documents=docs[ok], vocab=vocab, K=F2$Kstar,
                 prevalence=~ study_condition+age+bmi+gender,
                 data=covdf[ok,], max.em.its=75, init.type="Spectral", verbose=FALSE)
  timings$stm <- proc.time()[3] - t0
  cat(sprintf("STM done (%.2fs)\n", timings$stm))
}, error=function(e) { cat("STM FAILED:", conditionMessage(e), "\n") })

# ---- G5b recomputed: directional concordance + SNR ranking ----
fit <- F2$fit
d_dir <- fit$Bz["study_condition", ]; d_dir <- d_dir/sqrt(sum(d_dir^2))
loadings <- as.numeric(t(fit$Psi) %*% d_dir); names(loadings) <- colnames(Xg)
# SNR: bootstrap the loading to get its SE, rank by |loading|/SE
Bpsi <- 200L
lb <- matrix(NA_real_, Bpsi, N)
for (b in seq_len(Bpsi)) {
  idx <- sample.int(M, M, replace=TRUE)
  fb <- tryCatch(sgscatm(Xg[idx,], C[idx,], K=F2$Kstar, lambda=F2$lamstar, rotate=TRUE),
                 error=function(e) NULL)
  if (!is.null(fb)) {
    dd <- fb$Bz["study_condition",]; dd <- dd/sqrt(sum(dd^2))
    # sign-align to reference direction via Psi correlation
    lo <- as.numeric(t(fb$Psi) %*% dd)
    if (sum(lo*loadings) < 0) lo <- -lo
    lb[b,] <- lo
  }
}
lo_sd <- apply(lb, 2, sd, na.rm=TRUE)
snr   <- loadings / lo_sd
known <- D$known_crc; known_in <- intersect(known, colnames(Xg))

# directional: do known CRC genera load positive (CRC-ward) more than chance?
sign_known <- sign(loadings[known_in])
n_pos <- sum(sign_known > 0)
p_sign <- binom.test(n_pos, length(known_in), 0.5, alternative="greater")$p.value

# top-k by SNR overlap with known (hypergeometric)
for (topn in c(15L, 25L, 40L)) {
  top_snr <- names(sort(abs(snr), decreasing=TRUE))[seq_len(topn)]
  hit <- length(intersect(top_snr, known_in))
  ph  <- phyper(hit-1L, length(known_in), N-length(known_in), topn, lower.tail=FALSE)
  cat(sprintf("SNR top-%d: known hits=%d hypergeom p=%.4g | overlap: %s\n",
              topn, hit, ph, paste(intersect(top_snr, known_in), collapse=", ")))
}
cat(sprintf("\nDirectional: %d/%d known CRC genera load CRC-ward (+); sign-test p=%.4g\n",
            n_pos, length(known_in), p_sign))
cat("known genera loadings & SNR:\n")
print(data.frame(genus=known_in, loading=round(loadings[known_in],3),
                 snr=round(snr[known_in],2)), row.names=FALSE)

# concordance vs ALDEx2 (rank agreement on shared genera)
if (!is.null(aldex_p)) {
  shared <- intersect(names(aldex_p), colnames(Xg))
  rho <- suppressWarnings(cor(abs(snr[shared]), -log10(aldex_p[shared]+1e-8),
                              method="spearman", use="complete.obs"))
  cat(sprintf("\nSpearman(|SNR|, -log10 ALDEx2 p) over %d genera = %.3f\n", length(shared), rho))
  aldex_top <- names(sort(aldex_p))[seq_len(min(25L,length(aldex_p)))]
  snr_top   <- names(sort(abs(snr), decreasing=TRUE))[seq_len(25L)]
  cat(sprintf("sgscatm-SNR-top25 vs ALDEx2-top25 overlap = %d\n",
              length(intersect(aldex_top, snr_top))))
}

saveRDS(list(timings=timings, aldex_p=aldex_p, aldex_res=aldex_res,
             loadings=loadings, snr=snr, lo_sd=lo_sd,
             n_pos=n_pos, p_sign=p_sign, known_in=known_in),
        "output/phase2_competitors_fix.rds")
cat("\n=== G6 timings (s) ===\n"); print(unlist(timings))
cat("Saved output/phase2_competitors_fix.rds\n")
