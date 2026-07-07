#!/usr/bin/env Rscript
# ===================================================================
# BES Wave 25 (MII open-text) — sgscatm vs STM, strong vs weak scenario.
#
# Question: does sgscatm recover the covariate->content effect better/faster
# than STM when documents are selected into a STRONG-signal regime (politically
# committed respondents, |lr_std|>=1, where Brexit stance strongly structures
# content) vs a WEAK-signal regime (centrists, |lr_std|<=0.3)?
#
# Ground truth is MODEL-AGNOSTIC: the observed Leave-Remain word-frequency shift
# on a large held-out portion of each pool. Each method's recovery = correlation
# of its IMPLIED covariate->word shift (Delta topic-prevalence x topic-word) with
# that observed shift. Over R resamples we report recovery, timing, and stability.
# ===================================================================
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE); library(stm); library(Matrix)
})
set.seed(2026)
d <- readRDS("replication/data/bes_w25_dtm.rds")
W <- as.matrix(d$W); C <- d$C; vocab <- d$vocab
leave <- ifelse(C[,"leave"] > 0, 1L, 0L)
lr <- C[,"lr_std"]
K <- 6L; lam <- 3; n_tr <- 1500L; R <- 25L

pools <- list(strong = which(abs(lr) >= 1.0), weak = which(abs(lr) <= 0.3))
cat(sprintf("strong pool: %d docs (Leave %.1f%%) | weak pool: %d docs (Leave %.1f%%)\n",
    length(pools$strong), 100*mean(leave[pools$strong]),
    length(pools$weak),   100*mean(leave[pools$weak])))

delta_word <- function(Wm, lv) {
  fr <- Wm / pmax(rowSums(Wm), 1)
  colMeans(fr[lv==1,,drop=FALSE]) - colMeans(fr[lv==0,,drop=FALSE])
}

one_rep <- function(pool, seed) {
  set.seed(seed); pool <- sample(pool)
  tr <- pool[seq_len(n_tr)]; ho <- pool[-seq_len(n_tr)]
  Wtr <- W[tr,]; Ctr <- C[tr,]; lvtr <- leave[tr]
  keep <- colSums(Wtr) > 0
  Wtr <- Wtr[,keep]; voc <- vocab[keep]
  dw_obs <- delta_word(W[ho, keep, drop=FALSE], leave[ho])   # model-agnostic truth

  # --- sgscatm (fit + closed-form SE) ---
  t0 <- proc.time()[3]
  fs <- sgscatm(Wtr, Ctr, K=K, lambda=lam, rotate=TRUE)
  t_sg <- proc.time()[3]-t0
  t0 <- proc.time()[3]; vc <- sgscatm_vcov(fs, identified=TRUE); t_se <- proc.time()[3]-t0
  dprev_sg <- colMeans(fs$Pi[lvtr==1,]) - colMeans(fs$Pi[lvtr==0,])
  dw_sg <- as.numeric(dprev_sg %*% fs$Phi)

  # --- STM (matched covariates) ---
  docs <- lapply(seq_len(nrow(Wtr)), function(i){ idx<-which(Wtr[i,]>0)
    rbind(as.integer(idx), as.integer(Wtr[i,idx])) })
  ok <- vapply(docs, function(x) ncol(x)>0, logical(1))
  metatr <- data.frame(leave=lvtr, age=Ctr[,"age_std"], female=Ctr[,"female"],
                       educ=Ctr[,"educ_std"], lr=Ctr[,"lr_std"])
  t0 <- proc.time()[3]
  fstm <- stm(docs[ok], voc, K=K, prevalence=~leave+age+female+educ+lr,
              data=metatr[ok,], init.type="Spectral", max.em.its=50, verbose=FALSE)
  t_stm <- proc.time()[3]-t0
  beta <- exp(fstm$beta$logbeta[[1]])
  dth <- colMeans(fstm$theta[metatr$leave[ok]==1,]) -
         colMeans(fstm$theta[metatr$leave[ok]==0,])
  dw_stm <- as.numeric(dth %*% beta)

  data.frame(
    cor_sg  = cor(dw_sg, dw_obs),   cor_stm  = cor(dw_stm, dw_obs),
    rmse_sg = sqrt(mean((dw_sg-dw_obs)^2)), rmse_stm = sqrt(mean((dw_stm-dw_obs)^2)),
    t_sg = t_sg + t_se, t_stm = t_stm)
}

res <- list()
for (sc in names(pools)) {
  rows <- lapply(seq_len(R), function(r) {
    out <- tryCatch(one_rep(pools[[sc]], 7000L + r), error=function(e) NULL)
    if (!is.null(out)) out$scenario <- sc; out
  })
  res[[sc]] <- do.call(rbind, rows)
  cat(sprintf("[%s] done %d reps\n", sc, nrow(res[[sc]])))
}
df <- do.call(rbind, res)
saveRDS(df, "output/bes_strong_weak.rds")

# ---- summary ----
agg <- function(x) sprintf("%.3f (%.3f)", mean(x), sd(x))
cat("\n=================== SUMMARY (mean (sd) over reps) ===================\n")
for (sc in names(pools)) {
  s <- df[df$scenario==sc,]
  cat(sprintf("\n--- %s scenario ---\n", toupper(sc)))
  cat(sprintf("  recovery corr   sgscatm %s   STM %s\n", agg(s$cor_sg), agg(s$cor_stm)))
  cat(sprintf("  recovery RMSE   sgscatm %s   STM %s\n", agg(s$rmse_sg), agg(s$rmse_stm)))
  cat(sprintf("  wall-clock (s)  sgscatm %s   STM %s   speedup %.1fx\n",
              agg(s$t_sg), agg(s$t_stm), mean(s$t_stm)/mean(s$t_sg)))
  pt <- t.test(s$cor_sg, s$cor_stm, paired=TRUE)
  cat(sprintf("  paired t (corr sgscatm-STM): diff=%.3f  p=%.2g\n",
              mean(s$cor_sg-s$cor_stm), pt$p.value))
}
cat("\nSaved output/bes_strong_weak.rds\n")
