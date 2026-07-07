#!/usr/bin/env Rscript
# BES base application re-run (full data) + strong/weak figures/tables.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE); library(stm); library(Matrix); library(ggplot2)
})
dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)
dir.create("output/tables",  showWarnings=FALSE, recursive=TRUE)
set.seed(2026)
d <- readRDS("replication/data/bes_w25_dtm.rds")
W <- as.matrix(d$W); C <- d$C; vocab <- d$vocab
leave <- ifelse(C[,"leave"] > 0, 1L, 0L)
K <- 6L; lam <- 3

# ================= base sgscatm fit on ALL docs =================
t0 <- proc.time()[3]; fit <- sgscatm(W, C, K=K, lambda=lam, rotate=TRUE)
t_fit <- proc.time()[3]-t0
t0 <- proc.time()[3]; vc <- sgscatm_vcov(fit, identified=TRUE); t_se <- proc.time()[3]-t0
cat(sprintf("Full BES: %d docs x %d terms. sgscatm fit %.3fs + SE %.3fs\n",
            nrow(W), ncol(W), t_fit, t_se))

# STM on full data (timing at scale)
docs <- lapply(seq_len(nrow(W)), function(i){ idx<-which(W[i,]>0)
  rbind(as.integer(idx), as.integer(W[i,idx])) })
ok <- vapply(docs, function(x) ncol(x)>0, logical(1))
meta <- data.frame(leave=leave, age=C[,"age_std"], female=C[,"female"],
                   educ=C[,"educ_std"], lr=C[,"lr_std"])
t0 <- proc.time()[3]
fstm <- stm(docs[ok], vocab, K=K, prevalence=~leave+age+female+educ+lr,
            data=meta[ok,], init.type="Spectral", max.em.its=75, verbose=FALSE)
t_stm <- proc.time()[3]-t0
cat(sprintf("STM full fit %.2fs -> speedup %.0fx\n", t_stm, t_stm/(t_fit+t_se)))

# covariate significance (joint Wald per covariate)
Bz <- vc$B; SE <- vc$se; P <- nrow(Bz); Km1 <- ncol(Bz); vcm <- vc$vcov
rownames(Bz) <- colnames(C)
wald <- data.frame(covariate=colnames(C), stat=NA, df=NA, p=NA)
for (p in seq_len(P)) {
  idx <- (seq_len(Km1)-1L)*P + p; bp <- as.vector(Bz[p,]); Sp <- vcm[idx,idx,drop=FALSE]
  ei <- eigen(Sp, symmetric=TRUE); pos <- ei$values > 1e-10*max(ei$values)
  Si <- ei$vectors[,pos,drop=FALSE] %*% diag(1/ei$values[pos], sum(pos)) %*% t(ei$vectors[,pos,drop=FALSE])
  wald$stat[p] <- as.numeric(t(bp)%*%Si%*%bp); wald$df[p] <- sum(pos)
  wald$p[p] <- pchisq(wald$stat[p], sum(pos), lower.tail=FALSE)
}
cat("\n=== covariate joint Wald (full BES) ===\n"); print(wald, row.names=FALSE, digits=3)

# leave effect on topic prevalence + top terms
dprev <- colMeans(fit$Pi[leave==1,]) - colMeans(fit$Pi[leave==0,])
tt <- top_terms(fit, n=8L, vocab=vocab)
cat("\n=== topics: leave (Leave-Remain) prevalence effect + top terms ===\n")
ord <- order(dprev)
for (k in ord) cat(sprintf("  T%d  dprev=%+.4f  | %s\n", k, dprev[k],
                           paste(tt[k,], collapse=" ")))
saveRDS(list(fit=fit, vc=vc, wald=wald, dprev=dprev, tt=tt,
             t_fit=t_fit, t_se=t_se, t_stm=t_stm), "output/bes_base.rds")

# ================= figures =================
df <- readRDS("output/bes_strong_weak.rds")
long <- rbind(
  data.frame(scenario=df$scenario, method="sgscatm", cor=df$cor_sg, t=df$t_sg),
  data.frame(scenario=df$scenario, method="STM",     cor=df$cor_stm, t=df$t_stm))
long$scenario <- factor(long$scenario, levels=c("strong","weak"),
                        labels=c("strong (|lr|>=1)","weak (|lr|<=0.3)"))
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank()))

pA <- ggplot(long, aes(scenario, cor, fill=method)) +
  geom_boxplot(outlier.size=0.7, width=0.6, position=position_dodge(0.7)) +
  scale_fill_manual(values=c("sgscatm"="#534AB7","STM"="#B7743A")) +
  labs(x=NULL, y="recovery corr (implied vs observed Leave-Remain word shift)",
       title="BES: covariate-effect recovery, sgscatm vs STM", fill=NULL) +
  ylim(0,1)
ggsave("output/figures/bes_recovery_strong_weak.pdf", pA, width=6.5, height=4)

pB <- ggplot(long, aes(scenario, t, fill=method)) +
  geom_boxplot(outlier.size=0.7, width=0.6, position=position_dodge(0.7)) +
  scale_y_log10() + scale_fill_manual(values=c("sgscatm"="#534AB7","STM"="#B7743A")) +
  labs(x=NULL, y="wall-clock (s, log10)", title="BES: runtime per fit", fill=NULL)
ggsave("output/figures/bes_timing_strong_weak.pdf", pB, width=6, height=4)

# summary table (LaTeX)
fmt <- function(m,s) sprintf("%.3f (%.3f)", m, s)
mk <- function(sc){ s <- df[df$scenario==sc,]
  c(fmt(mean(s$cor_sg),sd(s$cor_sg)), fmt(mean(s$cor_stm),sd(s$cor_stm)),
    sprintf("%.1f$\\times$", mean(s$t_stm)/mean(s$t_sg)),
    formatC(t.test(s$cor_sg,s$cor_stm,paired=TRUE)$p.value, format="g", digits=2)) }
st <- mk("strong"); wk <- mk("weak")
lines <- c("\\begin{tabular}{lcccc}","\\toprule",
  "scenario & sgscatm corr & STM corr & speedup & paired $p$ \\\\","\\midrule",
  sprintf("strong (|lr|$\\geq$1) & %s & %s & %s & %s \\\\", st[1],st[2],st[3],st[4]),
  sprintf("weak (|lr|$\\leq$0.3) & %s & %s & %s & %s \\\\", wk[1],wk[2],wk[3],wk[4]),
  "\\bottomrule","\\end{tabular}")
writeLines(lines, "output/tables/bes_strong_weak.tex")
cat("\nFigures + table written.\n")
