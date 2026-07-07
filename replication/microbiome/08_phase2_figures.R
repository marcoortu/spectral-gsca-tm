#!/usr/bin/env Rscript
# Phase 2 figures + LaTeX tables.
suppressPackageStartupMessages({ devtools::load_all(".", quiet=TRUE); library(ggplot2) })
dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)
dir.create("output/tables",  showWarnings=FALSE, recursive=TRUE)
D  <- readRDS("output/phase2_data.rds")
F2 <- readRDS("output/phase2_fit.rds")
CP <- readRDS("output/phase2_competitors.rds")
CF <- readRDS("output/phase2_competitors_fix.rds")
fmt <- function(x,d=3) formatC(x, format="f", digits=d)
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank()))

# ---- Fig A: Bz forest (all covariates x components, 95% CI) ----
Bz <- F2$Bz; SE <- F2$SE; P <- nrow(Bz); Km1 <- ncol(Bz)
dff <- data.frame(cov=rep(rownames(Bz), Km1),
                  comp=factor(rep(seq_len(Km1), each=P)),
                  est=as.vector(Bz), se=as.vector(SE))
dff$lo <- dff$est-1.96*dff$se; dff$hi <- dff$est+1.96*dff$se
pA <- ggplot(dff, aes(est, interaction(comp,cov), colour=cov)) +
  geom_vline(xintercept=0, linetype=2, colour="grey60") +
  geom_point(size=2) + geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25) +
  guides(colour="none") +
  labs(x=expression(hat(B)[z]*" (standardized, 95% CI)"), y=NULL,
       title="sgscatm ILR path coefficients (Zeller CRC)")
ggsave("output/figures/phase2_Bz_forest.pdf", pA, width=6.5, height=5)

# ---- Fig B: G4 calibration analytical SE vs bootstrap SD ----
dcal <- data.frame(ana=as.vector(F2$SE), boot=as.vector(F2$boot_sd))
lim <- range(c(dcal$ana, dcal$boot))
pB <- ggplot(dcal, aes(boot, ana)) +
  geom_abline(slope=1, intercept=0, colour="grey60") +
  geom_abline(slope=1.25, intercept=0, linetype=3, colour="grey70") +
  geom_abline(slope=0.8,  intercept=0, linetype=3, colour="grey70") +
  geom_point(colour="#534AB7", size=2, alpha=0.8) +
  coord_equal(xlim=lim, ylim=lim) +
  labs(x="document-bootstrap SD (aligned)", y="analytical SE (sgscatm_vcov)",
       title=sprintf("G4 SE calibration (median ratio %.2f, %d%% within 25%%)",
                     median(F2$ratio), round(100*mean(abs(log(F2$ratio))<log(1.25)))))
ggsave("output/figures/phase2_se_calibration.pdf", pB, width=5, height=5)

# ---- Fig C: G6 timing bar (log scale) ----
tt <- CF$timings; tt <- tt[!is.na(unlist(tt))]
dtime <- data.frame(method=names(unlist(tt)), sec=as.numeric(unlist(tt)))
dtime$method <- sub("\\.elapsed$","",dtime$method)
dtime <- dtime[order(dtime$sec), ]
dtime$method <- factor(dtime$method, levels=dtime$method)
pC <- ggplot(dtime, aes(method, sec, fill=method=="sgscatm")) +
  geom_col(width=0.65) + scale_y_log10() +
  geom_text(aes(label=sprintf("%.2fs", sec)), vjust=-0.3, size=3.3) +
  scale_fill_manual(values=c("grey60","#B7434A"), guide="none") +
  labs(x=NULL, y="wall-clock (s, log10)", title="G6 runtime: sgscatm vs competitors")
ggsave("output/figures/phase2_timings.pdf", pC, width=6, height=4)

# ---- Fig D: known CRC taxa directional loadings ----
ki <- CF$known_in
dk <- data.frame(genus=ki, loading=CF$loadings[ki], snr=CF$snr[ki])
dk <- dk[order(dk$loading), ]; dk$genus <- factor(dk$genus, levels=dk$genus)
pD <- ggplot(dk, aes(loading, genus, fill=loading>0)) +
  geom_col(width=0.6) +
  geom_vline(xintercept=0, colour="grey50") +
  scale_fill_manual(values=c("#3b7dd8","#B7434A"), guide="none") +
  labs(x="loading on sgscatm disease (CRC) direction", y=NULL,
       title=sprintf("Known CRC genera: %d/%d load CRC-ward (sign p=%.3f)",
                     CF$n_pos, length(ki), CF$p_sign))
ggsave("output/figures/phase2_known_taxa.pdf", pD, width=6, height=3.2)

# ---- Table: G6 timings ----
sg <- dtime$sec[dtime$method=="sgscatm"]
lines <- c("\\begin{tabular}{lrr}","\\toprule",
  "method & wall-clock (s) & speedup vs sgscatm \\\\","\\midrule",
  apply(dtime, 1, function(r) sprintf("%s & %s & %s \\\\", r["method"],
        fmt(as.numeric(r["sec"]),3),
        ifelse(r["method"]=="sgscatm","--", paste0(fmt(as.numeric(r["sec"])/sg,1),"$\\times$")))),
  "\\bottomrule","\\end{tabular}")
writeLines(lines, "output/tables/phase2_timings.tex")

# ---- Table: G4 SE calibration ----
rat <- F2$ratio
lines4 <- c("\\begin{tabular}{lrrr}","\\toprule",
  "covariate & analytical SE & bootstrap SD & ratio \\\\","\\midrule",
  sapply(seq_len(nrow(Bz)), function(i) sprintf("%s & %s & %s & %s \\\\",
    rownames(Bz)[i], fmt(mean(F2$SE[i,])), fmt(mean(F2$boot_sd[i,])),
    fmt(median(rat[i,]),2))),
  "\\midrule",
  sprintf("median (all entries) & & & %s \\\\", fmt(median(rat),2)),
  sprintf("legacy \\texttt{ilr\\_se\\_analytical} & %s & & %s$\\times$ \\\\",
          fmt(median(F2$se_leg),2), fmt(median(F2$se_leg)/median(F2$SE),1)),
  "\\bottomrule","\\end{tabular}")
writeLines(lines4, "output/tables/phase2_se_calibration.tex")

# ---- Table: G1 delocalization ----
g1 <- F2$G1
lines1 <- c("\\begin{tabular}{llrrr}","\\toprule",
  "level & prev & N & $r_{\\mathrm{fit}}$ & $r_{\\mathrm{clr}}$ \\\\","\\midrule",
  apply(g1,1,function(r) sprintf("%s & %s & %d & %s & %s \\\\", r["level"],
    fmt(as.numeric(r["prev"]),2), as.integer(r["N"]),
    fmt(as.numeric(r["r_fit"]),2), fmt(as.numeric(r["r_clr"]),1))),
  "\\bottomrule","\\end{tabular}")
writeLines(lines1, "output/tables/phase2_delocalization.tex")

# ---- Table: G5a covariate significance ----
comp <- CP$comp
lines5 <- c("\\begin{tabular}{lrr}","\\toprule",
  "covariate & sgscatm Wald $p$ & PERMANOVA $p$ \\\\","\\midrule",
  apply(comp,1,function(r) sprintf("%s & %s & %s \\\\", r["covariate"],
    fmt(as.numeric(r["sgscatm_wald_p"]),4), fmt(as.numeric(r["permanova_p"]),3))),
  "\\bottomrule","\\end{tabular}")
writeLines(lines5, "output/tables/phase2_covariate_sig.tex")

cat("Phase 2 figures + tables written.\n")
cat(sprintf("G6 speedups: %s\n", paste(sprintf("%s=%.0fx", dtime$method, dtime$sec/sg), collapse=" ")))
