#!/usr/bin/env Rscript
# CRC gates on the frozen-Phi chain (sandwich primary). T4 delocalization (G1),
# G5b known genera on disease axis, T6 speed (G6). G4/G5a come from gate4_crc.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
D  <- readRDS("output/phase2_data.rds")
C  <- D$C; M <- nrow(C); K <- 5L
to_counts <- function(X) round(X / rowSums(X) * 1e6)
r_deloc <- function(X) {
  ch <- sgscatm_chain(to_counts(as.matrix(X)), C, K = K, refine = "frozen_phi")
  list(r = M * max(rowSums(ch$Z^2)) / (K-1L), ch = ch, N = ncol(X),
       sweeps = ch$sweeps)
}

## T4 / G1: delocalization on chain scores, genus (prev filters) + species
cat("== T4/G1 delocalization on chain z ==\n")
G1 <- data.frame(level=character(), prev=numeric(), N=integer(), r=numeric(),
                 sweeps=integer())
for (nm in c("Xg10","Xg20","Xg30","Xsp10","Xsp20")) {
  rr <- r_deloc(D[[nm]])
  lvl <- ifelse(grepl("g", nm), "genus", "species")
  pv  <- as.numeric(sub("[A-Za-z]+","",nm))/100
  G1 <- rbind(G1, data.frame(level=lvl, prev=pv, N=rr$N, r=rr$r, sweeps=rr$sweeps))
  cat(sprintf("  %-7s prev=%.2f N=%3d r=%.2f (sweeps=%d)\n", lvl, pv, rr$N, rr$r, rr$sweeps))
}

## fit at genus prev>=10 for the substantive gates
chg <- sgscatm_chain(to_counts(as.matrix(D$Xg10)), C, K = K, refine = "frozen_phi")
genera <- colnames(D$Xg10)

## G5b: known CRC genera on the disease axis (fitted composition shift, CRC-ctrl)
Xhat <- chg$Theta %*% chg$Phi                 # M x N fitted genus freq
crc  <- C[,"study_condition"] > 0
dshift <- colMeans(Xhat[crc,,drop=FALSE]) - colMeans(Xhat[!crc,,drop=FALSE])
names(dshift) <- genera
known <- D$known_crc; known_in <- intersect(known, genera)
n_pos <- sum(dshift[known_in] > 0)
p_sign <- binom.test(n_pos, length(known_in), 0.5, alternative="greater")$p.value
cat(sprintf("\n== G5b: %d/%d known CRC genera shift CRC-ward; sign-test p=%.4g ==\n",
            n_pos, length(known_in), p_sign))
print(round(dshift[known_in], 5))

## T6 / G6: speed, three paths + competitors
tp <- system.time(sgscatm_chain(to_counts(as.matrix(D$Xg10)), C, K=K, refine="frozen_phi"))[3]
ts <- system.time(vcov(chg))[3]
comp <- tryCatch(readRDS("output/phase2_competitors_fix.rds")$timings, error=function(e) NULL)
tb <- tryCatch(readRDS("output/gate4_crc.rds")$t_boot, error=function(e) NA)
cat("\n== T6/G6 speed (s) ==\n")
cat(sprintf("  chain point            : %.3f\n", tp))
cat(sprintf("  chain point+sandwich   : %.3f  (deliverable, fast path)\n", tp+ts))
cat(sprintf("  chain point+bootstrap  : %.1f  (conservative check)\n", tb))
if (!is.null(comp)) {
  for (nm in c("permanova","stm","aldex2")) if (!is.null(comp[[nm]]))
    cat(sprintf("  %-22s : %.3f  (%.0fx vs point+sandwich)\n", nm,
                comp[[nm]], comp[[nm]]/(tp+ts)))
}

saveRDS(list(G1=G1, dshift=dshift, known_in=known_in, n_pos=n_pos, p_sign=p_sign,
             t_point=tp, t_sandwich=ts, t_boot=tb, comp=comp, Bz=chg$Bz),
        "output/crc_chain.rds")
cat("\ncrc_chain DONE\n")
