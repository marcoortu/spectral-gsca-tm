#!/usr/bin/env Rscript
# A1-corrected feasible Part 2/3 (reuses the valid oracle Part 1 = diag_P1_oracle.rds).
# Point estimate is ALWAYS the chain B_hat at frozen@k* (A1); jackknife contributes
# only var_add=(B_A-B_B)^2/4. B_jk-centered coverage is a labelled SECONDARY column.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
sg <- getNamespace("sgscatm")
K <- 5L; P <- 4L; KSTAR <- 1L
set.seed(7)
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30, 0.05,-0.30,0.15,0.25), nrow = P, byrow = TRUE)
V <- ilr_contrast(K)
libfun <- function(mu) function(m) pmax(rnbinom(m, size = 3, mu = mu), 500L)

# A2 identity check (topic-term matrix is K x N)
d0 <- sim_dgp(M = 50L, N = 200L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
              alpha_beta = 0.05, doc_length = libfun(1e4), seed = 1L)
stopifnot(all(dim(d0$Beta) == c(K, ncol(d0$W))))
cat(sprintf("A2 OK: dim(Beta) = %d x %d = (K, N)\n", nrow(d0$Beta), ncol(d0$Beta)))

jackknife <- function(W, C, Vc, Z_full, Phi, seed) {
  set.seed(seed)
  A <- matrix(rbinom(length(W), as.vector(W), 0.5), nrow(W), ncol(W)); B <- W - A
  zfit <- function(Wp) { Wf <- Wp / pmax(rowSums(Wp), 1)
    sg$.sg_z_step(Z_full, Phi, Wf, Vc, 0, NULL, rep(1e-6, nrow(Wp)),
                  n_gn = 8L, dz_cap = 1)$Z }
  BA <- sg$.sg_b_step(zfit(A), C); BB <- sg$.sg_b_step(zfit(B), C)
  Bf <- sg$.sg_b_step(Z_full, C)
  list(B_jk = 2 * Bf - (BA + BB) / 2, var_add = (BA - BB)^2 / 4)
}
# aggregate aligned estimates + aligned SEs -> metric row
agg <- function(Bal, SE, tag) {
  R <- length(Bal); arr <- array(unlist(Bal), c(P, K-1L, R))
  Bbar <- apply(arr, c(1,2), mean); emp_sd <- apply(arr, c(1,2), sd)
  se_mean <- Reduce(`+`, SE) / R
  data.frame(tag = tag, bias_norm = sqrt(sum((Bbar - Bz0)^2)),
             se_sd = median(se_mean / emp_sd),
             true_cov = mean(mapply(function(b,s) mean(abs(b-Bz0) <= 1.96*s), Bal, SE)),
             biasrem_cov = mean(mapply(function(b,s) mean(abs(b-Bbar) <= 1.96*s), Bal, SE)))
}

feasible_cell <- function(M, R) {
  # storage: plain (center B_hat), jkinf (center B_hat + var_add), sec (center B_jk + var_add)
  bp <- sp <- bi <- si <- bs <- ss <- vector("list", 0)
  monodrop <- 0L
  for (r in seq_len(R)) {
    dat <- sim_dgp(M = M, N = 200L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                   alpha_beta = 0.05, doc_length = libfun(1e4), seed = 71000L + r)
    ch <- sgscatm_chain(dat$W, dat$C, K = K, refine = "frozen_phi", max_sweeps = KSTAR)
    if (!isTRUE(ch$monotone_ok)) monodrop <- monodrop + 1L
    C <- ch$C_centred
    Sig <- vcov(ch)
    jk <- jackknife(dat$W, C, V, ch$Z, ch$anchor_Phi, 71000L + r + 500L)
    Sig_infl <- Sig + diag(as.vector(jk$var_add))
    cvp <- perm_sign_coverage(ch$Bz,   Bz0, Sig,      V)      # plain, center B_hat
    cvi <- perm_sign_coverage(ch$Bz,   Bz0, Sig_infl, V)      # jk-inflated, center B_hat (PRIMARY)
    cvs <- perm_sign_coverage(jk$B_jk, Bz0, Sig_infl, V)      # SECONDARY: center B_jk
    bp <- c(bp, list(perm_sign_align(ch$Bz,Bz0,V)$B));   sp <- c(sp, list(cvp$se))
    bi <- c(bi, list(perm_sign_align(ch$Bz,Bz0,V)$B));   si <- c(si, list(cvi$se))
    bs <- c(bs, list(perm_sign_align(jk$B_jk,Bz0,V)$B)); ss <- c(ss, list(cvs$se))
  }
  out <- rbind(cbind(M=M, agg(bp, sp, "plain")),
               cbind(M=M, agg(bi, si, "jk_inflated")),
               cbind(M=M, agg(bs, ss, "SECONDARY_Bjk")))
  attr(out, "monodrop") <- monodrop
  out
}

cat(sprintf("== Feasible frozen@k*=%d, N=200, L=1e4 (A1-corrected) ==\n", KSTAR))
p2 <- list(); md <- 0L
for (M in c(1000L, 2000L, 5000L)) {
  R <- if (M == 2000L) 80L else 50L
  cell <- feasible_cell(M, R); md <- md + attr(cell, "monodrop")
  p2[[length(p2)+1L]] <- cell
  cat(sprintf("  done M=%d (monotone violations: %d/%d)\n", M, attr(cell,"monodrop"), R))
}
P2 <- do.call(rbind, p2)
saveRDS(P2, "output/diag_P2_feasible.rds")
cat(sprintf("total monotone violations: %d\n", md))
print(P2[, c("M","tag","bias_norm","se_sd","true_cov","biasrem_cov")], row.names=FALSE, digits=3)
cat("\ndiag_feasible DONE\n")
