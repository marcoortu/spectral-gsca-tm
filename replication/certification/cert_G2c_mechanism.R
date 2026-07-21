#!/usr/bin/env Rscript
# Decisive test of the G2c failure mechanism: is the sandwich under-coverage due
# to the anchor/orientation stage? Compare FEASIBLE (anchor) vs ORACLE (aligned
# to truth) orientation, same refinement + Lemma-17 sandwich, same cells.
# If oracle is calibrated and feasible is not -> anchor variance is the culprit.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
sg <- getNamespace("sgscatm")
K <- 5L; P <- 3L
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30), nrow = P, byrow = TRUE)

oracle_fit <- function(dat) {
  Wf <- dat$W / rowSums(dat$W); V <- dat$V; C <- scale(dat$C, TRUE, FALSE)
  fit <- sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE)
  gl <- sg$.sg_gl_align(fit$Z, scale(dat$Z_true, TRUE, FALSE))
  rf <- sg$.sg_refine(gl$Z, sg$.sg_phi_step(gl$Z, Wf, V), Wf, C, V,
                      mode = "joint", max_sweeps = 60L)
  list(Bz = sg$.sg_b_step(rf$Z, C), Z = rf$Z, C = C, V = V)
}

run_cell <- function(N, M, Lval, R = 20L) {
  fB <- oB <- vector("list", R); fSE <- oSE <- vector("list", R)
  fcov <- ocov <- numeric(R)
  for (rep in 1:R) {
    dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                   alpha_beta = 0.05, doc_length = as.integer(Lval),
                   seed = 63000L + rep)
    cf <- sgscatm_chain(dat$W, dat$C, K = K, refine = "joint", max_sweeps = 60L)
    fB[[rep]] <- perm_sign_align(cf$Bz, Bz0, dat$V)$B
    fcov[rep] <- mean(perm_sign_coverage(cf$Bz, Bz0, vcov(cf), dat$V)$covers)
    fSE[[rep]] <- perm_sign_coverage(cf$Bz, Bz0, vcov(cf), dat$V)$se

    orc <- oracle_fit(dat)
    Sig_o <- sg$.sg_sandwich(orc$Z, orc$C, orc$Bz)
    oB[[rep]] <- perm_sign_align(orc$Bz, Bz0, dat$V)$B
    ocov[rep] <- mean(perm_sign_coverage(orc$Bz, Bz0, Sig_o, dat$V)$covers)
    oSE[[rep]] <- perm_sign_coverage(orc$Bz, Bz0, Sig_o, dat$V)$se
  }
  arr <- function(L) array(unlist(L), c(P, K-1L, length(L)))
  fsd <- apply(arr(fB), c(1,2), sd); osd <- apply(arr(oB), c(1,2), sd)
  fse <- Reduce(`+`, fSE)/R; ose <- Reduce(`+`, oSE)/R
  cat(sprintf("N=%d M=%d L=%.0e | FEASIBLE cov=%.3f SE/SD=%.2f | ORACLE cov=%.3f SE/SD=%.2f\n",
              N, M, Lval, mean(fcov), median(fse/fsd), mean(ocov), median(ose/osd)))
}

run_cell(N = 500L, M = 2000L, Lval = 1e4)   # same as failing G2c cell
run_cell(N = 200L, M = 5000L, Lval = 1e4)   # cleaner anchors (small N, large M)
cat("mechanism test DONE\n")
