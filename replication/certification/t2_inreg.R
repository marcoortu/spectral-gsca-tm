#!/usr/bin/env Rscript
# T2-inreg-feasible (KEY): feasible frozen-Phi chain + analytic sandwich, large L
# (1e4) + genuinely clean anchors (N=100, exclusive topics alpha_beta=0.01), vs
# known Bz,0. Tests whether the FEASIBLE deliverable + sandwich is calibrated when
# BOTH Cor-18 gates close. PASS = coverage>=0.90 AND SE/SD in [0.9,1.6];
# under-coverage (SE/SD<0.9 & cov<0.90) = FAIL. Bootstrap on one cell = conservative check.
args <- commandArgs(trailingOnly = TRUE)
NC <- if ("--ncores" %in% args) as.integer(args[which(args=="--ncores")+1L]) else 8L
ROOT <- normalizePath(".")
suppressPackageStartupMessages(devtools::load_all(ROOT, quiet = TRUE))
source("replication/simulation/sim_dgp.R")
K <- 5L; P <- 4L
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30, 0.05,-0.30,0.15,0.25), nrow = P, byrow = TRUE)
V <- ilr_contrast(K)
libfun <- function(m) pmax(rnbinom(m, size = 3, mu = 1e4), 500L)

one_rep <- function(M, rep_id, do_boot = FALSE) {
  dat <- sim_dgp(M = M, N = 100L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                 alpha_beta = 0.01, doc_length = libfun, seed = 74000L + M + rep_id)
  ch <- sgscatm_chain(dat$W, dat$C, K = K, refine = "frozen_phi")
  al <- perm_sign_align(ch$Bz, Bz0, V)
  cv <- perm_sign_coverage(ch$Bz, Bz0, vcov(ch), V)
  # clean-anchor check: exclusivity of anchored Phi
  Phi <- ch$anchor_Phi
  excl <- min(apply(sweep(Phi,2L,pmax(colSums(Phi),1e-300),"/"),1L,max))
  out <- list(B = al$B, se = cv$se, covers = mean(cv$covers), excl = excl,
              se_boot = NA)
  if (do_boot) {
    bs <- chain_boot_se(ch, dat$W, dat$C, B = 80L, seed = 74000L + rep_id + 3L)
    out$se_boot <- median(bs$se / cv$se)   # boot/sandwich ratio
  }
  out
}

cl <- parallel::makePSOCKcluster(NC)
parallel::clusterExport(cl, c("ROOT","K","P","Bz0","V","libfun","one_rep"), envir=environment())
parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(devtools::load_all(ROOT, quiet=TRUE))
  source(file.path(ROOT,"replication/simulation/sim_dgp.R")); TRUE })

rows <- list()
for (M in c(1000L, 5000L)) {
  R <- 50L
  do_boot_ids <- if (M == 1000L) 1:10 else integer(0)   # bootstrap check on M=1000
  t0 <- proc.time()[3]
  res <- parallel::parLapply(cl, seq_len(R), function(r, Mv, dbi)
                             one_rep(Mv, r, do_boot = r %in% dbi),
                             Mv = M, dbi = do_boot_ids)
  Barr <- array(unlist(lapply(res,`[[`,"B")), c(P,K-1L,R))
  emp_sd <- apply(Barr, c(1,2), sd)
  se_mean <- Reduce(`+`, lapply(res,`[[`,"se"))/R
  boot_ratio <- mean(vapply(res,`[[`,numeric(1),"se_boot"), na.rm=TRUE)
  rows[[length(rows)+1L]] <- data.frame(
    M = M, R = R,
    coverage = mean(vapply(res,`[[`,numeric(1),"covers")),
    se_sd = median(se_mean/emp_sd),
    excl = mean(vapply(res,`[[`,numeric(1),"excl")),
    boot_over_sandwich = boot_ratio)
  cat(sprintf("M=%d cov=%.3f SE/SD=%.3f excl=%.3f boot/sandwich=%.2f [%.0fs]\n",
              M, rows[[length(rows)]]$coverage, rows[[length(rows)]]$se_sd,
              rows[[length(rows)]]$excl, boot_ratio, proc.time()[3]-t0))
}
parallel::stopCluster(cl)
T2 <- do.call(rbind, rows); saveRDS(T2, "output/t2_inreg.rds")
print(T2, row.names=FALSE, digits=3)
pass <- all(T2$coverage >= 0.90 & T2$se_sd >= 0.9 & T2$se_sd <= 1.6)
cat(sprintf("\nT2-inreg-feasible: %s\n", ifelse(pass,"PASS","FAIL")))
cat("t2_inreg DONE\n")
