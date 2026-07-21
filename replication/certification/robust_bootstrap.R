#!/usr/bin/env Rscript
# Phase B — robust bootstrap acceptance on the Gate-0 cell (N=200,L=1e4,M=2000).
# Primary robustification (pre-registered, no tuning): MAD-scale SE = 1.4826*MAD
# of the bootstrap B_hat*, vs the SD-scale. Reduced to R=20 outer (disclosed).
# Acceptance: MAD-scale SE/SD in [0.9,1.2] AND coverage in [0.90,0.97].
args <- commandArgs(trailingOnly = TRUE)
NC <- if ("--ncores" %in% args) as.integer(args[which(args=="--ncores")+1L]) else 8L
ROOT <- normalizePath(".")
suppressPackageStartupMessages(devtools::load_all(ROOT, quiet = TRUE))
source("replication/simulation/sim_dgp.R")
K <- 5L; P <- 4L
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30, 0.05,-0.30,0.15,0.25), nrow = P, byrow = TRUE)
V <- ilr_contrast(K)
libfun <- function(mu) function(m) pmax(rnbinom(m, size = 3, mu = mu), 500L)

one_rep <- function(rep_id) {
  dat <- sim_dgp(M = 2000L, N = 200L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                 alpha_beta = 0.05, doc_length = libfun(1e4), seed = 73000L + rep_id)
  ch <- sgscatm_chain(dat$W, dat$C, K = K, refine = "frozen_phi")
  bs <- chain_boot_se(ch, dat$W, dat$C, B = 100L, seed = 73000L + rep_id + 7L)
  al <- perm_sign_align(ch$Bz, Bz0, V)
  rot <- function(se) matrix(sqrt(diag(kronecker(al$Q, diag(P)) %*%
           diag(as.vector(se^2)) %*% t(kronecker(al$Q, diag(P))))), P, K-1L)
  se_sd <- rot(bs$se); se_mad <- rot(bs$se_mad)
  list(B = al$B, se_sd = se_sd, se_mad = se_mad,
       cov_sd  = mean(abs(al$B - Bz0) <= 1.96 * se_sd),
       cov_mad = mean(abs(al$B - Bz0) <= 1.96 * se_mad))
}

R <- 20L
cl <- parallel::makePSOCKcluster(NC)
parallel::clusterExport(cl, c("ROOT","K","P","Bz0","V","libfun","one_rep"), envir=environment())
parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(devtools::load_all(ROOT, quiet=TRUE))
  source(file.path(ROOT,"replication/simulation/sim_dgp.R")); TRUE })
cat(sprintf("== Phase B robust bootstrap: N=200 L=1e4 M=2000, R=%d x B=100 ==\n", R))
t0 <- proc.time()[3]
res <- parallel::parLapply(cl, seq_len(R), one_rep)
parallel::stopCluster(cl)

Barr <- array(unlist(lapply(res, `[[`, "B")), c(P, K-1L, R))
emp_sd <- apply(Barr, c(1,2), sd)
se_sd_mean  <- Reduce(`+`, lapply(res, `[[`, "se_sd"))  / R
se_mad_mean <- Reduce(`+`, lapply(res, `[[`, "se_mad")) / R
out <- data.frame(
  scale = c("SD","MAD"),
  se_sd_ratio = c(median(se_sd_mean/emp_sd), median(se_mad_mean/emp_sd)),
  coverage = c(mean(vapply(res,`[[`,numeric(1),"cov_sd")),
               mean(vapply(res,`[[`,numeric(1),"cov_mad"))))
cat(sprintf("[%.0fs]\n", proc.time()[3]-t0))
print(out, row.names=FALSE, digits=3)
pass_mad <- out$se_sd_ratio[2] >= 0.9 && out$se_sd_ratio[2] <= 1.2 &&
            out$coverage[2] >= 0.90 && out$coverage[2] <= 0.999
cat(sprintf("MAD-scale acceptance: %s (SE/SD=%.2f, cov=%.3f)\n",
            ifelse(pass_mad,"PASS","FAIL"), out$se_sd_ratio[2], out$coverage[2]))
saveRDS(list(out=out, emp_sd=emp_sd, R=R, pass_mad=pass_mad), "output/robust_bootstrap.rds")
cat("robust_bootstrap DONE\n")
