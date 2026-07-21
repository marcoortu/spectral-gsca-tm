#!/usr/bin/env Rscript
# Gate 0: full-chain bootstrap SE + in-regime coverage. STOP-IF-FAIL on 0a.
# Usage: Rscript gate0.R --ncores 8 [--full]  (--full also runs 0b)
args <- commandArgs(trailingOnly = TRUE)
NC   <- if ("--ncores" %in% args) as.integer(args[which(args=="--ncores")+1L]) else 8L
FULL <- "--full" %in% args
ROOT <- normalizePath(".")
suppressPackageStartupMessages(devtools::load_all(ROOT, quiet = TRUE))
source("replication/simulation/sim_dgp.R")
K <- 5L; P <- 4L
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30, 0.05,-0.30,0.15,0.25), nrow = P, byrow = TRUE)
V <- ilr_contrast(K)
libfun <- function(mu) function(m) pmax(rnbinom(m, size = 3, mu = mu), 500L)

one_rep <- function(M, Binner, rep_id) {
  dat <- sim_dgp(M = M, N = 200L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                 alpha_beta = 0.05, doc_length = libfun(1e4), seed = 72000L + rep_id)
  ch <- sgscatm_chain(dat$W, dat$C, K = K, refine = "frozen_phi")
  bs <- chain_boot_se(ch, dat$W, dat$C, B = Binner, seed = 72000L + rep_id + 7L)
  al <- perm_sign_align(ch$Bz, Bz0, V)                 # align point est to truth
  # rotate boot SE into the truth-aligned frame: SE is entrywise; align B* were
  # taken in the point-estimate frame, so rotate the point est by al$Q for coverage
  se <- bs$se                                          # in point-estimate frame
  # express coverage in the aligned frame: |al$B - Bz0| <= 1.96 * se_aligned;
  # se transforms by |Q| permutation/sign -> reorder/flip via al$Q
  se_al <- matrix(sqrt(diag(kronecker(al$Q, diag(P)) %*%
                    diag(as.vector(se^2)) %*% t(kronecker(al$Q, diag(P))))), P, K-1L)
  list(B = al$B, se = se_al,
       covers = mean(abs(al$B - Bz0) <= 1.96 * se_al))
}

run_cell <- function(M, R, Binner, cl) {
  res <- if (!is.null(cl)) {
    parallel::parLapply(cl, seq_len(R), function(r) one_rep(M, Binner, r))
  } else lapply(seq_len(R), function(r) one_rep(M, Binner, r))
  Barr <- array(unlist(lapply(res, `[[`, "B")), c(P, K-1L, R))
  emp_sd <- apply(Barr, c(1,2), sd)
  se_mean <- Reduce(`+`, lapply(res, `[[`, "se")) / R
  list(M = M, R = R, Binner = Binner,
       coverage = mean(vapply(res, `[[`, numeric(1), "covers")),
       se_sd = median(se_mean / emp_sd),
       bias_norm = sqrt(sum((apply(Barr, c(1,2), mean) - Bz0)^2)))
}

cl <- if (NC > 1L) {
  cl <- parallel::makePSOCKcluster(NC)
  parallel::clusterExport(cl, c("ROOT","K","P","Bz0","V","libfun","one_rep"), envir=environment())
  parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(devtools::load_all(ROOT, quiet=TRUE))
    source(file.path(ROOT,"replication/simulation/sim_dgp.R")); TRUE })
  cl
} else NULL

cat("== Gate 0a (STOP-IF-FAIL): N=200 L=1e4 M=2000, R=30 x B=100 ==\n")
t0 <- proc.time()[3]
g0a <- run_cell(2000L, 30L, 100L, cl)
cat(sprintf("  coverage=%.3f  SE/SD=%.3f  bias=%.4f  [%.0fs]\n",
            g0a$coverage, g0a$se_sd, g0a$bias_norm, proc.time()[3]-t0))
pass0a <- g0a$se_sd >= 0.85 && g0a$se_sd <= 1.25 &&
          g0a$coverage >= 0.88 && g0a$coverage <= 0.999
cat(sprintf("  Gate 0a: %s\n", ifelse(pass0a, "PASS", "FAIL")))
saveRDS(g0a, "output/gate0a.rds")

g0b <- NULL
if (pass0a && FULL) {
  cat("\n== Gate 0b: M in {2000,5000}, R=30 x B=100 (reduced knobs, disclosed) ==\n")
  g0b <- list()
  for (M in c(2000L, 5000L)) {
    t0 <- proc.time()[3]; cell <- run_cell(M, 30L, 100L, cl)
    cat(sprintf("  M=%d coverage=%.3f SE/SD=%.3f bias=%.4f [%.0fs]\n",
                M, cell$coverage, cell$se_sd, cell$bias_norm, proc.time()[3]-t0))
    g0b[[as.character(M)]] <- cell
  }
  saveRDS(g0b, "output/gate0b.rds")
}
if (!is.null(cl)) parallel::stopCluster(cl)
cat("\ngate0 DONE\n")
