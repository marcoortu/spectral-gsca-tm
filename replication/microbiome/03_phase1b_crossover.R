#!/usr/bin/env Rscript
# ===================================================================
# Phase 1b — Delocalization crossover (Proposition 16), done right.
#
# The main sweep measured the delocalization ratio on the STANDARDIZED
# eigenvector scores (fit$Z, unit-norm columns) -> scale-invariant, hence flat.
# Here we measure it on the NATURAL ILR scores and push covariate strength high
# enough (small residual noise) that compositions genuinely leave the centroid,
# so max_i||z_i|| grows and Assumption 5 is stressed.
#
# Records, per b_max cell: natural-scale delocalization ratio r_true, coverage
# of the (standardized) estimand, SE/SD calibration, and RMSE of recovering the
# true Bz0 (scale-recovered). Predicts: r_true up, recovery RMSE up, and at the
# extreme, coverage/calibration degrade.
# Output: output/phase1b_crossover.rds, figure phase1b_crossover.pdf
# ===================================================================
suppressPackageStartupMessages({ devtools::load_all(".", quiet = TRUE); library(ggplot2) })
source("replication/simulation/sim_dgp.R")
source("replication/simulation/sim_utils.R")

M <- 2000L; N <- 200L; K <- 5L; P <- 4L; Km1 <- K - 1L
lam <- 1; R <- 50L; sig_eps <- 0.12
libsize_fun <- function(m) pmax(rnbinom(m, size = 3, mu = 2e4), 500L)
align_to <- function(B, Ref) { sv <- svd(crossprod(B, Ref)); B %*% (sv$u %*% t(sv$v)) }

bmax_grid <- c(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
out <- vector("list", length(bmax_grid))

for (j in seq_along(bmax_grid)) {
  bm <- bmax_grid[j]
  set.seed(3000L + j)
  Bz0 <- matrix(runif(P * Km1, -bm, bm), P, Km1)
  Bstd <- vector("list", R); fits <- vector("list", R)
  rtrue <- numeric(R); rfit <- numeric(R); rmse_rec <- numeric(R); maxz <- numeric(R)
  for (r in seq_len(R)) {
    dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0, sigma_eps = sig_eps,
                   alpha_beta = 0.05, doc_length = libsize_fun, seed = j * 20000L + r)
    fit <- sgscatm(dat$W, dat$C, K = K, lambda = lam, rotate = FALSE)
    fits[[r]] <- fit; Bstd[[r]] <- sgscatm_vcov(fit, identified = FALSE)$B
    zt2 <- rowSums(dat$Z_true^2)
    rtrue[r] <- M * max(zt2) / Km1           # natural-scale delocalization ratio
    maxz[r]  <- sqrt(max(zt2))
    rfit[r]  <- M * max(rowSums(fit$Z^2)) / Km1
    Amap <- tryCatch(solve(crossprod(fit$Z), crossprod(fit$Z, dat$Z_true)),
                     error = function(e) NULL)
    rmse_rec[r] <- if (!is.null(Amap))
      sqrt(mean((solve(crossprod(dat$C), crossprod(dat$C, fit$Z %*% Amap)) - Bz0)^2)) else NA
  }
  Ref <- Bstd[[1]]
  for (p in 1:3) { Bal <- lapply(Bstd, align_to, Ref); Ref <- Reduce(`+`, Bal)/R }
  Bal <- lapply(Bstd, align_to, Ref)
  emp_sd <- apply(array(unlist(Bal), c(P, Km1, R)), c(1,2), sd)
  se_list <- vector("list", R); cov_arr <- array(NA, c(P, Km1, R))
  for (r in seq_len(R)) {
    Rr <- { sv <- svd(crossprod(Bstd[[r]], Ref)); sv$u %*% t(sv$v) }
    vc <- sgscatm_vcov(fits[[r]], rotation = Rr, identified = TRUE)
    se_list[[r]] <- vc$se
    cov_arr[,,r] <- (Ref >= Bal[[r]] - 1.96*vc$se) & (Ref <= Bal[[r]] + 1.96*vc$se)
  }
  se_mean <- Reduce(`+`, se_list)/R
  out[[j]] <- data.frame(b_max = bm, r_true = mean(rtrue), max_z = mean(maxz),
                         r_fit = mean(rfit), coverage = mean(cov_arr),
                         se_sd_ratio = median(se_mean/emp_sd),
                         rmse_recover = mean(rmse_rec, na.rm = TRUE))
  cat(sprintf("bmax=%.1f  r_true=%.1f max||z||=%.2f  cov=%.3f SE/SD=%.2f rmse=%.3f\n",
              bm, out[[j]]$r_true, out[[j]]$max_z, out[[j]]$coverage,
              out[[j]]$se_sd_ratio, out[[j]]$rmse_recover))
}

df <- do.call(rbind, out)
saveRDS(df, "output/phase1b_crossover.rds")
write.csv(df, "output/tables/phase1b_crossover.csv", row.names = FALSE)

# figure: recovery RMSE and coverage vs natural-scale delocalization ratio
theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()))
sc <- max(df$rmse_recover) / 1
p <- ggplot(df, aes(r_true)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey60") +
  geom_line(aes(y = coverage, colour = "coverage"), linewidth = 0.8) +
  geom_point(aes(y = coverage, colour = "coverage"), size = 2.5) +
  geom_line(aes(y = rmse_recover/sc*0.95 + 0.0, colour = "RMSE (scaled)"), linewidth = 0.8) +
  geom_point(aes(y = rmse_recover/sc*0.95, colour = "RMSE (scaled)"), size = 2.5) +
  scale_colour_manual(values = c("coverage" = "#534AB7", "RMSE (scaled)" = "#B7434A"),
                      name = NULL) +
  labs(x = expression("natural-scale delocalization ratio " * r[true]),
       y = "95% coverage  /  recovery RMSE (scaled)",
       title = "Delocalization crossover (Proposition 16)")
ggsave("output/figures/phase1b_crossover.pdf", p, width = 6.5, height = 4)
cat("Phase 1b done.\n")
