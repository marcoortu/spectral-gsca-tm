#!/usr/bin/env Rscript
# ===================================================================
#  K and lambda calibration for BES Wave 25 MII topic model
#
#  Diagnostics computed over K x lambda grid:
#    - Semantic coherence  (Mimno et al. 2011)
#    - Exclusivity / FREX  (Roberts et al. 2014)
#    - Eigenvalue gap      (normalised gap at position K-1)
#    - Covariate R²        (ILR variance explained by C)
#
#  Usage:
#    Rscript calibrate_K_lambda.R
# ===================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Source package
source("R/egscatm_fit.R")
source("R/ilr_contrast.R")
source("R/ilr_se.R")
source("R/refine_phi.R")
source("R/methods.R")
source("R/utils.R")

dir.create("output/bes", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# Load preprocessed data
# ---------------------------------------------------------------
cat("Loading preprocessed DTM...\n")
dat    <- readRDS("scripts/bes_case_study/bes_w25_dtm.rds")
W_sp   <- dat$W        # sparse M x N
C      <- dat$C        # M x 5
vocab  <- dat$vocab    # length-N character vector
M      <- nrow(W_sp)
N      <- ncol(W_sp)
P      <- ncol(C)
cat(sprintf("  M=%d  N=%d  P=%d\n", M, N, P))

W <- as.matrix(W_sp)

# ---------------------------------------------------------------
# Diagnostic functions
# ---------------------------------------------------------------

# Semantic coherence (Mimno 2011) for top-m words of each topic.
# Uses co-document frequency from the original DTM.
semantic_coherence <- function(Phi, W_bin, m = 10) {
  K   <- nrow(Phi)
  sc  <- numeric(K)
  for (k in seq_len(K)) {
    top <- order(Phi[k, ], decreasing = TRUE)[seq_len(m)]
    pairs <- combn(top, 2)
    co  <- colSums(W_bin[, pairs[1, ]] * W_bin[, pairs[2, ]])
    df1 <- colSums(W_bin[, pairs[1, ]])
    sc[k] <- mean(log((co + 1) / df1))
  }
  sc
}

# Exclusivity: for each topic word, fraction of its corpus weight
# that falls in this topic relative to all topics (FREX with w=0.5).
exclusivity <- function(Phi, m = 10) {
  K     <- nrow(Phi)
  # column-normalise Phi to get P(topic | word)
  col_s <- colSums(Phi)
  col_s[col_s == 0] <- 1
  Phi_n <- sweep(Phi, 2, col_s, "/")
  ex    <- numeric(K)
  for (k in seq_len(K)) {
    top  <- order(Phi[k, ], decreasing = TRUE)[seq_len(m)]
    freq <- rank(Phi[k, top])  / m
    excl <- rank(Phi_n[k, top]) / m
    ex[k] <- mean(1 / (0.5 / excl + 0.5 / freq))  # harmonic FREX
  }
  ex
}

# Normalised eigenvalue gap: (lambda_{K-1} - lambda_K) / lambda_1
# (uses one extra eigenvalue beyond K-1, so fit with K+1 topics)
eig_gap <- function(eigs) {
  n <- length(eigs)
  if (n < 2) return(NA_real_)
  gaps <- diff(-eigs)           # positive = drop from k to k+1
  gaps / eigs[1]
}

# Covariate R²: proportion of Z variance explained by C
# R² = tr(Z' H_C Z) / tr(Z' Z)  where H_C = C(C'C)^{-1}C'
covariate_r2 <- function(Z, C) {
  H  <- C %*% solve(crossprod(C), t(C))   # M x M hat matrix
  HZ <- H %*% Z
  sum(HZ * Z) / sum(Z * Z)
}

# ---------------------------------------------------------------
# Grid
# ---------------------------------------------------------------
K_grid      <- c(5, 6, 7, 8, 9, 10, 12)
lambda_grid <- c(0.5, 1, 2, 3, 5)

W_bin <- (W > 0) * 1.0   # binary document-term matrix for coherence

cat(sprintf("\nFitting %d models on %d x %d grid...\n",
            length(K_grid) * length(lambda_grid),
            length(K_grid), length(lambda_grid)))

results <- vector("list", length(K_grid) * length(lambda_grid))
idx <- 1L

for (K in K_grid) {
  for (lam in lambda_grid) {
    t0  <- proc.time()[3]
    fit <- egscatm(W, C, K = K, lambda = lam, rotate = TRUE)
    elapsed <- proc.time()[3] - t0

    sc   <- semantic_coherence(fit$Phi, W_bin, m = min(10L, N))
    ex   <- exclusivity(fit$Phi, m = min(10L, N))
    r2   <- covariate_r2(fit$Z, C)

    # eigenvalue gap at K-1 (last retained component)
    eigs <- fit$eigenvalues          # length K-1
    gap  <- if (length(eigs) >= 2) (eigs[length(eigs) - 1] - eigs[length(eigs)]) / eigs[1] else NA

    results[[idx]] <- data.frame(
      K           = K,
      lambda      = lam,
      coh_mean    = mean(sc),
      coh_min     = min(sc),
      excl_mean   = mean(ex),
      frex_mean   = mean(0.5 * scale(sc) + 0.5 * scale(ex)),
      eig_gap     = gap,
      cov_r2      = r2,
      elapsed_s   = elapsed
    )

    cat(sprintf("  K=%2d  lam=%.1f  coh=%.3f  excl=%.3f  R2=%.3f  gap=%.4f  [%.1fs]\n",
                K, lam, mean(sc), mean(ex), r2, ifelse(is.na(gap), 0, gap), elapsed))
    idx <- idx + 1L
  }
}

res <- bind_rows(results)
saveRDS(res, "output/bes/calibration_results.rds")

# ---------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------
cat("\n=== Calibration summary (sorted by FREX, lambda=1) ===\n")
summary_tab <- res %>%
  filter(lambda == 1) %>%
  arrange(desc(frex_mean)) %>%
  select(K, lambda, coh_mean, excl_mean, cov_r2, eig_gap)
print(summary_tab, digits = 3)

cat("\n=== Covariate R² by K and lambda ===\n")
r2_wide <- res %>%
  select(K, lambda, cov_r2) %>%
  pivot_wider(names_from = lambda, values_from = cov_r2,
              names_prefix = "lam=") %>%
  arrange(K)
print(r2_wide, digits = 3)

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()))

res$lambda_f <- factor(res$lambda)

# Coherence vs K
p_coh <- ggplot(res, aes(x = K, y = coh_mean, colour = lambda_f, group = lambda_f)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  labs(x = "K (topics)", y = "Mean semantic coherence",
       colour = expression(lambda)) +
  scale_x_continuous(breaks = K_grid)

# Exclusivity vs K
p_excl <- ggplot(res, aes(x = K, y = excl_mean, colour = lambda_f, group = lambda_f)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  labs(x = "K (topics)", y = "Mean exclusivity (FREX)",
       colour = expression(lambda)) +
  scale_x_continuous(breaks = K_grid)

# Covariate R² vs lambda, faceted by K
p_r2 <- ggplot(res, aes(x = lambda, y = cov_r2, colour = factor(K), group = K)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  labs(x = expression(lambda), y = expression("Covariate " * R^2),
       colour = "K") +
  scale_x_continuous(breaks = lambda_grid)

# Coherence-exclusivity trade-off (lambda=1)
p_tradeoff <- res %>%
  filter(lambda == 1) %>%
  ggplot(aes(x = coh_mean, y = excl_mean, label = K)) +
  geom_point(size = 3, colour = "#534AB7") +
  geom_text(nudge_y = 0.005, size = 3.5) +
  labs(x = "Mean semantic coherence",
       y = "Mean exclusivity (FREX)",
       title = expression(lambda * " = 1"))

ggsave("output/bes/calib_coherence.pdf",  p_coh,      width = 6, height = 3.5)
ggsave("output/bes/calib_exclusivity.pdf", p_excl,    width = 6, height = 3.5)
ggsave("output/bes/calib_r2.pdf",          p_r2,      width = 6, height = 3.5)
ggsave("output/bes/calib_tradeoff.pdf",    p_tradeoff, width = 5, height = 4)

cat("\nOutputs saved to output/bes/\n")
cat("  calib_coherence.pdf   — coherence vs K\n")
cat("  calib_exclusivity.pdf — exclusivity vs K\n")
cat("  calib_r2.pdf          — covariate R2 vs lambda\n")
cat("  calib_tradeoff.pdf    — coherence-exclusivity trade-off\n")
cat("  calibration_results.rds\n")
