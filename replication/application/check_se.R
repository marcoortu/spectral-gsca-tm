suppressPackageStartupMessages(library(Matrix))
source("R/sgscatm_fit.R"); source("R/ilr_contrast.R")
source("R/ilr_se.R"); source("R/refine_phi.R")
source("R/methods.R"); source("R/utils.R")

fit    <- readRDS("output/bes/fit_K7.rds")
se_res <- readRDS("output/bes/se_K7.rds")
dat    <- readRDS("scripts/bes_case_study/bes_w25_dtm.rds")
W      <- as.matrix(dat$W); C <- dat$C

cat("=== Z scale ===\n")
cat("Z column sd:", round(apply(fit$Z, 2, sd), 4), "\n")
cat("C column sd:", round(apply(C, 2, sd), 4), "\n")

cat("\n=== Bz (full precision) ===\n")
print(fit$Bz)

cat("\n=== SE (full precision) ===\n")
if (is.list(se_res) && "se" %in% names(se_res)) {
  print(se_res$se)
  cat("\n=== z-statistics (Bz / SE) ===\n")
  zstat <- fit$Bz / se_res$se
  print(round(zstat, 2))
  cat("\n=== Significant at |z|>1.96 ===\n")
  print(abs(zstat) > 1.96)
} else {
  cat("SE structure:", names(se_res), "\n")
  print(se_res)
}

cat("\n=== Marginal effects on Pi (Leave vs Remain) ===\n")
P <- ncol(C)
C_leave  <- matrix(0, 1, P); C_leave[1, 5]  <-  1 - mean(dat$df$euref)
C_remain <- matrix(0, 1, P); C_remain[1, 5] <- -mean(dat$df$euref)
Pi_l <- ilr_to_proportions(C_leave  %*% fit$Bz, fit$V)
Pi_r <- ilr_to_proportions(C_remain %*% fit$Bz, fit$V)
delta <- Pi_l - Pi_r
cat("Leave - Remain delta (pp):\n")
names(delta) <- paste0("T", seq_len(fit$K))
print(round(100 * delta, 2))
