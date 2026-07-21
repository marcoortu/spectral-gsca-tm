## =====================================================================
## 04_Lrobust.R — G5(c) barrier / L-robustness. Fixed M, b_max; vary
## document length L. The one-step residual bias must NOT grow as L
## shrinks (two-gate regime); weighted vs unweighted compared. Prereg G5(c).
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
source(file.path(getwd(), "replication/one_step/sweep_utils.R"))

LGRID <- c(50L, 100L, 200L, 400L, 1000L)
NREP  <- 40L
BASE  <- list(M = 2000L, N = 500L, K = 5L, P = 3L, b_max = 0.5,
              sigma_eps = 0.3, alpha_beta = 0.1)
WHICH <- c("proj", "onestep_uw", "onestep_mw")

cat("=== 04 L-robustness (M=2000, b_max=0.5, n_rep=", NREP, ") ===\n", sep = "")
res <- list()
for (L in LGRID) {
  t0 <- proc.time()[3]
  cell <- run_sweep_cell(c(BASE, list(doc_length = L)), NREP, WHICH,
                         anchor_fun = function(dat) dat$Beta)
  cell$L <- L
  res[[as.character(L)]] <- cell
  cat(sprintf("L=%d done in %.0fs\n", L, proc.time()[3] - t0))
  show <- cell[, c("estimator","rmse_norm","bias_norm","sd_med","se_sd","cov_Bz0")]
  show[-1] <- lapply(show[-1], round, 4)
  print(show, row.names = FALSE); cat("\n")
}
df <- do.call(rbind, res)
saveRDS(df, file.path(getwd(), "replication/one_step/out_04_Lrobust.rds"))
write.csv(df, file.path(getwd(), "replication/one_step/out_04_Lrobust.csv"),
          row.names = FALSE)
cat("saved out_04_Lrobust.{rds,csv}\n")
