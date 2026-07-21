## =====================================================================
## 03_Mscaling.R — G5(b) DECISIVE consistency gate. Fixed b_max, fixed L;
## RMSE(B_z0) of the one-step must decline toward 0 as M grows, while
## proj-only plateaus at the linearisation-bias floor. Prereg §4 G5(b).
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
source(file.path(getwd(), "replication/one_step/sweep_utils.R"))

MGRID <- c(500L, 1000L, 2000L, 4000L, 8000L)
NREP  <- 40L
BASE  <- list(N = 500L, K = 5L, P = 3L, b_max = 0.5,
              sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L)
WHICH <- c("proj", "onestep_uw", "onestep_mw")

cat("=== 03 M-scaling (b_max=0.5, L=200, n_rep=", NREP, ") ===\n", sep = "")
res <- list()
for (M in MGRID) {
  t0 <- proc.time()[3]
  cell <- run_sweep_cell(c(BASE, list(M = M)), NREP, WHICH,
                         anchor_fun = function(dat) dat$Beta)
  cell$M <- M
  res[[as.character(M)]] <- cell
  cat(sprintf("M=%d done in %.0fs\n", M, proc.time()[3] - t0))
  show <- cell[, c("estimator","rmse_norm","bias_norm","sd_med","se_sd","cov_Bz0")]
  show[-1] <- lapply(show[-1], round, 4)
  print(show, row.names = FALSE); cat("\n")
}
df <- do.call(rbind, res)
saveRDS(df, file.path(getwd(), "replication/one_step/out_03_Mscaling.rds"))
write.csv(df, file.path(getwd(), "replication/one_step/out_03_Mscaling.csv"),
          row.names = FALSE)
cat("saved out_03_Mscaling.{rds,csv}\n")
