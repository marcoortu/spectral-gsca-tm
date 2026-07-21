## =====================================================================
## 02_sweep_bmax.R — covariate-strength sweep. Evaluates G1, G2, G3, G4,
## G5(a) in one pass over all four estimators. Prereg §4, §6.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
source(file.path(getwd(), "replication/one_step/sweep_utils.R"))

BMAX  <- c(0.1, 0.25, 0.5, 0.75, 1.0, 1.5)
NREP  <- 60L
BASE  <- list(M = 2000L, N = 500L, K = 5L, P = 3L,
              sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L)
WHICH <- c("baseline_std", "proj", "onestep_uw", "onestep_mw")

cat("=== 02 b_max sweep (n_rep=", NREP, ") ===\n", sep = "")
res <- list()
for (bm in BMAX) {
  t0 <- proc.time()[3]
  cell <- run_sweep_cell(c(BASE, list(b_max = bm)), NREP, WHICH,
                         anchor_fun = function(dat) dat$Beta)
  cell$b_max <- bm
  res[[as.character(bm)]] <- cell
  cat(sprintf("b_max=%.2f done in %.0fs\n", bm, proc.time()[3] - t0))
  show <- cell[, c("estimator","rmse_norm","se_sd","cov_Bz0",
                   "cov_mean","bias_norm","se_med","sd_med")]
  show[-1] <- lapply(show[-1], round, 4)
  print(show, row.names = FALSE)
  cat("\n")
}
df <- do.call(rbind, res)
saveRDS(df, file.path(getwd(), "replication/one_step/out_02_bmax.rds"))
write.csv(df, file.path(getwd(), "replication/one_step/out_02_bmax.csv"),
          row.names = FALSE)
cat("saved out_02_bmax.{rds,csv}\n")
