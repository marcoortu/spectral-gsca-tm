## =====================================================================
## 05_anchor_G6.R — G6 Gate-2 (Delta_Phi) CLEAN probe. Scale-PRESERVING
## corruption of the topic-word anchor Phi (rows stay on the simplex, so
## the loading scale is unchanged; only topic-recovery error varies).
## Report RMSE(B_z0) and coverage vs delta_Phi, and the boundary. Prereg G6.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
source(file.path(getwd(), "replication/one_step/sweep_utils.R"))

DELTA <- c(0, 0.05, 0.10, 0.20, 0.40, 0.80)
NREP  <- 40L
BASE  <- list(M = 2000L, N = 500L, K = 5L, P = 3L, b_max = 0.5,
              sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L)
WHICH <- c("proj", "onestep_mw")

cat("=== 05 anchor Delta_Phi (M=2000,b_max=0.5,L=200,n_rep=", NREP, ") ===\n",
    sep = "")
# sanity: corruption preserves loading scale (row sums stay 1)
b0 <- sim_dgp(M = 10, N = 500, K = 5, P = 3, seed = 1)$Beta
for (d in DELTA)
  cat(sprintf("  delta=%.2f  anchor row-sum range=[%.3f,%.3f]  ||anchor||/||Beta||=%.3f\n",
      d, min(rowSums(corrupt_anchor(b0, d, seed = 1))),
      max(rowSums(corrupt_anchor(b0, d, seed = 1))),
      norm_F(corrupt_anchor(b0, d, seed = 1)) / norm_F(b0)))

res <- list()
for (d in DELTA) {
  t0 <- proc.time()[3]
  cell <- run_sweep_cell(c(BASE), NREP, WHICH,
                         anchor_fun = local({
                           dd <- d
                           function(dat) corrupt_anchor(dat$Beta, dd,
                                                        seed = 777L)
                         }))
  cell$delta <- d
  res[[as.character(d)]] <- cell
  cat(sprintf("delta=%.2f done in %.0fs\n", d, proc.time()[3] - t0))
  show <- cell[, c("estimator","rmse_norm","bias_norm","se_sd","cov_Bz0")]
  show[-1] <- lapply(show[-1], round, 4)
  print(show, row.names = FALSE); cat("\n")
}
df <- do.call(rbind, res)
saveRDS(df, file.path(getwd(), "replication/one_step/out_05_anchor.rds"))
write.csv(df, file.path(getwd(), "replication/one_step/out_05_anchor.csv"),
          row.names = FALSE)
cat("saved out_05_anchor.{rds,csv}\n")
