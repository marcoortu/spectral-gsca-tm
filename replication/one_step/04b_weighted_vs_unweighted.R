## =====================================================================
## 04b_weighted_vs_unweighted.R — genuinely exercise estimator (D).
## Under CONSTANT L the 1/L weight factors out of the GN solve, so
## onestep_mw == onestep_uw exactly. Only VARIABLE document lengths
## (heavy-tailed) can distinguish them: mw downweights short/noisy docs.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
source(file.path(getwd(), "replication/one_step/sweep_utils.R"))

NREP <- 60L
# heavy-tailed lengths: many short docs (where the multinomial weight bites)
len_fun <- function(m) pmax(rnbinom(m, size = 1.2, mu = 120), 10L)
BASE <- list(M = 3000L, N = 500L, K = 5L, P = 3L, b_max = 0.7,
             sigma_eps = 0.3, alpha_beta = 0.1, doc_length = len_fun)
WHICH <- c("proj", "onestep_uw", "onestep_mw", "onestep_wls")

cat("=== 04b weighted vs unweighted, VARIABLE L (n_rep=", NREP, ") ===\n",
    sep = "")
# show the length distribution once
set.seed(1); ll <- len_fun(3000)
cat(sprintf("doc lengths: min=%d med=%d mean=%.0f max=%d  (%.0f%% below 50)\n",
            min(ll), as.integer(median(ll)), mean(ll), max(ll),
            100 * mean(ll < 50)))

cell <- run_sweep_cell(BASE, NREP, WHICH, anchor_fun = function(dat) dat$Beta)
show <- cell[, c("estimator","rmse_norm","bias_norm","sd_med","se_sd","cov_Bz0")]
show[-1] <- lapply(show[-1], round, 4)
print(show, row.names = FALSE)
red_uw <- 1 - cell$rmse_norm[cell$estimator=="onestep_uw"] /
              cell$rmse_norm[cell$estimator=="proj"]
red_mw <- 1 - cell$rmse_norm[cell$estimator=="onestep_mw"] /
              cell$rmse_norm[cell$estimator=="proj"]
cat(sprintf("\nRMSE reduction vs proj:  unweighted=%.1f%%  weighted=%.1f%%\n",
            100*red_uw, 100*red_mw))
cat(sprintf("weighted better than unweighted? %s (mw=%.4f vs uw=%.4f)\n",
            cell$rmse_norm[cell$estimator=="onestep_mw"] <=
            cell$rmse_norm[cell$estimator=="onestep_uw"],
            cell$rmse_norm[cell$estimator=="onestep_mw"],
            cell$rmse_norm[cell$estimator=="onestep_uw"]))
saveRDS(cell, file.path(getwd(), "replication/one_step/out_04b_varL.rds"))
