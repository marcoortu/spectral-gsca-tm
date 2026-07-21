## =====================================================================
## sweep_utils.R  —  replicate harness + metric aggregation for the gates
## =====================================================================
## For a fixed DGP configuration, run n_rep replicates of all requested
## estimators, Procrustes-align each to B_z0, and aggregate:
##   rmse         mean over reps of RMSE(B_z0) after alignment
##   rmse_norm    rmse / ||B_z0||
##   se_sd        median over entries of  mean_SE / SD(across reps)
##   cov_Bz0      coverage of the TRUE B_z0 by per-rep 95% CIs
##   cov_mean     coverage of the across-rep MEAN by per-rep 95% CIs
##   bias_norm    ||mean_aligned - B_z0||_F
##   se_med       median per-entry SE (to detect collapse-to-1 / pinning)
## =====================================================================

# corrupt an anchor while PRESERVING the loading scale (G6). Rotate each
# topic-word row toward a random simplex direction by fraction delta, then
# renormalise to keep rows on the simplex (scale preserved).
corrupt_anchor <- function(Beta, delta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (delta <= 0) return(Beta)
  K <- nrow(Beta); N <- ncol(Beta)
  Bad <- matrix(rgamma(K * N, shape = 0.1), K, N); Bad <- Bad / rowSums(Bad)
  mix <- (1 - delta) * Beta + delta * Bad
  mix / rowSums(mix)
}

run_sweep_cell <- function(dgp_args, n_rep, which,
                           anchor_fun = function(dat) dat$Beta,
                           lambda = 1, conf = 0.95,
                           seed0 = 10000L, verbose = FALSE) {
  z <- qnorm(1 - (1 - conf) / 2)
  acc <- list()   # per estimator: list of aligned B, se, coverage arrays

  # ---- FIX the parameter of interest across replicates -------------
  # sim_dgp() draws a fresh random Bz0/Beta on every call unless supplied.
  # For meaningful SE-vs-SD, coverage, and bias we must hold the TRUE
  # Bz0 (and topics Beta) FIXED and regenerate only C, eps, W per rep.
  M0 <- dgp_args$M; N0 <- dgp_args$N; K0 <- dgp_args$K; P0 <- dgp_args$P
  bm <- if (is.null(dgp_args$b_max)) 0.5 else dgp_args$b_max
  set.seed(seed0)
  Bz0_ref  <- matrix(runif(P0 * (K0 - 1L), -bm, bm), P0, K0 - 1L)
  Beta_ref <- .rdirichlet_matrix(K0, N0,
                if (is.null(dgp_args$alpha_beta)) 0.1 else dgp_args$alpha_beta)
  dgp_fixed <- dgp_args
  dgp_fixed$Bz0 <- Bz0_ref; dgp_fixed$Beta <- Beta_ref
  dgp_fixed$b_max <- NULL   # ignored when Bz0 supplied

  for (r in seq_len(n_rep)) {
    args_r <- c(dgp_fixed, list(seed = seed0 + r))
    dat <- do.call(sim_dgp, args_r)
    anc <- anchor_fun(dat)
    est <- tryCatch(
      sg_all_estimators(dat$W, dat$C, K = dat$params$K, L = dat$doc_lengths,
                        lambda = lambda, V = dat$V, Phi_anchor = anc,
                        which = which),
      error = function(e) { if (verbose) message("rep ", r, ": ", e$message); NULL })
    if (is.null(est)) next
    for (nm in which) {
      ev <- eval_vs_Bz0(est[[nm]], dat$Bz0, conf = conf)
      acc[[nm]]$B   <- c(acc[[nm]]$B,   list(ev$B_aligned))
      acc[[nm]]$se  <- c(acc[[nm]]$se,  list(ev$se_aligned))
      acc[[nm]]$rmse <- c(acc[[nm]]$rmse, ev$rmse)
    }
    if (verbose && r %% 10L == 0L) message("  rep ", r, "/", n_rep)
  }

  P <- nrow(Bz0_ref); Km1 <- ncol(Bz0_ref); nB <- sqrt(mean(Bz0_ref^2))
  rows <- lapply(which, function(nm) {
    a <- acc[[nm]]; R <- length(a$B)
    if (R < 2L) return(NULL)
    Barr <- simplify2array(a$B)          # P x Km1 x R
    Sarr <- simplify2array(a$se)         # P x Km1 x R
    Bmean <- apply(Barr, c(1, 2), mean)
    SDm   <- apply(Barr, c(1, 2), sd)
    SEm   <- apply(Sarr, c(1, 2), mean)
    covB  <- mean(vapply(seq_len(R), function(j)
                mean(abs(Barr[, , j] - array(Bz0_ref, dim(Barr)[1:2])) <=
                     z * Sarr[, , j]), numeric(1)))
    covM  <- mean(vapply(seq_len(R), function(j)
                mean(abs(Barr[, , j] - Bmean) <= z * Sarr[, , j]), numeric(1)))
    data.frame(
      estimator = nm,
      rmse      = mean(a$rmse),
      rmse_norm = mean(a$rmse) / nB,
      se_sd     = median(SEm / pmax(SDm, 1e-12)),
      cov_Bz0   = covB,
      cov_mean  = covM,
      bias_norm = sqrt(mean((Bmean - Bz0_ref)^2)),
      se_med    = median(SEm),
      sd_med    = median(SDm),
      n_ok      = R,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}
