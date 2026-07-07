# ===================================================================
#  block4.R  —  Part B: four-method complementary comparison
# ===================================================================
#  Methods: sgscatm, sgscatm+refine (refine_phi k-means), stm, stm+warm.
#  Metrics: MSE(B_z) on the corrected standardized scale, MSE(Phi) under
#  optimal row permutation, wall time, EM iterations (STM only).
#  STM run in-process (R 4.5.1); reuses refine_phi() and std_Bz().
# ===================================================================

.SCEN4 <- list(
  high_sep = list(alpha_beta = 0.01, b_max = 0.30, sigma_eps = 0.3),
  low_sep  = list(alpha_beta = 1.00, b_max = 0.30, sigma_eps = 0.3),
  weak_sig = list(alpha_beta = 0.10, b_max = 0.10, sigma_eps = 0.3))

# exact optimal-row-permutation MSE(Phi) for small K
.perms <- function(x) if (length(x) <= 1L) list(x) else
  do.call(c, lapply(seq_along(x), function(i)
    lapply(.perms(x[-i]), function(s) c(x[i], s))))
.best_phi_mse <- function(Phi, Beta, K) {
  min(vapply(.perms(seq_len(K)), function(p) mean((Phi[p, ] - Beta)^2), numeric(1)))
}
.stm_phi <- function(fs, fmt, N) tryCatch({
  lb <- fs$beta$logbeta[[1L]]
  phi <- exp(lb - apply(lb, 1L, function(x){lse<-max(x); lse+log(sum(exp(x-lse)))}))
  P <- matrix(0, nrow(phi), N); P[, fmt$active] <- phi; P
}, error = function(e) NULL)

.stm_init <- function(fit, fmt, K) {
  Pi <- pmax(fit$Pi, 1e-10); Pi <- Pi / rowSums(Pi)
  eta <- log(Pi[, -K, drop = FALSE]) - log(Pi[, K])
  phi <- pmax(fit$Phi[, fmt$active, drop = FALSE], 1e-10); phi <- phi/rowSums(phi)
  list(mu = list(mu = matrix(colMeans(eta), ncol = 1L)), sigma = diag(K - 1L),
       beta = list(logbeta = list(log(phi))), eta = eta,
       convergence = list(bound = numeric(0), its = 0L, stopits = FALSE,
                          converged = FALSE, allow.neg.change = TRUE))
}

run_block4 <- function() {
  cat("\n====== BLOCK 4 : four-method comparison (Part B) ======\n")
  have_stm <- requireNamespace("stm", quietly = TRUE)
  if (have_stm) suppressPackageStartupMessages(library(stm))
  else cat("  NOTE: stm not available — STM methods = NA.\n")
  Vilr <- ilr_contrast(K_TOPICS)

  raw <- list();  ci <- 0L
  for (M in M_B4) {
    for (scn in SCEN_B4) {
      sc <- .SCEN4[[scn]]; ci <- ci + 1L
      cat(sprintf("  M=%d scenario=%s x %d reps\n", M, scn, N_REP_B4))
      met <- data.frame()
      prev <- as.formula(paste0("~ ", paste0("V", 1:P_COV, collapse = " + ")))
      for (r in seq_len(N_REP_B4)) {
        seed_r <- 40000L + ci * 1000L + r; set.seed(seed_r)
        Bz0r <- matrix(runif(P_COV*(K_TOPICS-1L), -sc$b_max, sc$b_max),
                       P_COV, K_TOPICS-1L)
        dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV, Bz0 = Bz0r,
                       sigma_eps = sc$sigma_eps, alpha_beta = sc$alpha_beta,
                       doc_length = DOC_LEN, seed = seed_r)
        B0std <- std_Bz(dat$Z_true, dat$C)
        fmt <- .to_stm_fmt(dat$W)
        cov_df <- as.data.frame(dat$C); colnames(cov_df) <- paste0("V", 1:P_COV)

        # sgscatm
        t0 <- proc.time()
        fit <- tryCatch(sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1,
                                rotate = TRUE), error = function(e) NULL)
        t_sg <- (proc.time() - t0)[3L]
        if (!is.null(fit)) {
          mse_b <- procrustes_align(std_Bz(fit$Z, dat$C), B0std)$mse
          mse_p <- .best_phi_mse(fit$Phi, dat$Beta, K_TOPICS)
          met <- rbind(met, data.frame(rep=r, method="sgscatm",
                       mse_Bz=mse_b, mse_phi=mse_p, time_s=t_sg, n_iter=NA_integer_))
          # sgscatm + refine
          t0 <- proc.time()
          fr <- tryCatch(refine_phi(fit, dat$W, method="kmeans", seed=seed_r),
                         error = function(e) NULL)
          t_rf <- t_sg + (proc.time()-t0)[3L]
          if (!is.null(fr)) {
            met <- rbind(met, data.frame(rep=r, method="sgscatm_ref",
                         mse_Bz=mse_b, mse_phi=.best_phi_mse(fr$Phi, dat$Beta, K_TOPICS),
                         time_s=t_rf, n_iter=NA_integer_))
          }
        }
        # stm
        if (have_stm) {
          t0 <- proc.time()
          fs <- tryCatch(stm(documents=fmt$docs, vocab=fmt$vocab, K=K_TOPICS,
                            prevalence=prev, data=cov_df, init.type="Spectral",
                            verbose=FALSE, max.em.its=STM_MAXIT), error=function(e) NULL)
          t_s <- (proc.time()-t0)[3L]
          if (!is.null(fs)) {
            Zs <- log(pmax(fs$theta,1e-10)) %*% Vilr
            Ph <- .stm_phi(fs, fmt, N_VOCAB)
            met <- rbind(met, data.frame(rep=r, method="stm",
                         mse_Bz=procrustes_align(std_Bz(Zs, dat$C), B0std)$mse,
                         mse_phi=if(!is.null(Ph)) .best_phi_mse(Ph, dat$Beta, K_TOPICS) else NA,
                         time_s=t_s, n_iter=as.integer(fs$convergence$its)))
          }
          # stm + warm (from sgscatm)
          if (!is.null(fit)) {
            init <- .stm_init(fit, fmt, K_TOPICS)
            t0 <- proc.time()
            fw <- tryCatch(stm(documents=fmt$docs, vocab=fmt$vocab, K=K_TOPICS,
                              prevalence=prev, data=cov_df, model=init,
                              verbose=FALSE, max.em.its=STM_MAXIT), error=function(e) NULL)
            t_w <- t_sg + (proc.time()-t0)[3L]
            if (!is.null(fw)) {
              Zw <- log(pmax(fw$theta,1e-10)) %*% Vilr
              Phw <- .stm_phi(fw, fmt, N_VOCAB)
              met <- rbind(met, data.frame(rep=r, method="stm_warm",
                           mse_Bz=procrustes_align(std_Bz(Zw, dat$C), B0std)$mse,
                           mse_phi=if(!is.null(Phw)) .best_phi_mse(Phw, dat$Beta, K_TOPICS) else NA,
                           time_s=t_w, n_iter=as.integer(fw$convergence$its)))
            }
          }
        }
        if (r %% 5L == 0L) cat(sprintf("    rep %d/%d\n", r, N_REP_B4))
      }
      raw[[ci]] <- list(M = M, scenario = scn, metrics = met)
    }
  }
  saveRDS(raw, file.path(DATA_DIR, "block4.rds"))
  summarise_block4(raw)
  invisible(raw)
}

summarise_block4 <- function(raw) {
  ord <- c("sgscatm","sgscatm_ref","stm","stm_warm")
  lab <- c(sgscatm="\\texttt{sgscatm}", sgscatm_ref="\\texttt{sgscatm}+refine",
           stm="\\texttt{stm}", stm_warm="\\texttt{stm}+warm")
  scl <- c(high_sep="High sep.", low_sep="Low sep.", weak_sig="Weak signal")

  df <- do.call(rbind, lapply(raw, function(x) {
    if (!nrow(x$metrics)) return(NULL)
    do.call(rbind, lapply(ord, function(m) {
      d <- x$metrics[x$metrics$method == m, ]
      if (!nrow(d)) return(NULL)
      data.frame(M=x$M, scenario=x$scenario, method=m,
                 MSE_Bz=mean(d$mse_Bz, na.rm=TRUE),
                 MSE_Phi=mean(d$mse_phi, na.rm=TRUE),
                 Time=mean(d$time_s, na.rm=TRUE),
                 Iter=mean(d$n_iter, na.rm=TRUE), n=nrow(d))
    }))
  }))
  write.csv(df, file.path(TAB_DIR, "block4.csv"), row.names = FALSE)
  cat("\n  Block 4 summary:\n"); print(round(df[, c("M","MSE_Bz","MSE_Phi","Time","Iter")], 4))

  body <- c()
  for (scn in SCEN_B4) {
    sub <- df[df$scenario == scn, ]; if (!nrow(sub)) next
    body <- c(body, sprintf("\\multicolumn{6}{l}{\\textit{Scenario: %s}} \\\\", scl[scn]),
              "\\midrule")
    for (M in M_B4) {
      first <- TRUE
      for (m in ord) {
        rr <- sub[sub$M == M & sub$method == m, ]; if (!nrow(rr)) next
        it <- if (m %in% c("stm","stm_warm") && is.finite(rr$Iter)) sprintf("%.0f", rr$Iter) else "---"
        body <- c(body, sprintf("%s & %s & %.4f & %.2e & %.2f & %s \\\\",
                  if (first) fmt_int(M) else "", lab[m],
                  rr$MSE_Bz, rr$MSE_Phi, rr$Time, it))
        first <- FALSE
      }
      body <- c(body, "\\addlinespace")
    }
    body <- c(body, "\\midrule")
  }
  cap <- sprintf(paste0(
    "Four-method comparison across corpus sizes and scenarios. ",
    "MSE$(\\hat{\\mathbf{B}}_z)$ on the corrected standardized scale; ",
    "MSE$(\\hat{\\boldsymbol{\\Phi}})$ under optimal row permutation; ",
    "Time in seconds (\\texttt{+warm}/\\texttt{+refine} include the ",
    "\\texttt{sgscatm} base); Iter = mean EM iterations (STM only). ",
    "$K=%d$, $P=%d$, $N=%d$, $\\lambda=1$, %d replicates per cell."),
    K_TOPICS, P_COV, N_VOCAB, N_REP_B4)
  if (REDUCED) cap <- paste(cap, "\\emph{Reduced replicate count (runtime).}")
  lines <- c("\\begin{table}[t]", "\\centering", "\\small",
    sprintf("\\caption{%s}", cap), "\\label{tab:block4}",
    "\\begin{tabular}{rlcccc}", "\\toprule",
    "$M$ & Method & MSE$(\\hat{\\mathbf{B}}_z)$ & MSE$(\\hat{\\boldsymbol{\\Phi}})$ & Time (s) & Iter \\\\",
    "\\midrule", head(body, -1L), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file.path(TAB_DIR, "block4.tex"))
  cat("  Wrote tables/block4.csv, tables/block4.tex\n")
  df
}
