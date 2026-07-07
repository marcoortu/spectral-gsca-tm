# ===================================================================
#  lambda_diag.R  —  Part 3: lambda-scaling diagnostic (Section 4)
# ===================================================================
#  Re-runs the Block-1 grid for lambda=1 (nominal) AND lambda=lambda_A
#  (data-driven) and asks whether the penalty is asymptotically inert.
#  The estimator here is the AUGMENTED, lambda-dependent standardized
#  fit (sqrt(M) * fit$Bz) so that any lambda effect is visible.
# ===================================================================

run_lambda_diag <- function() {
  cat("\n====== PART 3 : lambda diagnostic ======\n")
  B0     <- readRDS(file.path(DATA_DIR, "B0.rds"))
  B_star <- readRDS(file.path(DATA_DIR, "B_star.rds"))
  lambda_A_probe <- readRDS(file.path(DATA_DIR, "lambda_A.rds"))
  Km1 <- K_TOPICS - 1L

  proc_align <- function(B) {
    sv <- svd(crossprod(B, B_star)); B %*% (sv$u %*% t(sv$v))
  }
  # augmented standardized estimator at a given lambda
  aug_fit <- function(dat, lambda) {
    fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = lambda, rotate = TRUE)
    M <- nrow(dat$W)
    list(B = proc_align(sqrt(M) * fit$Bz),
         se = { sv <- svd(crossprod(sgscatm_vcov(fit)$B, B_star))
                sgscatm_vcov(fit, rotation = sv$u %*% t(sv$v))$se },
         weig = mean(sort(eigen(crossprod(fit$W_tilde), symmetric = TRUE,
                                only.values = TRUE)$values, decreasing = TRUE)[1:Km1]) / M,
         lam_A = lambda_A_rule(fit))
  }

  rows <- list(); ri <- 0L
  for (setting in c("lambda1", "lambdaA")) {
    for (mi in seq_along(M_B1)) {
      M <- M_B1[mi]; ri <- ri + 1L
      bias <- cover <- numeric(N_REP_LAMBDA)
      weig <- lamv <- numeric(N_REP_LAMBDA)
      for (r in seq_len(N_REP_LAMBDA)) {
        dat <- gen_clean(M, B0, SEED_B1 + mi * 10000L + r)
        lam <- if (setting == "lambda1") 1 else lambda_A_rule(
          sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE))
        af <- aug_fit(dat, lam)
        bias[r]  <- sqrt(sum((af$B - B_star)^2))          # Frobenius
        cover[r] <- mean(B_star >= af$B - ZQ * af$se & B_star <= af$B + ZQ * af$se)
        weig[r]  <- af$weig; lamv[r] <- lam
      }
      rows[[ri]] <- data.frame(
        setting = setting, M = M,
        lambda = mean(lamv),
        bias = mean(bias), coverage95 = mean(cover),
        top_word_eig = mean(weig),               # O(1)
        cov_block_eig = mean(lamv) / M)          # lambda / M
      cat(sprintf("  %s M=%4d  lambda=%.3g  bias=%.4f  cov=%.3f  topWeig=%.3g  covBlk=%.3g\n",
                  setting, M, mean(lamv), mean(bias), mean(cover),
                  mean(weig), mean(lamv) / M))
    }
  }
  df <- do.call(rbind, rows)
  write.csv(df, file.path(TAB_DIR, "lambda_diagnostic.csv"), row.names = FALSE)
  cat("  Wrote tables/lambda_diagnostic.csv\n")
  saveRDS(df, file.path(DATA_DIR, "lambda_diag.rds"))

  # --- verdict ----------------------------------------------------
  b1 <- df$bias[df$setting == "lambda1"]; bA <- df$bias[df$setting == "lambdaA"]
  Ms <- df$M[df$setting == "lambda1"]
  r_lo <- b1[Ms == min(Ms)] / bA[Ms == min(Ms)]
  r_hi <- b1[Ms == max(Ms)] / bA[Ms == max(Ms)]
  inert <- r_hi < r_lo && r_hi < 1.5
  cat(sprintf("\n  LAMBDA VERDICT: bias(1)/bias(A) at M=%d vs M=%d = %.2f , %.2f  ->  %s\n",
              min(Ms), max(Ms), r_lo, r_hi,
              if (inert) "PENALTY ASYMPTOTICALLY INERT (Section 4 lambda-free variance confirmed)"
              else "PENALTY NOT INERT (revise Section 4)"))
  invisible(list(df = df, inert = inert, r_lo = r_lo, r_hi = r_hi))
}
