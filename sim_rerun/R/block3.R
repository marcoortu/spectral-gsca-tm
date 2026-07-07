# ===================================================================
#  block3.R  —  Part B: sgscatm vs STM (corrected standardized scale)
# ===================================================================
#  Corrected scale: every method's document-topic proportions are mapped
#  to ILR scores, each score column is standardized to unit variance
#  (removing the raw sqrt(M) scale that produced the original metric bug),
#  regressed on C, and Procrustes-aligned to the identically-standardized
#  true coefficient.  Same treatment for sgscatm, STM, and ground truth.
#  STM is run in-process (R 4.5.1); lambda=1 to match the prior design.
# ===================================================================

# Standardized-scale path coefficients: scores -> unit-variance cols ->
# OLS on C.  Returns P x (K-1).
std_Bz <- function(Z, C) {
  Zc  <- scale(Z, center = TRUE, scale = FALSE)
  sds <- apply(Zc, 2L, sd);  sds[sds == 0] <- 1
  Zs  <- sweep(Zc, 2L, sds, "/")
  safe_solve(crossprod(C), crossprod(C, Zs))
}

.to_stm_fmt <- function(W) {
  active <- which(colSums(W) > 0L)
  Wa <- W[, active, drop = FALSE]
  docs <- lapply(seq_len(nrow(Wa)), function(i) {
    idx <- which(Wa[i, ] > 0L); rbind(as.integer(idx), as.integer(Wa[i, idx]))
  })
  list(docs = docs, vocab = as.character(active), active = active)
}

run_block3 <- function() {
  cat("\n====== BLOCK 3 : sgscatm vs STM (Part B) ======\n")
  have_stm <- requireNamespace("stm", quietly = TRUE)
  if (have_stm) suppressPackageStartupMessages(library(stm))
  else cat("  NOTE: stm not available — STM columns = NA.\n")
  Vilr <- ilr_contrast(K_TOPICS)

  raw <- list();  ci <- 0L
  for (M in M_B3) {
    for (snm in names(SIGNAL_B3)) {
      bmax <- SIGNAL_B3[[snm]]; ci <- ci + 1L
      cat(sprintf("  M=%d signal=%s (b_max=%.2f) x %d reps\n",
                  M, snm, bmax, N_REP_B3))
      met <- data.frame()
      for (r in seq_len(N_REP_B3)) {
        set.seed(30000L + ci * 1000L + r)
        Bz0r <- matrix(runif(P_COV * (K_TOPICS - 1L), -bmax, bmax),
                       P_COV, K_TOPICS - 1L)
        dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV, Bz0 = Bz0r,
                       sigma_eps = 0.3, alpha_beta = ALPHA_BETA,
                       doc_length = DOC_LEN, seed = 30000L + ci * 1000L + r)
        B0std <- std_Bz(dat$Z_true, dat$C)                 # normalized truth

        # --- sgscatm ---
        t0 <- proc.time()
        fit <- tryCatch(sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1,
                                rotate = TRUE), error = function(e) NULL)
        t_sg <- (proc.time() - t0)[3L]
        if (!is.null(fit)) {
          mse_sg <- procrustes_align(std_Bz(fit$Z, dat$C), B0std)$mse
          met <- rbind(met, data.frame(rep = r, method = "sgscatm",
                                        mse_Bz = mse_sg, time_s = t_sg))
        }
        # --- STM ---
        if (have_stm) {
          fmt <- .to_stm_fmt(dat$W)
          cov_df <- as.data.frame(dat$C); colnames(cov_df) <- paste0("V", 1:P_COV)
          prev <- as.formula(paste0("~ ", paste0("V", 1:P_COV, collapse = " + ")))
          t0 <- proc.time()
          fs <- tryCatch(stm(documents = fmt$docs, vocab = fmt$vocab,
                             K = K_TOPICS, prevalence = prev, data = cov_df,
                             init.type = "Spectral", verbose = FALSE,
                             max.em.its = STM_MAXIT), error = function(e) NULL)
          t_stm <- (proc.time() - t0)[3L]
          if (!is.null(fs)) {
            Zstm <- log(pmax(fs$theta, 1e-10)) %*% Vilr
            mse_stm <- procrustes_align(std_Bz(Zstm, dat$C), B0std)$mse
            met <- rbind(met, data.frame(rep = r, method = "STM",
                                          mse_Bz = mse_stm, time_s = t_stm))
          }
        }
        if (r %% 5L == 0L) cat(sprintf("    rep %d/%d\n", r, N_REP_B3))
      }
      raw[[ci]] <- list(M = M, signal = snm, b_max = bmax, metrics = met)
    }
  }

  saveRDS(raw, file.path(DATA_DIR, "block3.rds"))
  B3 <- summarise_block3(raw)
  invisible(B3)
}

summarise_block3 <- function(raw) {
  df <- do.call(rbind, lapply(raw, function(x) {
    eg <- x$metrics[x$metrics$method == "sgscatm", ]
    st <- x$metrics[x$metrics$method == "STM", ]
    data.frame(
      M = x$M, Signal = x$signal,
      MSE_Spectral = if (nrow(eg)) mean(eg$mse_Bz) else NA,
      MSE_STM      = if (nrow(st)) mean(st$mse_Bz) else NA,
      Time_Spectral= if (nrow(eg)) mean(eg$time_s) else NA,
      Time_STM     = if (nrow(st)) mean(st$time_s) else NA,
      Speedup      = if (nrow(eg) && nrow(st)) mean(st$time_s)/mean(eg$time_s) else NA)
  }))
  write.csv(df, file.path(TAB_DIR, "block3.csv"), row.names = FALSE)
  cat("\n  Block 3 summary:\n"); print(round(df[,-2],4))

  body <- apply(df, 1L, function(row) {
    sp <- as.numeric(row["Speedup"])
    sprintf("%s & %s & %.4f & %s & %.2f & %s & %s",
            fmt_int(row["M"]), row["Signal"],
            as.numeric(row["MSE_Spectral"]),
            ifelse(is.na(row["MSE_STM"]), "---", sprintf("%.4f", as.numeric(row["MSE_STM"]))),
            as.numeric(row["Time_Spectral"]),
            ifelse(is.na(row["Time_STM"]), "---", sprintf("%.2f", as.numeric(row["Time_STM"]))),
            ifelse(is.na(sp), "---", sprintf("%.1f$\\times$", sp)))
  })
  cap <- sprintf(paste0(
    "Structural comparison of \\texttt{sgscatm} and STM on the corrected ",
    "standardized scale (each method's ILR scores standardized to unit ",
    "variance before OLS; identical treatment for ground truth). ",
    "$K=%d$, $P=%d$, $N=%d$, $\\sigma_\\varepsilon=0.3$, $\\lambda=1$, ",
    "%d replicates."),
    K_TOPICS, P_COV, N_VOCAB, N_REP_B3)
  if (REDUCED) cap <- paste(cap, "\\emph{Reduced replicate count (runtime).}")
  write_booktabs(
    file.path(TAB_DIR, "block3.tex"),
    header_cells = c("$M$", "Signal", "MSE Spec.", "MSE STM",
                     "Time Spec. (s)", "Time STM (s)", "Speedup"),
    body_rows = body, caption = cap, label = "tab:block3", colspec = "clccccc")
  cat("  Wrote tables/block3.csv, tables/block3.tex\n")

  make_block3_figures(raw, df)
  df
}

make_block3_figures <- function(raw, df) {
  ## boxplot of MSE by signal x M (log y) --------------------------
  all <- do.call(rbind, lapply(raw, function(x) {
    if (!nrow(x$metrics)) return(NULL)
    x$metrics$M <- x$M; x$metrics$signal <- x$signal; x$metrics }))
  path1 <- file.path(IMG_DIR, "block3_mse_boxplot.pdf")
  open_pdf(path1, 6.5, 3.6)
  par(mar = c(5.5, 4.5, 1, 1))
  all$grp <- paste(all$signal, all$M, all$method, sep = "\n")
  ord <- with(expand.grid(method = c("sgscatm","STM"),
                          M = sort(unique(all$M)),
                          signal = names(SIGNAL_B3)),
              paste(signal, M, method, sep = "\n"))
  all$grp <- factor(all$grp, levels = ord)
  cols <- ifelse(grepl("sgscatm$", levels(all$grp)),
                 unname(OI["blue"]), unname(OI["vermilion"]))
  boxplot(mse_Bz ~ grp, data = all, log = "y", col = cols, las = 2,
          cex.axis = 0.6, xlab = "", ylab = expression("MSE of"~hat(bold(B))[z]))
  legend("topright", c("sgscatm","STM"), fill = c(unname(OI["blue"]),
         unname(OI["vermilion"])), bty = "n", cex = 0.8)
  dev.off()

  ## mean time vs M (log y) ---------------------------------------
  path2 <- file.path(IMG_DIR, "block3_timing.pdf")
  open_pdf(path2, 5.5, 3.6)
  par(mar = c(4.2, 4.5, 1, 1))
  ts <- aggregate(time_s ~ M + method, data = all, FUN = mean)
  Ms <- sort(unique(ts$M))
  yr <- range(ts$time_s, na.rm = TRUE)
  plot(NA, xlim = range(Ms), ylim = yr, log = "xy",
       xlab = expression(italic(M)~"(corpus size)"),
       ylab = "mean computation time (s)")
  for (m in c("sgscatm","STM")) {
    d <- ts[ts$method == m, ]; d <- d[order(d$M), ]
    col <- if (m == "sgscatm") unname(OI["blue"]) else unname(OI["vermilion"])
    lines(d$M, d$time_s, type = "b", pch = 19, col = col, lwd = 2)
  }
  legend("topleft", c("sgscatm","STM"), col = c(unname(OI["blue"]),
         unname(OI["vermilion"])), lty = 1, pch = 19, lwd = 2, bty = "n")
  dev.off()
  cat("  Wrote imgs/block3_mse_boxplot.pdf, imgs/block3_timing.pdf\n")
}
