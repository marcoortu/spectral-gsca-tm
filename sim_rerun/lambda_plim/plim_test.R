#!/usr/bin/env Rscript
# ===================================================================
#  plim_test.R  —  one-shot ESTIMAND test: does the lambda regime
#                  change the probability limit of B_z?
# ===================================================================
#  WHY.  The recommended data-driven lambda_A grows proportionally with M
#  (lambda_A/M held constant), whereas Section 4's theory assumed lambda
#  fixed (lambda/M -> 0). In principle these two regimes have DIFFERENT
#  probability limits: lambda fixed -> the covariate block vanishes ->
#  pure-word (LSA) subspace estimand; lambda_A proportional to M -> the
#  covariate block stays commensurate -> covariate-regularized estimand.
#  This test measures how far apart those two plims actually are, at one
#  fixed huge M where sampling noise is negligible, so we learn whether the
#  regime choice matters for the estimand — and hence whether Section 4's
#  fixed-lambda variance can be applied to lambda_A.
#
#  No replicate loops. One large pilot corpus (streamed in chunks; only the
#  N x N word Gram and P x N / P x P cross-moments are accumulated, never a
#  dense M x N matrix). One spectral fit per lambda + a half-split floor.
# ===================================================================

options(warn = 1)
source("sim_rerun/R/common.R")          # package solver, V, closure, varimax, procrustes
K <- 5L; P <- 4L; N <- 500L; L <- 200L; Km1 <- K - 1L
OUT <- file.path(ROOT, "sim_rerun", "lambda_plim")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- DGP: reuse the clean Block-1 B0/R0 if present ------------------
b0p <- file.path(DATA_DIR, "B0.rds"); r0p <- file.path(DATA_DIR, "R0.rds")
d_target <- c(1.00, 0.70, 0.50, 0.35); se2 <- 0.15
if (file.exists(b0p)) {
  B0 <- readRDS(b0p)
} else {
  set.seed(20260704L)
  R0 <- qr.Q(qr(matrix(rnorm(P * P), P, P)))
  B0 <- R0 %*% diag(sqrt(d_target - se2)); attr(B0, "sigma_eps") <- sqrt(se2)
  saveRDS(B0, b0p); saveRDS(R0, r0p)
}
sig_eps <- attr(B0, "sigma_eps")
cat("eig(Cov(z)) =", round(eigen(crossprod(B0) + se2 * diag(P),
    only.values = TRUE)$values, 3), "  (target 1.00,0.70,0.50,0.35)\n")
Vhel <- ilr_contrast(K)

# --- streaming accumulator of sufficient statistics ----------------
# Accumulates, on document term-FREQUENCIES (row-normalised as sgscatm),
#   SWW = sum f_i f_i'  (N x N), sw = sum f_i, CtC_raw, SCW_raw, sc, M.
new_acc <- function() list(SWW = matrix(0, N, N), sw = numeric(N),
  CtC = matrix(0, P, P), SCW = matrix(0, P, N), sc = numeric(P), M = 0L)
add_chunk <- function(acc, W, C) {
  rs <- rowSums(W); rs[rs == 0] <- 1; Fr <- W / rs         # frequencies
  acc$SWW <- acc$SWW + crossprod(Fr)
  acc$sw  <- acc$sw  + colSums(Fr)
  acc$CtC <- acc$CtC + crossprod(C)
  acc$SCW <- acc$SCW + crossprod(C, Fr)
  acc$sc  <- acc$sc  + colSums(C)
  acc$M   <- acc$M   + nrow(W)
  acc
}
# centred Gram / cross-moments from an accumulator
finalize <- function(a) {
  M <- a$M; wbar <- a$sw / M; cbar <- a$sc / M
  list(GW  = a$SWW - M * tcrossprod(wbar),          # W~'W~  (N x N)
       CtC = a$CtC - M * tcrossprod(cbar),          # C'C    (P x P)
       SCW = a$SCW - M * outer(cbar, wbar),         # C'W~   (P x N)
       M = M)
}

# --- pilot size ----------------------------------------------------
M_PILOT <- 2000000L; CHUNK <- 100000L
n_chunks <- M_PILOT / CHUNK
cat(sprintf("Streaming pilot: M=%s in %d chunks of %s ...\n",
            format(M_PILOT, big.mark = ","), n_chunks, format(CHUNK, big.mark = ",")))
t0 <- proc.time()
h1 <- new_acc(); h2 <- new_acc()
for (ch in seq_len(n_chunks)) {
  dat <- sim_dgp(M = CHUNK, N = N, K = K, P = P, Bz0 = B0, sigma_eps = sig_eps,
                 alpha_beta = 0.1, doc_length = L, V = Vhel, seed = 5000L + ch)
  if (ch <= n_chunks / 2) h1 <- add_chunk(h1, dat$W, dat$C)
  else                    h2 <- add_chunk(h2, dat$W, dat$C)
  rm(dat); if (ch %% 5L == 0L) { gc(FALSE)
    cat(sprintf("  chunk %d/%d  (%.1f min)\n", ch, n_chunks, (proc.time()-t0)[3]/60)) }
}
S1 <- finalize(h1); S2 <- finalize(h2)
Sf <- finalize(list(SWW = h1$SWW + h2$SWW, sw = h1$sw + h2$sw,
                    CtC = h1$CtC + h2$CtC, SCW = h1$SCW + h2$SCW,
                    sc = h1$sc + h2$sc, M = h1$M + h2$M))
cat(sprintf("Pilot accumulated: M_used = %s  (%.1f min)\n",
            format(Sf$M, big.mark = ","), (proc.time() - t0)[3] / 60))

# --- dual spectral fit from sufficient statistics ------------------
# Never forms a dense M x N matrix: uses the N x N word Gram + the
# (N+P)-dim dual Gram.  Reuses varimax() and ilr_contrast() (Helmert V).
dual_fit <- function(S, lambda) {
  GW <- S$GW; CtC <- S$CtC; SCW <- S$SCW; M <- S$M
  Rc  <- chol(CtC)                                  # Rc' Rc = C'C
  CtCi <- chol2inv(Rc)
  QcW <- backsolve(Rc, SCW, transpose = TRUE)       # R^{-T} C'W~ = Q_C' W~  (P x N)
  # dual Gram G = [[GW, sqrt(l) W~'Q_C],[sqrt(l) Q_C'W~, l I]]
  G <- rbind(cbind(GW,               sqrt(lambda) * t(QcW)),
             cbind(sqrt(lambda) * QcW, lambda * diag(P)))
  eg <- eigen(G, symmetric = TRUE)
  s  <- pmax(eg$values[seq_len(Km1)], .Machine$double.eps)
  Fm <- eg$vectors[, seq_len(Km1), drop = FALSE]     # (N+P) x (K-1)
  # B_raw[,k] = sqrt(M/s_k) (C'C)^{-1} [C'W~ | sqrt(l) C'Q_C] f_k ;  C'Q_C = Rc'
  CtH <- cbind(SCW, sqrt(lambda) * t(Rc))            # P x (N+P)
  Braw <- CtH %*% Fm                                  # P x (K-1) (pre-scale)
  Braw <- sweep(Braw, 2L, sqrt(M / s), "*")
  Braw <- CtCi %*% Braw
  # varimax on Psi = K * Z*'W~ ;  Z*'W~ = diag(1/sqrt s) F' [GW ; sqrt(l) Q_C'W~]
  HtW <- rbind(GW, sqrt(lambda) * QcW)               # (N+P) x N
  Psi <- K * (diag(1 / sqrt(s), Km1) %*% (t(Fm) %*% HtW))
  Rvm <- varimax(t(Psi), normalize = FALSE)$rotmat
  Braw %*% Rvm
}
word_topmean <- function(S) mean(eigen(S$GW, symmetric = TRUE,
  only.values = TRUE)$values[seq_len(Km1)]) / S$M
lambda_A_of  <- function(S) eigen(S$GW, symmetric = TRUE,
  only.values = TRUE)$values[Km1]

# --- lambda grid ---------------------------------------------------
lamA <- lambda_A_of(Sf)
cat(sprintf("lambda_A(pilot) = %.4g   (word signal ~ %.3g)\n", lamA, word_topmean(Sf)))
grid <- c(`0` = 0, `0.01A` = 0.01 * lamA, `0.1A` = 0.1 * lamA,
          `1` = 1, `A` = lamA, `3A` = 3 * lamA, `10A` = 10 * lamA)

fits <- lapply(grid, function(l) dual_fit(Sf, l))
B_ref <- fits[["A"]]                                  # reference = lambda_A plim

# --- SE reference at M=4000 from block1.csv ------------------------
SE4000 <- tryCatch({
  b1 <- read.csv(file.path(TAB_DIR, "block1.csv"))
  b1$mean_analytic_se[which.min(abs(b1$M - 4000))]
}, error = function(e) NA_real_)
if (is.na(SE4000)) { SE4000 <- 0.0062; cat("NOTE: block1.csv SE missing; using 0.0062\n") }

# --- distances -----------------------------------------------------
dist_row <- function(name, l, B) {
  pa <- procrustes_align(B, B_ref)                    # align to lambda_A plim
  d_abs <- sqrt(sum((pa$Bz_aligned - B_ref)^2))
  d_entry <- d_abs / sqrt(P * Km1)
  data.frame(lambda_label = name, lambda = l,
             cov_block_eig = l / Sf$M, word_eig = word_topmean(Sf),
             d_abs = d_abs, d_entry = d_entry)
}
rows <- do.call(rbind, Map(dist_row, names(grid), grid, fits))

# --- noise floor (half-split) for lambda in {1, lambda_A} ----------
floor_entry <- function(lam1, lam2) {
  Bh1 <- dual_fit(S1, lam1); Bh2 <- dual_fit(S2, lam2)
  pa <- procrustes_align(Bh2, Bh1)
  sqrt(sum((pa$Bz_aligned - Bh1)^2)) / sqrt(P * Km1)
}
lamA1 <- lambda_A_of(S1); lamA2 <- lambda_A_of(S2)   # per-half commensurate lambda_A
floor_1 <- floor_entry(1, 1)
floor_A <- floor_entry(lamA1, lamA2)
f_mean  <- mean(c(floor_1, floor_A))

# --- finalize table columns ---------------------------------------
rows$d_entry_over_floor  <- rows$d_entry / f_mean
rows$d_entry_over_SE4000 <- rows$d_entry / SE4000
out <- rows[, c("lambda_label","lambda","cov_block_eig","word_eig",
                "d_abs","d_entry","d_entry_over_floor","d_entry_over_SE4000")]
floor_block <- data.frame(
  lambda_label = c("floor_entry(1)", "floor_entry(lambda_A)"),
  lambda = c(1, lamA), cov_block_eig = c(1 / Sf$M, lamA / Sf$M),
  word_eig = word_topmean(Sf), d_abs = NA,
  d_entry = c(floor_1, floor_A), d_entry_over_floor = NA, d_entry_over_SE4000 = NA)
csv_path <- file.path(OUT, "plim_distances.csv")
write.csv(rbind(out, floor_block), csv_path, row.names = FALSE)

cat("\n--- per-lambda distances to the lambda_A plim ---\n")
print(format(out, digits = 3))
cat(sprintf("\nnoise floor: floor_entry(1)=%.5f  floor_entry(A)=%.5f  mean=%.5f\n",
            floor_1, floor_A, f_mean))

# --- figure --------------------------------------------------------
fig_path <- file.path(OUT, "plim_sensitivity.pdf")
if (capabilities("cairo")) cairo_pdf(fig_path, 6.2, 4.2) else pdf(fig_path, 6.2, 4.2)
par(mar = c(4.4, 4.6, 1, 1))
pos <- out[out$cov_block_eig > 0, ]                   # x on log scale
plot(pos$cov_block_eig, pos$d_entry, log = "x", type = "b", pch = 19,
     col = unname(OI["blue"]), lwd = 2,
     xlab = expression("covariate-block eigenvalue"~lambda/M~"(log)"),
     ylab = "per-entry distance to  " , ylim = c(0, max(out$d_entry) * 1.1))
title(ylab = expression(hat(B)[z]^{(lambda[A])}~"plim"), line = 2.3)
abline(h = f_mean, lty = 2, col = "grey40", lwd = 1.5)
abline(h = SE4000, lty = 3, col = unname(OI["vermilion"]), lwd = 1.5)
lam1x <- out$cov_block_eig[out$lambda_label == "1"]
lamAx <- out$cov_block_eig[out$lambda_label == "A"]
abline(v = lam1x, col = "grey70"); abline(v = lamAx, col = "grey70")
text(lam1x, max(out$d_entry), expression(lambda == 1), pos = 4, cex = .85)
text(lamAx, max(out$d_entry) * 0.9, expression(lambda[A]), pos = 2, cex = .85)
legend("topright", bty = "n", cex = .85,
       legend = c("d_entry(lambda)", "noise floor", "SE(M=4000)"),
       col = c(unname(OI["blue"]), "grey40", unname(OI["vermilion"])),
       lty = c(1, 2, 3), pch = c(19, NA, NA), lwd = 2)
dev.off()

# --- VERDICT -------------------------------------------------------
h <- out$d_entry[out$lambda_label == "1"]             # fixed-lambda plim vs lambda_A
f <- f_mean; sref <- SE4000
verdict <- {
  if (h < max(2 * f, 0.2 * sref)) "ESTIMAND lambda-INSENSITIVE"
  else if (h < 1.0 * sref)        "ESTIMAND lambda-DEPENDENT (mild)"
  else                            "ESTIMAND lambda-DEPENDENT (strong)"
}
cat("\n=====================================================\n")
cat("            ESTIMAND VERDICT\n")
cat("=====================================================\n")
cat(sprintf("  d_entry(lambda=1) = %.5f\n", h))
cat(sprintf("  noise floor       = %.5f\n", f))
cat(sprintf("  SE(M=4000)        = %.5f\n", sref))
cat(sprintf("  h/floor  = %.2f\n", h / f))
cat(sprintf("  h/SE4000 = %.2f\n", h / sref))
cat(sprintf("  -> %s\n", verdict))
cat("\n  full d_entry(lambda) curve (does the estimand move with cov-block strength?):\n")
for (i in seq_len(nrow(out)))
  cat(sprintf("    lambda=%-10s cov_block=%.3g   d_entry=%.5f  (%.2f x floor, %.2f x SE4000)\n",
      out$lambda_label[i], out$cov_block_eig[i], out$d_entry[i],
      out$d_entry_over_floor[i], out$d_entry_over_SE4000[i]))
cat("=====================================================\n")
cat(sprintf("\nCSV:    %s\nFigure: %s\nPilot M used: %s\n",
            normalizePath(csv_path), normalizePath(fig_path),
            format(Sf$M, big.mark = ",")))
