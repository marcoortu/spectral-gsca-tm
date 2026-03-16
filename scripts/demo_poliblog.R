# =============================================================================
# sgscatm demo: Political Blog Posts (poliblog5k, stm package)
#
# Dataset: 5000 US political blog posts from the 2008 presidential campaign.
#   - Covariates: rating (Conservative / Liberal), day (day of campaign)
#   - Pre-processed vocabulary (stemmed, stopwords removed): 2632 terms
#
# We sample 1000 documents and compare three configurations:
#   (A) sgscatm  — default (varimax rotation)
#   (B) sgscatm  — + M-step refinement of Phi
#   (C) sgscatm  — + M-step refinement + temperature sharpening (tau = 0.5)
# =============================================================================

# ---- 0. source package functions --------------------------------------------
# Locate the package root: works when called via Rscript and interactively.
.pkg_root <- local({
  args  <- commandArgs(trailingOnly = FALSE)
  farg  <- grep("--file=", args, value = TRUE)
  if (length(farg) > 0L)
    dirname(dirname(normalizePath(sub("--file=", "", farg[1L]))))
  else
    getwd()   # interactive: set working directory to package root
})
for (.f in list.files(file.path(.pkg_root, "R"), pattern = "[.]R$",
                      full.names = TRUE)) source(.f)
rm(.pkg_root, .f)

# ---- 1. load data -----------------------------------------------------------
library(stm)
data(poliblog5k)

vocab <- poliblog5k.voc
V     <- length(vocab)
M_all <- length(poliblog5k.docs)
cat(sprintf("Full corpus: %d documents, %d terms\n", M_all, V))

# ---- 2. sample 1000 documents -----------------------------------------------
set.seed(2024)
N_sample <- 1000L
idx      <- sample.int(M_all, N_sample)
docs_sub <- poliblog5k.docs[idx]
meta_sub <- poliblog5k.meta[idx, ]
cat(sprintf("Sample: %d documents\n", N_sample))
cat("Rating distribution:\n"); print(table(meta_sub$rating))

# ---- 3. dense DTM -----------------------------------------------------------
stm_to_dtm <- function(docs, V) {
  W <- matrix(0.0, length(docs), V)
  for (i in seq_along(docs)) { d <- docs[[i]]; W[i, d[1L,]] <- d[2L,] }
  W
}
cat("Building DTM... ")
W <- stm_to_dtm(docs_sub, V)
cat(sprintf("done. %d x %d,  sparsity: %.1f%%\n", nrow(W), ncol(W),
            100 * mean(W == 0)))

# ---- 4. covariate matrix ----------------------------------------------------
C <- scale(cbind(rating = as.integer(meta_sub$rating == "Liberal"),
                 day    = meta_sub$day),
           center = TRUE, scale = FALSE)

# ---- 5. fit -----------------------------------------------------------------
K <- 7L; lambda <- 3.0

cat(sprintf("\nFitting sgscatm: K=%d, lambda=%.1f, rotate=TRUE (default) ...\n",
            K, lambda))
t0  <- proc.time()
fit <- sgscatm(W = W, C = C, K = K, lambda = lambda, scale_W = TRUE)
t_fit <- (proc.time() - t0)["elapsed"]
cat(sprintf("  Done in %.2f s\n", t_fit))

# ---- 6. M-step refinement of Phi — k-means (default) -----------------------
cat("Applying refine_phi (method=kmeans, temp=1) ...\n")
t0    <- proc.time()
fit_r <- refine_phi(fit, W, method = "kmeans", smooth = 1e-4,
                    temp = 1.0, seed = 42L)
t_ref <- (proc.time() - t0)["elapsed"]
cat(sprintf("  Done in %.2f s\n", t_ref))

# ---- 7. k-means + temperature sharpening (Option 1+2) ----------------------
fit_rt <- refine_phi(fit, W, method = "kmeans", smooth = 1e-4,
                     temp = 0.5, seed = 42L)

# ---- 8. extract topic-word distributions for metrics -----------------------
phi_base <- topic_word_dist(fit)          # softmax(V*Psi)
phi_ref  <- topic_word_dist(fit_r)        # M-step
phi_rt   <- topic_word_dist(fit_rt)       # M-step + temp=0.5

# ---- 9. model header --------------------------------------------------------
sep <- paste(rep("=", 70), collapse = "")
cat("\n", sep, "\n", sep = "")
cat(sprintf("  sgscatm  |  K=%d  |  M=%d  |  N=%d  |  lambda=%.1f\n",
            K, nrow(W), ncol(W), lambda))
cat(sep, "\n"); print(fit_r)

cat("\nTop K-1 eigenvalues of S_z (invariant under rotation and refinement):\n")
print(round(fit$eigenvalues, 3))
cat(sprintf("Top %d eigenvalues: %.1f%% of tr(S_z)\n\n",
            K-1L, 100*sum(fit$eigenvalues)/sum(fit$eigenvalues_all)))

# ---- 10. topic diversity comparison -----------------------------------------
cosine_div <- function(A) {
  A <- A / sqrt(rowSums(A^2) + 1e-12)
  S <- tcrossprod(A); 1 - mean(S[upper.tri(S)])
}
excl_mean <- function(phi, n=10L) {
  phi  <- pmax(phi, 1e-12)
  excl <- sweep(phi, 2L, colSums(phi), "/")
  mean(sapply(seq_len(nrow(phi)), function(k)
    mean(excl[k, order(phi[k,], decreasing=TRUE)[seq_len(n)]])))
}

cat("QUALITY METRICS (averaged over topics)\n")
cat(sprintf("  %-35s %10s %10s %10s\n",
            "Metric", "base", "refined", "ref+temp"))
cat(paste(rep("-", 68), collapse=""), "\n")
cat(sprintf("  %-35s %10.4f %10.4f %10.4f\n", "Topic Diversity (higher better)",
            cosine_div(phi_base), cosine_div(phi_ref), cosine_div(phi_rt)))
cat(sprintf("  %-35s %10.4f %10.4f %10.4f\n", "Mean Exclusivity top-10 (higher)",
            excl_mean(phi_base), excl_mean(phi_ref), excl_mean(phi_rt)))
cat("\n")

# ---- 11. top-10 keywords per topic ------------------------------------------
n_top <- 10L
tt_b  <- top_terms(fit,    n = n_top, vocab = vocab)
tt_r  <- top_terms(fit_r,  n = n_top, vocab = vocab)
tt_rt <- top_terms(fit_rt, n = n_top, vocab = vocab)

cat(sep, "\n")
cat(sprintf("  TOP %d KEYWORDS PER TOPIC\n", n_top))
cat(sep, "\n\n")

for (k in seq_len(K)) {
  cat(sprintf("Topic %d  (mean prevalence: %.1f%%)\n",
              k, round(colMeans(fit$Pi)[k] * 100, 1)))
  cat(sprintf("  base    : %s\n", paste(tt_b[k,],  collapse = " | ")))
  cat(sprintf("  refined : %s\n", paste(tt_r[k,],  collapse = " | ")))
  cat(sprintf("  ref+tmp : %s\n\n", paste(tt_rt[k,], collapse = " | ")))
}

# ---- 12. topic prevalence by political leaning (refined fit) ----------------
cat(sep, "\n")
cat("  TOPIC PREVALENCE BY POLITICAL LEANING  (refined fit)\n")
cat(sep, "\n\n")

for (k in seq_len(K)) {
  pi_cons <- mean(fit_r$Pi[meta_sub$rating == "Conservative", k])
  pi_lib  <- mean(fit_r$Pi[meta_sub$rating == "Liberal",       k])
  lean    <- if (pi_lib > pi_cons) "Liberal >" else "Conservative >"
  cat(sprintf("Topic %d  |  Conservative: %.3f  Liberal: %.3f  [%s]\n",
              k, pi_cons, pi_lib, lean))
}
cat("\n")

# ---- 13. path coefficients --------------------------------------------------
cat("Path coefficients Bz  [P x (K-1)]  (unchanged by refinement)\n")
Bz_df <- as.data.frame(round(fit$Bz, 4))
rownames(Bz_df) <- c("rating (Liberal)", "day")
colnames(Bz_df) <- paste0("rILR", seq_len(K-1L))
print(Bz_df); cat("\n")

# ---- 14. timing recap -------------------------------------------------------
cat(sep, "\n")
cat("  TIMING SUMMARY\n")
cat(sep, "\n")
cat(sprintf("  sgscatm fit (rotate=TRUE) : %.2f s\n", t_fit))
cat(sprintf("  refine_phi (M-step)       : %.2f s\n", t_ref))
cat(sprintf("  Total                     : %.2f s\n", t_fit + t_ref))
cat(sprintf("  Documents: %d  |  Terms: %d  |  Topics: %d\n", N_sample, V, K))
