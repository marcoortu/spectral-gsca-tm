# =============================================================================
# egscatm vs STM — Benchmark on FULL CORPUS (poliblog5k, 5000 documents)
#
# Covariates: rating (Conservative/Liberal), day (centred)
# Models  :  (1) egscatm  — default (rotate=TRUE)
#            (2) egscatm  — + M-step refinement (kmeans, temp=0.5)
#            (3) STM variational EM (reference)
# Focus   : TIMING on full corpus (5000 documents)
# =============================================================================

.pkg_root <- local({
  args  <- commandArgs(trailingOnly = FALSE)
  farg  <- grep("--file=", args, value = TRUE)
  if (length(farg) > 0L)
    dirname(dirname(normalizePath(sub("--file=", "", farg[1L]))))
  else
    getwd()
})
for (.f in list.files(file.path(.pkg_root, "R"), pattern = "[.]R$",
                      full.names = TRUE)) source(.f)
rm(.pkg_root, .f)
library(stm)
set.seed(2024)

# ---- 1. full corpus ---------------------------------------------------------
data(poliblog5k)
K <- 7L

out  <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta,
                      lower.thresh = 2L, verbose = FALSE)
docs  <- out$documents
vocab <- out$vocab
meta  <- out$meta
V     <- length(vocab)
M     <- length(docs)
cat(sprintf("Full corpus: %d documents x %d terms\n\n", M, V))

stm_to_dtm <- function(docs, V) {
  W <- matrix(0.0, length(docs), V)
  for (i in seq_along(docs)) { d <- docs[[i]]; W[i, d[1L,]] <- d[2L,] }
  W
}
cat("Building DTM... ")
W <- stm_to_dtm(docs, V)
cat(sprintf("done. sparsity: %.1f%%\n\n", 100 * mean(W == 0)))

C <- scale(cbind(rating = as.integer(meta$rating == "Liberal"),
                 day    = meta$day),
           center = TRUE, scale = FALSE)

# ---- 2. train/test split (80/20) -------------------------------------------
idx_test  <- sample.int(M, round(M * 0.2))
idx_train <- setdiff(seq_len(M), idx_test)
W_train <- W[idx_train, ]; W_test  <- W[idx_test, ]
C_train <- C[idx_train, ]; meta_train <- meta[idx_train, ]
docs_train <- docs[idx_train]
cat(sprintf("Train: %d  |  Test: %d\n\n", length(idx_train), length(idx_test)))

# =============================================================================
# FIT 1 — egscatm (rotate=TRUE, default)
# =============================================================================
cat(">>> egscatm  (rotate=TRUE, no refinement) ...\n")
t0    <- proc.time()
fit_eg <- egscatm(W_train, C_train, K = K, lambda = 3.0, scale_W = TRUE)
t_eg  <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s\n\n", t_eg))

# =============================================================================
# FIT 2 — egscatm + refine_phi (kmeans, temp=0.5)
# =============================================================================
cat(">>> egscatm  + refine_phi (kmeans, temp=0.5) ...\n")
t0      <- proc.time()
fit_rt  <- refine_phi(fit_eg, W_train, method = "kmeans",
                      smooth = 1e-4, temp = 0.5, seed = 2024L)
t_ref   <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s  (overhead)\n\n", t_ref))

# =============================================================================
# FIT 3 — STM (reference)
# =============================================================================
cat(">>> STM  (K=7, Spectral init, max.em.its=75) ...\n")
t0 <- proc.time()
fit_stm <- stm(documents = docs_train, vocab = vocab, K = K,
               prevalence = ~ rating + day, data = meta_train,
               init.type = "Spectral", max.em.its = 75L, verbose = FALSE)
t_stm <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s  (%d EM iter)\n\n", t_stm, fit_stm$convergence$its))

# =============================================================================
# METRICS
# =============================================================================
phi_eg  <- topic_word_dist(fit_eg)
phi_rt  <- topic_word_dist(fit_rt)
phi_stm <- exp(fit_stm$beta$logbeta[[1]])

sem_coherence <- function(phi, W_bin, n_top = 10L) {
  sapply(seq_len(nrow(phi)), function(k) {
    top <- order(phi[k,], decreasing = TRUE)[seq_len(n_top)]
    s   <- 0
    for (m in 2:n_top)
      for (l in 1:(m - 1L)) {
        D_ml <- sum(W_bin[, top[m]] & W_bin[, top[l]])
        s    <- s + log((D_ml + 1) / sum(W_bin[, top[l]]))
      }
    s
  })
}
frex_score <- function(phi, n_top = 10L, w = 0.7) {
  phi  <- pmax(phi, 1e-12)
  excl <- sweep(phi, 2L, colSums(phi), "/")
  sc   <- 1 / (w / excl + (1 - w) / phi)
  sapply(seq_len(nrow(phi)), function(k)
    mean(sc[k, order(sc[k,], decreasing = TRUE)[seq_len(n_top)]]))
}
excl_mean <- function(phi, n_top = 10L) {
  phi  <- pmax(phi, 1e-12)
  excl <- sweep(phi, 2L, colSums(phi), "/")
  sapply(seq_len(nrow(phi)), function(k)
    mean(excl[k, order(phi[k,], decreasing = TRUE)[seq_len(n_top)]]))
}
cosine_div <- function(A) {
  A <- A / sqrt(rowSums(A^2) + 1e-12)
  S <- tcrossprod(A); 1 - mean(S[upper.tri(S)])
}
held_out_ll <- function(phi_mat, W_test_raw) {
  W_t  <- W_test_raw / (rowSums(W_test_raw) + 1e-10)
  PtP  <- tcrossprod(phi_mat)
  th   <- tcrossprod(W_t, phi_mat) %*% solve(PtP + diag(1e-6, K))
  th   <- pmax(th, 0); th <- th / rowSums(th)
  sum(W_test_raw * log(pmax(th %*% phi_mat, 1e-300))) / sum(W_test_raw)
}

W_bin   <- W_train > 0
sc_eg   <- sem_coherence(phi_eg,  W_bin)
sc_rt   <- sem_coherence(phi_rt,  W_bin)
sc_stm  <- sem_coherence(phi_stm, W_bin)

frex_eg  <- frex_score(phi_eg);  frex_rt  <- frex_score(phi_rt)
frex_stm <- frex_score(phi_stm)

excl_eg  <- excl_mean(phi_eg);   excl_rt  <- excl_mean(phi_rt)
excl_stm <- excl_mean(phi_stm)

div_eg   <- cosine_div(phi_eg);  div_rt   <- cosine_div(phi_rt)
div_stm  <- cosine_div(phi_stm)

ll_eg    <- held_out_ll(phi_eg,  W_test)
ll_rt    <- held_out_ll(phi_rt,  W_test)
ll_stm   <- held_out_ll(phi_stm, W_test)

# =============================================================================
# OUTPUT
# =============================================================================
sep <- paste(rep("=", 72), collapse = "")
cat("\n", sep, "\n", sep = "")
cat("  BENCHMARK  egscatm / egscatm+refine+temp0.5 / STM\n")
cat(sprintf("  poliblog5k FULL  |  K=%d  |  M=%d  |  V=%d\n", K, M, V))
cat(sep, "\n\n")

# ---- timing -----------------------------------------------------------------
cat("TIMING\n")
cat(sprintf("  %-42s %8.2f s\n",  "egscatm (rotate=TRUE)",        t_eg))
cat(sprintf("  %-42s %8.2f s  (overhead +%.2f s)\n",
            "egscatm + refine_phi (kmeans, temp=0.5)",
            t_eg + t_ref, t_ref))
cat(sprintf("  %-42s %8.2f s  (%d EM iter)\n",
            "STM (variational EM)", t_stm, fit_stm$convergence$its))
cat(sprintf("  Speedup  egscatm / STM               : %.1fx\n",
            t_stm / t_eg))
cat(sprintf("  Speedup  egscatm+refine / STM        : %.1fx\n\n",
            t_stm / (t_eg + t_ref)))

# ---- aggregate metrics ------------------------------------------------------
cat("AGGREGATE METRICS (averaged over all topics)\n")
cat(sprintf("  %-36s %10s %10s %10s\n",
            "Metric", "eg-base", "eg-ref+t", "STM"))
cat(paste(rep("-", 70), collapse = ""), "\n")

rows <- list(
  list("Sem. Coherence (higher better)", "%10.2f",
       mean(sc_eg),   mean(sc_rt),   mean(sc_stm)),
  list("Exclusivity top-10 (higher)",   "%10.4f",
       mean(excl_eg), mean(excl_rt), mean(excl_stm)),
  list("FREX score (higher better)",    "%10.6f",
       mean(frex_eg), mean(frex_rt), mean(frex_stm)),
  list("Topic Diversity (higher better)","%10.4f",
       div_eg,        div_rt,        div_stm),
  list("Held-out log-lik/token (higher)","%10.4f",
       ll_eg,         ll_rt,         ll_stm)
)
for (r in rows) {
  fmt <- sprintf("  %%-36s %s %s %s\n", r[[2]], r[[2]], r[[2]])
  cat(sprintf(fmt, r[[1]], r[[3]], r[[4]], r[[5]]))
}
cat("\n")

# ---- gap vs STM + gap closed ------------------------------------------------
cat("GAP vs STM  (eg-ref+temp0.5 - STM)  and % gap closed by refinement\n")
eg_vals  <- c(mean(sc_eg),  mean(excl_eg),  mean(frex_eg),  div_eg,  ll_eg)
rt_vals  <- c(mean(sc_rt),  mean(excl_rt),  mean(frex_rt),  div_rt,  ll_rt)
stm_vals <- c(mean(sc_stm), mean(excl_stm), mean(frex_stm), div_stm, ll_stm)
lbl      <- c("Sem.Coherence","Exclusivity","FREX","Diversity","Held-out ll")

cat(sprintf("  %-20s %10s %10s %10s %10s\n",
            "Metric", "eg-base", "eg-ref+t", "STM", "gap closed"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_along(lbl)) {
  pct <- if (stm_vals[i] != eg_vals[i])
    100 * (rt_vals[i] - eg_vals[i]) / (stm_vals[i] - eg_vals[i]) else NA
  winner <- if (rt_vals[i] >= stm_vals[i]) " [eg>=STM]" else ""
  cat(sprintf("  %-20s %10.4f %10.4f %10.4f %9.0f%%%s\n",
              lbl[i], eg_vals[i], rt_vals[i], stm_vals[i], pct, winner))
}
cat("\n")

# ---- top-10 keywords per topic (eg-ref+temp vs STM) ------------------------
match_topics <- function(phi_a, phi_b) {
  na <- phi_a / sqrt(rowSums(phi_a^2) + 1e-12)
  nb <- phi_b / sqrt(rowSums(phi_b^2) + 1e-12)
  S  <- na %*% t(nb); assigned <- logical(K); p <- integer(K)
  for (k in seq_len(K)) {
    best <- which.max(S[k,] * !assigned); p[k] <- best; assigned[best] <- TRUE
  }
  p
}
perm_stm     <- match_topics(phi_rt, phi_stm)
phi_stm_al   <- phi_stm[perm_stm, ]

tt_eg  <- top_terms(fit_eg,  n = 10L, vocab = vocab)
tt_rt  <- top_terms(fit_rt,  n = 10L, vocab = vocab)
stm_top <- apply(phi_stm_al, 1L, function(p)
  paste(vocab[order(p, decreasing = TRUE)[1:10]], collapse = " | "))

cat(sep, "\n")
cat("  TOP-10 KEYWORDS PER TOPIC\n")
cat(sep, "\n\n")

prev_eg  <- colMeans(fit_eg$Pi)
prev_stm <- colMeans(fit_stm$theta)[perm_stm]

for (k in seq_len(K)) {
  cat(sprintf("Topic %d  |  prev_eg=%.3f  prev_stm=%.3f\n",
              k, prev_eg[k], prev_stm[k]))
  cat(sprintf("  eg-base : %s\n", paste(tt_eg[k,],  collapse = " | ")))
  cat(sprintf("  eg-rt   : %s\n", paste(tt_rt[k,],  collapse = " | ")))
  cat(sprintf("  STM     : %s\n", stm_top[k]))
  cos_sim <- sum(phi_rt[k,] * phi_stm_al[k,]) /
    (sqrt(sum(phi_rt[k,]^2)) * sqrt(sum(phi_stm_al[k,]^2)))
  cat(sprintf("  cos-sim (eg-rt, stm): %.3f\n\n", cos_sim))
}

# ---- final summary ----------------------------------------------------------
cat(sep, "\n")
cat("  FINAL SUMMARY\n")
cat(sep, "\n")
cat(sprintf("  Total documents               : %d\n", M))
cat(sprintf("  Vocabulary size               : %d\n", V))
cat(sprintf("  K                             : %d\n", K))
cat(sprintf("\n  Time egscatm (fit only)       : %.2f s\n", t_eg))
cat(sprintf("  Time refine_phi overhead      : %.2f s\n", t_ref))
cat(sprintf("  Time egscatm+refine TOTAL     : %.2f s\n", t_eg + t_ref))
cat(sprintf("  Time STM                      : %.2f s  (%d iter)\n",
            t_stm, fit_stm$convergence$its))
cat(sprintf("\n  Speedup egscatm / STM         : %.1fx\n", t_stm / t_eg))
cat(sprintf("  Speedup egscatm+refine / STM  : %.1fx\n\n", t_stm / (t_eg + t_ref)))
