# =============================================================================
# egscatm vs STM — Direct comparison
#
# Dataset : poliblog5k (1000 sampled documents, US political blogs 2008)
# Covariates: rating (Conservative/Liberal), day (centred)
# Models  :  (1) egscatm  — default (rotate=TRUE)
#            (2) egscatm  — + M-step refinement of Phi
#            (3) egscatm  — + refinement + temperature sharpening (tau=0.5)
#            (4) STM variational EM (reference)
# Metrics : timing, semantic coherence, exclusivity, FREX,
#           topic diversity, held-out log-lik
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

# ---- 1. data ----------------------------------------------------------------
data(poliblog5k)
N_sample <- 1000L; K <- 7L

idx      <- sample.int(5000L, N_sample)
out      <- prepDocuments(poliblog5k.docs[idx], poliblog5k.voc,
                          poliblog5k.meta[idx,],
                          lower.thresh = 2L, verbose = FALSE)
docs  <- out$documents; vocab <- out$vocab; meta <- out$meta
V_sub <- length(vocab);  M    <- length(docs)
cat(sprintf("Corpus: %d doc x %d terms\n\n", M, V_sub))

stm_to_dtm <- function(docs, V) {
  W <- matrix(0.0, length(docs), V)
  for (i in seq_along(docs)) { d <- docs[[i]]; W[i,d[1L,]] <- d[2L,] }
  W
}
W <- stm_to_dtm(docs, V_sub)
C <- scale(cbind(rating = as.integer(meta$rating == "Liberal"),
                 day    = meta$day), center = TRUE, scale = FALSE)

# ---- 2. train/test split (80/20) -------------------------------------------
idx_test  <- sample.int(M, round(M * 0.2))
idx_train <- setdiff(seq_len(M), idx_test)
docs_train <- docs[idx_train]; docs_test  <- docs[idx_test]
W_train <- W[idx_train,];      W_test    <- W[idx_test,]
C_train <- C[idx_train,];      meta_train<- meta[idx_train,]
cat(sprintf("Train: %d  |  Test: %d\n\n", length(idx_train), length(idx_test)))

# =============================================================================
# FIT 1 — egscatm (rotate=TRUE, default)
# =============================================================================
cat(">>> egscatm  (rotate=TRUE, no refinement) ...\n")
t0 <- proc.time()
fit_eg <- egscatm(W_train, C_train, K = K, lambda = 3.0, scale_W = TRUE)
t_eg   <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s\n\n", t_eg))

# =============================================================================
# FIT 2 — egscatm + refine_phi  (k-means, temp=1)
# =============================================================================
cat(">>> egscatm  + refine_phi (kmeans, temp=1) ...\n")
t0      <- proc.time()
fit_ref <- refine_phi(fit_eg, W_train, method = "kmeans",
                      smooth = 1e-4, temp = 1.0, seed = 2024L)
t_ref   <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s  (overhead over base fit)\n\n", t_ref))

# =============================================================================
# FIT 3 — egscatm + refine_phi + temperature sharpening (tau=0.5)
# =============================================================================
cat(">>> egscatm  + refine_phi + temp=0.5 ...\n")
fit_rt <- refine_phi(fit_eg, W_train, method = "kmeans",
                     smooth = 1e-4, temp = 0.5, seed = 2024L)
cat("    Done (same k-means assignment, only temp differs)\n\n")

# =============================================================================
# FIT 4 — STM (reference)
# =============================================================================
cat(">>> STM  (K=7, Spectral init, max.em.its=75) ...\n")
t0 <- proc.time()
fit_stm <- stm(documents = docs_train, vocab = vocab, K = K,
               prevalence = ~ rating + day, data = meta_train,
               init.type = "Spectral", max.em.its = 75L, verbose = FALSE)
t_stm <- (proc.time() - t0)["elapsed"]
cat(sprintf("    %.2f s  (%d EM iter)\n\n", t_stm, fit_stm$convergence$its))

# =============================================================================
# TOPIC-WORD MATRICES
# =============================================================================
phi_eg  <- topic_word_dist(fit_eg)    # softmax(V*Psi)
phi_ref <- topic_word_dist(fit_ref)   # M-step
phi_rt  <- topic_word_dist(fit_rt)    # M-step + temp=0.5
phi_stm <- exp(fit_stm$beta$logbeta[[1]])

# =============================================================================
# METRICS
# =============================================================================

sem_coherence <- function(phi, W_bin, n_top = 10L) {
  sapply(seq_len(nrow(phi)), function(k) {
    top <- order(phi[k,], decreasing = TRUE)[seq_len(n_top)]
    s   <- 0
    for (m in 2:n_top)
      for (l in 1:(m-1L)) {
        D_ml <- sum(W_bin[,top[m]] & W_bin[,top[l]])
        s    <- s + log((D_ml + 1) / sum(W_bin[,top[l]]))
      }
    s
  })
}
frex_score <- function(phi, n_top = 10L, w = 0.7) {
  phi  <- pmax(phi, 1e-12)
  excl <- sweep(phi, 2L, colSums(phi), "/")
  sc   <- 1 / (w / excl + (1-w) / phi)
  sapply(seq_len(nrow(phi)), function(k)
    mean(sc[k, order(sc[k,], decreasing=TRUE)[seq_len(n_top)]]))
}
excl_mean <- function(phi, n_top = 10L) {
  phi  <- pmax(phi, 1e-12)
  excl <- sweep(phi, 2L, colSums(phi), "/")
  sapply(seq_len(nrow(phi)), function(k)
    mean(excl[k, order(phi[k,], decreasing=TRUE)[seq_len(n_top)]]))
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

W_bin  <- W_train > 0
sc_eg  <- sem_coherence(phi_eg,  W_bin)
sc_ref <- sem_coherence(phi_ref, W_bin)
sc_rt  <- sem_coherence(phi_rt,  W_bin)
sc_stm <- sem_coherence(phi_stm, W_bin)

frex_eg  <- frex_score(phi_eg);  frex_ref <- frex_score(phi_ref)
frex_rt  <- frex_score(phi_rt);  frex_stm <- frex_score(phi_stm)

excl_eg  <- excl_mean(phi_eg);   excl_ref <- excl_mean(phi_ref)
excl_rt  <- excl_mean(phi_rt);   excl_stm <- excl_mean(phi_stm)

div_eg   <- cosine_div(phi_eg);  div_ref  <- cosine_div(phi_ref)
div_rt   <- cosine_div(phi_rt);  div_stm  <- cosine_div(phi_stm)

ll_eg    <- held_out_ll(phi_eg,  W_test); ll_ref <- held_out_ll(phi_ref, W_test)
ll_rt    <- held_out_ll(phi_rt,  W_test); ll_stm <- held_out_ll(phi_stm, W_test)

# =============================================================================
# OUTPUT
# =============================================================================
sep <- paste(rep("=", 72), collapse = "")
cat("\n", sep, "\n", sep = "")
cat("  COMPARISON  egscatm / +refine / +refine+temp / STM\n")
cat("  poliblog5k  |  K=7  |  train N=", length(idx_train), "\n", sep = "")
cat(sep, "\n\n")

# ---- timing -----------------------------------------------------------------
cat("TIMING\n")
cat(sprintf("  %-38s %8.2f s\n",  "egscatm (rotate=TRUE)",      t_eg))
cat(sprintf("  %-38s %8.2f s  (+%.3f s overhead)\n",
            "egscatm + refine_phi",
            t_eg + t_ref, t_ref))
cat(sprintf("  %-38s %8.2f s  (same as above)\n",
            "egscatm + refine_phi + temp=0.5",
            t_eg + t_ref))
cat(sprintf("  %-38s %8.2f s  (%d EM iter)\n",
            "STM (variational EM)", t_stm, fit_stm$convergence$its))
cat(sprintf("  Speedup  egscatm / STM          : %.1fx\n",
            t_stm / t_eg))
cat(sprintf("  Speedup  egscatm+refine / STM   : %.1fx\n\n",
            t_stm / (t_eg + t_ref)))

# ---- aggregate metrics ------------------------------------------------------
cat("AGGREGATE METRICS (averaged over all topics)\n")
cat(sprintf("  %-34s %9s %9s %9s %9s\n",
            "Metric", "eg-base", "eg-ref", "eg-ref+t", "STM"))
cat(paste(rep("-", 75), collapse=""), "\n")

rows <- list(
  list("Sem. Coherence (higher better)", "%9.2f",
       mean(sc_eg),   mean(sc_ref),   mean(sc_rt),   mean(sc_stm)),
  list("Exclusivity top-10 (higher)",    "%9.4f",
       mean(excl_eg), mean(excl_ref), mean(excl_rt), mean(excl_stm)),
  list("FREX score (higher better)",     "%9.6f",
       mean(frex_eg), mean(frex_ref), mean(frex_rt), mean(frex_stm)),
  list("Topic Diversity (higher better)","%9.4f",
       div_eg,        div_ref,        div_rt,        div_stm),
  list("Held-out log-lik/token (higher)","%9.4f",
       ll_eg,         ll_ref,         ll_rt,         ll_stm)
)
for (r in rows) {
  fmt <- sprintf("  %%-34s %s %s %s %s\n", r[[2]], r[[2]], r[[2]], r[[2]])
  cat(sprintf(fmt, r[[1]], r[[3]], r[[4]], r[[5]], r[[6]]))
}
cat("\n")

# ---- improvement from refinement over baseline ------------------------------
cat("IMPROVEMENT FROM REFINEMENT  (eg-ref - eg-base)\n")
deltas <- c(
  `Sem.Coherence` = mean(sc_ref)   - mean(sc_eg),
  Exclusivity     = mean(excl_ref) - mean(excl_eg),
  FREX            = mean(frex_ref) - mean(frex_eg),
  Diversity       = div_ref        - div_eg,
  `Held-out ll`   = ll_ref         - ll_eg
)
for (nm in names(deltas)) {
  flag <- if (deltas[nm] > 0) "+" else "-"
  cat(sprintf("  %-20s : %+.4f  [%s]\n",
              nm, deltas[nm],
              if (deltas[nm] > 0) "improved" else "worsened"))
}
cat("\n")

cat("IMPROVEMENT FROM TEMP=0.5  (eg-ref+temp - eg-ref)\n")
deltas_t <- c(
  `Sem.Coherence` = mean(sc_rt)   - mean(sc_ref),
  Exclusivity     = mean(excl_rt) - mean(excl_ref),
  FREX            = mean(frex_rt) - mean(frex_ref),
  Diversity       = div_rt        - div_ref,
  `Held-out ll`   = ll_rt         - ll_ref
)
for (nm in names(deltas_t)) {
  cat(sprintf("  %-20s : %+.4f  [%s]\n",
              nm, deltas_t[nm],
              if (deltas_t[nm] > 0) "improved" else "worsened"))
}
cat("\n")

# ---- residual gap vs STM ----------------------------------------------------
cat("RESIDUAL GAP vs STM  (eg-ref+temp - STM)\n")
gap <- c(
  `Sem.Coherence` = mean(sc_rt)   - mean(sc_stm),
  Exclusivity     = mean(excl_rt) - mean(excl_stm),
  FREX            = mean(frex_rt) - mean(frex_stm),
  Diversity       = div_rt        - div_stm,
  `Held-out ll`   = ll_rt         - ll_stm
)
for (nm in names(gap)) {
  cat(sprintf("  %-20s : %+.4f  [%s]\n",
              nm, gap[nm],
              if (gap[nm] >= 0) "egscatm >= STM" else "STM better"))
}
cat("\n")

# ---- per-topic metrics (ref+temp vs STM) ------------------------------------
match_topics <- function(phi_a, phi_b) {
  na <- phi_a / sqrt(rowSums(phi_a^2)+1e-12)
  nb <- phi_b / sqrt(rowSums(phi_b^2)+1e-12)
  S  <- na %*% t(nb); assigned <- logical(K); p <- integer(K)
  for (k in seq_len(K)) {
    best <- which.max(S[k,]*!assigned); p[k] <- best; assigned[best] <- TRUE
  }
  p
}
perm_stm <- match_topics(phi_rt, phi_stm)

cat("PER-TOPIC METRICS  (eg-ref+temp vs STM aligned)\n")
cat(sprintf("  %-7s %10s %10s %10s %10s\n",
            "Topic","SC-eg","SC-stm","FREX-eg","FREX-stm"))
cat(paste(rep("-", 52), collapse=""), "\n")
for (k in seq_len(K)) {
  cat(sprintf("  T%-6d %10.2f %10.2f %10.6f %10.6f\n",
              k, sc_rt[k], sc_stm[perm_stm[k]],
              frex_rt[k], frex_stm[perm_stm[k]]))
}
cat("\n")

# ---- top-10 keywords per topic (all 4 models) ------------------------------
cat(sep, "\n")
cat("  TOP-10 KEYWORDS PER TOPIC\n")
cat(sep, "\n\n")

tt_eg  <- top_terms(fit_eg,  n=10L, vocab=vocab)
tt_ref <- top_terms(fit_ref, n=10L, vocab=vocab)
tt_rt2 <- top_terms(fit_rt,  n=10L, vocab=vocab)
phi_stm_al <- phi_stm[perm_stm,]
stm_top <- apply(phi_stm_al, 1L, function(p)
  paste(vocab[order(p, decreasing=TRUE)[1:10]], collapse=" | "))

prev_eg <- colMeans(fit_eg$Pi)
prev_stm <- colMeans(fit_stm$theta)[perm_stm]

for (k in seq_len(K)) {
  cat(sprintf("Topic %d  |  prev_eg=%.3f  prev_stm=%.3f\n",
              k, prev_eg[k], prev_stm[k]))
  cat(sprintf("  eg-base : %s\n", paste(tt_eg[k,],  collapse=" | ")))
  cat(sprintf("  eg-ref  : %s\n", paste(tt_ref[k,], collapse=" | ")))
  cat(sprintf("  eg-rt   : %s\n", paste(tt_rt2[k,], collapse=" | ")))
  cat(sprintf("  STM     : %s\n", stm_top[k]))
  cos_rt_stm <- sum(phi_rt[k,]*phi_stm_al[k,]) /
    (sqrt(sum(phi_rt[k,]^2))*sqrt(sum(phi_stm_al[k,]^2)))
  cat(sprintf("  cos-sim (eg-rt, stm): %.3f\n\n", cos_rt_stm))
}

# ---- rating effect ----------------------------------------------------------
cat(sep, "\n")
cat("  RATING EFFECT  (Liberal - Conservative, delta prevalence)\n")
cat(sep, "\n\n")
delta_eg  <- colMeans(fit_eg$Pi[meta_train$rating=="Liberal",]) -
             colMeans(fit_eg$Pi[meta_train$rating=="Conservative",])
delta_stm <- colMeans(fit_stm$theta[meta_train$rating=="Liberal",]) -
             colMeans(fit_stm$theta[meta_train$rating=="Conservative",])
cat(sprintf("  %-7s %12s %12s\n", "Topic", "egscatm", "STM"))
cat(paste(rep("-", 34), collapse=""), "\n")
for (k in seq_len(K))
  cat(sprintf("  T%-6d %+12.4f %+12.4f\n",
              k, delta_eg[k], delta_stm[perm_stm[k]]))
cat("\n")

# ---- final summary ----------------------------------------------------------
cat(sep, "\n")
cat("  FINAL SUMMARY\n")
cat(sep, "\n")
cat(sprintf("  Speedup egscatm / STM             : %.1fx\n", t_stm/t_eg))
cat(sprintf("  Speedup egscatm+refine / STM      : %.1fx\n", t_stm/(t_eg+t_ref)))
cat(sprintf("  refine_phi overhead               : %.3f s\n\n", t_ref))
cat("  Improvement eg-base -> eg-ref+temp (vs STM):\n")
metrics_lbl <- c("Sem.Coherence","Exclusivity","FREX","Diversity","Held-out ll")
eg_vals  <- c(mean(sc_eg),  mean(excl_eg),  mean(frex_eg),  div_eg,  ll_eg)
rt_vals  <- c(mean(sc_rt),  mean(excl_rt),  mean(frex_rt),  div_rt,  ll_rt)
stm_vals <- c(mean(sc_stm), mean(excl_stm), mean(frex_stm), div_stm, ll_stm)
for (i in seq_along(metrics_lbl)) {
  pct_gap_closed <- if (stm_vals[i] != eg_vals[i])
    100 * (rt_vals[i] - eg_vals[i]) / (stm_vals[i] - eg_vals[i]) else NA
  cat(sprintf("    %-20s : eg=%.4f  rt=%.4f  stm=%.4f  gap closed: %.0f%%\n",
              metrics_lbl[i], eg_vals[i], rt_vals[i], stm_vals[i],
              pct_gap_closed))
}
cat("\n")
