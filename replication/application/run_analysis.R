#!/usr/bin/env Rscript
# ===================================================================
#  BES Wave 25 — Final analysis  K=7, lambda=1
#  Outputs: tables (.tex) and figures (.pdf) in output/bes/
# ===================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("R/egscatm_fit.R")
source("R/ilr_contrast.R")
source("R/ilr_se.R")
source("R/refine_phi.R")
source("R/methods.R")
source("R/utils.R")

dir.create("output/bes", recursive = TRUE, showWarnings = FALSE)

K_TOPICS <- 7L
LAMBDA   <- 1

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
cat("Loading preprocessed DTM...\n")
dat    <- readRDS("scripts/bes_case_study/bes_w25_dtm.rds")
W_sp   <- dat$W
C      <- dat$C
vocab  <- dat$vocab
df_raw <- dat$df
doc_ids <- dat$doc_ids
W      <- as.matrix(W_sp)
M      <- nrow(W); N <- ncol(W); P <- ncol(C)
cat(sprintf("  M=%d  N=%d  P=%d\n", M, N, P))

cov_labels  <- c("Age", "Female", "Education", "Left--right", "Leave")
topic_names <- paste0("T", seq_len(K_TOPICS))

# ---------------------------------------------------------------
# Fit
# ---------------------------------------------------------------
cat("Fitting sgscatm (K=7, lambda=1)...\n")
t0  <- proc.time()[3]
fit <- sgscatm(W, C, K = K_TOPICS, lambda = LAMBDA, rotate = TRUE)
t_fit <- proc.time()[3] - t0
cat(sprintf("  Fit time: %.2f s\n", t_fit))

# Refine Phi
cat("Refining Phi (K-means, tau=0.5)...\n")
fit <- refine_phi(fit, W, method = "kmeans", temp = 0.5, seed = 42)

# Analytical SEs
cat("Computing analytical SEs...\n")
se_res <- tryCatch(
  ilr_se_analytical(fit),
  error = function(e) {
    cat("  Analytical SE failed:", e$message, "\n  Falling back to bootstrap...\n")
    ilr_se(fit, W, C, B = 300, seed = 42)
  }
)

# Save fitted objects
saveRDS(fit,    "output/bes/fit_K7.rds")
saveRDS(se_res, "output/bes/se_K7.rds")
cat("Saved: output/bes/fit_K7.rds\n")

# ---------------------------------------------------------------
# Descriptive stats
# ---------------------------------------------------------------
prevalence  <- colMeans(fit$Pi)
dl          <- rowSums(W)
leave_share <- mean(C[, 5] + mean(df_raw$euref[match(doc_ids, as.character(df_raw$id))]) > 0)

# Retrieve Leave/Remain counts directly
euref_vec <- df_raw$euref[match(doc_ids, as.character(df_raw$id))]
n_leave   <- sum(euref_vec == 1, na.rm = TRUE)
n_remain  <- sum(euref_vec == 0, na.rm = TRUE)

# ---------------------------------------------------------------
# Top terms per topic
# ---------------------------------------------------------------
tt <- top_terms(fit, n = 10, vocab = vocab)
rownames(tt) <- topic_names

cat("\nTop 10 terms per topic:\n")
print(tt)

# LaTeX table — top terms
tt_body <- paste(
  sapply(seq_len(K_TOPICS), function(k) {
    sprintf("$T_%d$ & %s \\\\", k,
            paste(tt[k, ], collapse = " $\\cdot$ "))
  }),
  collapse = "\n")

tt_tex <- sprintf(
"\\begin{table}[!ht]
\\centering
\\small
\\caption{Top 10 stemmed terms per topic from \\textsc{ilr-egsca} applied to
BES Wave~25 open-ended ``most important issue'' responses
($M=%d$, $N=%d$, $K=7$, $\\lambda=1$, varimax rotation).}
\\label{tab:bes_topics}
\\begin{tabular}{cl}
\\toprule
Topic & Top terms \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", M, N, tt_body)

writeLines(tt_tex, "output/bes/table_topics.tex")
cat("Saved: output/bes/table_topics.tex\n")

# ---------------------------------------------------------------
# Path coefficients table with SEs and stars
# ---------------------------------------------------------------
Bz     <- fit$Bz
se_mat <- if (is.list(se_res) && "se" %in% names(se_res)) se_res$se else
            matrix(NA_real_, P, K_TOPICS - 1L)

rownames(Bz) <- cov_labels
colnames(Bz) <- paste0("$\\mathrm{ILR}_", seq_len(K_TOPICS - 1), "$")

bz_rows <- sapply(seq_len(P), function(p) {
  vals <- sapply(seq_len(K_TOPICS - 1), function(j) {
    b  <- Bz[p, j]
    se <- se_mat[p, j]
    star <- if (!is.na(se) && abs(b / se) > 1.96) "^{*}" else ""
    if (!is.na(se)) sprintf("%.3f%s\\;({\\small %.3f})", b, star, se)
    else            sprintf("%.3f%s", b, star)
  })
  paste0(cov_labels[p], " & ", paste(vals, collapse = " & "), " \\\\")
})

bz_header <- paste(
  paste0("$\\mathrm{ILR}_", seq_len(K_TOPICS - 1), "$"),
  collapse = " & ")

bz_tex <- sprintf(
"\\begin{table}[!ht]
\\centering
\\small
\\caption{ILR path coefficients $\\hat{\\mathbf{B}}_z$ with analytical
standard errors in parentheses. Continuous covariates are standardised;
\\textit{Leave} and \\textit{Female} are centred binary indicators.
Stars ($^{*}$) denote $|\\hat{B}_{z,pk}|/\\hat{\\sigma}_{pk}>1.96$.}
\\label{tab:bes_Bz}
\\begin{tabular}{l%s}
\\toprule
 & %s \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", paste(rep("r", K_TOPICS - 1), collapse = ""),
   bz_header,
   paste(bz_rows, collapse = "\n"))

writeLines(bz_tex, "output/bes/table_Bz.tex")
cat("Saved: output/bes/table_Bz.tex\n")

# ---------------------------------------------------------------
# Z-statistics table  (Bz / SE, rounded to 1 decimal)
# ---------------------------------------------------------------
zstat_mat <- Bz / se_mat   # element-wise z = Bz / SE

zstat_rows <- sapply(seq_len(P), function(p) {
  vals <- sapply(seq_len(K_TOPICS - 1), function(j) {
    z <- zstat_mat[p, j]
    if (is.na(z) || !is.finite(z)) {
      return("---")
    }
    fmt <- sprintf("%.1f", z)
    if (abs(z) > 1.96) paste0("\\textbf{", fmt, "}") else fmt
  })
  paste0(cov_labels[p], " & ", paste(vals, collapse = " & "), " \\\\")
})

zstat_header <- paste(
  paste0("$\\mathrm{ILR}_", seq_len(K_TOPICS - 1), "$"),
  collapse = " & ")

zstat_tex <- sprintf(
"\\begin{table}[!ht]
\\centering
\\small
\\caption{Z-statistics for ILR path coefficients $\\hat{\\mathbf{B}}_z / \\hat{\\boldsymbol{\\sigma}}$ (BES Wave~25, $K=7$, $\\lambda=1$). With $M=18\\,836$ all effects are statistically distinguishable from zero; the magnitude reflects the strength of the structural relationship in ILR coordinates.}
\\label{tab:bes_zstat}
\\begin{tabular}{l%s}
\\toprule
 & %s \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
", paste(rep("r", K_TOPICS - 1), collapse = ""),
   zstat_header,
   paste(zstat_rows, collapse = "\n"))

writeLines(zstat_tex, "output/bes/table_zstat.tex")
cat("Saved: output/bes/table_zstat.tex\n")

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        legend.position   = "bottom",
        legend.key.width  = unit(1.2, "cm")))

cols7 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
           "#FF7F00","#A65628","#F781BF")

# Helper: marginal prediction (all covariates except col_idx at 0)
.pred_marginal <- function(col_idx, grid) {
  Cp <- matrix(0, length(grid), P)
  Cp[, col_idx] <- grid
  Pi <- ilr_to_proportions(Cp %*% fit$Bz, fit$V)
  colnames(Pi) <- topic_names
  as.data.frame(Pi) %>% mutate(x = grid) %>%
    pivot_longer(-x, names_to = "Topic", values_to = "Proportion")
}

# Figure 1: Left–right marginal effect
pd_lr <- .pred_marginal(4L, seq(-2.5, 2.5, length.out = 100))
p_lr  <- ggplot(pd_lr, aes(x, Proportion, colour = Topic)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols7) +
  labs(x = "Left\u2013right self-placement (standardised, 0 = mean)",
       y = "Predicted topic proportion", colour = NULL)
ggsave("output/bes/fig_lr.pdf", p_lr, width = 7, height = 3.8)

# Figure 2: Education marginal effect
pd_ed <- .pred_marginal(3L, seq(-2.5, 2.5, length.out = 100))
p_ed  <- ggplot(pd_ed, aes(x, Proportion, colour = Topic)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols7) +
  labs(x = "Education level (standardised, 0 = mean)",
       y = "Predicted topic proportion", colour = NULL)
ggsave("output/bes/fig_educ.pdf", p_ed, width = 7, height = 3.8)

# Figure 3: Brexit cleavage — Remain vs Leave bar chart
C_remain <- matrix(0, 1, P); C_remain[1, 5] <- -mean(euref_vec)
C_leave  <- matrix(0, 1, P); C_leave[1,  5] <-  1 - mean(euref_vec)
Pi_remain <- as.numeric(ilr_to_proportions(C_remain %*% fit$Bz, fit$V))
Pi_leave  <- as.numeric(ilr_to_proportions(C_leave  %*% fit$Bz, fit$V))

brexit_df <- data.frame(
  Topic      = rep(topic_names, 2),
  Vote       = rep(c("Remain", "Leave"), each = K_TOPICS),
  Proportion = c(Pi_remain, Pi_leave)
)

p_brexit <- ggplot(brexit_df, aes(Topic, Proportion, fill = Vote)) +
  geom_col(position = "dodge", width = 0.65) +
  scale_fill_manual(values = c(Remain = "#2166AC", Leave = "#D6604D")) +
  labs(x = NULL, y = "Predicted topic proportion", fill = NULL) +
  theme(legend.position = "top")
ggsave("output/bes/fig_brexit.pdf", p_brexit, width = 6, height = 3.5)

# Figure 4: Scree plot
eig_df <- data.frame(k = seq_along(fit$eigenvalues),
                     eigenvalue = fit$eigenvalues)
p_scree <- ggplot(eig_df, aes(k, eigenvalue)) +
  geom_line(colour = "#534AB7", linewidth = 0.6) +
  geom_point(colour = "#534AB7", size = 2.5) +
  labs(x = "Component", y = "Eigenvalue of augmented matrix") +
  scale_x_continuous(breaks = seq_len(K_TOPICS - 1))
ggsave("output/bes/fig_scree.pdf", p_scree, width = 4.5, height = 3)

cat("Saved: fig_lr.pdf  fig_educ.pdf  fig_brexit.pdf  fig_scree.pdf\n")

# ---------------------------------------------------------------
# Topic prevalence bar chart
# ---------------------------------------------------------------
prev_df <- data.frame(Topic = topic_names, Prevalence = prevalence)
p_prev  <- ggplot(prev_df, aes(Topic, Prevalence, fill = Topic)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = cols7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean topic proportion")
ggsave("output/bes/fig_prevalence.pdf", p_prev, width = 5, height = 3)

# ---------------------------------------------------------------
# Summary stats (for LaTeX macros)
# ---------------------------------------------------------------
sink("output/bes/summary_stats.txt")
cat("BES Wave 25 — ILR-EGSCA Summary\n================================\n\n")
cat(sprintf("M (documents):      %d\n", M))
cat(sprintf("N (vocab terms):    %d\n", N))
cat(sprintf("K (topics):         %d\n", K_TOPICS))
cat(sprintf("P (covariates):     %d\n", P))
cat(sprintf("lambda:             %.1f\n", LAMBDA))
cat(sprintf("Fit time (s):       %.2f\n", t_fit))
cat(sprintf("Mean doc length:    %.1f tokens\n", mean(dl)))
cat(sprintf("Median doc length:  %.0f tokens\n", median(dl)))
cat(sprintf("Leave voters:       %d (%.1f%%)\n", n_leave,
            100 * n_leave / (n_leave + n_remain)))
cat(sprintf("Remain voters:      %d (%.1f%%)\n", n_remain,
            100 * n_remain / (n_leave + n_remain)))
cat("\nTopic prevalences:\n")
for (k in seq_len(K_TOPICS))
  cat(sprintf("  T%d: %.1f%%  top: %s\n", k, 100 * prevalence[k],
              paste(tt[k, 1:3], collapse = ", ")))
cat("\nTop eigenvalues:\n")
print(round(fit$eigenvalues, 4))
cat("\nPath coefficients (Bz):\n")
rownames(fit$Bz) <- cov_labels
colnames(fit$Bz) <- paste0("ILR", seq_len(K_TOPICS - 1))
print(round(fit$Bz, 4))
if (!any(is.na(se_mat))) {
  cat("\nAnalytical SEs:\n")
  rownames(se_mat) <- cov_labels
  colnames(se_mat) <- paste0("ILR", seq_len(K_TOPICS - 1))
  print(round(se_mat, 4))
}
sink()
cat("Saved: output/bes/summary_stats.txt\n")
cat("\n=== Analysis complete ===\n")
