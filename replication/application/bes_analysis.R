#!/usr/bin/env Rscript
# ===================================================================
#  British Election Study — Case Study for Biometrika Paper
# ===================================================================
#
#  Analysis of open-ended responses from the BES Internet Panel
#  using egscatm. Focus: "Most Important Issue" (MII) question.
#
#  Input:  BES SPSS file (open-ended responses)
#  Output: output/bes/  — tables, figures, fitted objects
#
#  Usage:
#    Rscript bes_analysis.R --explore    # Step 1: explore data
#    Rscript bes_analysis.R --fit        # Step 2: fit models
#    Rscript bes_analysis.R --all        # Everything
# ===================================================================

args <- commandArgs(trailingOnly = TRUE)
DO_EXPLORE <- any(c("--explore", "--all") %in% args) || length(args) == 0
DO_FIT     <- any(c("--fit", "--all") %in% args) || length(args) == 0

suppressPackageStartupMessages({
  library(haven)       # read SPSS
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tidytext)    # tokenisation
  library(SnowballC)   # stemming
  library(Matrix)      # sparse matrices
  library(ggplot2)
  library(xtable)      # LaTeX tables
})

# Source egscatm
# devtools::load_all(".")
source("R/sgscatm_fit.R")
source("R/ilr_contrast.R")
source("R/ilr_se.R")
source("R/refine_phi.R")
source("R/methods.R")
source("R/utils.R")

dir.create("output/bes", recursive = TRUE, showWarnings = FALSE)


# ===================================================================
#  STEP 1: EXPLORE AND IDENTIFY VARIABLES
# ===================================================================

if (DO_EXPLORE) {
  cat("=== Step 1: Exploring BES data (Wave 25) ===\n")

  STRINGS_FILE <- "scripts/bes_case_study/BES2024_W30Strings_v30.1.sav.zip"
  PANEL_FILE   <- "scripts/bes_case_study/BES_teaching_long_v30.1.sav-1.zip"
  WAVE         <- 25L

  cat("  Reading Strings file (metadata only)...\n")
  strings_raw <- read_sav(unz(STRINGS_FILE, "BES2024_W30Strings_v30.1.sav"),
                          n_max = 0)
  cat("  Strings columns:", ncol(strings_raw), "\n")
  mii_cols <- grep("^MII_text", names(strings_raw), value = TRUE)
  cat("  MII text columns:", length(mii_cols), "\n")
  for (col in mii_cols) cat(sprintf("    %s\n", col))

  cat("\n  Reading panel (wave", WAVE, ")...\n")
  panel_w <- read_sav(unz(PANEL_FILE, "BES_teaching_long_v30.1.sav")) %>%
    filter(wave == WAVE)
  cat("  Panel wave", WAVE, "rows:", nrow(panel_w), "\n")

  cov_check <- c("age","gender","p_edlevel","lr_scale",
                 "euRefVote","euRefVoteAfter")
  cat("\n  Covariate availability in wave", WAVE, ":\n")
  for (col in cov_check) {
    if (col %in% names(panel_w)) {
      n <- sum(!is.na(panel_w[[col]]))
      cat(sprintf("    %-20s  n_valid=%d\n", col, n))
    } else {
      cat(sprintf("    %-20s  NOT FOUND\n", col))
    }
  }

  # Sample MII responses for wave 25
  strings_sample <- read_sav(unz(STRINGS_FILE, "BES2024_W30Strings_v30.1.sav"),
                             n_max = 500) %>%
    select(id, MII_textW25)
  cat("\n  Sample MII_textW25 responses:\n")
  resp <- strings_sample$MII_textW25
  resp <- resp[!is.na(resp) & nchar(trimws(resp)) > 0]
  for (r in head(resp, 10)) cat(sprintf("    %s\n", r))

  sink("output/bes/exploration_summary.txt")
  cat("BES Wave 25 Exploration\n=======================\n\n")
  cat("Strings file:", STRINGS_FILE, "\n")
  cat("Panel file:  ", PANEL_FILE, "\n")
  cat("Wave:", WAVE, "| Panel rows:", nrow(panel_w), "\n\n")
  cat("Panel columns:\n")
  cat(paste(names(panel_w), collapse = "\n"), "\n")
  sink()

  cat("  Exploration summary saved to output/bes/exploration_summary.txt\n\n")
}


# ===================================================================
#  STEP 2: PREPROCESS AND FIT
# ===================================================================

if (DO_FIT) {
  cat("=== Step 2: Preprocessing and fitting ===\n")

  # ---------------------------------------------------------------
  # CONFIGURATION
  # ---------------------------------------------------------------
  # BES2024 W30 Strings — open-ended MII text (one col per wave)
  STRINGS_FILE <- "scripts/bes_case_study/BES2024_W30Strings_v30.1.sav.zip"
  # BES Teaching version — numeric panel (long format: id x wave)
  PANEL_FILE   <- "scripts/bes_case_study/BES_teaching_long_v30.1.sav-1.zip"

  WAVE <- 25L   # Wave to analyse

  # Column names — Teaching Version / Strings file
  TEXT_COL   <- paste0("MII_textW", WAVE)  # "MII_textW25"
  ID_COL     <- "id"
  AGE_COL    <- "age"            # continuous, 18–97
  GENDER_COL <- "gender"         # 1=Men, 2=Female
  EDUC_COL   <- "p_edlevel"      # ordinal 0–5 (No qual … Postgrad)
  LR_COL     <- "lr_scale"       # composite 0–10 (0=Left, 10=Right)
  EUREF_COL  <- "euRefVoteAfter" # 0=Remain, 1=Leave, 2/9999=excl.

  # Model parameters
  K_TOPICS <- 10L
  LAMBDA   <- 3
  MIN_DOC_LENGTH <- 3L   # minimum tokens per response (MII answers are short)
  MIN_TERM_FREQ  <- 10L  # minimum corpus frequency for a term
  MAX_DOC_PROP   <- 0.5  # maximum document proportion for a term
  # ---------------------------------------------------------------

  cat("  Reading BES Strings file...\n")
  strings_raw <- read_sav(unz(STRINGS_FILE,
                              "BES2024_W30Strings_v30.1.sav"))

  cat("  Reading BES Teaching panel (long)...\n")
  panel_raw <- read_sav(unz(PANEL_FILE,
                            "BES_teaching_long_v30.1.sav")) %>%
    filter(wave == WAVE) %>%
    select(all_of(c(ID_COL, AGE_COL, GENDER_COL, EDUC_COL,
                    LR_COL, EUREF_COL)))

  cat(sprintf("  Panel wave %d: %d respondents\n", WAVE, nrow(panel_raw)))

  cat("  Joining on id...\n")
  bes_raw <- strings_raw %>%
    select(all_of(c(ID_COL, TEXT_COL))) %>%
    inner_join(panel_raw, by = ID_COL)

  # --- Extract and clean text ---
  cat("  Extracting text and covariates...\n")

  # Build analysis dataframe
  bes_df <- bes_raw %>%
    transmute(
      id    = as.character(.data[[ID_COL]]),
      text  = as.character(.data[[TEXT_COL]]),
      age   = as.numeric(.data[[AGE_COL]]),
      gender = as.numeric(.data[[GENDER_COL]]),
      educ  = as.numeric(.data[[EDUC_COL]]),
      lr    = as.numeric(.data[[LR_COL]]),
      euref = as.numeric(.data[[EUREF_COL]])
    ) %>%
    filter(
      !is.na(text), nchar(trimws(text)) > 0,
      !is.na(age), age >= 18,
      !is.na(gender), gender %in% c(1, 2),
      !is.na(educ), educ >= 0,
      !is.na(lr), lr >= 0,
      euref %in% c(0, 1)   # keep only Remain (0) and Leave (1)
    ) %>%
    mutate(
      text = str_to_lower(text),
      text = str_replace_all(text, "[^a-z\\s]", " "),
      text = str_squish(text)
    )

  cat(sprintf("  %d responses with complete covariates\n", nrow(bes_df)))

  # --- Tokenise and build DTM ---
  cat("  Tokenising...\n")

  # Custom stop words (add BES-specific ones)
  custom_stops <- c(
    stop_words$word,               # tidytext default
    "dont", "didnt", "doesnt", "cant", "wont", "isnt", "arent",
    "thats", "theres", "theyre", "youre", "ive", "theyd",
    "im", "id", "hes", "shes", "its", "weve", "theyve",
    "think", "important", "issue", "country", "today",  # question stems
    "people", "thing", "lot", "really", "much", "just",
    "also", "need", "get", "make", "going", "one", "well"
  )

  tokens <- bes_df %>%
    select(id, text) %>%
    unnest_tokens(word, text) %>%
    filter(
      !word %in% custom_stops,
      nchar(word) >= 3L
    ) %>%
    mutate(word = wordStem(word, language = "english"))

  # Term frequencies
  term_freq <- tokens %>%
    count(word) %>%
    arrange(desc(n))

  # Document frequencies
  doc_freq <- tokens %>%
    distinct(id, word) %>%
    count(word, name = "n_docs")

  n_docs_total <- n_distinct(tokens$id)

  # Vocabulary selection
  vocab <- term_freq %>%
    inner_join(doc_freq, by = "word") %>%
    filter(
      n >= MIN_TERM_FREQ,
      n_docs / n_docs_total <= MAX_DOC_PROP
    ) %>%
    pull(word)

  cat(sprintf("  Vocabulary: %d terms (from %d unique tokens)\n",
              length(vocab), nrow(term_freq)))

  # Build DTM
  dtm_tidy <- tokens %>%
    filter(word %in% vocab) %>%
    count(id, word)

  # Convert to matrix
  doc_ids <- unique(dtm_tidy$id)
  dtm_sparse <- dtm_tidy %>%
    mutate(
      i = match(id, doc_ids),
      j = match(word, vocab)
    ) %>%
    {sparseMatrix(i = .$i, j = .$j, x = .$n,
                  dims = c(length(doc_ids), length(vocab)),
                  dimnames = list(doc_ids, vocab))}

  # Filter short documents
  doc_lengths <- rowSums(dtm_sparse)
  keep_docs <- doc_lengths >= MIN_DOC_LENGTH
  dtm_sparse <- dtm_sparse[keep_docs, ]
  doc_ids <- doc_ids[keep_docs]

  W <- as.matrix(dtm_sparse)
  cat(sprintf("  DTM: %d documents x %d terms\n", nrow(W), ncol(W)))

  # --- Build covariate matrix ---
  # Block 1 — demographic (structural): age, female, educ
  # Block 2 — attitudinal (political framing): lr, leave
  # All centred; continuous vars standardised.
  cov_df <- bes_df %>%
    filter(id %in% doc_ids) %>%
    arrange(match(id, doc_ids)) %>%  # ensure same order as W
    transmute(
      age_std  = (age  - mean(age))  / sd(age),
      female   = as.numeric(gender == 2) - mean(gender == 2),
      educ_std = (educ - mean(educ)) / sd(educ),
      lr_std   = (lr   - mean(lr))   / sd(lr),
      leave    = euref - mean(euref) # euref: 0=Remain, 1=Leave
    )

  C <- as.matrix(cov_df)
  P <- ncol(C)
  cat(sprintf("  Covariates: %d (age, female, educ, lr, leave)\n", P))

  # Verify alignment
  stopifnot(nrow(W) == nrow(C))

  # --- Fit sgscatm ---
  cat("  Fitting sgscatm...\n")
  t0 <- proc.time()
  fit <- sgscatm(W, C, K = K_TOPICS, lambda = LAMBDA, rotate = TRUE)
  t_fit <- (proc.time() - t0)[3]
  cat(sprintf("  Fit time: %.2f seconds\n", t_fit))

  # --- Refine Phi ---
  cat("  Refining topic-word distributions...\n")
  fit <- refine_phi(fit, W, method = "kmeans", temp = 0.5, seed = 42)

  # --- Analytical SEs ---
  cat("  Computing analytical standard errors...\n")
  se_res <- tryCatch(
    ilr_se_analytical(fit),
    error = function(e) {
      cat("  Analytical SE failed:", e$message, "\n")
      cat("  Falling back to bootstrap...\n")
      ilr_se(fit, W, C, B = 200, seed = 42)
    }
  )

  # --- Save fitted objects ---
  saveRDS(fit, "output/bes/fit_egscatm.rds")
  saveRDS(se_res, "output/bes/se_results.rds")
  saveRDS(list(W = W, C = C, vocab = vocab, doc_ids = doc_ids,
               cov_df = cov_df, cov_names = names(cov_df)),
          "output/bes/data_processed.rds")

  cat("  Saved: fit_egscatm.rds, se_results.rds, data_processed.rds\n")


  # =================================================================
  #  RESULTS FOR THE PAPER
  # =================================================================

  cat("\n=== Generating paper outputs ===\n")

  # --- Table: Top terms per topic ---
  tt <- top_terms(fit, n = 10, vocab = vocab)
  cat("\n  Top 10 terms per topic:\n")
  print(tt)

  # LaTeX table
  tt_tex <- paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\small\n",
    "\\caption{Top 10 terms per topic from egscatm applied to ",
    "BES open-ended ``most important issue'' responses. ",
    "$K{=}", K_TOPICS, "$, $\\lambda{=}", LAMBDA,
    "$, with varimax rotation and $K$-means refinement ",
    "($\\tau{=}0.5$).}\n",
    "\\label{tab:bes_topics}\n",
    "\\begin{tabular}{cl}\n",
    "\\toprule\n",
    "Topic & Top terms \\\\\n",
    "\\midrule\n"
  )
  for (k in seq_len(K_TOPICS)) {
    terms_str <- paste(tt[k, ], collapse = " $\\cdot$ ")
    tt_tex <- paste0(tt_tex,
      sprintf("T%d & %s \\\\\n", k, terms_str))
  }
  tt_tex <- paste0(tt_tex,
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  writeLines(tt_tex, "output/bes/table_topics.tex")

  # --- Table: Topic prevalence ---
  prevalence <- colMeans(fit$Pi)
  prev_df <- data.frame(
    Topic = paste0("T", seq_len(K_TOPICS)),
    Prevalence = sprintf("%.1f\\%%", 100 * prevalence)
  )

  # --- Table: Path coefficients with SEs ---
  Bz <- fit$Bz
  se_mat <- if (is.list(se_res) && "se" %in% names(se_res)) {
    se_res$se
  } else {
    matrix(NA, nrow(Bz), ncol(Bz))
  }

  cov_labels <- c("Age", "Female", "Education", "Left--right", "Leave")

  bz_tex <- paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\caption{ILR path coefficients $\\hat{\\mathbf{B}}_z$ with ",
    "analytical standard errors (in parentheses). ",
    "Covariates are standardised. Stars denote significance at ",
    "$|\\hat{B}_{z,pk}|/\\hat\\sigma_{pk} > 1.96$.}\n",
    "\\label{tab:bes_Bz}\n",
    "\\begin{tabular}{l", paste(rep("c", K_TOPICS - 1), collapse = ""), "}\n",
    "\\toprule\n",
    " & ", paste(paste0("ILR", seq_len(K_TOPICS - 1)), collapse = " & "),
    " \\\\\n",
    "\\midrule\n"
  )
  for (p in seq_len(P)) {
    vals <- character(K_TOPICS - 1)
    for (j in seq_len(K_TOPICS - 1)) {
      b_val <- Bz[p, j]
      se_val <- se_mat[p, j]
      star <- if (!is.na(se_val) && abs(b_val / se_val) > 1.96) "*" else ""
      vals[j] <- sprintf("%.3f%s", b_val, star)
      if (!is.na(se_val)) {
        vals[j] <- paste0(vals[j], sprintf(" (%.3f)", se_val))
      }
    }
    bz_tex <- paste0(bz_tex,
      cov_labels[p], " & ", paste(vals, collapse = " & "), " \\\\\n")
  }
  bz_tex <- paste0(bz_tex,
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  writeLines(bz_tex, "output/bes/table_Bz.tex")

  # --- Figures: marginal effect of each covariate (others held at 0) ---
  cat("  Generating figures...\n")

  theme_set(theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom"))

  .marginal_fig <- function(col_idx, grid, x_label, file) {
    C_pred <- matrix(0, length(grid), P)
    C_pred[, col_idx] <- grid
    Pi_pred <- ilr_to_proportions(C_pred %*% fit$Bz, fit$V)
    colnames(Pi_pred) <- paste0("T", seq_len(K_TOPICS))
    pd <- as.data.frame(Pi_pred) %>%
      mutate(x = grid) %>%
      pivot_longer(-x, names_to = "Topic", values_to = "Proportion")
    gg <- ggplot(pd, aes(x = x, y = Proportion, colour = Topic)) +
      geom_line(linewidth = 0.7) +
      labs(x = x_label, y = "Predicted topic proportion", colour = NULL) +
      scale_colour_brewer(palette = "Set3")
    ggsave(file, gg, width = 7, height = 4)
    invisible(gg)
  }

  # col 4: lr_std  |  col 5: leave  |  col 3: educ_std
  .marginal_fig(4L, seq(-2, 2, length.out = 50),
                "Left\u2013right self-placement (standardised)",
                "output/bes/fig_prevalence_lr.pdf")

  .marginal_fig(3L, seq(-2, 2, length.out = 50),
                "Education level (standardised)",
                "output/bes/fig_prevalence_educ.pdf")

  # Brexit cleavage: Leave vs Remain (centred binary, ±0.5 span)
  .marginal_fig(5L, seq(-0.5, 0.5, length.out = 2),
                "Remain \u2013 Leave",
                "output/bes/fig_prevalence_leave.pdf")

  # --- Figure: Eigenvalue scree plot ---
  eig_df <- data.frame(
    k = seq_along(fit$eigenvalues),
    eigenvalue = fit$eigenvalues
  )

  p_scree <- ggplot(eig_df, aes(x = k, y = eigenvalue)) +
    geom_point(size = 2, colour = "#534AB7") +
    geom_line(colour = "#534AB7", linewidth = 0.5) +
    labs(x = "Component", y = "Eigenvalue of augmented matrix") +
    scale_x_continuous(breaks = seq_len(K_TOPICS - 1))

  ggsave("output/bes/fig_screeplot.pdf", p_scree,
         width = 5, height = 3)

  # --- Summary statistics for paper text ---
  sink("output/bes/summary_stats.txt")
  cat("BES Case Study — Summary Statistics\n")
  cat("====================================\n\n")
  cat(sprintf("Documents (after filtering): M = %d\n", nrow(W)))
  cat(sprintf("Vocabulary size: N = %d\n", ncol(W)))
  cat(sprintf("Topics: K = %d\n", K_TOPICS))
  cat(sprintf("Covariates: P = %d\n", P))
  cat(sprintf("Lambda: %.1f\n", LAMBDA))
  cat(sprintf("Fitting time: %.2f seconds\n", t_fit))
  cat(sprintf("Mean document length: %.1f tokens\n", mean(rowSums(W))))
  cat(sprintf("Median document length: %.0f tokens\n", median(rowSums(W))))
  cat("\nTopic prevalences:\n")
  for (k in seq_len(K_TOPICS)) {
    cat(sprintf("  T%d: %.1f%%\n", k, 100 * prevalence[k]))
  }
  cat("\nTop eigenvalues:\n")
  print(round(fit$eigenvalues, 2))
  cat("\nPath coefficients (Bz):\n")
  rownames(Bz) <- cov_labels
  colnames(Bz) <- paste0("ILR", seq_len(K_TOPICS - 1))
  print(round(Bz, 4))
  if (!is.null(se_mat) && all(!is.na(se_mat))) {
    cat("\nAnalytical SEs:\n")
    rownames(se_mat) <- cov_labels
    colnames(se_mat) <- paste0("ILR", seq_len(K_TOPICS - 1))
    print(round(se_mat, 4))
  }
  sink()

  cat("\n  All outputs saved to output/bes/\n")
  cat("  Files:\n")
  cat("    output/bes/table_topics.tex    — Top terms per topic\n")
  cat("    output/bes/table_Bz.tex        — Path coefficients\n")
  cat("    output/bes/fig_prevalence_lr.pdf\n")
  cat("    output/bes/fig_prevalence_educ.pdf\n")
  cat("    output/bes/fig_prevalence_leave.pdf\n")
  cat("    output/bes/fig_screeplot.pdf\n")
  cat("    output/bes/summary_stats.txt\n")
  cat("\n  To include in LaTeX:\n")
  cat("    \\input{output/bes/table_topics}\n")
  cat("    \\input{output/bes/table_Bz}\n")
}

cat("\n=== BES analysis complete ===\n")
