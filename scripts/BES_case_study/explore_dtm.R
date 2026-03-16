suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidytext)
  library(SnowballC); library(Matrix)
})

df <- read.csv("scripts/bes_case_study/bes_w25_filtered.csv",
               stringsAsFactors = FALSE)
cat(sprintf("Input rows: %d\n", nrow(df)))

df$text <- str_squish(str_replace_all(str_to_lower(df$text), "[^a-z\\s]", " "))

custom_stops <- c(
  stop_words$word,
  "dont","didnt","doesnt","cant","wont","isnt","arent",
  "thats","theres","theyre","youre","ive","theyd","im",
  "hes","shes","its","weve","theyve",
  "think","important","issue","country","today",
  "people","thing","lot","really","much","just",
  "also","need","get","make","going","one","well")

tokens <- df %>%
  unnest_tokens(word, text) %>%
  filter(!word %in% custom_stops, nchar(word) >= 3) %>%
  mutate(word = wordStem(word, language = "english"))

term_freq <- tokens %>% count(word) %>% arrange(desc(n))
doc_freq  <- tokens %>% distinct(id, word) %>% count(word, name = "n_docs")
n_total   <- n_distinct(tokens$id)

# Survey the effect of MIN_TERM_FREQ and MIN_DOC_LENGTH on corpus size
cat("\n=== Sensitivity to MIN_TERM_FREQ (MIN_DOC_LENGTH = 1) ===\n")
for (mtf in c(3, 5, 10, 20)) {
  v <- term_freq %>% inner_join(doc_freq, by = "word") %>%
    filter(n >= mtf, n_docs / n_total <= 0.5) %>% pull(word)
  d <- tokens %>% filter(word %in% v) %>% count(id) %>% filter(n >= 1)
  cat(sprintf("  MIN_TERM_FREQ=%2d  vocab=%3d terms  docs=%5d\n",
              mtf, length(v), nrow(d)))
}

cat("\n=== Sensitivity to MIN_DOC_LENGTH (MIN_TERM_FREQ = 5) ===\n")
vocab5 <- term_freq %>% inner_join(doc_freq, by = "word") %>%
  filter(n >= 5, n_docs / n_total <= 0.5) %>% pull(word)
dtm5 <- tokens %>% filter(word %in% vocab5) %>% count(id, word)
doc_ids5 <- unique(dtm5$id)
W5 <- dtm5 %>%
  mutate(i = match(id, doc_ids5), j = match(word, vocab5)) %>%
  { sparseMatrix(i = .$i, j = .$j, x = .$n,
      dims = c(length(doc_ids5), length(vocab5))) }
rownames(W5) <- doc_ids5

for (mdl in c(1, 2, 3, 5)) {
  keep <- rowSums(W5) >= mdl
  dl   <- rowSums(W5)[keep]
  cat(sprintf("  MIN_DOC_LENGTH=%d  docs=%5d  mean_len=%.1f  med=%.0f\n",
              mdl, sum(keep), mean(dl), median(dl)))
}

# Chosen parameters: MIN_TERM_FREQ=5, MIN_DOC_LENGTH=1
cat("\n=== Final DTM (MIN_TERM_FREQ=5, MIN_DOC_LENGTH=1) ===\n")
keep <- rowSums(W5) >= 1
W    <- W5[keep, ]
doc_ids <- rownames(W)[keep]

dl <- rowSums(W)
cat(sprintf("DTM:     %d docs x %d terms\n", nrow(W), ncol(W)))
cat(sprintf("Mean len: %.1f  |  Median: %.0f\n", mean(dl), median(dl)))
cat(sprintf("Quantiles 25/50/75/90/99: %s\n",
  paste(round(quantile(dl, c(.25,.5,.75,.9,.99))), collapse = " / ")))

euref_vec <- df$euref[match(doc_ids, as.character(df$id))]
cat(sprintf("Leave share: %.1f%%  (%d Leave / %d Remain)\n",
  100 * mean(euref_vec == 1, na.rm = TRUE),
  sum(euref_vec == 1, na.rm = TRUE), sum(euref_vec == 0, na.rm = TRUE)))

cat("\nTop 30 vocab terms:\n")
top30 <- term_freq %>% filter(word %in% vocab5) %>% head(30)
print(as.data.frame(top30), row.names = FALSE)

# Build covariate matrix
cov_df <- df %>%
  filter(as.character(id) %in% doc_ids) %>%
  arrange(match(as.character(id), doc_ids)) %>%
  transmute(
    age_std  = (age  - mean(age))  / sd(age),
    female   = as.numeric(gender == 2) - mean(gender == 2),
    educ_std = (educ - mean(educ)) / sd(educ),
    lr_std   = (lr   - mean(lr))   / sd(lr),
    leave    = euref - mean(euref)
  )
C <- as.matrix(cov_df)
stopifnot(nrow(W) == nrow(C))
cat(sprintf("\nCovariate matrix C: %d x %d\n", nrow(C), ncol(C)))

saveRDS(list(W = W, C = C, vocab = vocab5, doc_ids = doc_ids, df = df),
        "scripts/bes_case_study/bes_w25_dtm.rds")
cat("Saved: scripts/bes_case_study/bes_w25_dtm.rds\n")
