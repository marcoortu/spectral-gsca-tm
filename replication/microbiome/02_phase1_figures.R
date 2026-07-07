#!/usr/bin/env Rscript
# Phase 1 figures + LaTeX tables from output/phase1_cells.rds
suppressPackageStartupMessages({ library(ggplot2) })
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

df <- readRDS("output/phase1_cells.rds")
sweep <- df[df$M == 2000 & df$N == 200 & df$K == 5, ]
sweep <- sweep[order(sweep$b_max), ]

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank()))

# --- Fig 1: coverage vs b_max, with nominal + delocalization ratio ---
p_cov <- ggplot(sweep, aes(b_max, coverage)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
  geom_line(colour = "#534AB7", linewidth = 0.8) +
  geom_point(aes(size = r_deloc), colour = "#534AB7") +
  scale_size_continuous(name = "deloc. ratio r", range = c(2, 6)) +
  labs(x = expression("covariate strength " * b[max]),
       y = "analytical 95% CI coverage",
       title = "Coverage vs covariate strength (crossover)") +
  ylim(min(0.8, min(sweep$coverage)), 1.0)
ggsave("output/figures/phase1_coverage_vs_bmax.pdf", p_cov, width = 6, height = 4)

# --- Fig 2: SE/SD calibration ratio vs b_max ---
p_cal <- ggplot(sweep, aes(b_max, se_sd_ratio)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = 0.8, ymax = 1.25), fill = "grey85", alpha = 0.4) +
  geom_line(colour = "#B7434A", linewidth = 0.8) + geom_point(size = 3, colour = "#B7434A") +
  labs(x = expression("covariate strength " * b[max]),
       y = "analytical SE / empirical SD",
       title = "SE calibration vs covariate strength")
ggsave("output/figures/phase1_calibration_vs_bmax.pdf", p_cal, width = 6, height = 4)

# --- Fig 3: delocalization ratio vs b_max ---
p_r <- ggplot(sweep, aes(b_max, r_deloc)) +
  geom_line(colour = "#2E7D32", linewidth = 0.8) +
  geom_point(size = 3, colour = "#2E7D32") +
  labs(x = expression("covariate strength " * b[max]),
       y = expression("delocalization ratio r"),
       title = "Delocalization boundary")
ggsave("output/figures/phase1_deloc_vs_bmax.pdf", p_r, width = 6, height = 4)

# --- Fig 4: recovery RMSE vs b_max (scale-recovered) ---
p_rmse <- ggplot(sweep, aes(b_max, rmse_recover)) +
  geom_line(colour = "#00695C", linewidth = 0.8) +
  geom_point(size = 3, colour = "#00695C") +
  labs(x = expression("covariate strength " * b[max]),
       y = expression("RMSE(" * B[z] * ") after scale-recovery"),
       title = "Bz recovery error vs covariate strength")
ggsave("output/figures/phase1_rmse_vs_bmax.pdf", p_rmse, width = 6, height = 4)

# --- LaTeX table: b_max sweep ---
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
tab <- sweep[, c("b_max", "r_deloc", "coverage", "se_sd_ratio",
                 "cc_min_mean", "rmse_recover", "boot_se_ratio")]
lines <- c(
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  "$b_{\\max}$ & $r$ & coverage & SE/SD & $\\min\\rho_{cc}$ & RMSE & boot SE/SD \\\\",
  "\\midrule",
  apply(tab, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
        fmt(r["b_max"],2), fmt(r["r_deloc"],2), fmt(r["coverage"]),
        fmt(r["se_sd_ratio"]), fmt(r["cc_min_mean"]), fmt(r["rmse_recover"]),
        ifelse(is.na(r["boot_se_ratio"]), "--", fmt(r["boot_se_ratio"],2)))),
  "\\bottomrule", "\\end{tabular}")
writeLines(lines, "output/tables/phase1_bmax_sweep.tex")

# --- LaTeX table: dimension robustness (favorable b_max=0.25) ---
grd <- df[df$b_max %in% c(0.25, 1.00) & !(df$M == 2000 & df$N == 200 & df$K == 5), ]
grd <- grd[order(grd$b_max, grd$K, grd$N, grd$M), ]
lines2 <- c(
  "\\begin{tabular}{rrrrrrr}",
  "\\toprule",
  "$b_{\\max}$ & M & N & K & $r$ & coverage & SE/SD \\\\",
  "\\midrule",
  apply(grd, 1, function(r) sprintf("%s & %d & %d & %d & %s & %s & %s \\\\",
        fmt(r["b_max"],2), as.integer(r["M"]), as.integer(r["N"]),
        as.integer(r["K"]), fmt(r["r_deloc"],2), fmt(r["coverage"]),
        fmt(r["se_sd_ratio"]))),
  "\\bottomrule", "\\end{tabular}")
writeLines(lines2, "output/tables/phase1_dim_robustness.tex")

cat("Figures + tables written.\n")
print(sweep[, c("b_max","r_deloc","coverage","se_sd_ratio","cc_min_mean","rmse_recover")])
