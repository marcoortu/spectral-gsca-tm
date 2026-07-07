#!/usr/bin/env Rscript
# ===================================================================
#  Deviance probe — tables and figures from results/
#
#  Writes: results/summary_dev.csv, results/tables_dev.md,
#          results/p1_pathology_overlay.png, results/p2_k_curves.png,
#          results/p2_coverage.png (only if P2b ran)
#
#  Usage: Rscript replication/deviance_probe/03_report.R
# ===================================================================

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "deviance_probe", "results")
stopifnot(dir.exists(RES_DIR))

p1  <- readRDS(file.path(RES_DIR, "p1_results.rds"))
p2a <- readRDS(file.path(RES_DIR, "p2a_results.rds"))
p2b_f <- file.path(RES_DIR, "p2b_results.rds")
p2b <- if (file.exists(p2b_f)) readRDS(p2b_f) else NULL
p3  <- readRDS(file.path(RES_DIR, "p3_results.rds"))
p3x_f <- file.path(RES_DIR, "p3_m5000_results.rds")
p3x <- if (file.exists(p3x_f)) readRDS(p3x_f) else NULL
# LS references on identical data
f3  <- readRDS(file.path(ROOT, "replication", "feasibility", "results",
                         "f3_results.rds"))
f1c <- readRDS(file.path(ROOT, "replication", "feasibility", "results",
                         "f1c_results.rds"))
b2  <- readRDS(file.path(ROOT, "replication", "audit_block1_stm",
                         "results", "b2_results.rds"))
b2x_f <- file.path(ROOT, "replication", "audit_block1_stm", "results",
                   "b2_m5000_results.rds")
b2x <- if (file.exists(b2x_f)) readRDS(b2x_f) else NULL

ok <- function(x) Filter(function(r) is.null(r$error), x$results)
md <- character(0); add <- function(...) md <<- c(md, ...)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
msd <- function(x, d = 4) sprintf("%s (%s)", fmt(mean(x), d), fmt(sd(x), d))
REGIMES <- c("weak", "strong")
sum_rows <- list()
srow <- function(...) sum_rows[[length(sum_rows) + 1L]] <<- data.frame(...)

ok1 <- ok(p1); oka <- ok(p2a); ok3f <- ok(f3); okp3 <- ok(p3)

# ===================================================================
#  P1 — pathology probe
# ===================================================================
add("## Table P1 — truth-start pathology probe: deviance vs LS on identical data",
    "",
    paste("| regime | criterion | mse sweep 1/10/50/100 | norm ratio @100 |",
          "gauge @100 | tilt (perp) @100 | monotone |"),
    "|---|---|---|---|---|---|---|")
idx4 <- c(1, 10, 50, 100)
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok1)
  for (nm in c("dev_truth", "ls_truth")) {
    g <- function(f, s) mean(vapply(rr, function(x) x[[nm]][[f]][s],
                                    numeric(1)))
    add(sprintf("| %s | %s | %s | %s | %s | %s | %d/%d |", rg,
      if (nm == "dev_truth") "deviance" else "LS (constrained)",
      paste(fmt(vapply(idx4, function(s) g("mse_path", s), numeric(1)), 4),
            collapse = " / "),
      fmt(g("nr_path", 100), 2), fmt(g("gauge_path", 100), 1),
      fmt(g("perp_path", 100), 1),
      sum(vapply(rr, function(x) x[[nm]]$monotone_ok, logical(1))),
      length(rr)))
    srow(block = "P1", regime = rg, metric = paste0(nm, "_mse100"),
         value = g("mse_path", 100))
  }
  mse0 <- mean(vapply(rr, function(x) x$mse0, numeric(1)))
  add(sprintf("| %s | (mse at sweep 0 = truth) | %s | | | | |", rg,
              fmt(mse0, 5)))
}
add("", paste("Gate P1 (dev mse@100 <= 2x mse@0): evaluated per",
              "replicate in the run log; ratios are reported there and",
              "in REPORT_DEV.md."), "")

# ===================================================================
#  P2a — k curves, deviance vs LS, oracle start
# ===================================================================
KS <- c(1L, 3L, 5L, 10L, 20L, 50L, 100L)
add("## Table P2a — oracle-start k curves: deviance vs LS (identical data)",
    "",
    paste("| regime | criterion |",
          paste(sprintf("k=%d", KS), collapse = " | "),
          "| rule stop (med) | rule MSE |"),
    paste0("|", paste(rep("---", length(KS) + 4), collapse = "|"), "|"))
for (rg in REGIMES) {
  ra <- Filter(function(x) x$regime == rg, oka)
  rf <- Filter(function(x) x$regime == rg, ok3f)
  mk_a <- vapply(KS, function(k)
    mean(vapply(ra, function(x) x$mse_path[k], numeric(1))), numeric(1))
  mk_f <- vapply(KS, function(k)
    mean(vapply(rf, function(x) x$mse_path[k], numeric(1))), numeric(1))
  rs_a <- vapply(ra, function(x) x$rule_stop, numeric(1))
  mr_a <- mean(vapply(ra, function(x) x$mse_path[x$rule_stop], numeric(1)))
  rs_f <- vapply(rf, function(x) x$rule_stop, numeric(1))
  mr_f <- mean(vapply(rf, function(x)
    x$mse_path[min(x$rule_stop, length(x$mse_path))], numeric(1)))
  add(sprintf("| %s | deviance | %s | %d | %s |", rg,
              paste(fmt(mk_a, 5), collapse = " | "),
              as.integer(median(rs_a)), fmt(mr_a, 5)))
  add(sprintf("| %s | LS (F3)  | %s | %d | %s |", rg,
              paste(fmt(mk_f, 5), collapse = " | "),
              as.integer(median(rs_f)), fmt(mr_f, 5)))
  srow(block = "P2a", regime = rg, metric = "dev_rule_mse", value = mr_a)
  srow(block = "P2a", regime = rg, metric = "ls_rule_mse", value = mr_f)
}
add("")

# ===================================================================
#  P2b — Block 1 grid (if it ran)
# ===================================================================
if (!is.null(p2b)) {
  okb <- ok(p2b)
  dfb <- do.call(rbind, lapply(okb, function(x) data.frame(
    M = x$M, mse = x$mse, cov = x$coverage, rn = mean(x$rn_covers))))
  add("## Table P2b — Block 1 grid, deviance-refined + sandwich",
      "", "| M | RMSE | cov entry | cov rownorm |", "|---|---|---|---|")
  for (M in sort(unique(dfb$M))) {
    d <- dfb[dfb$M == M, ]
    add(sprintf("| %d | %s | %s | %s |", M, fmt(sqrt(mean(d$mse))),
                fmt(mean(d$cov), 3), fmt(mean(d$rn), 3)))
  }
  sl <- lm(log(sqrt(mse)) ~ log(M), data = dfb)
  ci <- confint(sl)[2, ]
  add("", sprintf("RMSE slope: %.3f [%.3f, %.3f].", coef(sl)[2],
                  ci[1], ci[2]), "")
}

# ===================================================================
#  P3 — feasible chain under deviance
# ===================================================================
add("## Table P3 — anchor polish + feasible deviance chain (M = 1000)",
    "",
    paste("| regime | anchor TV | polished TV | chain mse paper |",
          "chain mse perm | norm ratio | coverage | time (s) |"),
    "|---|---|---|---|---|---|---|---|")
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg && x$M == 1000L, okp3)
  g <- function(f) vapply(rr, function(x) x[[f]], numeric(1))
  add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", rg,
              msd(g("tv_anchor"), 3), msd(g("tv_polish"), 3),
              msd(g("mse_paper"), 4), msd(g("mse_perm"), 4),
              fmt(mean(g("norm_ratio")), 2), fmt(mean(g("coverage")), 3),
              fmt(mean(g("time_s")), 1)))
  srow(block = "P3", regime = rg, metric = "chain_mse_paper",
       value = mean(g("mse_paper")))
  srow(block = "P3", regime = rg, metric = "tv_polish",
       value = mean(g("tv_polish")))
  # references on the same data
  rrb <- Filter(function(x) is.null(x$error) && is.null(x$stm_error) &&
                  x$regime == rg && x$M == 1000L, b2$results)
  rrf <- Filter(function(x) is.null(x$error) && x$regime == rg,
                ok(f1c))
  add(sprintf("| %s | (refs) STM %s | LS-V4+jk %s | oracle-LS rule %s | | | | |",
    rg,
    fmt(mean(vapply(rrb, function(x) x$stm_gamma$mse_paper, numeric(1))), 4),
    fmt(mean(vapply(rrf, function(x) x$V4_jk$paper, numeric(1))), 4),
    fmt(mean(vapply(rrf, function(x) x$oracle_rule$paper, numeric(1))), 4)))
}
if (!is.null(p3x)) {
  rr <- ok(p3x)
  g <- function(f) vapply(rr, function(x) x[[f]], numeric(1))
  stm5 <- if (!is.null(b2x)) mean(vapply(
    Filter(function(x) is.null(x$error) && is.null(x$stm_error),
           b2x$results),
    function(x) x$stm_gamma$mse_paper, numeric(1))) else NA
  add(sprintf("| strong M=5000 | %s | %s | %s | %s | %s | %s | %s |",
              msd(g("tv_anchor"), 3), msd(g("tv_polish"), 3),
              msd(g("mse_paper"), 4), msd(g("mse_perm"), 4),
              fmt(mean(g("norm_ratio")), 2), fmt(mean(g("coverage")), 3),
              fmt(mean(g("time_s")), 1)))
  add(sprintf("| strong M=5000 | (ref) STM %s | | | | | | |", fmt(stm5, 4)))
  srow(block = "P3", regime = "strong_M5000", metric = "chain_mse_paper",
       value = mean(g("mse_paper")))
}
add("")

writeLines(md, file.path(RES_DIR, "tables_dev.md"))
write.csv(do.call(rbind, sum_rows),
          file.path(RES_DIR, "summary_dev.csv"), row.names = FALSE)

# ===================================================================
#  Figures
# ===================================================================
png(file.path(RES_DIR, "p1_pathology_overlay.png"), width = 1250,
    height = 580, res = 130)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok1)
  mp <- function(nm) rowMeans(vapply(rr, function(x)
    x[[nm]]$mse_path[1:100], numeric(100)))
  d <- mp("dev_truth"); l <- mp("ls_truth"); o <- mp("dev_oracle")
  plot(1:100, l, log = "y", type = "l", col = "#D85A30", lwd = 1.6,
       ylim = range(c(d, l, o)) * c(0.8, 1.2),
       xlab = "sweep", ylab = "mean mse_Bz (log)",
       main = sprintf("P1 pathology overlay — %s", rg))
  lines(1:100, d, col = "#534AB7", lwd = 1.6)
  lines(1:100, o, col = "#1D9E75", lwd = 1.2, lty = 3)
  legend("bottomright",
         c("LS from truth (D1 reference)", "deviance from truth",
           "deviance from oracle-GL start"),
         col = c("#D85A30", "#534AB7", "#1D9E75"), lwd = c(1.6, 1.6, 1.2),
         lty = c(1, 1, 3), bty = "n", cex = 0.75)
}
par(op); dev.off()

png(file.path(RES_DIR, "p2_k_curves.png"), width = 1250, height = 580,
    res = 130)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (rg in REGIMES) {
  ra <- Filter(function(x) x$regime == rg, oka)
  rf <- Filter(function(x) x$regime == rg, ok3f)
  da <- rowMeans(vapply(ra, function(x) x$mse_path[1:100], numeric(100)))
  df <- rowMeans(vapply(rf, function(x) x$mse_path[1:100], numeric(100)))
  plot(1:100, df, log = "y", type = "l", col = "#D85A30", lwd = 1.6,
       ylim = range(c(da, df)) * c(0.8, 1.2),
       xlab = "sweep k", ylab = "mean mse_Bz (log)",
       main = sprintf("P2a k curves (oracle start) — %s", rg))
  lines(1:100, da, col = "#534AB7", lwd = 1.6)
  abline(v = median(vapply(ra, function(x) x$rule_stop, numeric(1))),
         lty = 2, col = "#534AB7")
  abline(v = median(vapply(rf, function(x) x$rule_stop, numeric(1))),
         lty = 2, col = "#D85A30")
  legend("right", c("LS (feasibility F3)", "deviance",
                    "median rule stops (dashed)"),
         col = c("#D85A30", "#534AB7", "grey40"), lwd = c(1.6, 1.6, NA),
         lty = c(1, 1, 2), bty = "n", cex = 0.75)
}
par(op); dev.off()

if (!is.null(p2b)) {
  okb <- ok(p2b)
  dfb <- do.call(rbind, lapply(okb, function(x) data.frame(
    M = x$M, cov = x$coverage, rn = mean(x$rn_covers))))
  Ms <- sort(unique(dfb$M))
  png(file.path(RES_DIR, "p2_coverage.png"), width = 950, height = 650,
      res = 130)
  plot(Ms, vapply(Ms, function(M) mean(dfb$cov[dfb$M == M]), numeric(1)),
       log = "x", type = "b", pch = 17, col = "#534AB7",
       ylim = c(0.5, 1.02), xlab = "M", ylab = "coverage",
       main = "P2b: deviance-refined sandwich coverage")
  lines(Ms, vapply(Ms, function(M) mean(dfb$rn[dfb$M == M]), numeric(1)),
        type = "b", pch = 15, col = "#1D9E75")
  abline(h = 0.95, lty = 2, col = "grey40")
  legend("bottomleft", c("entrywise", "row norms", "nominal"),
         col = c("#534AB7", "#1D9E75", "grey40"), pch = c(17, 15, NA),
         lty = c(1, 1, 2), bty = "n", cex = 0.8)
  dev.off()
}

cat("Written: tables_dev.md, summary_dev.csv, p1_pathology_overlay.png,",
    "p2_k_curves.png", if (!is.null(p2b)) ", p2_coverage.png" else "",
    "\n")
