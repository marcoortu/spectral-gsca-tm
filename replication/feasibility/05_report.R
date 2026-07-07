#!/usr/bin/env Rscript
# ===================================================================
#  Feasibility round — tables and figures from results/
#
#  Writes: results/summary_feas.csv, results/tables_feas.md,
#          results/f1_alpha_boundary.png, results/f2_slope.png,
#          results/f2_coverage.png, results/f3_k_curves.png
#
#  Usage: Rscript replication/feasibility/05_report.R
# ===================================================================

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "feasibility", "results")
stopifnot(dir.exists(RES_DIR))

f3  <- readRDS(file.path(RES_DIR, "f3_results.rds"))
f1c <- readRDS(file.path(RES_DIR, "f1c_results.rds"))
f1d_f <- file.path(RES_DIR, "f1d_results.rds")
f1d <- if (file.exists(f1d_f)) readRDS(f1d_f) else NULL
f2  <- readRDS(file.path(RES_DIR, "f2_results.rds"))
f2l <- readRDS(file.path(RES_DIR, "f2_lgrid_results.rds"))
f4  <- readRDS(file.path(RES_DIR, "f4_results.rds"))
# STM columns from the audit (identical seeds/data)
b2  <- readRDS(file.path(ROOT, "replication", "audit_block1_stm",
                         "results", "b2_results.rds"))

ok <- function(x) Filter(function(r) is.null(r$error), x$results)
md <- character(0); add <- function(...) md <<- c(md, ...)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
fmt_e <- function(x, d = 2) formatC(x, digits = d, format = "e")
msd <- function(x, d = 4) sprintf("%s (%s)", fmt(mean(x), d), fmt(sd(x), d))
REGIMES <- c("weak", "strong")
sum_rows <- list()
srow <- function(...) sum_rows[[length(sum_rows) + 1L]] <<- data.frame(...)

# ===================================================================
#  F3 — k curves and criterion pathology
# ===================================================================
ok3 <- ok(f3)
KS <- c(3L, 5L, 10L, 20L, 50L, 100L)
add("## Table F3 — MSE(B_hat) by sweep count k (oracle start, unconstrained)",
    "",
    paste("| regime |", paste(sprintf("k=%d", KS), collapse = " | "),
          "| rule sweeps (med) | rule MSE |"),
    paste0("|", paste(rep("---", length(KS) + 3), collapse = "|"), "|"))
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok3)
  mk <- vapply(KS, function(k)
    mean(vapply(rr, function(x) x$mse_path[k], numeric(1))), numeric(1))
  rs <- vapply(rr, function(x) x$rule_stop, numeric(1))
  mrule <- mean(vapply(rr, function(x) x$mse_path[x$rule_stop], numeric(1)))
  add(sprintf("| %s | %s | %d | %s |", rg,
              paste(fmt(mk, 5), collapse = " | "),
              as.integer(median(rs)), fmt(mrule, 5)))
  srow(block = "F3", regime = rg, metric = "rule_sweeps_median",
       value = median(rs))
  srow(block = "F3", regime = rg, metric = "rule_mse", value = mrule)
}
add("")

# slow-replicate diagnosis (strong regime)
rr <- Filter(function(x) x$regime == "strong", ok3)
sw10 <- vapply(rr, function(x) {
  fin <- x$mse_path[length(x$mse_path)]
  min(which(x$mse_path <= 1.1 * fin + 1e-12))
}, numeric(1))
dg <- data.frame(
  sw10 = sw10,
  relgap = vapply(rr, function(x) x$relgap, numeric(1)),
  rho = vapply(rr, function(x) x$rho_glZ, numeric(1)),
  sat = vapply(rr, function(x) x$sat_true, numeric(1)),
  nB = vapply(rr, function(x) x$nBz0, numeric(1)))
add("### F3 diagnosis (strong): sweeps to reach 1.1x final MSE vs replicate traits",
    "",
    sprintf("- slow tercile mean(relgap, rho_GL, nBz0): %s vs fast tercile: %s",
            paste(fmt(colMeans(dg[dg$sw10 >= quantile(dg$sw10, 2/3),
                                  c("relgap", "rho", "nB")]), 3),
                  collapse = ", "),
            paste(fmt(colMeans(dg[dg$sw10 <= quantile(dg$sw10, 1/3),
                                  c("relgap", "rho", "nB")]), 3),
                  collapse = ", ")),
    sprintf("- correlations of sweeps-to-1.1x with (relgap, rho_GL, sat, nBz0): %s",
            paste(fmt(c(cor(dg$sw10, dg$relgap), cor(dg$sw10, dg$rho),
                        suppressWarnings(cor(dg$sw10, dg$sat)),
                        cor(dg$sw10, dg$nB)), 2), collapse = ", ")),
    "")

# criterion pathology table (constrained descent from the truth)
pt <- Filter(function(x) !is.null(x$path_truth), ok3)
if (length(pt)) {
  add("### Criterion pathology — simplex-constrained descent FROM THE TRUTH",
      "",
      "| regime | F(truth) | F sweep 1/10/50/100 | mse sweep 1/10/50/100 | monotone |",
      "|---|---|---|---|---|")
  for (rg in REGIMES) {
    rp <- Filter(function(x) x$regime == rg, pt)
    gF <- function(s) mean(vapply(rp, function(x)
      x$path_truth$F_path[s], numeric(1)))
    gm <- function(s) mean(vapply(rp, function(x)
      x$path_truth$mse_path[s], numeric(1)))
    add(sprintf("| %s | %s | %s | %s | %d/%d |", rg,
      fmt(mean(vapply(rp, function(x) x$path_truth$F_truth, numeric(1))), 4),
      paste(fmt(vapply(c(1, 10, 50, 100), gF, numeric(1)), 4),
            collapse = " / "),
      paste(fmt(vapply(c(1, 10, 50, 100), gm, numeric(1)), 5),
            collapse = " / "),
      sum(vapply(rp, function(x) x$path_truth$monotone_ok, logical(1))),
      length(rp)))
  }
  add("", paste("F decreases monotonically while mse_Bz rises by two",
      "orders: the exact frequency-LS criterion (constrained or not)",
      "does not identify B at finite L — the estimator must be k-step,",
      "not argmin."), "")
}

# ===================================================================
#  F1c — feasible variants vs oracle (+ STM columns from audit B3)
# ===================================================================
okc <- ok(f1c)
VAR_ORDER <- c("V1", "V2", "V2_jk", "V3", "V3_jk", "V4", "V4_jk",
               "oracle_k5", "oracle_rule")
add("## Table F1c — feasible variants vs oracle reference (M = 1000, 20 reps)",
    "",
    paste("| regime | variant | mse paper | mse permutation | norm ratio |",
          "sweeps (med) | time (s) |"),
    "|---|---|---|---|---|---|---|")
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, okc)
  for (nm in VAR_ORDER) {
    g <- function(f) vapply(rr, function(x) {
      v <- x[[nm]][[f]]; if (is.null(v)) NA_real_ else v
    }, numeric(1))
    add(sprintf("| %s | %s | %s | %s | %s | %s | %s |", rg, nm,
      msd(g("paper"), 5),
      if (grepl("oracle", nm)) "—" else msd(g("perm"), 5),
      fmt(mean(g("norm_ratio")), 2),
      if (all(is.na(g("sweeps")))) "—" else
        sprintf("%d", as.integer(median(g("sweeps"), na.rm = TRUE))),
      if (all(is.na(g("time_s")))) "—" else fmt(mean(g("time_s"),
                                                     na.rm = TRUE), 1)))
    srow(block = "F1c", regime = rg, metric = paste0(nm, "_mse_paper"),
         value = mean(g("paper")))
  }
  # STM columns (audit B3, same seeds)
  rrb <- Filter(function(x) is.null(x$error) && is.null(x$stm_error) &&
                  x$regime == rg && x$M == 1000L, b2$results)
  gS <- function(m, f) vapply(rrb, function(x) x[[m]][[f]], numeric(1))
  add(sprintf("| %s | STM (native gamma; audit B3) | %s | — | %s | — | %s |",
              rg, msd(gS("stm_gamma", "mse_paper"), 5),
              fmt(mean(gS("stm_gamma", "norm_ratio")), 2),
              fmt(mean(gS("stm_gamma", "time_s")), 1)))
  add(sprintf("| %s | anchored Phi TV (context) | %s | | | | |",
              rg, msd(vapply(rr, function(x) x$anchor_tv, numeric(1)), 3)))
}
add("")

# ===================================================================
#  F1d — identification boundary
# ===================================================================
if (!is.null(f1d)) {
  okd <- ok(f1d)
  add("## Table F1d — identification boundary (strong regime, M = 1000)",
      "",
      paste("| alpha_beta | true exclusivity | anchor TV | V2 mse paper |",
            "V2 mse perm |"),
      "|---|---|---|---|---|")
  for (a in sort(unique(vapply(okd, function(x) x$alpha, numeric(1))))) {
    rr <- Filter(function(x) x$alpha == a, okd)
    g <- function(f) vapply(rr, function(x) x[[f]], numeric(1))
    add(sprintf("| %.2f | %s | %s | %s | %s |", a,
                fmt(mean(g("excl_true")), 3), msd(g("anchor_tv"), 3),
                msd(g("mse_paper"), 4), msd(g("mse_perm"), 4)))
    srow(block = "F1d", regime = sprintf("alpha_%.2f", a),
         metric = "V2_mse_paper", value = mean(g("mse_paper")))
  }
  add("")
}

# ===================================================================
#  F2 — jackknife: slope, coverage, L-grid
# ===================================================================
ok2 <- ok(f2)
df2 <- do.call(rbind, lapply(ok2, function(x) data.frame(
  M = x$M, mse_full = x$mse_full, mse_jk = x$mse_jk,
  cov_full = x$cov_full, cov_jk = x$cov_jk, cov_jk_infl = x$cov_jk_infl,
  rn_full = mean(x$rn_full), rn_jk = mean(x$rn_jk),
  rn_jk_infl = mean(x$rn_jk_infl), sweeps = x$sweeps)))
add("## Table F2 — split-document jackknife (oracle start, Block 1 grid)",
    "",
    paste("| M | RMSE uncorrected | RMSE jackknifed | cov entry (unc) |",
          "cov entry (jk) | cov entry (jk, inflated) | cov rownorm (unc) |",
          "cov rownorm (jk) | cov rownorm (jk, infl) |"),
    "|---|---|---|---|---|---|---|---|---|")
for (M in sort(unique(df2$M))) {
  d <- df2[df2$M == M, ]
  add(sprintf("| %d | %s | %s | %s | %s | %s | %s | %s | %s |", M,
    fmt(sqrt(mean(d$mse_full))), fmt(sqrt(mean(d$mse_jk))),
    fmt(mean(d$cov_full), 3), fmt(mean(d$cov_jk), 3),
    fmt(mean(d$cov_jk_infl), 3), fmt(mean(d$rn_full), 3),
    fmt(mean(d$rn_jk), 3), fmt(mean(d$rn_jk_infl), 3)))
  srow(block = "F2", regime = paste0("M", M), metric = "rmse_jk",
       value = sqrt(mean(d$mse_jk)))
  srow(block = "F2", regime = paste0("M", M), metric = "cov_entry_jk",
       value = mean(d$cov_jk))
}
sl_f <- lm(log(sqrt(mse_full)) ~ log(M), data = df2)
sl_j <- lm(log(sqrt(mse_jk)) ~ log(M), data = df2)
ci_f <- confint(sl_f)[2, ]; ci_j <- confint(sl_j)[2, ]
add("", sprintf(paste("Log-log RMSE slopes: uncorrected **%.3f** [%.3f, %.3f];",
                      "jackknifed **%.3f** [%.3f, %.3f] (target -0.5)."),
                coef(sl_f)[2], ci_f[1], ci_f[2],
                coef(sl_j)[2], ci_j[1], ci_j[2]), "")
srow(block = "F2", regime = "all", metric = "slope_jk", value = coef(sl_j)[2])
srow(block = "F2", regime = "all", metric = "slope_jk_lo", value = ci_j[1])
srow(block = "F2", regime = "all", metric = "slope_jk_hi", value = ci_j[2])

# bias share after correction (needs B_tilde matrices)
Bz0 <- matrix(c(0.40, -0.20, 0.10, 0.30, -0.15, 0.35, -0.25, 0.05,
                0.20, 0.10, 0.40, -0.30), 3, 4, byrow = TRUE)
add("### F2 bias/variance decomposition (aligned scale)",
    "", "| M | bias^2 share (uncorrected) | bias^2 share (jackknifed) |",
    "|---|---|---|")
for (M in sort(unique(df2$M))) {
  rr <- Filter(function(x) x$M == M, ok2)
  shr <- function(f) {
    arr <- simplify2array(lapply(rr, `[[`, f))
    b2v <- mean((apply(arr, c(1, 2), mean) - Bz0)^2)
    b2v / (b2v + mean(apply(arr, c(1, 2), var)))
  }
  add(sprintf("| %d | %.2f | %.2f |", M, shr("B_tilde_full"),
              shr("B_tilde_jk")))
}
add("")

okl <- ok(f2l)
add("## Table F2-L — jackknife vs L (weak regime, M = 1000, oracle start)",
    "", "| L | mse uncorrected | mse jackknifed | reduction |",
    "|---|---|---|---|")
for (dl in sort(unique(vapply(okl, function(x) x$dl, numeric(1))))) {
  rr <- Filter(function(x) x$dl == dl, okl)
  mu <- mean(vapply(rr, function(x) x$mse_full, numeric(1)))
  mj <- mean(vapply(rr, function(x) x$mse_jk, numeric(1)))
  add(sprintf("| %d | %s | %s | %.2fx |", dl, fmt(mu, 5), fmt(mj, 5),
              mu / mj))
  srow(block = "F2L", regime = paste0("L", dl), metric = "reduction",
       value = mu / mj)
}
add("")

# ===================================================================
#  F4(ii) — feasible chain on the Block 1 grid
# ===================================================================
ok4 <- ok(f4)
df4 <- do.call(rbind, lapply(ok4, function(x) data.frame(
  M = x$M, mse_full = x$mse_full, mse_jk = x$mse_jk,
  mse_perm = x$mse_perm_jk, cov = x$cov_jk, cov_infl = x$cov_jk_infl,
  rn = mean(x$rn_jk), rn_infl = mean(x$rn_jk_infl))))
add(sprintf("## Table F4(ii) — feasible chain (init %s) on the Block 1 grid",
            f4$init_variant),
    "",
    paste("| M | RMSE (jk) | RMSE perm (jk) | cov entry | cov entry (infl) |",
          "cov rownorm | cov rownorm (infl) |"),
    "|---|---|---|---|---|---|---|")
for (M in sort(unique(df4$M))) {
  d <- df4[df4$M == M, ]
  add(sprintf("| %d | %s | %s | %s | %s | %s | %s |", M,
    fmt(sqrt(mean(d$mse_jk))), fmt(sqrt(mean(d$mse_perm))),
    fmt(mean(d$cov), 3), fmt(mean(d$cov_infl), 3),
    fmt(mean(d$rn), 3), fmt(mean(d$rn_infl), 3)))
  srow(block = "F4ii", regime = paste0("M", M), metric = "rmse_jk",
       value = sqrt(mean(d$mse_jk)))
}
sl4 <- lm(log(sqrt(mse_jk)) ~ log(M), data = df4)
ci4 <- confint(sl4)[2, ]
add("", sprintf("Feasible-chain RMSE slope: **%.3f** [%.3f, %.3f].",
                coef(sl4)[2], ci4[1], ci4[2]), "")

writeLines(md, file.path(RES_DIR, "tables_feas.md"))
write.csv(do.call(rbind, sum_rows),
          file.path(RES_DIR, "summary_feas.csv"), row.names = FALSE)

# ===================================================================
#  Figures
# ===================================================================
png(file.path(RES_DIR, "f3_k_curves.png"), width = 1200, height = 560,
    res = 130)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok3)
  mpath <- rowMeans(vapply(rr, function(x) x$mse_path[1:100],
                           numeric(100)))
  plot(1:100, mpath, log = "y", type = "l", col = "#534AB7", lwd = 1.5,
       xlab = "sweep k", ylab = "mean mse_Bz",
       main = sprintf("F3: oracle-start k curve — %s", rg))
  rs <- median(vapply(rr, function(x) x$rule_stop, numeric(1)))
  abline(v = rs, lty = 2, col = "#D85A30")
  # criterion pathology overlay (constrained from truth)
  rp <- Filter(function(x) !is.null(x$path_truth), rr)
  if (length(rp)) {
    tp <- rowMeans(vapply(rp, function(x) x$path_truth$mse_path[1:100],
                          numeric(100)))
    lines(1:100, tp, col = "#1D9E75", lwd = 1.5, lty = 3)
  }
  legend("right", c("oracle start (unconstrained)", "median rule stop",
                    "constrained from TRUTH (pathology)"),
         col = c("#534AB7", "#D85A30", "#1D9E75"), lty = c(1, 2, 3),
         bty = "n", cex = 0.75)
}
par(op); dev.off()

if (!is.null(f1d)) {
  okd <- ok(f1d)
  als <- sort(unique(vapply(okd, function(x) x$alpha, numeric(1))))
  gg <- function(f) vapply(als, function(a)
    mean(vapply(Filter(function(x) x$alpha == a, okd), function(x)
      x[[f]], numeric(1))), numeric(1))
  png(file.path(RES_DIR, "f1_alpha_boundary.png"), width = 950,
      height = 650, res = 130)
  plot(als, gg("mse_paper"), log = "xy", type = "b", pch = 19,
       col = "#534AB7", xlab = "alpha_beta (Dirichlet concentration)",
       ylab = "value (log)", ylim = range(c(gg("mse_paper"),
                                            gg("anchor_tv"),
                                            gg("excl_true"))),
       main = "F1d: identification boundary (strong, M = 1000)")
  lines(als, gg("anchor_tv"), type = "b", pch = 17, col = "#D85A30")
  lines(als, gg("excl_true"), type = "b", pch = 15, col = "#1D9E75")
  legend("bottomright", c("V2 mse_Bz (paper)", "anchor TV error",
                          "true exclusivity"),
         col = c("#534AB7", "#D85A30", "#1D9E75"), pch = c(19, 17, 15),
         bty = "n", cex = 0.8)
  dev.off()
}

png(file.path(RES_DIR, "f2_slope.png"), width = 950, height = 650,
    res = 130)
Ms <- sort(unique(df2$M))
r_f <- vapply(Ms, function(M) sqrt(mean(df2$mse_full[df2$M == M])),
              numeric(1))
r_j <- vapply(Ms, function(M) sqrt(mean(df2$mse_jk[df2$M == M])),
              numeric(1))
plot(Ms, r_f, log = "xy", type = "b", pch = 19, col = "#D85A30",
     ylim = range(c(r_f, r_j)) * c(0.7, 1.3), xlab = "M",
     ylab = "RMSE of Bz",
     main = "F2: RMSE vs M, uncorrected vs jackknifed")
lines(Ms, r_j, type = "b", pch = 17, col = "#534AB7")
lines(Ms, r_j[2] * sqrt(1000 / Ms), lty = 2, col = "grey50")
legend("bottomleft", c("uncorrected", "jackknifed", "M^{-1/2} reference"),
       col = c("#D85A30", "#534AB7", "grey50"), pch = c(19, 17, NA),
       lty = c(1, 1, 2), bty = "n", cex = 0.85)
dev.off()

png(file.path(RES_DIR, "f2_coverage.png"), width = 950, height = 650,
    res = 130)
c_f <- vapply(Ms, function(M) mean(df2$cov_full[df2$M == M]), numeric(1))
c_j <- vapply(Ms, function(M) mean(df2$cov_jk[df2$M == M]), numeric(1))
c_i <- vapply(Ms, function(M) mean(df2$cov_jk_infl[df2$M == M]),
              numeric(1))
rnj <- vapply(Ms, function(M) mean(df2$rn_jk[df2$M == M]), numeric(1))
plot(Ms, c_f, log = "x", type = "b", pch = 19, col = "#D85A30",
     ylim = c(0.5, 1.02), xlab = "M", ylab = "coverage (nominal 95%)",
     main = "F2: coverage, uncorrected vs jackknifed")
lines(Ms, c_j, type = "b", pch = 17, col = "#534AB7")
lines(Ms, c_i, type = "b", pch = 18, col = "#7A6FD0")
lines(Ms, rnj, type = "b", pch = 15, col = "#1D9E75")
abline(h = 0.95, lty = 2, col = "grey40")
legend("bottomleft", c("entrywise, uncorrected", "entrywise, jk",
                       "entrywise, jk + inflated SE", "row norms, jk",
                       "nominal"),
       col = c("#D85A30", "#534AB7", "#7A6FD0", "#1D9E75", "grey40"),
       pch = c(19, 17, 18, 15, NA), lty = c(1, 1, 1, 1, 2), bty = "n",
       cex = 0.8)
dev.off()

cat("Written: tables_feas.md, summary_feas.csv, f3_k_curves.png,",
    "f1_alpha_boundary.png, f2_slope.png, f2_coverage.png\n")
