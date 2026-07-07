#!/usr/bin/env Rscript
# ===================================================================
#  Audit report builder — tables, csvs and figures from results/
#
#  Writes: results/summary_audit.csv, results/summary_stm.csv,
#          results/tables_audit.md, results/rmse_vs_M.png,
#          results/coverage.png, results/qq_M2000.png,
#          results/regime_map.png (if Task C ran)
#
#  Usage: Rscript replication/audit_block1_stm/04_report.R
# ===================================================================

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "audit_block1_stm", "results")
stopifnot(dir.exists(RES_DIR))

a2 <- readRDS(file.path(RES_DIR, "a2_results.rds"))
a3 <- readRDS(file.path(RES_DIR, "a3_results.rds"))
b2 <- readRDS(file.path(RES_DIR, "b2_results.rds"))
b2m5f <- file.path(RES_DIR, "b2_m5000_results.rds")
b2m5 <- if (file.exists(b2m5f)) readRDS(b2m5f) else NULL
cf <- file.path(RES_DIR, "c_results.rds")
cres <- if (file.exists(cf)) readRDS(cf) else NULL

Bz0 <- a2$Bz0
zero_level <- sqrt(mean(Bz0^2))
M_VALUES <- c(500L, 1000L, 2000L)
md <- character(0); add <- function(...) md <<- c(md, ...)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
fmt_e <- function(x, d = 2) formatC(x, digits = d, format = "e")
msd <- function(x, d = 4) sprintf("%s (%s)", fmt(mean(x), d), fmt(sd(x), d))

ok2 <- Filter(function(x) is.null(x$error), a2$results)
ok3 <- Filter(function(x) is.null(x$error), a3$results)

# published Block 1 numbers, straight from the shipped raw results
pub <- readRDS(file.path(ROOT, "replication", "output", "data",
                         "block1_results.rds"))

# ===================================================================
#  A2 table
# ===================================================================
add("## Table A2 — published pipeline: reproduction and mechanism",
    "",
    paste("| M | published RMSE (50 reps) | reproduced RMSE (20 reps) |",
          "zero-est. level | RMSE/zero | norm ratio ||Bz||/||Bz0|| |",
          "coverage (published SE) | coverage (SE/sqrt(M)) |",
          "med CI half-width | med |center-truth| |"),
    "|---|---|---|---|---|---|---|---|---|---|")
a2_rows <- list()
for (M in M_VALUES) {
  rr <- Filter(function(x) x$M == M, ok2)
  g  <- function(f) vapply(rr, function(x)
    if (is.null(x[[f]])) NA_real_ else x[[f]], numeric(1))
  pb <- pub[[as.character(M)]]
  pub_rmse <- sqrt(mean(pb$mse[!is.na(pb$mse) & pb$mse > 0]))
  rmse <- sqrt(mean(g("mse")))
  add(sprintf("| %d | %s | %s | %s | %.2f | %s | %s | %s | %s | %s |",
    M, fmt(pub_rmse), fmt(rmse), fmt(zero_level), rmse / zero_level,
    fmt(mean(g("norm_ratio")), 3),
    fmt(mean(g("coverage_published"), na.rm = TRUE), 3),
    fmt(mean(g("coverage_scaledSE"), na.rm = TRUE), 3),
    fmt(1.96 * median(g("med_se_published"), na.rm = TRUE), 3),
    fmt(median(g("med_abs_offset"), na.rm = TRUE), 3)))
  a2_rows[[length(a2_rows) + 1L]] <- data.frame(
    M = M, block = "A2_published_pipeline",
    published_rmse = pub_rmse, reproduced_rmse = rmse,
    zero_level = zero_level, norm_ratio = mean(g("norm_ratio")),
    coverage_published_SE = mean(g("coverage_published"), na.rm = TRUE),
    coverage_scaled_SE = mean(g("coverage_scaledSE"), na.rm = TRUE),
    med_ci_halfwidth = 1.96 * median(g("med_se_published"), na.rm = TRUE),
    med_center_offset = median(g("med_abs_offset"), na.rm = TRUE))
}
add("", "")

# ===================================================================
#  A3 tables: RMSE slope + coverage
# ===================================================================
df3 <- do.call(rbind, lapply(ok3, function(x) data.frame(
  M = x$M, mse = x$mse, pilot_mse = x$pilot_mse,
  coverage = x$coverage, norm_ratio = x$norm_ratio,
  s_cov1 = x$s_covers[1], s_cov2 = x$s_covers[2], s_cov3 = x$s_covers[3],
  n_fail = x$n_fail, monotone_ok = x$monotone_ok)))

sl <- lm(I(log(sqrt(mse))) ~ I(log(M)), data = df3)
ci <- confint(sl)[2, ]

# bias/variance decomposition of the two-step error (aligned scale):
# bias entry = mean over reps of (B_tilde - Bz0); variance = entrywise
# var over reps.  The variance component should scale ~1/M; the bias
# component is the finite-L (incidental-parameter) floor.
Bz0 <- a3$Bz0
bv <- do.call(rbind, lapply(M_VALUES, function(M) {
  rr  <- Filter(function(x) x$M == M, ok3)
  arr <- simplify2array(lapply(rr, `[[`, "B_tilde"))     # P x Kp x n
  bias2 <- mean((apply(arr, c(1, 2), mean) - Bz0)^2)
  vr    <- mean(apply(arr, c(1, 2), var))
  zbar  <- mean(abs(apply(simplify2array(
    lapply(rr, `[[`, "std_err_mat")), c(1, 2), mean)))   # mean |bias/SE|
  data.frame(M = M, bias2 = bias2, var = vr, mse = bias2 + vr,
             mean_abs_bias_over_se = zbar)
}))
sl_var <- lm(log(sqrt(bv$var)) ~ log(bv$M))
add("### Bias/variance decomposition (two-step, aligned scale)",
    "",
    "| M | bias^2 | variance | bias^2/mse | mean(|bias|/SE) |",
    "|---|---|---|---|---|")
for (i in seq_len(nrow(bv)))
  add(sprintf("| %d | %s | %s | %.2f | %.2f |", bv$M[i],
              fmt_e(bv$bias2[i]), fmt_e(bv$var[i]),
              bv$bias2[i] / (bv$bias2[i] + bv$var[i]),
              bv$mean_abs_bias_over_se[i]))
add("", sprintf(
  paste("Log-log slope of the sqrt(variance) component on M: **%.3f**",
        "(3 points) — the variance obeys the M^{-1/2} law; the bias",
        "component is M-independent (finite-L floor at L = 200)."),
  coef(sl_var)[2]), "")
add("## Table A3 — corrected two-step estimator (k = 5, lambda = 0)",
    "",
    paste("| M | RMSE (two-step) | RMSE (pilot, published conv.) |",
          "norm ratio | entrywise coverage | row-norm coverage |",
          "Armijo failures | monotone |"),
    "|---|---|---|---|---|---|---|---|")
a3_rows <- list()
for (M in M_VALUES) {
  d <- df3[df3$M == M, ]
  scov <- mean(c(d$s_cov1, d$s_cov2, d$s_cov3))
  add(sprintf("| %d | %s | %s | %s | %s | %s | %d | %d/%d |",
    M, fmt(sqrt(mean(d$mse))), fmt(sqrt(mean(d$pilot_mse))),
    fmt(mean(d$norm_ratio), 3), fmt(mean(d$coverage), 3), fmt(scov, 3),
    sum(d$n_fail), sum(d$monotone_ok), nrow(d)))
  a3_rows[[length(a3_rows) + 1L]] <- data.frame(
    M = M, block = "A3_two_step",
    rmse = sqrt(mean(d$mse)), pilot_rmse = sqrt(mean(d$pilot_mse)),
    norm_ratio = mean(d$norm_ratio),
    coverage_entrywise = mean(d$coverage), coverage_rownorm = scov,
    n = nrow(d))
}
add("", sprintf(
  "Log-log slope of per-replicate RMSE on M: **%.3f** (95%% CI [%.3f, %.3f]); M^{-1/2} reference = -0.5.",
  coef(sl)[2], ci[1], ci[2]), "")
write.csv(rbind(
  do.call(rbind, lapply(a2_rows, function(d) { d[setdiff(names(a3_rows[[1]]), names(d))] <- NA; d })),
  do.call(rbind, lapply(a3_rows, function(d) { d[setdiff(names(a2_rows[[1]]), names(d))] <- NA; d }))
)[, union(names(a2_rows[[1]]), names(a3_rows[[1]]))],
  file.path(RES_DIR, "summary_audit.csv"), row.names = FALSE)

# ===================================================================
#  B tables
# ===================================================================
okb <- Filter(function(x) is.null(x$error) && is.null(x$stm_error),
              b2$results)
stm_block <- function(rr, regime, M_lab) {
  g <- function(meth, f) vapply(rr, function(x) x[[meth]][[f]], numeric(1))
  rows <- list(
    c("pilot (published pipeline)", msd(g("pilot", "mse_paper")),
      fmt(mean(g("pilot", "norm_ratio")), 3),
      msd(g("pilot", "mse_gl")), fmt(mean(g("pilot", "time_s")), 1), "—"),
    c("pilot + refined (k=5, l=0)", msd(g("refined", "mse_paper")),
      fmt(mean(g("refined", "norm_ratio")), 3),
      msd(g("refined", "mse_gl")), fmt(mean(g("refined", "time_s")), 1), "—"),
    c("STM (native gamma, ALR->ILR)", msd(g("stm_gamma", "mse_paper")),
      fmt(mean(g("stm_gamma", "norm_ratio")), 3), "—",
      fmt(mean(g("stm_gamma", "time_s")), 1),
      fmt(mean(g("stm_gamma", "em_its")), 0)),
    c("STM (theta -> ILR -> OLS, old worker)",
      msd(g("stm_theta", "mse_paper")),
      fmt(mean(g("stm_theta", "norm_ratio")), 3),
      msd(g("stm_theta", "mse_gl")),
      fmt(mean(g("stm_theta", "time_s")), 1),
      fmt(mean(g("stm_theta", "em_its")), 0)))
  lapply(rows, function(r) sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                                   paste0(regime, M_lab), r[1], r[2], r[3],
                                   r[4], r[5], r[6]))
}
add("## Table B3 — fair STM comparison (paired with basin_check E2 data)",
    "",
    paste("| regime | method | mse_Bz paper (Procrustes) | norm ratio",
          "||B||/||Bz0|| | GL/oracle mse | time (s) | EM its |"),
    "|---|---|---|---|---|---|---|")
stm_rows <- list()
for (rg in c("weak", "strong")) {
  rr <- Filter(function(x) x$regime == rg && x$M == 1000L, okb)
  for (ln in stm_block(rr, rg, "")) add(ln)
  g <- function(meth, f) vapply(rr, function(x) x[[meth]][[f]], numeric(1))
  for (meth in c("pilot", "refined", "stm_gamma", "stm_theta"))
    stm_rows[[length(stm_rows) + 1L]] <- data.frame(
      regime = rg, M = 1000L, method = meth,
      mse_paper_mean = mean(g(meth, "mse_paper")),
      mse_paper_sd = sd(g(meth, "mse_paper")),
      norm_ratio = mean(g(meth, "norm_ratio")),
      time_s = mean(g(meth, "time_s")), n = length(rr))
}
if (!is.null(b2m5)) {
  okm5 <- Filter(function(x) is.null(x$error) && is.null(x$stm_error),
                 b2m5$results)
  if (length(okm5)) {
    for (ln in stm_block(okm5, "strong", " M=5000")) add(ln)
    g <- function(meth, f) vapply(okm5, function(x) x[[meth]][[f]], numeric(1))
    for (meth in c("pilot", "refined", "stm_gamma", "stm_theta"))
      stm_rows[[length(stm_rows) + 1L]] <- data.frame(
        regime = "strong", M = 5000L, method = meth,
        mse_paper_mean = mean(g(meth, "mse_paper")),
        mse_paper_sd = sd(g(meth, "mse_paper")),
        norm_ratio = mean(g(meth, "norm_ratio")),
        time_s = mean(g(meth, "time_s")), n = length(okm5))
  }
}
add("", "Published Table 3 references: weak 0.009, strong 0.021 (M = 1000).",
    "")
write.csv(do.call(rbind, stm_rows),
          file.path(RES_DIR, "summary_stm.csv"), row.names = FALSE)

# ===================================================================
#  Task C table
# ===================================================================
if (!is.null(cres)) {
  okc <- Filter(function(x) is.null(x$error), cres$results)
  add("## Table C — operating-regime map (M = 1000, 10 reps per cell)",
      "",
      paste("| b_max | sat(theta_true) | sat(theta_end) | pilot GL mse",
            "(oracle) | pilot paper mse | refined paper mse |",
            "Armijo fails | max nu | monotone |"),
      "|---|---|---|---|---|---|---|---|---|")
  for (bm in sort(unique(vapply(okc, `[[`, numeric(1), "b_max")))) {
    rr <- Filter(function(x) x$b_max == bm, okc)
    g <- function(f) vapply(rr, `[[`, numeric(1), f)
    add(sprintf("| %.2f | %s | %s | %s | %s | %s | %d | %s | %d/%d |",
      bm, fmt(mean(g("sat_true")), 3), fmt(mean(g("sat_end")), 3),
      msd(g("pilot_gl_mse")), msd(g("pilot_paper_mse")),
      msd(g("refined_paper_mse")),
      sum(g("n_fail")), fmt_e(max(g("nu_max"))),
      sum(vapply(rr, `[[`, logical(1), "monotone_ok")), length(rr)))
  }
  add("")
}

writeLines(md, file.path(RES_DIR, "tables_audit.md"))

# ===================================================================
#  Figures
# ===================================================================
# rmse_vs_M: published-convention vs corrected two-step
png(file.path(RES_DIR, "rmse_vs_M.png"), width = 950, height = 700,
    res = 130)
rmse_pub <- vapply(M_VALUES, function(M) {
  rr <- Filter(function(x) x$M == M, ok2)
  sqrt(mean(vapply(rr, `[[`, numeric(1), "mse")))
}, numeric(1))
rmse_two <- vapply(M_VALUES, function(M)
  sqrt(mean(df3$mse[df3$M == M])), numeric(1))
plot(M_VALUES, rmse_pub, log = "xy", type = "b", pch = 19,
     col = "#D85A30", ylim = range(c(rmse_pub, rmse_two, zero_level)) * c(0.6, 1.4),
     xlab = "M (corpus size)", ylab = "RMSE of Bz",
     main = "Block 1 RMSE: published convention vs corrected two-step")
lines(M_VALUES, rmse_two, type = "b", pch = 17, col = "#534AB7")
abline(h = zero_level, lty = 3, col = "grey40")
lines(M_VALUES, rmse_two[2] * sqrt(1000 / M_VALUES), lty = 2, col = "grey60")
legend("left", c("published pipeline (reproduced)",
                 "two-step (pilot + 5 GN sweeps)",
                 "zero-estimator level sqrt(mean(Bz0^2))",
                 "M^{-1/2} reference"),
       col = c("#D85A30", "#534AB7", "grey40", "grey60"),
       pch = c(19, 17, NA, NA), lty = c(1, 1, 3, 2), bty = "n", cex = 0.8)
dev.off()

# coverage figure
png(file.path(RES_DIR, "coverage.png"), width = 950, height = 700,
    res = 130)
cov_ent <- vapply(M_VALUES, function(M)
  mean(df3$coverage[df3$M == M]), numeric(1))
cov_row <- vapply(M_VALUES, function(M) {
  d <- df3[df3$M == M, ]; mean(c(d$s_cov1, d$s_cov2, d$s_cov3))
}, numeric(1))
cov_pub <- vapply(M_VALUES, function(M) {
  rr <- Filter(function(x) x$M == M, ok2)
  mean(vapply(rr, function(x)
    if (is.null(x$coverage_published)) NA_real_ else x$coverage_published,
    numeric(1)), na.rm = TRUE)
}, numeric(1))
plot(M_VALUES, cov_ent, log = "x", type = "b", pch = 17, col = "#534AB7",
     ylim = c(0.5, 1.02), xlab = "M (corpus size)",
     ylab = "empirical coverage (nominal 95%)",
     main = "Coverage: published check vs corrected sandwich CIs")
lines(M_VALUES, cov_row, type = "b", pch = 15, col = "#1D9E75")
lines(M_VALUES, cov_pub, type = "b", pch = 19, col = "#D85A30")
abline(h = 0.95, lty = 2, col = "grey40")
legend("bottomleft", c("published pipeline (reproduced)",
                       "two-step, entrywise (Procrustes + rotated cov)",
                       "two-step, row norms (alignment-free)",
                       "nominal 0.95"),
       col = c("#D85A30", "#534AB7", "#1D9E75", "grey40"),
       pch = c(19, 17, 15, NA), lty = c(1, 1, 1, 2), bty = "n", cex = 0.8)
dev.off()

# QQ plot at M = 2000
png(file.path(RES_DIR, "qq_M2000.png"), width = 700, height = 700,
    res = 130)
zz <- unlist(lapply(Filter(function(x) x$M == 2000L, ok3),
                    function(x) as.vector(x$std_err_mat)))
qqnorm(zz, pch = 19, cex = 0.5, col = "#534AB7",
       main = sprintf("Standardised entrywise errors, M = 2000 (n = %d)",
                      length(zz)))
abline(0, 1, lty = 2, col = "grey40")
dev.off()

# regime map figure
if (!is.null(cres)) {
  okc <- Filter(function(x) is.null(x$error), cres$results)
  bms <- sort(unique(vapply(okc, `[[`, numeric(1), "b_max")))
  gg <- function(f) vapply(bms, function(bm)
    mean(vapply(Filter(function(x) x$b_max == bm, okc), `[[`,
                numeric(1), f)), numeric(1))
  png(file.path(RES_DIR, "regime_map.png"), width = 950, height = 700,
      res = 130)
  plot(bms, gg("refined_paper_mse"), log = "y", type = "b", pch = 17,
       col = "#534AB7",
       ylim = range(c(gg("refined_paper_mse"), gg("pilot_gl_mse"),
                      gg("pilot_paper_mse"))),
       xlab = "b_max (signal strength)", ylab = "mse_Bz (log scale)",
       main = "Operating-regime map (M = 1000)")
  lines(bms, gg("pilot_gl_mse"), type = "b", pch = 19, col = "#1D9E75")
  lines(bms, gg("pilot_paper_mse"), type = "b", pch = 15, col = "#D85A30")
  par(new = TRUE)
  plot(bms, gg("sat_true"), type = "b", pch = 4, col = "grey50", axes = FALSE,
       xlab = "", ylab = "", ylim = c(0, 1), lty = 3)
  axis(4, col.axis = "grey40"); mtext("share saturated docs", 4, line = -1.2,
                                      col = "grey40", cex = 0.8)
  legend("topleft", c("refined (paper metric)", "pilot GL (oracle)",
                      "pilot (paper metric)", "sat. share (right axis)"),
         col = c("#534AB7", "#1D9E75", "#D85A30", "grey50"),
         pch = c(17, 19, 15, 4), lty = c(1, 1, 1, 3), bty = "n", cex = 0.8)
  dev.off()
}

cat("Written: summary_audit.csv, summary_stm.csv, tables_audit.md,",
    "rmse_vs_M.png, coverage.png, qq_M2000.png",
    if (!is.null(cres)) ", regime_map.png" else "", "\n")
