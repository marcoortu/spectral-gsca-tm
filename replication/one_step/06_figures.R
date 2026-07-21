## =====================================================================
## 06_figures.R — build all figures/tables from the saved sweep outputs.
## Run after 02–05 complete.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
od <- file.path(getwd(), "replication/one_step")
cols <- c(baseline_std = "#b2182b", proj = "#2166ac",
          onestep_uw = "#66a61e", onestep_mw = "#000000",
          onestep_wls = "#e6ab02")
pch  <- c(baseline_std = 1, proj = 16, onestep_uw = 17, onestep_mw = 15,
          onestep_wls = 18)

lineplot <- function(df, xvar, yvar, ylab, main, log = "", legpos = "topright",
                     ests = unique(df$estimator), hline = NULL) {
  xs <- sort(unique(df[[xvar]]))
  yl <- range(df[[yvar]][is.finite(df[[yvar]])], hline, na.rm = TRUE)
  plot(NA, xlim = range(xs), ylim = yl, log = log, xlab = xvar, ylab = ylab,
       main = main, cex.main = 0.95)
  if (!is.null(hline)) abline(h = hline, col = "grey60", lty = 3)
  for (e in ests) {
    s <- df[df$estimator == e, ]; s <- s[order(s[[xvar]]), ]
    lines(s[[xvar]], s[[yvar]], col = cols[e], lwd = 2)
    points(s[[xvar]], s[[yvar]], col = cols[e], pch = pch[e])
  }
  legend(legpos, legend = ests, col = cols[ests], pch = pch[ests],
         lwd = 2, bty = "n", cex = 0.8)
}

## ---- Fig 1: b_max sweep (G1,G2,G5a + G3,G4) -------------------------
if (file.exists(file.path(od, "out_02_bmax.rds"))) {
  d2 <- readRDS(file.path(od, "out_02_bmax.rds"))
  pdf(file.path(od, "figures/fig1_bmax.pdf"), width = 10, height = 7)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  lineplot(d2, "b_max", "rmse_norm", "RMSE(B_z0)/||B_z0||",
           "G1/G2/G5a: recovery vs covariate strength", log = "y", hline = 1)
  lineplot(d2, "b_max", "se_sd", "median SE / SD",
           "G3: SE calibration", legpos = "topright", hline = c(0.8, 1.25),
           ests = c("proj","onestep_uw","onestep_mw"))
  lineplot(d2, "b_max", "cov_Bz0", "coverage of B_z0",
           "G4: coverage of TRUE B_z0 (bias-driven decline)",
           legpos = "bottomleft", hline = 0.95,
           ests = c("proj","onestep_uw","onestep_mw"))
  lineplot(d2, "b_max", "cov_mean", "coverage of across-rep mean",
           "G3/G4: coverage of MEAN (SE stays calibrated)",
           legpos = "bottomleft", hline = 0.95,
           ests = c("proj","onestep_uw","onestep_mw"))
  dev.off(); cat("wrote fig1_bmax.pdf\n")
}

## ---- Fig 2: M-scaling (G5b DECISIVE) --------------------------------
if (file.exists(file.path(od, "out_03_Mscaling.rds"))) {
  d3 <- readRDS(file.path(od, "out_03_Mscaling.rds"))
  pdf(file.path(od, "figures/fig2_Mscaling.pdf"), width = 10, height = 4.5)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  lineplot(d3, "M", "rmse_norm", "RMSE(B_z0)/||B_z0||",
           "G5b: consistency — RMSE vs M", log = "xy", legpos = "bottomleft")
  lineplot(d3, "M", "bias_norm", "||mean_est - B_z0||",
           "G5b: systematic bias vs M (proj plateaus)", log = "xy",
           legpos = "bottomleft")
  dev.off(); cat("wrote fig2_Mscaling.pdf\n")
}

## ---- Fig 3: L-robustness (G5c) --------------------------------------
if (file.exists(file.path(od, "out_04_Lrobust.rds"))) {
  d4 <- readRDS(file.path(od, "out_04_Lrobust.rds"))
  pdf(file.path(od, "figures/fig3_Lrobust.pdf"), width = 10, height = 4.5)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  lineplot(d4, "L", "bias_norm", "||mean_est - B_z0||",
           "G5c: one-step bias vs document length", log = "x",
           legpos = "topright")
  lineplot(d4, "L", "rmse_norm", "RMSE(B_z0)/||B_z0||",
           "G5c: RMSE vs document length", log = "x", legpos = "topright")
  dev.off(); cat("wrote fig3_Lrobust.pdf\n")
}

## ---- Fig 4: anchor G6 -----------------------------------------------
if (file.exists(file.path(od, "out_05_anchor.rds"))) {
  d5 <- readRDS(file.path(od, "out_05_anchor.rds"))
  pdf(file.path(od, "figures/fig4_anchor.pdf"), width = 10, height = 4.5)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  lineplot(d5, "delta", "rmse_norm", "RMSE(B_z0)/||B_z0||",
           "G6: RMSE vs anchor corruption delta_Phi", legpos = "topleft",
           ests = intersect(c("proj","onestep_mw"), unique(d5$estimator)))
  lineplot(d5, "delta", "cov_Bz0", "coverage of B_z0",
           "G6: coverage vs delta_Phi", legpos = "bottomleft", hline = 0.95,
           ests = intersect(c("proj","onestep_mw"), unique(d5$estimator)))
  dev.off(); cat("wrote fig4_anchor.pdf\n")
}

## ---- Fig 5: variable-L weighted vs unweighted (G5c') ----------------
if (file.exists(file.path(od, "out_04b_varL.rds"))) {
  d6 <- readRDS(file.path(od, "out_04b_varL.rds"))
  d6 <- d6[order(d6$rmse_norm), ]
  pdf(file.path(od, "figures/fig5_weighting.pdf"), width = 7, height = 5)
  par(mar = c(7, 4, 3, 1))
  bp <- barplot(d6$rmse_norm, names.arg = d6$estimator, las = 2,
                col = cols[d6$estimator], ylab = "RMSE(B_z0)/||B_z0||",
                main = "Variable L (30% docs <50 words): weighting matters\nonly in the final pooled regression")
  text(bp, d6$rmse_norm, round(d6$rmse_norm, 3), pos = 3, cex = 0.8)
  dev.off(); cat("wrote fig5_weighting.pdf\n")
}
cat("figures done\n")
