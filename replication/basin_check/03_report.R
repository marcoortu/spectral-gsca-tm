#!/usr/bin/env Rscript
# ===================================================================
#  Basin-condition verification — tables and figures from results/
# ===================================================================
#
#  Reads results/*.rds produced by 02_run_experiments.R and writes:
#    results/summary.csv        — E2 accuracy table (flat)
#    results/tables.md          — all tables in markdown (E1-E5, audit)
#    results/spectrum_tail.png  — spectrum tail figure (+ .pdf)
#    results/e5_biasfloor.png   — E5 figure (if E5 was run)
#
#  Usage (from the package root):
#    Rscript replication/basin_check/03_report.R
# ===================================================================

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "basin_check", "results")
stopifnot(dir.exists(RES_DIR))

e12 <- readRDS(file.path(RES_DIR, "e12_results.rds"))
e34 <- readRDS(file.path(RES_DIR, "e34_results.rds"))
e5f <- file.path(RES_DIR, "e5_results.rds")
e5  <- if (file.exists(e5f)) readRDS(e5f) else NULL

REGIMES <- c("weak", "strong")
STM_REF <- c(weak = 0.009, strong = 0.021)   # published Table 3, M = 1000

ok12 <- Filter(function(x) is.null(x$error), e12$results)
ok34 <- Filter(function(x) is.null(x$error), e34$results)
md <- character(0)
add <- function(...) md <<- c(md, ...)

fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
fmt_e <- function(x, d = 2) formatC(x, digits = d, format = "e")
msd <- function(x, d = 4)
  sprintf("%s (%s)", fmt(mean(x), d), fmt(sd(x), d))

# ===================================================================
#  Pilot alignment table (feeds the audit section)
# ===================================================================
add("## Table A — Pilot accuracy under the three alignments",
    "",
    paste("| regime | reps | mse_Bz paper (Procrustes on fit$Bz) |",
          "mse_Bz GL-aligned Z | mse_Bz OP-aligned Z |",
          "rho_pil GL (Z / Phi) | rho_pil OP |"),
    "|---|---|---|---|---|---|---|")
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok12)
  p  <- lapply(rr, `[[`, "pilot")
  rho_gl <- t(vapply(p, `[[`, numeric(3), "rho_gl"))
  rho_op <- t(vapply(p, `[[`, numeric(3), "rho_op"))
  add(sprintf("| %s | %d | %s | %s | %s | %s (%s / %s) | %s |",
    rg, length(rr),
    msd(vapply(p, `[[`, numeric(1), "paper_mse")),
    msd(vapply(p, `[[`, numeric(1), "gl_mse")),
    msd(vapply(p, `[[`, numeric(1), "op_mse")),
    fmt(mean(rho_gl[, 1]), 2), fmt(mean(rho_gl[, 2]), 2),
    fmt(mean(rho_gl[, 3]), 2), fmt(mean(rho_op[, 1]), 2)))
}
add("", "rho values are means over replicates of",
    "sqrt(||Z_pil - Z_true||_F^2 + ||Phi_pil - Phi_true||_F^2).", "")

# ===================================================================
#  E1 — operational basin check
# ===================================================================
add("## Table E1 — Endpoint coincidence (refine from pilot vs from truth)",
    "",
    paste("| regime | lambda | same basin (strict) | same basin (gauge-aware) |",
          "med dZ_rel | med dFit_rel | med dEtaPerp_rel | med dF_rel |",
          "med sweeps (pilot) | converged p/t |"),
    "|---|---|---|---|---|---|---|---|---|---|")
for (rg in REGIMES) for (lam in c("lambda0", "lambda1")) {
  rr <- Filter(function(x) x$regime == rg, ok12)
  rf <- lapply(rr, function(x) x$refits[[lam]])
  add(sprintf("| %s | %s | %d/%d | %d/%d | %s | %s | %s | %s | %d | %d/%d |",
    rg, sub("lambda", "", lam),
    sum(vapply(rf, `[[`, logical(1), "same_basin")), length(rf),
    sum(vapply(rf, `[[`, logical(1), "same_basin_gauge")), length(rf),
    fmt_e(median(vapply(rf, `[[`, numeric(1), "dZ_rel"))),
    fmt_e(median(vapply(rf, `[[`, numeric(1), "dFit_rel"))),
    fmt_e(median(vapply(rf, `[[`, numeric(1), "dEta_perp_rel"))),
    fmt_e(median(vapply(rf, `[[`, numeric(1), "dF_rel"))),
    as.integer(median(vapply(rf, `[[`, numeric(1), "sweeps_pilot"))),
    sum(vapply(rf, `[[`, logical(1), "converged_pilot")),
    sum(vapply(rf, `[[`, logical(1), "converged_truth"))))
}
add("", paste("Strict criterion: dZ_rel < 1e-3 AND dF_rel < 1e-8 (pre-registered).",
    "Gauge-aware: dFit_rel < 1e-3 AND dF_rel < 1e-8, where dFit compares the",
    "gauge-invariant fitted matrices Theta(Z)Phi."), "")

mono_all <- all(vapply(ok12, function(x)
  all(vapply(x$refits, `[[`, logical(1), "monotone_ok")), logical(1)))
add(sprintf("Monotone F decrease held in **all** refinement runs: %s.",
            if (mono_all) "yes" else "**NO — see raw traces**"), "")

# ===================================================================
#  E2 — accuracy table
# ===================================================================
add("## Table E2 — mse_Bz, mean (sd) over replicates",
    "",
    "| regime | estimator | paper metric (Procrustes) | GL / direct metric |",
    "|---|---|---|---|")
sum_rows <- list()
for (rg in REGIMES) {
  rr <- Filter(function(x) x$regime == rg, ok12)
  p  <- lapply(rr, `[[`, "pilot")
  grab <- function(lam, f1, f2) list(
    paper = vapply(rr, function(x) x$refits[[lam]][[f1]], numeric(1)),
    gl    = vapply(rr, function(x) x$refits[[lam]][[f2]], numeric(1)))
  ests <- list(
    "pilot (paper alignment)" = list(
      paper = vapply(p, `[[`, numeric(1), "paper_mse"),
      gl    = vapply(p, `[[`, numeric(1), "gl_mse")),
    "pilot + refined (lambda = 0)" = grab("lambda0", "refined_paper", "refined_gl"),
    "pilot + refined (lambda = 1)" = grab("lambda1", "refined_paper", "refined_gl"),
    "refined from truth (lambda = 0, oracle floor)" =
      grab("lambda0", "truthref_paper", "truthref_gl"),
    "refined from truth (lambda = 1)" =
      grab("lambda1", "truthref_paper", "truthref_gl"))
  for (nm in names(ests)) {
    add(sprintf("| %s | %s | %s | %s |", rg, nm,
                msd(ests[[nm]]$paper), msd(ests[[nm]]$gl)))
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      regime = rg, estimator = nm,
      mse_paper_mean = mean(ests[[nm]]$paper),
      mse_paper_sd   = sd(ests[[nm]]$paper),
      mse_gl_mean    = mean(ests[[nm]]$gl),
      mse_gl_sd      = sd(ests[[nm]]$gl),
      n = length(ests[[nm]]$paper))
  }
  add(sprintf("| %s | STM (published Table 3, M = 1000) | %.4f | — |",
              rg, STM_REF[[rg]]))
  sum_rows[[length(sum_rows) + 1L]] <- data.frame(
    regime = rg, estimator = "STM (published)",
    mse_paper_mean = STM_REF[[rg]], mse_paper_sd = NA,
    mse_gl_mean = NA, mse_gl_sd = NA, n = NA)
}
add("", paste("Pilot GL column: OLS of the GL-aligned scores on C, direct",
    "entry-wise MSE (the GL alignment already absorbs rotation and scale).",
    "Refined estimates live in model coordinates, so their two columns",
    "differ only by the final Procrustes rotation."), "")
write.csv(do.call(rbind, sum_rows),
          file.path(RES_DIR, "summary.csv"), row.names = FALSE)

# ===================================================================
#  E3 — Hessian analysis
# ===================================================================
add("## Table E3 — Hessian spectrum diagnostics (5 reps per regime)",
    "",
    paste("| regime | lambda | anchor | lambda_max | gamma",
          "(perp at l=0 / raw min at l=1) | max |QtHQ| rel |",
          "max ||H q||/lmax | cosines > 0.99 (of 20) | sep. ratio |"),
    "|---|---|---|---|---|---|---|---|")
for (rg in REGIMES) for (lam in c("lambda0", "lambda1")) {
  rr <- Filter(function(x) x$regime == rg, ok34)
  for (anchor in c("at_truth", "at_star")) {
    aa <- lapply(rr, function(x) x[[lam]][[anchor]])
    lmax <- vapply(aa, `[[`, numeric(1), "lambda_max")
    gam  <- vapply(aa, function(a)
      if (!is.null(a$gamma_perp)) a$gamma_perp else a$gamma_raw, numeric(1))
    qth  <- vapply(aa, function(a) max(abs(a$QtHQ_eigs)), numeric(1)) / lmax
    hqt  <- vapply(aa, function(a) max(a$HQt_relnorms), numeric(1))
    ncos <- vapply(aa, function(a)
      sum(a$principal_cosines > 0.99), numeric(1))
    sep  <- gam / (qth * lmax)
    add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %.1f | %s |",
      rg, sub("lambda", "", lam),
      if (anchor == "at_truth") "truth eta0" else "M-est. eta*",
      fmt(mean(lmax), 1), fmt_e(mean(gam)), fmt_e(mean(qth)),
      fmt_e(mean(hqt)), mean(ncos), fmt(mean(sep), 1)))
  }
}
add("", paste("QtHQ = 20 x 20 compression of H onto the gauge basis; its",
    "eigenvalues measure curvature along the gauge orbit directly (ARPACK",
    "cannot resolve the 20-fold degenerate cluster; principal cosines of",
    "the raw Lanczos vectors undercount the null space for that reason).",
    "sep. ratio = gamma / max|gauge eigenvalue|: how cleanly the",
    "non-gauge curvature separates from the gauge block."), "")

# ===================================================================
#  E4 — Kantorovich diagnostic
# ===================================================================
add("## Table E4 — Kantorovich diagnostic (same 5 reps)",
    "",
    paste("| regime | lambda | L_H | rho_perp (from eta0) |",
          "rho_perp (from eta*) | gamma(eta0) | gamma(eta*) |",
          "r = 2 L_H rho_perp/gamma (eta0) | r (eta*) |"),
    "|---|---|---|---|---|---|---|---|---|")
for (rg in REGIMES) for (lam in c("lambda0", "lambda1")) {
  rr <- Filter(function(x) x$regime == rg, ok34)
  g  <- function(f) vapply(rr, function(x) x[[lam]][[f]], numeric(1))
  rp <- vapply(rr, `[[`, numeric(1), "rho_perp")
  add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    rg, sub("lambda", "", lam),
    msd(g("L_H"), 3), fmt(mean(rp), 2), fmt(mean(g("rho_star_perp")), 2),
    fmt_e(mean(g("gamma_truth"))), fmt_e(mean(g("gamma_star"))),
    fmt_e(mean(g("r_truth"))), fmt_e(mean(g("r_star")))))
}
add("", "")

# ===================================================================
#  E5 — finite-L bias floor
# ===================================================================
if (!is.null(e5)) {
  ok5 <- Filter(function(x) is.null(x$error), e5$results)
  add("## Table E5 — refined-from-truth (lambda = 0) vs document length",
      "(weak regime)", "",
      "| doc_length | mse_Bz paper | mse_Bz direct | mse ratio vs 1/L ratio |",
      "|---|---|---|---|")
  dls <- sort(unique(vapply(ok5, `[[`, numeric(1), "doc_length")))
  m_by_dl <- vapply(dls, function(dl) {
    v <- vapply(Filter(function(x) x$doc_length == dl, ok5),
                `[[`, numeric(1), "mse_paper")
    mean(v)
  }, numeric(1))
  for (i in seq_along(dls)) {
    rr5 <- Filter(function(x) x$doc_length == dls[i], ok5)
    ref <- which(dls == 200)
    add(sprintf("| %d | %s | %s | %s |",
      dls[i],
      msd(vapply(rr5, `[[`, numeric(1), "mse_paper")),
      msd(vapply(rr5, `[[`, numeric(1), "mse_gl")),
      if (length(ref) == 1)
        sprintf("%.2f vs %.2f", m_by_dl[i] / m_by_dl[ref],
                (200 / dls[i])) else "—"))
  }
  add("")

  png(file.path(RES_DIR, "e5_biasfloor.png"), width = 900, height = 650,
      res = 130)
  plot(dls, m_by_dl, log = "xy", type = "b", pch = 19,
       xlab = "document length L", ylab = "mse_Bz (paper metric)",
       main = "E5: refined-from-truth error vs L (weak regime)")
  lines(dls, m_by_dl[dls == 200] * 200 / dls, lty = 2, col = "grey40")
  legend("topright", c("refined from truth", "1/L reference"),
         lty = c(1, 2), pch = c(19, NA), col = c("black", "grey40"),
         bty = "n")
  dev.off()
}

# ===================================================================
#  Spectrum tail figure (replicate 1 of each regime, at eta*)
# ===================================================================
for (dev_fun in c("png", "pdf")) {
  f <- file.path(RES_DIR, paste0("spectrum_tail.", dev_fun))
  if (dev_fun == "png") png(f, width = 1300, height = 620, res = 130)
  else pdf(f, width = 10, height = 4.8)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  for (rg in REGIMES) {
    x1 <- Filter(function(x) x$regime == rg, ok34)[[1]]
    a0 <- x1$lambda0$at_star
    a1 <- x1$lambda1$at_star
    vals <- list(
      "raw H, lambda=0"       = a0$raw_smallest,
      "deflated P H P, l=0"   = a0$deflated_smallest,
      "gauge block Qt'HQt"    = sort(abs(a0$QtHQ_eigs)),
      "raw H, lambda=1"       = a1$raw_smallest)
    cols <- c("#534AB7", "#1D9E75", "#D85A30", "#C49A00")
    ylim <- range(unlist(lapply(vals, function(v) pmax(abs(v), 1e-13))))
    plot(NA, xlim = c(1, max(lengths(vals))), ylim = ylim, log = "y",
         xlab = "eigenvalue index (ascending |value|)",
         ylab = "|eigenvalue|",
         main = sprintf("Spectrum tail at eta* — %s (rep 1)", rg))
    for (j in seq_along(vals)) {
      v <- sort(pmax(abs(vals[[j]]), 1e-13))
      points(seq_along(v), v, col = cols[j], pch = c(19, 17, 15, 18)[j],
             cex = 0.7)
      lines(seq_along(v), v, col = cols[j], lwd = 0.8)
    }
    legend("bottomright", names(vals), col = cols,
           pch = c(19, 17, 15, 18), cex = 0.7, bty = "n")
  }
  par(op)
  dev.off()
}

writeLines(md, file.path(RES_DIR, "tables.md"))
cat("Written:\n  ", file.path(RES_DIR, "summary.csv"), "\n  ",
    file.path(RES_DIR, "tables.md"), "\n  ",
    file.path(RES_DIR, "spectrum_tail.png"), " (+ .pdf)\n")
if (!is.null(e5)) cat("  ", file.path(RES_DIR, "e5_biasfloor.png"), "\n")
