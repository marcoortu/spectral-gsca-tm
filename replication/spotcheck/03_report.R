#!/usr/bin/env Rscript
# ===================================================================
#  Spot-check — tables, gate evaluation and figures from results/
#
#  Writes: results/summary_spotcheck.csv, results/tables_spotcheck.md,
#          results/a_bias_extrapolation.png, results/b_fd_vs_formula.png
#
#  Usage: Rscript replication/spotcheck/03_report.R
# ===================================================================

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "spotcheck", "results")
stopifnot(dir.exists(RES_DIR))

scA  <- readRDS(file.path(RES_DIR, "sca_results.rds"))
scAb <- readRDS(file.path(RES_DIR, "sca_brute_results.rds"))
scB1 <- readRDS(file.path(RES_DIR, "scb1_results.rds"))
scB2 <- readRDS(file.path(RES_DIR, "scb2_results.rds"))
scC  <- readRDS(file.path(RES_DIR, "scc_results.rds"))

ok <- function(x) Filter(function(r) is.null(r$error), x$results)
okA <- ok(scA); okAb <- ok(scAb); okB1 <- ok(scB1); okC <- ok(scC)
okB2 <- ok(scB2)
md <- character(0); add <- function(...) md <<- c(md, ...)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
sum_rows <- list()
srow <- function(...) sum_rows[[length(sum_rows) + 1L]] <<- data.frame(...)

# ===================================================================
#  SC-A: componentwise gate per test point
# ===================================================================
pts <- unique(vapply(okA, function(r) r$point, character(1)))
Km1 <- length(okA[[1]]$bhat)
gateA_pass <- TRUE
add("## SC-A — bias field b(z): componentwise verdicts",
    "",
    paste("| point | comp | b closed | intercept 1/sqrt(L) (se) |",
          "intercept 1/L (se) | Richardson 400/800 (se) | bhat_800 (se) |",
          "gate sqrtL | gate 1/L | gate Rich | gate bhat_800 |"),
    "|---|---|---|---|---|---|---|---|---|---|---|")
extrap <- list()
for (pt in pts) {
  rr <- Filter(function(r) r$point == pt, okA)
  Ls <- vapply(rr, function(r) r$L, numeric(1))
  ord <- order(Ls); rr <- rr[ord]; Ls <- Ls[ord]
  bh <- t(vapply(rr, function(r) r$bhat, numeric(Km1)))
  bs <- t(vapply(rr, function(r) r$bse, numeric(Km1)))
  bc <- rr[[1]]$b_closed
  for (cc in seq_len(Km1)) {
    # registered fit (b + c/sqrt(L)) and the diagnosis fit (b + c/L):
    # multinomial third cumulants are O(L^-2), so the remainder of
    # L * E[zhat - z - m] is O(1/L), not O(1/sqrt(L)) — the 1/L fit is
    # the theoretically correct extrapolation (draft correction).
    fit  <- lm(bh[, cc] ~ I(1 / sqrt(Ls)), weights = 1 / bs[, cc]^2)
    fitL <- lm(bh[, cc] ~ I(1 / Ls), weights = 1 / bs[, cc]^2)
    ic <- coef(fit)[1];  ic_se <- sqrt(vcov(fit)[1, 1])
    il <- coef(fitL)[1]; il_se <- sqrt(vcov(fitL)[1, 1])
    b8 <- bh[Ls == 800, cc]; b8se <- bs[Ls == 800, cc]
    gi <- abs(ic - bc[cc]) <= pmax(0.05 * abs(bc[cc]), 3 * ic_se)
    gl <- abs(il - bc[cc]) <= pmax(0.05 * abs(bc[cc]), 3 * il_se)
    g8 <- abs(b8 - bc[cc]) <= pmax(0.05 * abs(bc[cc]), 3 * b8se)
    # Richardson on the CLEAN cells only (L in {400, 800}; certificate
    # exclusions <= 0.7% there, vs up to 12% at L = 100 for large-norm
    # points, which conditions the sample and corrupts full-grid fits):
    # with an O(1/L) remainder, 2*bhat_800 - bhat_400 removes it.
    b4 <- bh[Ls == 400, cc]; b4se <- bs[Ls == 400, cc]
    ir <- 2 * b8 - b4
    ir_se <- sqrt(4 * b8se^2 + b4se^2)
    gr <- abs(ir - bc[cc]) <= pmax(0.05 * abs(bc[cc]), 3 * ir_se)
    if (!gi) gateA_pass <- FALSE
    add(sprintf("| %s | %d | %s | %s (%s) | %s (%s) | %s (%s) | %s (%s) | %s | %s | %s | %s |",
                pt, cc, fmt(bc[cc]), fmt(ic), fmt(ic_se),
                fmt(il), fmt(il_se), fmt(ir), fmt(ir_se),
                fmt(b8), fmt(b8se),
                if (gi) "PASS" else "FAIL", if (gl) "PASS" else "FAIL",
                if (gr) "PASS" else "FAIL", if (g8) "PASS" else "FAIL"))
    srow(block = "SCA", point = pt, comp = cc, b_closed = bc[cc],
         intercept = ic, intercept_se = ic_se, intercept_invL = il,
         intercept_invL_se = il_se, richardson = ir,
         richardson_se = ir_se, bhat800 = b8,
         bhat800_se = b8se, gate_intercept = gi, gate_invL = gl,
         gate_richardson = gr, gate_b800 = g8)
    extrap[[paste(pt, cc)]] <- list(pt = pt, cc = cc, Ls = Ls,
                                    bh = bh[, cc], bs = bs[, cc],
                                    bc = bc[cc], ic = ic, il = il)
  }
}
nf <- sum(vapply(okA, function(r) r$n_fail, numeric(1)))
add("", sprintf(paste("GN certificate failures across all SC-A cells: %d",
                      "of %s solves."), nf,
                format(sum(vapply(okA, function(r) r$R_used, numeric(1))) +
                         nf, big.mark = ",")), "")

# localization columns (b1 alone / -Hi a2/2 alone), reported per point
add("### SC-A localization (means over components of |bhat800 - x| / se)",
    "", "| point | vs b | vs b1 only | vs -0.5 Hi a2 only |",
    "|---|---|---|---|")
for (pt in pts) {
  r8 <- Filter(function(r) r$point == pt && r$L == 800, okA)[[1]]
  zs <- function(x) mean(abs(r8$bhat - x) / r8$bse)
  add(sprintf("| %s | %.1f | %.1f | %.1f |", pt, zs(r8$b_closed),
              zs(r8$b1_closed), zs(r8$a2half)))
}
add("")

# brute cell
Rtot <- sum(vapply(okAb, function(c) c$R_used, numeric(1)))
L_br <- okAb[[1]]$L
sr <- Reduce(`+`, lapply(okAb, function(c) c$sum_raw))
sc <- Reduce(`+`, lapply(okAb, function(c) c$sum_cv))
ssr <- Reduce(`+`, lapply(okAb, function(c) c$ss_raw))
ssc <- Reduce(`+`, lapply(okAb, function(c) c$ss_cv))
m_raw <- L_br * sr / Rtot
m_cv  <- L_br * sc / Rtot
se_raw <- L_br * sqrt(pmax(ssr / Rtot - (sr / Rtot)^2, 0) / Rtot)
se_cv  <- L_br * sqrt(pmax(ssc / Rtot - (sc / Rtot)^2, 0) / Rtot)
add("### SC-A brute-force cross-check (K=3, N=30, L=50, no CV)",
    "",
    sprintf("- R = %s certified solves.", format(Rtot, big.mark = ",")),
    sprintf("- raw  L*bias: %s (se %s)",
            paste(fmt(m_raw), collapse = ", "),
            paste(fmt(se_raw), collapse = ", ")),
    sprintf("- CV   L*bias: %s (se %s)",
            paste(fmt(m_cv), collapse = ", "),
            paste(fmt(se_cv), collapse = ", ")),
    sprintf("- closed-form b: %s",
            paste(fmt(okAb[[1]]$b_closed), collapse = ", ")),
    "")
srow(block = "SCA_brute", point = "brute", comp = NA,
     b_closed = okAb[[1]]$b_closed[1], intercept = m_cv[1],
     intercept_se = se_cv[1], bhat800 = m_raw[1], bhat800_se = se_raw[1],
     gate_intercept = NA, gate_b800 = NA)

# ===================================================================
#  SC-B route 1
# ===================================================================
add("## SC-B route 1 — ||E[G0]||_F under the Block-3 z law",
    "", "| regime | ||G0||_F | MC se | ratio | gate B1 |",
    "|---|---|---|---|---|")
gateB1 <- TRUE
for (r in okB1) {
  rat <- r$fro / r$fro_se
  if (rat <= 10) gateB1 <- FALSE
  add(sprintf("| %s | %.5f | %.5f | %.0f | %s |", r$regime, r$fro,
              r$fro_se, rat, if (rat > 10) "PASS" else "FAIL"))
  srow(block = "SCB1", point = r$regime, comp = NA, b_closed = NA,
       intercept = r$fro, intercept_se = r$fro_se, bhat800 = rat,
       bhat800_se = NA, gate_intercept = rat > 10, gate_b800 = NA)
}
rS <- okB1[[which(vapply(okB1, function(x)
  identical(x$regime, "strong"), logical(1)))]]
cn <- sqrt(colSums(rS$Gbar^2))
topc <- order(cn, decreasing = TRUE)[1:5]
corlp <- cor(log(cn + 1e-300), log(rS$pbar))
add("",
    sprintf("Strong regime: top-5 columns by ||G0[, j]||: %s; their word marginals p_j: %s (median p_j over all words = %.5f) — the mass sits on **%s**-frequency columns (correlation of log ||G0[, j]|| with log p_j: %.2f).",
            paste(topc, collapse = ", "),
            paste(sprintf("%.5f", rS$pbar[topc]), collapse = ", "),
            median(rS$pbar),
            if (corlp > 0) "HIGH" else "low", corlp),
    if (corlp > 0)
      paste("This CONTRADICTS the draft's step-criterion-lemma claim that",
            "the action concentrates on low-frequency columns: for",
            "unit-Frobenius directions confined to the 50 lowest-p",
            "columns, <G0, U> is ~300x smaller than for generic",
            "directions (see route 2, dirs 5-6).") else "",
    "")

# ===================================================================
#  SC-B route 2 — FD vs formula
# ===================================================================
g0_by_L <- setNames(lapply(scB2$g0, function(g) g$Gbar),
                    vapply(scB2$g0, function(g) as.character(g$L),
                           character(1)))
# Directions and the base draws are deterministic: rebuild them to
# (a) evaluate the formula and (b) sharpen the FD estimator with the
# analytic control variate c_r = -2 <e_r, theta(z_r)' U> — the leading
# linear-in-e component of the per-replicate directional derivative,
# which has mean exactly 0 and absorbs the O(L^-1/2) noise that made
# the raw CRN estimator vacuous (se ~ |formula|).
if (!exists("bc_theta")) {
  source(file.path(ROOT, "R", "ilr_contrast.R"))
  source(file.path(ROOT, "replication", "simulation", "sim_dgp.R"))
  source(file.path(ROOT, "replication", "basin_check", "01_functions.R"))
}
K_W <- 5L; N_W <- 500L
R_REP <- 20000L
sc_phi0 <- function(seed, K = K_W, N = N_W, alpha = 0.1) {
  set.seed(seed); .rdirichlet_matrix(K, N, alpha)
}
sc_draw_z <- function(S, b_max, seed) {
  set.seed(seed)
  t(vapply(seq_len(S), function(s) {
    B <- matrix(runif(3L * (K_W - 1L), -b_max, b_max), 3L, K_W - 1L)
    as.vector(crossprod(B, rnorm(3L))) + rnorm(K_W - 1L, 0, 0.3)
  }, numeric(K_W - 1L)))
}
set.seed(78100L)
U <- lapply(1:4, function(i) {
  u <- matrix(rnorm(K_W * N_W), K_W, N_W); u / sqrt(sum(u^2))
})
low <- order(rS$pbar)[1:50]
for (i in 1:2) {
  u <- matrix(0, K_W, N_W)
  u[, low] <- rnorm(K_W * 50L)
  U[[4L + i]] <- u / sqrt(sum(u^2))
}
cvs <- list()
for (L in c(200, 400)) {
  Phi0 <- sc_phi0(77001L)
  V <- ilr_contrast(K_W)
  Zs <- sc_draw_z(R_REP, 0.50, 78200L + L)
  Th <- bc_theta(Zs, V)
  P0 <- Th %*% Phi0
  set.seed(78300L + L)
  n <- matrix(0L, R_REP, N_W)
  for (r in seq_len(R_REP)) n[r, ] <- rmultinom(1L, L, P0[r, ])
  E <- n / L - P0
  cvs[[as.character(L)]] <- lapply(1:6, function(d)
    -2 * rowSums(E * (Th %*% U[[d]])))
}

cells <- list()
for (L in c(200, 400)) for (d in 1:6) for (h in c(1e-3, 1e-4)) {
  rp <- Filter(function(r) r$L == L && r$dir == d && r$h == h &&
                 r$sgn == 1, okB2)
  rm <- Filter(function(r) r$L == L && r$dir == d && r$h == h &&
                 r$sgn == -1, okB2)
  if (!length(rp) || !length(rm)) next
  okk <- rp[[1]]$ok & rm[[1]]$ok
  dr <- (rp[[1]]$loss[okk] - rm[[1]]$loss[okk]) / (2 * h)
  dcv <- dr - cvs[[as.character(L)]][[d]][okk]
  cells[[length(cells) + 1L]] <- data.frame(
    L = L, dir = d, h = h, Dhat = mean(dcv),
    se = sd(dcv) / sqrt(sum(okk)),
    Dhat_raw = mean(dr), se_raw = sd(dr) / sqrt(sum(okk)),
    n_used = sum(okk))
}
cells <- do.call(rbind, cells)

gateB2 <- TRUE
add("## SC-B route 2 — finite differences (CRN + control variate) vs -(2/L) <G0, U>, strong regime",
    "",
    paste("| L | dir | h | Dhat CV (se) | Dhat raw (se) | formula |",
          "ratio (CV) | gate |"),
    "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  form <- -(2 / ce$L) * sum(g0_by_L[[as.character(ce$L)]] * U[[ce$dir]])
  tol <- pmax(0.10 * abs(form), 3 * ce$se)
  g <- abs(ce$Dhat - form) <= tol
  if (!g) gateB2 <- FALSE
  add(sprintf("| %d | %d%s | %.0e | %.3e (%.1e) | %.3e (%.1e) | %.3e | %.3f | %s |",
              ce$L, ce$dir, if (ce$dir >= 5) " (low-freq)" else "",
              ce$h, ce$Dhat, ce$se, ce$Dhat_raw, ce$se_raw, form,
              ce$Dhat / form, if (g) "PASS" else "FAIL"))
  srow(block = "SCB2", point = sprintf("L%d_d%d_h%g", ce$L, ce$dir, ce$h),
       comp = NA, b_closed = form, intercept = ce$Dhat,
       intercept_se = ce$se, bhat800 = ce$Dhat / form, bhat800_se = NA,
       gate_intercept = g, gate_b800 = NA)
}
# 1/L scaling: Dhat(L=200)/Dhat(L=400) per (dir, h)
add("", "### 1/L scaling (Dhat_200 / Dhat_400, target 2)",
    "", "| dir | h | ratio |", "|---|---|---|")
for (d in 1:6) for (h in c(1e-3, 1e-4)) {
  c2 <- cells[cells$L == 200 & cells$dir == d & cells$h == h, ]
  c4 <- cells[cells$L == 400 & cells$dir == d & cells$h == h, ]
  if (nrow(c2) && nrow(c4))
    add(sprintf("| %d | %.0e | %.2f |", d, h, c2$Dhat / c4$Dhat))
}
add("")

# ===================================================================
#  SC-C — cosine vs LS pathology displacement
# ===================================================================
gm <- -rS$Gbar
cosv <- vapply(okC, function(r)
  sum(r$dPhi * gm) / sqrt(sum(r$dPhi^2) * sum(gm^2)), numeric(1))
add("## SC-C — cosine(-G0_strong, Phi_sweep10 - Phi0), LS pathology reps",
    "",
    sprintf("cosines: %s (mean %.3f) — prediction cosine > 0: %s",
            paste(fmt(cosv, 3), collapse = ", "), mean(cosv),
            if (all(cosv > 0)) "PASS" else "MIXED/FAIL"),
    "")
for (i in seq_along(cosv))
  srow(block = "SCC", point = sprintf("rep%d", i), comp = NA,
       b_closed = NA, intercept = cosv[i], intercept_se = NA,
       bhat800 = NA, bhat800_se = NA, gate_intercept = cosv[i] > 0,
       gate_b800 = NA)

writeLines(md, file.path(RES_DIR, "tables_spotcheck.md"))
all_cols <- Reduce(union, lapply(sum_rows, names))
write.csv(do.call(rbind, lapply(sum_rows, function(d) {
  d[setdiff(all_cols, names(d))] <- NA
  d[all_cols]
})), file.path(RES_DIR, "summary_spotcheck.csv"), row.names = FALSE)

# ===================================================================
#  Figures
# ===================================================================
png(file.path(RES_DIR, "a_bias_extrapolation.png"), width = 1300,
    height = 900, res = 130)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
show <- c("A_z0 1", "A_u3_n1.0 4", "A_u4_n1.0 1", "A_u5_n1.5 1",
          "A_u6_n1.5 2", "B_u3_n1.0 4")
for (key in show) {
  e <- extrap[[key]]
  if (is.null(e)) next
  x <- 1 / sqrt(e$Ls)
  plot(x, e$bh, pch = 19, col = "#534AB7",
       ylim = range(c(e$bh - 3 * e$bs, e$bh + 3 * e$bs, e$bc)),
       xlim = c(0, max(x) * 1.05),
       xlab = expression(1 / sqrt(L)), ylab = "L x CV bias estimate",
       main = sprintf("%s, comp %d", e$pt, e$cc), cex.main = 0.9)
  arrows(x, e$bh - 2 * e$bs, x, e$bh + 2 * e$bs, angle = 90, code = 3,
         length = 0.03, col = "#534AB7")
  fit <- lm(e$bh ~ x, weights = 1 / e$bs^2)
  abline(fit, col = "grey50", lty = 2)
  points(0, e$bc, pch = 17, col = "#D85A30", cex = 1.4)
  points(0, coef(fit)[1], pch = 1, col = "grey30", cex = 1.4)
  if (key == show[1])
    legend("topleft", c("bhat_L (2 se)", "closed-form b", "intercept"),
           pch = c(19, 17, 1), col = c("#534AB7", "#D85A30", "grey30"),
           bty = "n", cex = 0.75)
}
par(op); dev.off()

png(file.path(RES_DIR, "b_fd_vs_formula.png"), width = 800, height = 750,
    res = 130)
form_v <- vapply(seq_len(nrow(cells)), function(i)
  -(2 / cells$L[i]) * sum(g0_by_L[[as.character(cells$L[i])]] *
                            U[[cells$dir[i]]]), numeric(1))
plot(form_v, cells$Dhat, pch = ifelse(cells$L == 200, 19, 17),
     col = ifelse(cells$dir >= 5, "#D85A30", "#534AB7"),
     xlab = "-(2/L) <G0, U> (closed form)",
     ylab = "central FD of profiled criterion (CRN)",
     main = "SC-B: finite differences vs formula")
abline(0, 1, lty = 2, col = "grey40")
arrows(form_v, cells$Dhat - 2 * cells$se, form_v, cells$Dhat + 2 * cells$se,
       angle = 90, code = 3, length = 0.03,
       col = adjustcolor("grey30", 0.5))
legend("topleft", c("L = 200", "L = 400", "Gaussian dir", "low-freq dir",
                    "identity"),
       pch = c(19, 17, 15, 15, NA), lty = c(NA, NA, NA, NA, 2),
       col = c("grey30", "grey30", "#534AB7", "#D85A30", "grey40"),
       bty = "n", cex = 0.8)
dev.off()

cat(sprintf("Gate A (intercepts all pass): %s | Gate B1: %s | Gate B2: %s | SC-C all>0: %s\n",
            gateA_pass, gateB1, gateB2, all(cosv > 0)))
cat("Written: tables_spotcheck.md, summary_spotcheck.csv,",
    "a_bias_extrapolation.png, b_fd_vs_formula.png\n")
