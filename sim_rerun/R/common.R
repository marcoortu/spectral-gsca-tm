# ===================================================================
#  common.R  —  shared setup for the sim_rerun replication package
# ===================================================================
#  Sources the sgscatm package solver (does NOT modify package source),
#  defines the figure style (Okabe-Ito palette, theme_minimal), and a
#  handful of small helpers used across blocks.
#
#  GLOBAL RULE honoured here: every block saves per-replicate raw
#  results to sim_rerun/data/<block>.rds BEFORE aggregating, and every
#  figure is drawn from the SAME numeric object that populates the
#  corresponding table.
# ===================================================================

options(warn = 1)                       # print warnings as they occur
options(stringsAsFactors = FALSE)

# --- Paths ---------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

# Scripts are invoked from the package root; allow override via env var.
ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
if (!dir.exists(file.path(ROOT, "R")))
  stop("common.R: cannot locate package R/ under ROOT=", ROOT,
       " — run from the package root or set SGSCATM_ROOT.")

SIM_DIR    <- file.path(ROOT, "sim_rerun")
DATA_DIR   <- file.path(SIM_DIR, "data")
TAB_DIR    <- file.path(SIM_DIR, "tables")
IMG_DIR    <- file.path(SIM_DIR, "imgs")
for (d in c(DATA_DIR, TAB_DIR, IMG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Reuse the sgscatm package solver (source, do NOT modify) ------
# We deliberately source individual package files rather than the old
# egscatm_fit.R alias to avoid duplicate definitions.
.pkg_files <- c("R/sgscatm_fit.R", "R/ilr_contrast.R", "R/utils.R",
                "R/ilr_se.R", "R/refine_phi.R", "R/methods.R", "R/vcov.R")
for (f in .pkg_files) source(file.path(ROOT, f))

# Simulation infrastructure shipped with the replication package
source(file.path(ROOT, "replication/simulation/sim_dgp.R"))
source(file.path(ROOT, "replication/simulation/sim_utils.R"))

.print_reused <- function() {
  cat("\n--- Reused sgscatm package functions ------------------------\n")
  reused <- c(
    "sgscatm()          spectral solver (truncated SVD + augmented eigen)",
    "ilr_contrast()     Helmert ILR contrast matrix V (V'V=I, V'1=0)",
    "ilr_to_proportions() closure map  Pi = softmax(V z)",
    "proportions_to_ilr() inverse closure map (log-ratio)",
    "refine_phi()       one-step k-means M-step (Block 4)",
    "varimax()          orthogonal rotation, invoked inside sgscatm()",
    "sim_dgp()          structural-topic-model data generator",
    "procrustes_align() SVD-based O(K-1) alignment to ground truth",
    "eval_linearisation() exact-vs-linear closure error + Prop. bound"
  )
  cat(paste0("  * ", reused), sep = "\n")
  cat("-------------------------------------------------------------\n\n")
}

# --- Okabe-Ito colourblind-safe palette ----------------------------
OI <- c(black    = "#000000", orange   = "#E69F00", skyblue = "#56B4E9",
        green    = "#009E73", yellow   = "#F0E442", blue    = "#0072B2",
        vermilion= "#D55E00", purple   = "#CC79A7", grey    = "#999999")

# --- ggplot theme (used where ggplot2 is convenient) ---------------
.have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (.have_ggplot) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_rerun <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_blank(),          # captions live in LaTeX
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  theme_set(theme_rerun)
}

# --- PDF device helper: embed fonts via cairo_pdf where available --
open_pdf <- function(path, width, height) {
  ok <- capabilities("cairo")
  if (isTRUE(ok)) cairo_pdf(path, width = width, height = height)
  else            pdf(path, width = width, height = height)
  invisible(path)
}

# --- Numerics: guarded inverse / solve -----------------------------
safe_solve <- function(A, b = NULL, tol = 1e-12) {
  A <- as.matrix(A)
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(ch)) {
    if (is.null(b)) return(chol2inv(ch))
    return(backsolve(ch, backsolve(ch, b, transpose = TRUE)))
  }
  # fall back to pivoted QR
  qr_A <- qr(A, tol = tol)
  if (is.null(b)) return(qr.solve(qr_A, diag(nrow(A))))
  qr.solve(qr_A, b)
}

# --- LaTeX helpers -------------------------------------------------
# Minimal booktabs table writer.
write_booktabs <- function(path, header_cells, body_rows, caption, label,
                           colspec, note = NULL, small = FALSE) {
  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    if (small) "\\small" else NULL,
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{%s}", colspec),
    "\\toprule",
    paste0(paste(header_cells, collapse = " & "), " \\\\"),
    "\\midrule",
    paste0(body_rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    if (!is.null(note)) sprintf("\\\\[2pt]{\\footnotesize %s}", note) else NULL,
    "\\end{table}"
  )
  writeLines(lines, path)
  invisible(path)
}

# Format an integer with a thousands separator for LaTeX (\, ).
fmt_int <- function(x) formatC(as.integer(x), big.mark = "\\,", format = "d")

cat(sprintf("common.R loaded. ROOT = %s\n", ROOT))
cat(sprintf("cairo PDF available: %s\n", capabilities("cairo")))
