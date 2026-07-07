# ===================================================================
#  config.R  —  shared design constants (v2: clean P=K-1 DGP)
# ===================================================================
REDUCED <- nzchar(Sys.getenv("SGSCATM_REDUCED"))

# Fixed dimensions -- CLEAN design removes the P<K-1 degeneracy.
K_TOPICS <- 5L
P_COV    <- 4L                 # P = K-1 : no rank deficiency
N_VOCAB  <- 500L
DOC_LEN  <- 200L
ALPHA_BETA <- 0.1
CONF     <- 0.95
ZQ       <- qnorm(1 - (1 - CONF) / 2)

# Clean DGP: well-separated score variances (all gaps > 0.1)
SCORE_D  <- c(1.00, 0.70, 0.50, 0.35)   # eig(Cov(z))
SIGMA_EPS2 <- 0.15                        # < min(SCORE_D)

# ---- Part 3 lambda rule ------------------------------------------
# DATA-DRIVEN, truth-free: lambda_A = (K-1)-th eigenvalue of the word
# Gram W_tilde W_tilde^T  (== sigma_{K-1}(W_tilde W_tilde^T)); makes the
# covariate-augmentation block commensurate with the smallest retained
# word-signal eigenvalue.  This exact rule is used by the paper and the
# BES application.
lambda_A_rule <- function(fit) {
  # (K-1)-th eigenvalue of the WORD Gram W~ W~^T (== (K-1)-th of crossprod(W~)),
  # raw O(M) scale, independent of lambda and of the truth.
  ev <- eigen(crossprod(fit$W_tilde), symmetric = TRUE, only.values = TRUE)$values
  ev[K_TOPICS - 1L]
}

# ---- Part 2 (Block 1) grid ---------------------------------------
M_B1      <- c(500L, 1000L, 2000L, 4000L)
N_REP_B1  <- if (REDUCED) 40L else 50L
SEED_B1   <- 100000L
M_PILOT   <- 50000L
JK_M      <- M_B1                  # delete-block jackknife at all M
JK_REPS   <- 6L
JK_BLOCKS <- 100L

# ---- Part 3 diagnostic -------------------------------------------
N_REP_LAMBDA <- 30L

# ---- Part 4 ------------------------------------------------------
M_B2      <- 2000L
K_B2      <- c(3L, 5L, 10L)
SIG_B2    <- c(0.05, 0.10, 0.20, 0.30, 0.50, 0.80)
N_REP_B2  <- 20L

M_B3      <- c(1000L, 5000L)
SIGNAL_B3 <- c(weak = 0.15, strong = 0.50)
N_REP_B3  <- if (REDUCED) 12L else 15L

M_B4      <- c(500L, 2000L)
SCEN_B4   <- c("high_sep", "low_sep", "weak_sig")
N_REP_B4  <- if (REDUCED) 12L else 15L
STM_MAXIT <- 75L

if (REDUCED) cat("WARNING: SGSCATM_REDUCED set — Block1=40, Blocks3-4=12 reps.\n")
