# sgscatm

**Structural Topic Modeling via ILR-Spectral-GSCA**

`sgscatm` fits structural topic models that operate natively in the Aitchison geometry of the simplex. Document-topic proportions are represented in isometric log-ratio (ILR) coordinates and linked to covariates through a Generalized Structured Component Analysis (GSCA) formulation. Estimation reduces to a single eigenvalue decomposition of a covariate-augmented document similarity matrix — yielding a non-iterative, globally optimal solution.

## Requirements

- R >= 4.0
- Imports: `Matrix`, `methods`
- Suggests (optional): `testthat >= 3.0`, `knitr`, `rmarkdown`

## Installation

### From GitHub (development version)

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("marcoortu/spectral-gsca-tm")
```

### From source

Clone the repository and build/install from the project root:

```bash
git clone https://github.com/marcoortu/spectral-gsca-tm.git
cd spectral-gsca-tm
R CMD build .
R CMD INSTALL sgscatm_0.1.0.tar.gz
```

Or directly from R:

```r
# from the project root directory
devtools::install(".")
```

## Build & check

```bash
# Build source package
R CMD build .

# Full CRAN check
R CMD check sgscatm_0.1.0.tar.gz

# Or with devtools
devtools::check()
```

## Quick start

```r
library(sgscatm)

set.seed(42)
M <- 200   # documents
N <- 100   # vocabulary size
K <- 4     # number of topics
P <- 2     # number of covariates

# Simulate a document-term matrix and covariate matrix
W <- matrix(rpois(M * N, 5), M, N)
C <- scale(matrix(rnorm(M * P), M, P))

# Fit the model
fit <- sgscatm(W, C, K = K, lambda = 1)

# Print a brief summary
print(fit)
summary(fit)
```

## Main output

The fitted object is of class `"sgscatm"` and contains:

| Component | Dimensions | Description |
|-----------|-----------|-------------|
| `Pi` | M × K | Document-topic proportions (rows sum to 1) |
| `Phi` | K × N | Topic-term loading matrix |
| `Bz` | P × (K-1) | ILR path coefficients (covariate effects) |
| `Z` | M × (K-1) | ILR topic scores per document |
| `eigenvalues` | K-1 | Top eigenvalues of the augmented similarity matrix |

## Key functions

```r
# Top terms for each topic
vocab <- paste0("word", seq_len(N))
top_terms(fit, n = 10, vocab = vocab)

# Path coefficients (covariate → topic)
coef(fit)

# Fitted topic proportions
fitted(fit)

# Predict topic proportions for new documents
newW <- matrix(rpois(10 * N, 5), 10, N)
predict(fit, newW)
```

## Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `W` | — | M × N document-term matrix (non-negative) |
| `C` | — | M × P covariate matrix |
| `K` | — | Number of topics (≥ 2) |
| `lambda` | `1` | Regularisation weight balancing reconstruction and covariate structure |
| `r` | `min(M,N,100)` | Rank of the truncated SVD |
| `scale_W` | `TRUE` | Row-normalise W to term frequencies before fitting |
| `rotate` | `TRUE` | Apply varimax rotation to improve topic interpretability |

## Vignettes

```r
vignette("introduction", package = "sgscatm")
vignette("poliblog5k",   package = "sgscatm")
```

## Simulation study

The repository includes a full simulation study (`scripts/simulation/run_simulation.R`) with four blocks:

| Block | Purpose | Theoretical backing |
|-------|---------|---------------------|
| 1 | Consistency & asymptotic normality of B̂z | Thm 12 & 14 |
| 2 | Linearisation error bound | Prop 15 |
| 3 | Structural comparison with STM | empirical |
| **4** | **Four-method comparison: quality & speed** | **empirical** |

### Block 4 — Four-method comparison

Block 4 benchmarks four estimators on the same synthetic corpora across two corpus sizes (M = 500, 2 000 in quick mode; up to M = 20 000 in full mode) and three scenarios:

| Scenario | Description |
|----------|-------------|
| `high_sep` | Very distinct topics (`alpha_beta = 0.01`), moderate covariate signal |
| `low_sep` | Overlapping topics (`alpha_beta = 1.0`), moderate covariate signal |
| `weak_sig` | Medium topics, weak covariate signal (`b_max = 0.10`) |

**Methods compared:**

| Method | Description |
|--------|-------------|
| `sgscatm` | Spectral baseline — closed-form, zero EM iterations |
| `sgscatm+refine` | Spectral + one k-means M-step to re-estimate Φ |
| `stm` | Variational EM with default Spectral initialisation |
| `stm+warm` | Variational EM warm-started from the `sgscatm` solution |

**Metrics:** MSE(B̂z) after Procrustes alignment, MSE(Φ̂) under optimal row permutation, wall time, and EM iterations to convergence.

**Selected results (quick mode, 20 replicates, K = 5, P = 3, N = 500):**

| Scenario | M | Method | MSE(B̂z) | MSE(Φ̂) | Time (s) | Iter |
|----------|---|--------|---------|--------|---------|------|
| High sep. | 500 | sgscatm | 0.0222 | 1.85e-02 | **0.1** | — |
| | | sgscatm+refine | 0.0222 | **2.08e-04** | **0.1** | — |
| | | stm | 0.0101 | 4.43e-05 | 1.2 | 30 |
| | | **stm+warm** | **0.0015** | 1.67e-05 | 1.7 | 35 |
| Low sep. | 500 | stm+warm | 0.0052 | 7.73e-06 | **0.8** | **10** |
| Weak sig. | 500 | **sgscatm** | **0.0013** | 1.04e-03 | **0.1** | — |
| | | stm | 0.0087 | 1.68e-05 | 0.9 | 18 |

**Key findings:**
- `stm+warm` achieves the lowest MSE(B̂z) in most scenarios; the sgscatm warm start gives STM a ~4× accuracy advantage over vanilla STM in high-separation settings.
- `sgscatm+refine` reduces MSE(Φ̂) by ~100× over the baseline at negligible additional cost.
- `sgscatm` is the clear winner under **weak signal** (MSE(B̂z) = 0.0013 vs 0.0087 for STM), because the spectral closed-form solution is robust where variational EM struggles.
- `stm+warm` halves EM iterations in low-separation settings (10 vs 22), giving a measurable speed advantage on large corpora.

To run the simulation:

```bash
# Quick mode (M = 500 / 2 000, 20 replicates)
Rscript scripts/simulation/run_simulation.R --block 4 --quick

# Full mode (M = 500 / 1 000 / 5 000 / 20 000, 50 replicates)
Rscript scripts/simulation/run_simulation.R --block 4
```

## License

GPL-3
