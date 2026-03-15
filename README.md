# egscatm

**Structural Topic Modeling via ILR-EGSCA**

`egscatm` fits structural topic models that operate natively in the Aitchison geometry of the simplex. Document-topic proportions are represented in isometric log-ratio (ILR) coordinates and linked to covariates through a Generalized Structured Component Analysis (GSCA) formulation. Estimation reduces to a single eigenvalue decomposition of a covariate-augmented document similarity matrix — yielding a non-iterative, globally optimal solution.

## Requirements

- R >= 4.0
- Imports: `Matrix`, `methods`
- Suggests (optional): `testthat >= 3.0`, `knitr`, `rmarkdown`

## Installation

### From GitHub (development version)

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("marcoortu/irl-egsca-tm")
```

### From source

Clone the repository and build/install from the project root:

```bash
git clone https://github.com/marcoortu/irl-egsca-tm.git
cd irl-egsca-tm
R CMD build .
R CMD INSTALL egscatm_0.1.0.tar.gz
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
R CMD check egscatm_0.1.0.tar.gz

# Or with devtools
devtools::check()
```

## Quick start

```r
library(egscatm)

set.seed(42)
M <- 200   # documents
N <- 100   # vocabulary size
K <- 4     # number of topics
P <- 2     # number of covariates

# Simulate a document-term matrix and covariate matrix
W <- matrix(rpois(M * N, 5), M, N)
C <- scale(matrix(rnorm(M * P), M, P))

# Fit the model
fit <- egscatm(W, C, K = K, lambda = 1)

# Print a brief summary
print(fit)
summary(fit)
```

## Main output

The fitted object is of class `"egscatm"` and contains:

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
vignette("introduction", package = "egscatm")
vignette("poliblog5k",   package = "egscatm")
```

## License

GPL-3
