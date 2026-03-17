# Replication Package

This folder contains all code and data required to reproduce the results in the companion paper.

## Structure

```
replication/
├── simulation/          ← scripts for the simulation study (Sections 4–5)
│   ├── run_simulation.R         ← main entry point
│   ├── sim_dgp.R                ← data-generating process
│   ├── sim_utils.R              ← utility functions
│   ├── block3_stm_worker.R      ← parallel worker for Block 3
│   ├── block4_comparison_worker.R ← parallel worker for Block 4
│   ├── compare_stm.R            ← STM comparison (poliblog)
│   ├── compare_stm_5k.R         ← STM comparison (5k corpus)
│   └── demo_poliblog.R          ← demo on poliblog data
├── application/         ← scripts for the BES real-data application (Section 6)
│   ├── run_analysis.R           ← main entry point for BES analysis
│   ├── bes_analysis.R           ← full analysis pipeline
│   ├── calibrate_K_lambda.R     ← hyperparameter calibration (K, lambda)
│   ├── check_se.R               ← standard error diagnostics
│   └── explore_dtm.R            ← exploratory analysis of the DTM
├── data/                ← pre-processed data (not raw — see below for raw data)
│   ├── bes_w25_dtm.rds          ← document-term matrix (Wave 25)
│   └── bes_w25_filtered.csv     ← filtered corpus metadata
└── output/              ← generated figures and tables (reproducible from scripts)
```

## Requirements

- R >= 4.0
- Package `sgscatm` (install from the root of this repository or from the `.tar.gz`)
- Additional packages: `stm`, `quanteda`, `ggplot2`, `dplyr`, `parallel`, `haven`

Install all dependencies:

```r
install.packages(c("stm", "quanteda", "ggplot2", "dplyr", "parallel", "haven"))
# install sgscatm from root:
devtools::install("..")
# or from the tarball:
install.packages("../sgscatm_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Raw data

The raw BES (British Election Study) survey data is not included due to file size. It can be downloaded from:

- [British Election Study](https://www.britishelectionstudy.com/data-objects/panel-study-data/) — Wave 25 / Wave 30

Place the downloaded files in `application/` before running the scripts.

## Reproducing the simulation study

```bash
# All blocks, quick mode (M = 500/2000, 20 replicates)
Rscript simulation/run_simulation.R --quick

# Block 4 only, full mode (50 replicates)
Rscript simulation/run_simulation.R --block 4

# Output: replication/output/figures/ and replication/output/tables/
```

## Reproducing the BES application

```bash
# Step 1: calibrate K and lambda
Rscript application/calibrate_K_lambda.R

# Step 2: run the full analysis
Rscript application/run_analysis.R

# Output: replication/output/bes/
```
