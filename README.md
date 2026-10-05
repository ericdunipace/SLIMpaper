# Interpretable Model Summaries Using the Wasserstein Distance

Code to reproduce the analyses and figures in Dunipace, E. and Trippa, L. (2020). *Interpretable Model Summaries Using the Wasserstein Distance.* <https://arxiv.org/abs/2012.09999>

This repository is an R package (`SLIMpaper`) that holds the helper functions used by the analyses, plus the scripts that run the simulations and data analyses and draw the paper's figures. The methods themselves are in the [`WpProj`](https://github.com/ericdunipace/WpProj) package.

## Setup

1. **Install the package.** Clone the repository and install it, which also installs `WpProj` and the other dependencies:

   ```r
   # install.packages("devtools")
   devtools::install("SLIMpaper")   # from the folder containing the clone
   # or: devtools::install_github("ericdunipace/SLIMpaper")
   ```

2. **Run everything from the repository root.** All scripts read and write paths relative to it (`Output/`, `inst/figure/`, ...).

3. **Optional: MOSEK.** The scripts default to the commercial [MOSEK](https://www.mosek.com/) solver through `Rmosek` (free academic licenses are available). Without it:
   - in `vignettes/Simulation.R`, set `solver <- "cone"` to use the free ECOS solver instead;
   - the Lagrangian binary program timings and the ovarian cluster step (`ovar_cluster.R`) call MOSEK directly and need it.

4. **Python (binomial simulation only).** The neural network in the binomial simulation runs in Python through `reticulate`. If `python.path` is `NULL` in `vignettes/Simulation.R`, the script installs Python 3.10.16 and creates a virtual environment named `SLIM` with `numpy 2.2.1`, `torch 2.5.1` and `scipy 1.15.1`.

Intermediate results (`Output*/` folders) are not included in the repository because of their size. Every figure can be regenerated from the scripts below, except the GBM figures (see below).

## Which script makes which figure

| Figure(s) | Script(s), in order | Figure files |
|---|---|---|
| Toy example: selection order and ridge plots | `R/Simulations/single_indiv.R` | `inst/figure/toy_eg/` |
| Nonlinear approximation | `inst/figure/figure_code/nonlinear_approx.R` | `inst/figure/nonlinear_approx.pdf` |
| Simulations (Gaussian and binomial) | `vignettes/Simulation.R` → `vignettes/CombineSimulation.R` | `inst/figure/simulation/` |
| Timing comparison (Supplement) | `vignettes/SimulationTimings.R` → `vignettes/CombineTiming.R` | `inst/figure/timing/` |
| Ovarian cancer analysis (Supplement) | `R/DataAnalysis/Ovarian/`: `ovar_estimate.R` → `ovar_cluster.R` → `ovar_interp.R` | `inst/figure/applied/Ovar/` |
| Glioblastoma (GBM) analysis | `R/DataClean/` → `R/DataAnalysis/GBM/`: `gbm.R` → `gbm_cluster.R` → `gbm_interp.R` | `inst/figure/applied/GBM/` |

## Re-running each analysis

### Toy example and nonlinear approximation

Each is a single script that writes its figures directly:

```r
source("R/Simulations/single_indiv.R")
source("inst/figure/figure_code/nonlinear_approx.R")
```

### Simulations

1. **Run the simulations** with `vignettes/Simulation.R`. Settings are at the top of the file:
   - `families` chooses which simulations run. The script currently sets `families <- "gaussian"` just after listing both; change it to `c("gaussian", "binomial")` (or `"binomial"`) to run the binomial simulation.
   - Each family runs 100 replicates for each predictor correlation (0, 0.5, 0.9), in parallel on up to 8 cores. Seeds come from `data/seed_array.rda` (made by `data-raw/seeds.R`), so results are reproducible.
   - Gaussian runs use n = 1,024. Binomial runs fit a neural network with n = 131,072 and take much longer.
   - Results are written to `Output/<family>/mcp.net/none/exact/Corr_<correlation>/<n>/21/`.

2. **Make the figures** by running `vignettes/CombineSimulation.R` from top to bottom. It reads only result files dated after `date` (set at the top of the file), so set `date` to just before your simulation run. It writes the three simulation figures to `inst/figure/simulation/`.

### Timing comparison

1. `vignettes/SimulationTimings.R` times the binary program, the Lagrangian binary program (`R/lbp.R`, which requires MOSEK) and the relaxed binary program across sample sizes, numbers of covariates and numbers of posterior samples. Results go to `Output_timing/`.
2. `vignettes/CombineTiming.R` reads every file in `Output_timing/` and writes the three figures in `inst/figure/timing/`.

### Ovarian cancer analysis

The data are in the package (`data(ovar)`). They were built by `data-raw/ovar.R` from the TCGA data in the Bioconductor package `curatedOvarianData`, so you don't need to download anything.

1. `ovar_estimate.R` fits the Cox model and saves it to `Output/Ovar/recurrence_cox.rds`.
2. `ovar_cluster.R` fits the interpretable models. It is written as a SLURM array job (one task per penalty value, read from `SLURM_ARRAY_TASK_ID`) and needs MOSEK. Each task saves `Output/Ovar/global_model_<task>.RDS` and `local_model_<task>.RDS`.
3. Combine the task results into `Output/Ovar/global_model_cluster.RDS` and `local_model_cluster.RDS`. The code for this step is the commented-out `combine.cluster()` block near the top of `ovar_interp.R` (around line 147); uncomment it when running from scratch.
4. `ovar_interp.R` evaluates the models and writes the three figures in `inst/figure/applied/Ovar/`.

### Glioblastoma (GBM) analysis

**The GBM data are not included**: they come from a clinical database that can't be shared publicly. The scripts expect the data files in `../Data/GBM/` (a `Data` folder next to the repository), so these figures can only be regenerated with access to that data.

1. `R/DataClean/clean_gbm_df.R` and `clean_gbm_pubmed.R` prepare the data.
2. `gbm.R` fits the BART survival model, saves its predictions to `Output/GBM/`, and draws `bs_time_gbm.pdf` and `intbs_gbm.pdf`.
3. `gbm_cluster.R` fits the interpretable models as a SLURM array job and saves one file per task to `Output/GBM/cluster/`.
4. `gbm_interp.R` combines the task results and draws `heat_total.pdf`, `rank_total.pdf` and `oddsratio_gbm.pdf`.

## Note on file names

Binomial result files created before October 2026 have `1024` in their file names even though they were run with n = 131,072. The folder name (`.../131072/21/`) is correct, and the combine script uses the folder name, so this doesn't affect the figures.
