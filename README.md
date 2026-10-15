# Multi-fidelity Decision Analysis

R code for **A Decision Analysis Framework for High-fidelity and Low-fidelity Systems with Applications in Manufacturing Processes**. The repository covers the illustrative example (Section 5), curing process simulation (Section 6), and injection molding (Section 7).

The analyses use Gaussian process models, leave-one-out refitting of the scaling parameter, and predictive simulation to obtain distributions of selected optimal inputs.

## Setup

Install the required R packages:

```r
install.packages(c("laGP", "plgp", "lhs", "mvtnorm", "Matrix",
                   "MASS", "knitr", "ggplot2", "tidyr", "gridExtra"))
```

Set the working directory to the repository root. Use a fresh R session for each main analysis and load the following packages before running a calibration script:

```r
library(laGP)
library(plgp)
library(lhs)
library(mvtnorm)
```

Simulation settings can be changed in the scripts. Full simulation batches can require substantial computation and storage.

## Section 5: Illustrative example

| Script | Purpose |
| --- | --- |
| `simulation_generatedata.R` | Generates 50 dataset pairs, each with 200 low-fidelity and 50 high-fidelity observations. |
| `simulation_GPonly.R` | Fits separate low-fidelity-only and high-fidelity-only GPs to all 50 datasets and saves the optimization results. |
| `simulation_calibration.R` | Runs the multi-fidelity analysis on all 50 dataset pairs and saves the optimization results. |
| `simulation_histogram.R` | Plots distributions of minima locations from a loaded result workspace. |
| `simulation_boxplot.R` | Reads the four supplied summary CSVs to plot RMSE and signed errors of the median selected inputs. |

### RMSE and bias plots

The repository includes `x1_mse.csv`, `x2_mse.csv`, `x1_med.csv`, and `x2_med.csv`. Each contains 50 rows, with columns ordered as low fidelity, high fidelity, and multi-fidelity calibration. The MSE files are converted to RMSE by taking square roots; the median-location files are compared with the true high-fidelity optimum, `(-0.8, 0.4)`.

```r
source("simulation_boxplot.R")
```

These plots require no saved workspaces or model fitting.

### Model fitting and distributions of minima locations

The 50 generated dataset pairs and saved `.RData` files are not included. To create new results, run:

```r
dir.create("Calidata", showWarnings = FALSE)
source("simulation_generatedata.R")
source("simulation_GPonly.R")
```

Then run `simulation_calibration.R` in a fresh session after loading the packages above. Each method produces 10,000 selected input locations per dataset, using 500 candidates per decision replication.

The analysis scripts save workspaces in the repository root. For dataset index `i`, the filenames are:

```text
GP_2025618M_p5_noise02_i_u100rep100lhd500.RData        # Low fidelity
GP_2025618F_p5_noise02_i_u100rep100lhd500.RData        # High fidelity
cali_2025618F50M200_p5_noise02_i_u100rep100lhd500.RData # Multi-fidelity
```

Load the workspace for the desired method and dataset before plotting. For example:

```r
load("GP_2025618F_p5_noise02_1_u100rep100lhd500.RData")
source("simulation_histogram.R")
```

Replace the filename to display another result. Each histogram figure represents one method and one dataset, not all 50 datasets pooled together. The supplied summary CSVs cannot reconstruct these histograms. The repository does not include the aggregation step used to export those four CSVs.

The final inline plotting block in `simulation_calibration.R` may fail if selected `x1` values fall outside its histogram breaks. Results are saved before that block; use `simulation_histogram.R` to plot a saved workspace. Its `xlim` settings can be changed to `c(-1, 1)` to display the full domain.

## Section 6: Curing process simulation

`Curing_generatedata_nested.R` constructs nested datasets by subsampling the original simulation records in `L_Q44_complete.csv`, which is not included. It first selects 200 distinct input settings and then 50 settings from within this set, using repeated random subsampling with a distance-based space-filling criterion. Both responses come from the same original dataset; no polynomial models are used to generate them.

The two subsampled datasets are supplied:

- `Curing_44DoC_subsample200.csv`: 200 low-fidelity observations.
- `Curing_44u_subsample50_nested.csv`: 50 high-fidelity observations whose input locations are included in the low-fidelity design.

Both files contain columns `T1`, `T2`, `t1`, `t2`, and `y`. The high-fidelity input settings match the first 50 rows of the low-fidelity file in the same order.

To run calibration and plot the distributions of selected inputs, load the packages listed above and run:

```r
source("Curing_calibration.R")
```

The script saves `Results_NestedData_Calibration.RData`. The original data are only needed to regenerate the subsamples; they are not required to run the analysis with the supplied CSV files.

## Section 7: Injection molding

`InjectionMolding_calibration.R` performs calibration for four responses, selects inputs by minimizing their summed squared predictive draws, and produces marginal and pairwise plots.

It requires `sim.csv` and `exp.csv`, which are not included in this anonymized repository. The script cannot run without these data. Before execution, update its `setwd()` and `load()` paths to match the local working directory.

With the required data available and packages loaded, run:

```r
source("InjectionMolding_calibration.R")
```

The script saves `NEW_InjectionMolding_cali_y1234.RData` and `Figure7_pairwise_xopt.png`; the latter corresponds to Figure 8 in the revised manuscript.

## Data availability

The public files support plotting the Section 5 summaries, generating new illustrative-example results, and running the Section 6 analysis with the supplied nested datasets. The original curing data and injection molding data are not distributed with this anonymized version. Please contact corresponding author.
