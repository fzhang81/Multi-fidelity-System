# Multi-fidelity Decision Analysis

This repository contains the R code accompanying **A Decision Analysis Framework for High-fidelity and Low-fidelity Systems with Applications in Manufacturing Processes**. The code is organized around the three numerical studies in the manuscript: the illustrative example in Section 5, the curing process case study in Section 6, and the injection molding case study in Section 7.

The analyses combine low- and high-fidelity information through Gaussian process (GP) models, estimate the scaling parameter and model discrepancy, and use predictive simulation to obtain distributions of selected optimal inputs. The files retain the term `calibration` used in the manuscript and in the original implementation.

## What can be run with the public files?

Not all raw data and intermediate results are included. The following distinction is important when reproducing the analyses.

| Analysis | Files available in this repository | How to proceed |
| --- | --- | --- |
| Section 5: RMSE and bias plots across 50 datasets | Four summary CSV files and the plotting script | Run `simulation_boxplot.R`. No model fitting or simulation workspaces are needed. |
| Section 5: LF-only, HF-only, and multi-fidelity distributions of minima locations | Data-generation, analysis, and plotting scripts | Generate the 50 paired datasets, run the relevant analysis, and load one resulting `.RData` file before running `simulation_histogram.R`. The datasets and saved workspaces are not included. |
| Section 6: curing calibration and distributions of minima locations | Both generated curing datasets and the analysis script | Run `Curing_calibration.R` using the supplied CSV files. The original simulation data are not needed for this step. |
| Section 6: regeneration of the curing datasets | Data-generation script only | Obtain the original curing data from the corresponding author before running `Curing_generatedata_nested.R`. |
| Section 7: injection molding calibration and plots | Analysis and plotting code only | Obtain the injection molding datasets from the corresponding author and update the local paths before running the script. |

## Repository contents

```text
.
|-- simulation_generatedata.R
|-- simulation_GPonly.R
|-- simulation_calibration.R
|-- simulation_histogram.R
|-- simulation_boxplot.R
|-- x1_med.csv
|-- x2_med.csv
|-- x1_mse.csv
|-- x2_mse.csv
|-- Curing_generatedata_nested.R
|-- Curing_calibration.R
|-- Curing_44DoCS_noise0_lhd200.csv
|-- Curing_33uF_noise01_lhd50_nested.csv
`-- InjectionMolding_calibration.R
```

Files beginning with `simulation_` correspond to Section 5. Files beginning with `Curing_` correspond to Section 6. `InjectionMolding_calibration.R` corresponds to Section 7.

The `Calidata/` directory and the `.RData` workspaces described below are generated locally; they are not part of the public file set.

## Software and working directory

The scripts are written in R and can be run from an R console or RStudio. Install the packages loaded or called by the scripts:

```r
packages <- c(
  "laGP", "plgp", "lhs", "mvtnorm", "Matrix",
  "MASS", "knitr", "ggplot2", "tidyr", "gridExtra"
)
missing_packages <- setdiff(packages, rownames(installed.packages()))
if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}
```

Set the working directory to the repository root, where the R scripts and supplied CSV files are located:

```r
setwd("/path/to/Multi-fidelity-System-main")
```

Use a fresh R session for each main analysis. Before running any GP calibration script, load the following packages:

```r
library(laGP)
library(plgp)
library(lhs)
library(mvtnorm)
```

This setup is necessary because `Curing_calibration.R` does not contain its own package-loading statements, and the calibration scripts call `rmvnorm()` without explicitly loading `mvtnorm`. Other packages are loaded within the scripts that use them.

The repository does not include an R version record or a package-version lockfile. To document the environment used for a new run, save the session information after loading the required packages:

```r
capture.output(sessionInfo(), file = "sessionInfo.txt")
```

## Section 5: Illustrative example

### Recreate the RMSE and bias plots from the supplied summaries

This is the shortest route to the 50-dataset comparison plots:

```r
source("simulation_boxplot.R")
```

The script reads all four summary files directly from the repository root. Each file has 50 rows, one per dataset, and three columns. The column order is the same in all four files:

| CSV column | Method |
| --- | --- |
| `V1` | Low-fidelity GP only |
| `V2` | High-fidelity GP only |
| `V3` | Multi-fidelity calibration |

`x1_mse.csv` and `x2_mse.csv` contain the saved dataset-level MSE summaries for the two input coordinates. The plotting script takes their square roots to obtain RMSE values. `x1_med.csv` and `x2_med.csv` contain the corresponding median selected input locations. The bias plots display the signed errors of those medians relative to the high-fidelity optimum, `(-0.8, 0.4)`. Thus, each boxplot summarizes 50 dataset-level values; the bias panels do not display a single average bias estimate.

**Preprocessing in the distributed plotting script.** Before taking square roots, the script replaces high-fidelity MSE values greater than `0.3` with the median of the corresponding high-fidelity MSE column. This affects 5 entries in `x1_mse.csv` and 11 entries in `x2_mse.csv`. The CSV files themselves are not modified. The bias panels also use `scale_y_continuous(limits = c(-0.8, 0.8))`; values outside those limits are excluded from the plotted statistical summaries. The distributed median files contain one such high-fidelity value for each coordinate. These operations should be taken into account when interpreting the figures. To examine the complete, unmodified summaries, remove the two MSE-replacement statements and the fixed bias-axis limits before plotting.

The script draws the RMSE and bias figures sequentially. To save them in a two-page PDF rather than only displaying them in the graphics window:

```r
pdf("simulation_rmse_bias.pdf", width = 10, height = 5)
source("simulation_boxplot.R")
dev.off()
```

The supplied CSVs are sufficient for these plots, but they do not contain the individual selected optimal inputs needed for the distributions of minima locations.

### Generate the 50 paired datasets

Create the output directory before running the generator:

```r
dir.create("Calidata", showWarnings = FALSE)
source("simulation_generatedata.R")
```

The script generates 50 low-/high-fidelity dataset pairs. Each pair contains 200 low-fidelity observations and 50 high-fidelity observations at separately generated Latin hypercube designs in `[-1, 1]^2`. The illustrative-example designs are not constructed to be nested.

The functions implemented in the distributed generator are

```text
y_L(x) = (x1 + 0.6)^2 + (x2 - 0.2)^2 + 5
y_H(x) = (x1 + 0.8)^2 + (x2 - 0.4)^2 + 5
Y_H(x) = y_H(x) + epsilon,  epsilon ~ N(0, 0.2^2)
```

No noise is added to the low-fidelity response. These expressions describe the current code, including the common offset of `5` and the high-fidelity noise standard deviation of `0.2`.

For dataset index `i = 1, ..., 50`, the generator writes:

```text
Calidata/2025618M_p5_noise02_i.csv    # Low fidelity:  x1, x2, ys
Calidata/2025618F_p5_noise02_i.csv    # High fidelity: x1, x2, yp
```

Here, `M` denotes low-fidelity model data and `F` denotes high-fidelity data. The numeric prefixes are retained as filename identifiers.

### Run the single-fidelity baselines

After generating the datasets, run:

```r
# This script does not set its own seed.
# The following is an example seed for a reproducible new run.
set.seed(123)
source("simulation_GPonly.R")
```

The script processes both `M` and `F` datasets, fitting a separate GP to each dataset without combining the two fidelities. For each fit, it performs 10,000 decision replications, each using 500 Latin hypercube candidate points. It samples the GP response jointly over the candidate set, squares the sampled responses, and selects the candidate with the smallest squared response.

The script saves 100 workspaces in the repository root, one for each dataset and fidelity:

```text
GP_2025618M_p5_noise02_i_u100rep100lhd500.RData
GP_2025618F_p5_noise02_i_u100rep100lhd500.RData
```

The suffix is a filename convention; the single-fidelity script itself does not sample a calibration parameter.

### Run the multi-fidelity analysis

In a fresh session, set the working directory and load the GP packages listed above, then run:

```r
source("simulation_calibration.R")
```

For each of the 50 dataset pairs, the script fits a low-fidelity GP and estimates the scaling parameter using a discrepancy GP fitted to high-fidelity residuals. It then performs leave-one-out refitting of the high-fidelity observations while retaining the fitted low-fidelity emulator. A density estimate is constructed from the resulting scaling-parameter estimates, and 100 values are selected from the density grid with probability weights proportional to the estimated density, without replacement.

For each selected scaling value, the code performs 100 decision replications with 500 candidate inputs. Joint predictive draws are used to select the candidate minimizing the squared response. This produces 10,000 selected input locations per dataset.

The resulting workspaces are saved as:

```text
cali_2025618F50M200_p5_noise02_i_u100rep100lhd500.RData
```

The script also contains a final inline plotting block. Its first histogram uses breaks ending at `0`, which can cause a plotting error when selected `x1` values exceed `0`. The workspaces are saved before that block. For plotting a completed run, use `simulation_histogram.R` as described next; a plotting error at the end does not require rerunning the saved analyses.

The analysis scripts save individual workspaces, but the repository does not include a separate aggregation script that reconstructs the four summary CSVs from those workspaces. Use the supplied summary CSVs for the RMSE and bias plots.

### Plot distributions of minima locations

`simulation_histogram.R` is a plotting-only script. First load the workspace for the method and dataset to be displayed, then run the script. For example, after generating the workspaces for dataset 1:

```r
# Choose ONE workspace at a time:

# Low-fidelity GP only
load("GP_2025618M_p5_noise02_1_u100rep100lhd500.RData")
source("simulation_histogram.R")

# High-fidelity GP only
load("GP_2025618F_p5_noise02_1_u100rep100lhd500.RData")
source("simulation_histogram.R")

# Multi-fidelity calibration
load("cali_2025618F50M200_p5_noise02_1_u100rep100lhd500.RData")
source("simulation_histogram.R")
```

Replace the dataset index to display another run. Dataset 1 is an example, not an identification of the run used for the manuscript figures.

The plotting script uses `x1_recov`, `x2_recov`, `med_x1`, and `med_x2` from the loaded workspace. It produces two marginal histograms with density curves and median markers, corresponding to the three types of distribution plots in Section 5. Each figure describes one method on one dataset, rather than pooling the 50 datasets.

The default display ranges are `[-1, 0]` for `x1` and `[0, 1]` for `x2`, although the histogram breaks cover `[-1, 1]`. Set both `xlim` arguments to `c(-1, 1)` to view the full input domain.

## Section 6: Curing process simulation

### Run the analysis with the supplied nested datasets

The repository includes both generated datasets required by `Curing_calibration.R`:

| File | Rows | Response |
| --- | --- | --- |
| `Curing_44DoCS_noise0_lhd200.csv` | 200 | Low-fidelity Degree of Cure (DoC) |
| `Curing_33uF_noise01_lhd50_nested.csv` | 50 | High-fidelity deformation response with added noise |

Both files have columns `T1`, `T2`, `t1`, `t2`, and `y`. The 50 high-fidelity input locations coincide with the first 50 low-fidelity input locations. The CSVs contain generated responses from polynomial models fitted to the original simulation data, rather than the original finite element outputs themselves.

In a fresh session, set the working directory, load the GP packages, and run:

```r
source("Curing_calibration.R")
```

The script estimates the scaling parameter, performs leave-one-out refitting, and obtains 100 density-weighted scaling values. For each value, it carries out 100 decision replications over 200 candidate inputs. As implemented, selection minimizes the squared sampled scalar response. The selected inputs are transformed back to the scale of the supplied data.

The script saves `Results_NestedData_Calibration.RData` and draws four marginal distributions for `T1`, `T2`, `t1`, and `t2`, corresponding to the curing distribution plot in Section 6. The estimated scaling parameter and leave-one-out summaries are printed during the run and retained in the workspace.

### Regenerate the nested datasets from the original data

`Curing_generatedata_nested.R` documents the construction of the two supplied CSV files. This is an optional upstream step and is not needed to analyze the public datasets.

The script requires `L_Q44_complete.csv` and `L_Q33_complete.csv`, which are not included. The former must contain `T1`, `T2`, `t1`, `t2`, and `DoC`; the latter must contain those input variables together with `DoC` and the deformation response `u`. In these raw data, the column named `u` is a response variable, not the scaling parameter subsequently estimated by the calibration script.

The generator retains original observations with `DoC > 0.96`, fits separate second-order polynomial models with stepwise AIC selection for the two responses, and constructs a 200-point nested Latin hypercube design with a 50-point high-fidelity subset. The design is mapped to coordinate-wise bounds obtained from the retained original inputs. It adds Gaussian noise with standard deviation `0.1` to the high-fidelity predictions and no noise to the low-fidelity predictions.

After obtaining the original data and placing them in the working directory, run:

```r
source("Curing_generatedata_nested.R")
```

This overwrites the two generated curing CSV files. The `DoC > 0.96` filter is applied to the original training records; the script does not apply an additional DoC filter to the generated responses or to the later optimization candidates.

## Section 7: Injection molding

`InjectionMolding_calibration.R` contains the calibration, decision analysis, marginal plotting, and pairwise analysis for the injection molding case study. It requires two files that are not distributed with the repository:

```text
sim.csv    # Low-fidelity simulation data
exp.csv    # High-fidelity experimental data
```

The script selects columns by position. In both files, the first four columns must contain mold temperature, injection speed, packing pressure, and packing time, in that order. Columns 5-8 must contain the four responses, with names `Horizontal_Left`, `Horizontal_Right`, `Vertical_Up`, and `Vertical_Down`, in the same order in both datasets.

### Configure the paths and run

Before running the script, obtain the datasets from the corresponding author and place them in the repository root. Update the two machine-specific path statements in the script:

```r
# Remove or comment out this statement, because the working directory
# should already be set to your local repository:
# setwd("~/Downloads/calibration")

# Replace the absolute load() path later in the script with:
load("NEW_InjectionMolding_cali_y1234.RData")
```

The `load()` statement reloads a workspace saved immediately before it; it does not refer to a file that must be obtained separately.

Then, in a fresh session with the GP packages loaded, run:

```r
source("InjectionMolding_calibration.R")
```

The script fits a separate multi-fidelity model for each of the four responses. At each decision replication, it draws the four response vectors independently conditional on their fitted models and selects the candidate minimizing their sum of squares. There are 100 sets of sampled scaling values, 100 decision replications per set, and 200 candidates per replication.

The outputs include `NEW_InjectionMolding_cali_y1234.RData`, four marginal distributions of selected process inputs, and pairwise visualizations. The final pairwise scatterplot matrix displays Spearman rank correlations in its upper panels and is written to `Figure7_pairwise_xopt.png`. This filename is retained from the script; the pairwise plot corresponds to Figure 8 in the revised manuscript, while the marginal distributions correspond to Figure 7. The script also evaluates the pairwise correlation matrix at the end; run that final expression interactively to display it in the console.

The workspace is saved before the recovered coordinates and plotting objects are created. When reopening it for plotting, execute the code beginning with the construction of `post_x` after the `load()` statement, not only the final plotting commands.

Without `sim.csv` and `exp.csv`, the analysis cannot run. Updating the paths alone does not resolve the missing-data requirement.

## Computational and reproducibility notes

The descriptions above refer to the scripts as distributed. In the scalar-response analyses, the implemented decision objective is the **squared predictive draw**; in injection molding, it is the **sum of four squared predictive draws**. Candidate points are drawn within coordinate-wise ranges of the observed inputs: each dataset's ranges for a single-fidelity fit, and the combined low-/high-fidelity ranges for a multi-fidelity fit.

Scaling values are generated by density-weighted sampling from leave-one-out estimates, not by an MCMC sampler. The scripts retain full covariance matrices when generating predictive response draws over a candidate set.

The scripts contain the following random-seed settings:

| Script | Seed |
| --- | --- |
| `simulation_generatedata.R` | `2025619` |
| `simulation_GPonly.R` | Not set within the script |
| `simulation_calibration.R` | `20256201` |
| `Curing_generatedata_nested.R` | `10092026` |
| `Curing_calibration.R` | `20250734` |
| `InjectionMolding_calibration.R` | `20250815` |

Full simulation batches involve repeated GP fitting and multivariate predictive sampling. They can require substantial computation and storage, particularly because `save.image()` saves the entire R workspace, including stored candidate sets. Existing files with the same output names may be overwritten. Use separate fresh sessions and retain copies of results that need to be preserved.

The provided summary CSVs allow the Section 5 comparison plots to be recreated without repeating those computations. Rerunning the simulation pipeline produces new workspaces; it does not establish an exact match to a particular manuscript figure without the corresponding original workspace, random-number state, and software environment.

## Data access and correspondence

The original curing simulation data and the injection molding datasets are not included in this repository. For questions about the code or access to those data, please contact the manuscript's corresponding author:

**Qiong Zhang**  
School of Mathematical and Statistical Sciences, Clemson University  
Email: `qiongz@clemson.edu`
