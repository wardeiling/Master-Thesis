# From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes

This repository contains all materials associated with the manuscript:
**"From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes"**
It includes code, supplementary documentation, and simulation results to ensure full transparency and reproducibility of the study.

## Ethics Assessment

This simulation study was approved by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. The approval is based on the documents submitted by the researchers as required by the Ethics Committee and filed under FETC number 24-2003. The approval is valid through 31 May 2025. The approval pertains to ethical considerations, data management, and privacy issues (including GDPR compliance).

## Study Design

This simulation study evaluates the generalizability of disaggregation methods—commonly applied in multilevel linear models (MLMs)—to *generalized* multilevel models (GLMMs) and *generalized estimating equations* (GEEs) in the context of binary predictors and/or outcomes.

We address two questions:

1. Can disaggregation methods (uncentered, centering-within-cluster, Mundlak's contextual model) reliably recover within-person and contextual effects in GLMMs with binary predictors and/or outcomes?
2. Do GEEs require explicit disaggregation to correctly estimate within-person effects in the presence of contextual effects?

We simulate data under four data-generating mechanisms (DGMs) that vary the scale of the predictor and outcome variables (binary or continuous). Across DGMs, the following parameters are held constant: the within-cluster SD of the continuous predictor, the fixed intercept, the within-cluster effect, and the level-1 residual SD (for DGMs with a continuous outcome). The table below summarizes the manipulated design factors:

| Factor                               | Levels    |
| ------------------------------------ | --------- |
| Sample size *(N)*                    | 100, 200  |
| Number of time points *(T)*          | 5, 10, 20 |
| Between-cluster SD in continuous *X* | 0, 1, 3   |
| SD in latent *Z* for binary *X*      | 0, 1, 3   |
| Contextual effect                    | 0, 1, 3   |
| Random intercept residual SD         | 1, 3      |

Each dataset is analyzed using 12 strategies: all combinations of 3 disaggregation methods and 4 estimation approaches (GLMM and GEE with independence, exchangeable, and AR(1) correlation structures). Model performance is evaluated in terms of estimation bias for the within-person effect ($\beta_1$) and the contextual effect ($\gamma_{01}$).

## Repository Structure

**`renv.lock`**

Contains dependency information for full reproducibility with `renv`.

### `scripts/`

Contains all core scripts for running and analyzing the simulation study.

* **`main-simulation-function-future-simul-part1.R`**
  Modularized script that runs simulations for DGMs 2–4.

* **`main-simulation-function-future-simul-part2.R`**
  Modularized script that runs simulations for DGM 1.

* **`results-plotting.R`**
  Produces plots used in the result section of the manuscript.

* **`helper-functions/`** (subfolder with modular components):

  * `data-generation-centeredX.R`: Data-generating mechanisms based on the hybrid model.
  * `data-generation-mundlak.R`: Data-generating mechanisms based on Mundlak’s contextual model (the model used in the main manuscript).
  * `model-fitting.R`: Model-fitting procedures for both GLMMs and GEEs.
  * `result-formatting.R`: Function to clean and format model-fitting output.

### `docs/`

Contains supporting materials that provide additional context and in-depth explanations of specific aspects of the study.

* **`data-exploration.qmd`** / **`.html`**
  Allows users to explore all four data-generating mechanisms (DGMs) considered in the study. [View HTML](https://wardeiling.github.io/multilevel-vs-gee-binary/data-exploration.html)

* **`supplementary_materials.qmd`** / **`.html`**
  Supplementary materials accompanying the manuscript [View HTML](https://wardeiling.github.io/multilevel-vs-gee-binary/supplementary_materials.html)
  1. A comparison between the hybrid and Mundlak's contextual model.
  2. Discussion of boundary/extreme estimates in GEEs and how they were handled.

### `output/`

Contains raw and processed simulation outputs, organized into subfolders corresponding to different simulation runs:

* **`April10_fullsimulation/`**: Part 1 of the simulations, covering DGMs 2–4.
* **`April17_fullsimulation_contxy/`**: Part 2 of the simulations, covering DGM 1.
* **`April18_fullsimulation_combined/figures/`**: Final figures used in result section of the manuscript.

Each folder includes:

* `i.RDS`: Raw output for each scenario *i*.
* `settings.RDS`: Simulation settings used for that part.
* `log.txt`: Logs containing warnings and errors during simulation.
* `summary-results-bias.RDS` and `.csv`: Summary files quantifying bias in the estimates.

### `renv/`

Contains internal `renv` files storing the project-specific package environment.

## Reproducing Results

1. Run `scripts/main-simulation-function-future-simul-part1.R` to simulate DGMs 2–4. Output is saved to `April10_fullsimulation/`.
2. Run `scripts/main-simulation-function-future-simul-part2.R` for DGM 1. Output is saved to `April17_fullsimulation_contxy/`.
3. Run `scripts/results-plotting.R` to process and visualize the simulation results used in the manuscript.

## Reproducibility: Step-by-Step Guide

This repository uses the [`renv`](https://rstudio.github.io/renv/) package to create a reproducible R environment. To replicate the computational setup and rerun the analyses:

### Step 1: Setup R and RStudio

1. Install **R version 4.2.2** from CRAN ([download link](https://cran.rstudio.com/bin/windows/base/old/4.2.2/R-4.2.2-win.exe))
2. Install **RStudio** (latest stable release)

### Step 2: Clone or Download the Repository

Clone the repository via GitHub or download the ZIP file and unzip it locally.

### Step 3: Restore the Project Environment via `renv`

1. Open **`Master-Thesis.Rproj`** with RStudio.
2. Run the following in the R console:

   ```r
   renv::restore()
   ```

This restores the exact package versions as specified in the `renv.lock` file, ensuring a consistent and reproducible computational environment.

### Step 4: Run the Simulation Scripts

1. Execute **`scripts/main-simulation-function-future-simul-part1.R`**

   * Runs simulations for DGMs 2–4
   * Outputs will be saved in `simulation_results/April10_fullsimulation/`

2. Execute **`scripts/main-simulation-function-future-simul-part2.R`**

   * Runs simulations for DGM 1
   * Outputs will be saved in `simulation_results/April17_fullsimulation_contxy/`

### Step 5: Reproduce the Figures

Run **`scripts/results-plotting.R`** to generate the plots used in the manuscript.

This script automatically collects and merges the simulation outputs from both parts and creates boxplots for each design condition.

---
