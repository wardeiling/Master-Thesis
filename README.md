# From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes

This repository contains all materials associated with the manuscript:
**"From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes"**
It includes code, supplementary documentation, and simulation results to ensure full transparency and reproducibility of the study.

## Ethics Assessment

This simulation study was approved by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. The approval is based on the documents sent by the researchers as requested in the form of the Ethics committee and filed under FETC number 24-2003. The approval is valid through 31 May 2025. The approval of the Ethical Review Board concerns ethical aspects, as well as data management and privacy issues (including the GDPR).

## Study Design

This simulation study evaluates the generalizability of disaggregation methods—commonly applied in multilevel linear models (MLMs)—to *generalized* multilevel models (GLMMs) and *generalized estimating equations* (GEEs) when dealing with binary predictors and/or binary outcomes.

We address two primary questions:

1. Can disaggregation methods (UC, CWC, MuCo) reliably recover within-person and contextual effects in GLMMs with binary predictors and/or outcomes?
2. Do GEEs require explicit disaggregation to correctly estimate within-person effects, especially when contextual effects are present?

We simulate data across four data-generating models (DGMs) that vary in the scale of the predictor and outcome variables (continuous or binary). Across DGMs, we keep constant the within-cluser standard deviation (SD) of the continuous predictor, the fixed intercept, the within-cluster effect and the level 1 residual SD (for DGMs with continuous outcome). For each of these DGMs, we systematically vary: 

CREATE MARKDOWN TABLE

* Sample size *(N = 100, 200)*
* Number of time points *(T = 5, 10, 20)*
* Between-cluster SD in continuous predictor (0, 1, 3)
* SD in Z (the latent trait underlying between-person variability in binary X): (0, 1, 3)
* Contextual effect: (0, 1, 3)
* Random intercept residual SD: (1, 3)

Each dataset is analyzed using 12 strategies: combinations of 3 disaggregation methods (uncentered, centering-within-clusters and mundlak's contextual model) and 4 estimation approaches (GLMM, and GEE with independence, exchangeable, and AR(1) correlation structures).

Model performance is assessed via estimation bias in fixed effects (within-person: β₁; contextual: γ₀₁).

---

## Repository Structure

**`renv.lock`**

Contains information on the requirements of all the dependencies of R-packages used in the simulation study.

### `scripts/`

Contains all core scripts for running and analyzing the simulation study.

* **`main-simulation-function-future-simul-part1.R`**
  Modularized main script that executes the simulation across various design conditions for DGM 2, 3 and 4
* **`main-simulation-function-future-simul-part2.R`**
  Modularized main script that executes the simulation across various design conditions for DGM 1.

* **`results-plotting.R`**
  Scripts used to produce the plots featured in the Results section of the manuscript.

* **`helper-functions/`** (subfolder with modular components):

  * `data-generation-centeredX.R`: Data-generating mechanisms based on the hybrid model.
  * `data-generation-mundlak.R`: Data-generating mechanisms based on Mundlak’s contextual model (the model used in the main manuscript).
  * `model-fitting.R`: Model-fitting procedures for both GLMMs and GEEs.
  * `result-formatting.R`: Functions to clean and format simulation output.

### `docs/`

Contains interactive and rendered documents to explore the simulation designs.

* **`data-exploration.qmd`**: Allows users to explore all four data-generating mechanisms (DGMs) considered in the study.
* **`data-exploration.html`**: Rendered HTML version for [direct inspection](https://wardeiling.github.io/multilevel-vs-gee-binary/data-exploration.html).

Supplementary materials accompanying the manuscript.

* **`supplementary_materials.qmd`**: Quarto document with:

  1. A comparison between the hybrid and Mundlak's contextual model.
  2. Discussion of boundary/extreme estimates in GEEs and how they were handled.
* **`supplementary_materials.html`**: Rendered HTML version for [direct inspection](https://wardeiling.github.io/multilevel-vs-gee-binary/supplementary_materials.html).

### `simulation_results/`

Contains raw and processed simulation outputs, organized into subfolders corresponding to different simulation runs:

* **`April10_fullsimulation/`**: Part 1, covering DGMs 2–4.
* **`April17_fullsimulation_contxy/`**: Part 2, covering DGM 1.
* **`April18_fullsimulation_combined/figures/`**: Final figures used in the manuscript.

Each run folder contains:

* `i.RDS`: Raw output for each design/scenario `i`.
* `settings.RDS`: Simulation settings used for that batch.
* `log.txt`: Logs containing warnings and errors during simulation.
* `summary-results-bias.RDS` & `.csv`: Summary files quantifying bias in the estimates.

### **`renv/`**

Contains documents that save the settings of the `renv` environment.

## Reproducibility via `renv`

This repository uses the [`renv`](https://rstudio.github.io/renv/) package to create a reproducible R environment. To replicate the computational setup:

1. Download `R` version 4.2.2 from CRAN ([link](https://cran.rstudio.com/bin/windows/base/old/4.4.2/R-4.4.2-win.exe)), install in RStudio and set as R version.
2. Clone or download the entire repository.
3. Open the project in RStudio.
4. Run:

   ```r
   renv::restore()
   ```

This restores all package versions as specified in the `renv.lock` file, ensuring consistent results across systems and over time.

After the computational setup is replicated, we can run the main simulation and reproduce the output as follows

1. Run **`scripts/main-simulation-function-future-simul-part1.R`**, which should automatically retrieve the helper functions and produce output in the folder `April10_fullsimulation/`
2. Run **`scripts/main-simulation-function-future-simul-part2.R`**, which should automatically retrieve the helper functions and produce output in the folder `April17_fullsimulation_contXY/`

Now we can use the output to reproduce the figures shown in the result section of the manuscript as follows

1. Run **`post-processing/results-plotting.R`**,  which should retrieve the simulation output, merge them together and process it for creating the boxplots.

---
