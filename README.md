# From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes

This repository contains all materials associated with the manuscript:
**"From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes"**
It includes code, supplementary documentation, and simulation results to ensure full transparency and reproducibility of the study.

## Ethics Assessment

This simulation study was approved by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. The approval is based on the documents sent by the researchers as requested in the form of the Ethics committee and filed under FETC number 24-2003. The approval is valid through 31 May 2025. The approval of the Ethical Review Board concerns ethical aspects, as well as data management and privacy issues (including the GDPR).

## Study Design


## Repository Structure

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
