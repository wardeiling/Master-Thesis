# From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes

This repository contains all materials associated with the manuscript:
**"From Multilevel Modeling to GEE: Revisiting the Within- and Between-Person Debate with Binary Predictors and Outcomes"**
It includes code, supplementary documentation, and simulation results to ensure full transparency and reproducibility of the study.

## Repository Structure

### `scripts/`

Contains all core scripts for running and analyzing the simulation study.

* **`main-simulation-function-future-simul.R`**
  Modularized main script that executes the simulation across various design conditions.

* **`results-plotting.R`**
  Scripts used to produce the plots featured in the Results section of the manuscript.

* **`helper-functions/`** (subfolder with modular components):

  * `data-generation-centeredX.R`: Data-generating mechanisms based on the hybrid model.
  * `data-generation-mundlak.R`: Data-generating mechanisms based on Mundlak’s contextual model (the model used in the main manuscript).
  * `model-fitting.R`: Model-fitting procedures for both GLMMs and GEEs.
  * `result-formatting.R`: Functions to clean and format simulation output.

### `data_exploration/`

Contains interactive and rendered documents to explore the simulation designs.

* **`data-exploration.qmd`**: Allows users to explore all four data-generating mechanisms (DGMs) considered in the study.
* **`data-exploration.html`**: Rendered HTML version for direct inspection.

### `supplementary_materials/`

Supplementary materials accompanying the manuscript.

* **`supplementary_materials.qmd`**: Quarto document with:

  1. A comparison between the hybrid and Mundlak's contextual model.
  2. Discussion of boundary/extreme estimates in GEEs and how they were handled.
* **`supplementary_materials.html`**: Rendered HTML version for direct inspection.

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

---
