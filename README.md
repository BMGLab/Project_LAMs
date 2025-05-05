# --- Sample README.md for Bioinformatics Project ---

"""
# Project Title: [Descriptive title of your analysis/project]

## Overview
This project investigates [brief description of the biological context and goals].

## Project Structure
```
your_project/
├── data/            # raw and processed data
├── scripts/         # analysis scripts (numbered)
├── figures/         # generated plots and visuals
├── results/         # intermediate result files
├── notebooks/       # exploratory Jupyter/RMarkdown notebooks
├── docs/            # documentation and manuscript text
├── README.md        # this file
├── environment.yml  # reproducible environment
```

## Setup Instructions
```bash
conda env create -f environment.yml
conda activate EMAdown
```

## Pipeline Summary
1. **01_Preprocess_LuCa_and_AUCell_macs.ipynb** – Load and QC the input data. Then, run AUCell for mac subtypes.
2. **02_macs_predictions.R** – Create HieFIT cell type prediction model. 
3. **03_inspect_ct_composition.R** – Inspect cell type compositions across samples.
4. **04_PB-and-DESeq2.ipynb** – Perform DEG on pseudo-bulked profiles.
5. **05_PB-DEG_inspect_results.R** - Inspect DEG expressions across macrophages.
6. **06_create_signature_matrix.R** - Create cell type signature matrix for deconvolution with CiberSortX.
7. **07_prepare_bulk-data.R** - Fetch TCGA lung cancer bulk datasets and prepare.
8. **08_deconvolute_bulk-data.R** - Estimate macrophage composition of bulk RNAseq data.
9. **09_run_survival_analysis.R** - Correlate macrophage compositons with prognosis using Kaplan-Meier survival analysis.

to run sccoda scripts: conda activate sccoda-gpu


## Figures
All key plots are saved in `figures/`. The naming follows the manuscript figure order (e.g., `Figure1_PCA.png`, `Figure2_Volcano.pdf`).

## Reproducibility
- All scripts use fixed random seeds.
- Software versions are logged via `sessionInfo()` or `pip freeze`.

## Citation
If you use this project, please cite: [Your paper or bioRxiv link]
"""