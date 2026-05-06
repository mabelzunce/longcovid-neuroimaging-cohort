# PET intersubject network analysis

This folder contains a reproducible analysis script for group-level PET metabolic covariance networks.

## What it does

- Reads `data/processed_pet/PET_CERMEP_CEUNIM_harmonized.csv`
- Uses groups:
  - **Long COVID**: `Grupo == "LC"`
  - **Controls**: `Grupo == "Control"` and `SITE == "CERMEP"`
- Computes ROI-by-ROI connectivity matrices (correlations across subjects)
- Compares whole-network and frontal-involved network metrics
- Runs permutation tests for group differences
- Exports matrices, stats tables, and paper-ready plots

## Script

- `pet_intersubject_network_analysis.py`

Outputs are written to:
- `data/results/pet_network_analysis/`
- `data/plots/pet_network_analysis/`
