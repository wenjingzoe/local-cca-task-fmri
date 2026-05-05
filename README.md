# Spatial Local CCA for Task-related fMRI

This repository provides MATLAB code to reproduce the simulation studies and
figures in the paper:

**Local Canonical Correlation Analysis for Voxel-wise Detection of
Task-Related Dependence in fMRI**

## Overview

The code implements:
- Construction of block-design stimuli and Toeplitz design matrices
- Local linear drift removal
- Spatially weighted second-moment estimation
- Raw local CCA estimators and Wald-type test statistics
- Monte Carlo simulations for finite-sample evaluation (including QQ plots)

The repository is organized to allow direct reproduction of the results
reported in the simulation sections of the paper.

## Repository Structure

- `Simulation_singlevoxel/`  
  Scripts for single-voxel Monte Carlo experiments, including verification of
  asymptotic normality and Wald statistics.

- `Simulation_2d/`  
  Scripts for two-dimensional spatial simulations, including activation
  detection and neighborhood-size sensitivity analysis.

- `OpenNeuro_Emotional_ds000171/`  
  Scripts for the emotional task fMRI application based on OpenNeuro dataset
  ds000171.

- `block_stimulus.m`  
  Generate deterministic block-design stimulus sequences.

- `build_local_linear_smoother.m`  
  Construct the local-linear smoother used for drift removal.

- `generate_ar1.m`  
  Generate AR(1) noise sequences.

- `compute_TPR_FPR.m`  
  Compute detection performance metrics.

- `estimate_h_FGLS.m`  
  Benchmark function for semiparametric estimation.

## Data availability

The real fMRI dataset used in this study is publicly available from OpenNeuro:

https://openneuro.org/datasets/ds000171

Due to data size and licensing considerations, the full raw dataset is not
required to be downloaded from this repository. Instead, users are encouraged
to obtain the data directly from OpenNeuro.

For convenience, this repository may include preprocessed data used in the
analysis.

## Requirements

- MATLAB R2023a (or later)

Toolboxes that may be required:
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox
- SPM (for real fMRI analysis only)

