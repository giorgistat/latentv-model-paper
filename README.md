# Latent Variable Models for Antibody Response Data

This repository contains reproducible R code for the manuscript:

**"A flexible class of latent variable models for the analysis of antibody response data"**  
by Emanuele Giorgi and Jonas Wallin


## Repository Structure

```
├── 00_helper_functions.R       # Shared utility functions
├── 01_fit_ama_by_age_group.R   # Age-stratified AMA1 analysis (Section 6.1)
├── 02_fit_msp_by_age_group.R   # Age-stratified MSP1 analysis (Section 6.2)
├── 03_joint_fit_ama.R          # Joint AMA1 model across all ages
├── 04_joint_fit_msp.R          # Joint MSP1 model across all ages
├── 05_simulation_study1.R      # Simulation study on the performance of the L2 estimator (Section 5.1)
├── 06_simulation_study2.R      # Simulation study on the impact of model misspecification (Section 5.2)
├── generate_simulated_data.R   # Generate synthetic datasets
└── README.md                   # This file
```

## Installation

### Required R Packages

```r
install.packages(c(
  "dplyr",
  "tibble",
  "ggplot2",
  "scales",
  "purrr",
  "tidyr",
  "patchwork"
))
```

## Quick Start

### 1. Generate Simulated Data

First, create synthetic datasets that mimic the statistical properties of the original malaria serology data:

```r
source("generate_simulated_data.R")
```

This creates:
- `simulated_ama1_data.rds` (AMA1 antibody data)
- `simulated_msp1_data.rds` (MSP1 antibody data)

### 2. Run Analyses

#### Age-Stratified Analysis

Fit models separately to each age group (Figures 9 and 11 in manuscript):

```r
source("01_fit_ama_by_age_group.R")  # AMA1
source("02_fit_msp_by_age_group.R")  # MSP1
```

#### Joint Models

Fit models jointly across all ages (Tables 3 and 5, Figures 10 and 12):

```r
source("03_joint_fit_ama.R")         # AMA1 with change-point
source("04_joint_fit_msp.R")         # MSP1 with power law
```

### 3. Create Results Directory

Before running analyses, create a directory for outputs:

```r
dir.create("results", showWarnings = FALSE)
```

## Model Specifications

### AMA1 Model (Section 6.1)

**Age < τ (mechanistic)**:
- T ~ mixture of Dirac(0) and Beta(α₂, 1)
- Mixing weight: π(age) = 1 - exp(-λ · age)

**Age ≥ τ (regression)**:
- T ~ Beta(α(age), β(age))
- logit(μ(age)) = η₀ + η₁ · log(age)
- α(age) = μ(age) · φ, β(age) = (1 - μ(age)) · φ

**Conditional model (all ages)**:
- Y | T ~ N(μ(T), σ²(T))
- μ(T) = (1-T)μ₀ + T·μ₁
- σ²(T) = (1-T)σ₀² + T·σ₁²

### MSP1 Model (Section 6.2)

**Latent distribution**:
- T | age ~ Beta(α(age), β(age))
- α(age) = exp(α₀ + γ · log(age))
- β(age) = exp(β₀ + δ(age) · log(age))
- δ(age) = δ₁ (if age ≤ ζ), δ₁ + δ₂ (if age > ζ)

**Conditional model**: Same as AMA1

## Output Files

### Parameter Estimates

- `results/ama1_alpha_beta_by_age_group.csv`
- `results/msp1_alpha_beta_by_age_group.csv`
- `results/ama1_joint_model_parameters.csv`
- `results/msp1_joint_model_parameters.csv`

### Model Fits

- `results/joint_model_ama_fit.rds`
- `results/joint_model_msp_fit.rds`

### Predictions

- `results/msp1_expected_antibody_by_age.csv`

## Key Functions

### From `00_helper_functions.R`

| Function | Description |
|----------|-------------|
| `dnorm_latent()` | Conditional normal density Y\|T |
| `log_sum_exp()` | Numerically stable log-sum-exp |
| `make_bracket_labels()` | Format age group labels |
| `simulate_serology_data()` | Generate synthetic datasets |
| `histogram_envelope()` | Model validation via simulation |
| `compute_BIC()`, `compute_AIC()` | Model comparison |

## Important Notes

### Data Privacy

The simulated datasets preserve statistical properties of the original data while protecting participant privacy. Parameters are loosely based on published estimates.


### Bootstrap Inference

The original analysis included parametric bootstrap for uncertainty quantification (1000 replicates). This is **not included** in the repository code to reduce complexity, but can be implemented by:

1. Fitting the model to get parameter estimates θ̂
2. Simulating B datasets from the fitted model
3. Re-fitting to each bootstrap sample
4. Computing quantiles of bootstrap parameter distributions


## Contact

For questions or issues:

- Emanuele Giorgi: e.giorgi@bham.ac.uk
- Jonas Wallin: jonas.wallin@stat.lu.se


**Repository**: https://github.com/giorgistat/latentv-model-paper

**Last updated**: 27 December 2024
