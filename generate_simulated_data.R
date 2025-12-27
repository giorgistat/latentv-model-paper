################################################################################
# Generate Simulated Malaria Serology Data
#
# This script generates synthetic antibody response data similar to the 
# malaria serology data analyzed in the manuscript. The simulated data
# preserve the key statistical properties while protecting the privacy of
# the original study participants.
#
# Generates two datasets:
#   1. AMA1 (Apical Membrane Antigen 1)
#   2. MSP1 (Merozoite Surface Protein 1)
#
# Each dataset contains:
#   - age: Individual age in years
#   - *_log: Log-transformed antibody concentration
#   - T: True latent immune state (for validation only)
#
# The simulation follows the models described in Sections 5.1 and 5.2
################################################################################

# Clear workspace
rm(list = ls())
set.seed(20231215)  # For reproducibility

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

cat("Generating simulated malaria serology data...\n\n")

################################################################################
# 1. Simulate Population Age Distribution
################################################################################

cat("Step 1: Simulating population age distribution\n")

# Sample size (same as original data-set after filtering age > 1)
n_total <- 17503

# Age distribution: weighted sampling from realistic age structure
# 8 breaks define 7 age intervals
age_breaks <- c(1, 6, 10, 15, 20, 30, 45, 100)
n_intervals <- length(age_breaks) - 1

# Proportions for each age interval (must sum to 1.0)
age_proportions <- c(0.35, 0.25, 0.20, 0.12, 0.05, 0.02, 0.01)

# Verify proportions sum to 1
if (abs(sum(age_proportions) - 1.0) > 1e-10) {
  stop("age_proportions must sum to 1.0")
}

# Convert proportions to counts
age_counts <- round(age_proportions * n_total)

# Adjust for rounding (ensure exact total)
diff_n <- n_total - sum(age_counts)
if (diff_n != 0) {
  # Add/subtract difference to largest group
  age_counts[which.max(age_counts)] <- age_counts[which.max(age_counts)] + diff_n
}

cat(sprintf("  - Target sample size: %d\n", n_total))
cat(sprintf("  - Age counts by interval: %s\n", 
            paste(age_counts, collapse = ", ")))
cat(sprintf("  - Total: %d\n", sum(age_counts)))

# Sample ages uniformly within each bracket
ages <- numeric(n_total)
idx_start <- 1

for (i in seq_len(n_intervals)) {
  n_i <- age_counts[i]
  
  if (n_i > 0) {
    idx_end <- idx_start + n_i - 1
    ages[idx_start:idx_end] <- runif(n_i, age_breaks[i], age_breaks[i + 1])
    idx_start <- idx_end + 1
  }
}

cat(sprintf("  - Generated %d ages (range: %.1f to %.1f years)\n", 
            length(ages), min(ages), max(ages)))
cat(sprintf("  - Mean age: %.1f years (SD: %.1f)\n", mean(ages), sd(ages)))

################################################################################
# 2. Simulate AMA1 Data
################################################################################

cat("\nStep 2: Simulating AMA1 antibody data\n")

# True parameter values (loosely based on Table 3 in manuscript)
# These can be adjusted to match desired data characteristics

# Conditional distribution parameters
mu0_ama <- -3.19
mu1_ama <- 0.75
sigma0_ama <- 0.75
sigma1_ama <- 0.09

# Change point
tau_ama <- 21

# Parameters for age < tau (mechanistic model)
lambda_ama <- 0.15  # Seroconversion rate
alpha2_ama <- 1.50  # Beta shape for seropositive component
beta2_ama <- 1.00

# Parameters for age >= tau (regression model)
phi_ama <- 4.5
eta1_ama <- -0.14

# Compute eta0 from continuity constraint
t_grid <- seq(0.01, 0.99, length.out = 100)

# Component 1: approximated as point mass at 0
ET1 <- 0.001

# Component 2: Beta(alpha2, beta2)
prior_t2 <- dbeta(t_grid, alpha2_ama, beta2_ama)
prior_t2 <- prior_t2 / sum(prior_t2)
ET2 <- sum(t_grid * prior_t2)

# Mean at tau (from below)
pi_tau <- 1 - exp(-lambda_ama * tau_ama)
m_tau <- (1 - pi_tau) * ET1 + pi_tau * ET2

# Continuity constraint
eta0_ama <- qlogis(m_tau) - eta1_ama * log(tau_ama)

cat("  AMA1 simulation parameters:\n")
cat(sprintf("    - Change point: %.1f years\n", tau_ama))
cat(sprintf("    - Seroconversion rate: λ = %.3f\n", lambda_ama))
cat(sprintf("    - Mean at boundary: %.3f\n", m_tau))

# Simulate latent variables T
T_ama <- numeric(n_total)
for (i in seq_len(n_total)) {
  a <- ages[i]
  
  if (a < tau_ama) {
    # Mechanistic model: mixture of point mass and Beta
    pi_a <- 1 - exp(-lambda_ama * a)
    
    if (runif(1) < (1 - pi_a)) {
      T_ama[i] <- 0.001  # Approximation of point mass at 0
    } else {
      T_ama[i] <- rbeta(1, alpha2_ama, beta2_ama)
    }
    
  } else {
    # Regression model: Beta distribution
    eta_a <- eta0_ama + eta1_ama * log(a)
    mu_a <- plogis(eta_a)
    
    alpha_a <- mu_a * phi_ama
    beta_a <- (1 - mu_a) * phi_ama
    
    T_ama[i] <- rbeta(1, alpha_a, beta_a)
  }
}

# Simulate observed antibody levels Y | T
mu_ama <- T_ama * mu1_ama + (1 - T_ama) * mu0_ama
sigma_ama <- sqrt(T_ama * sigma1_ama^2 + (1 - T_ama) * sigma0_ama^2)
Y_ama <- rnorm(n_total, mean = mu_ama, sd = sigma_ama)

# Create data frame
data_ama <- tibble(
  age = ages,
  T = T_ama,
  ama_log = Y_ama
)

cat(sprintf("  - Simulated %d observations\n", nrow(data_ama)))
cat(sprintf("  - Mean log-antibody: %.2f (SD: %.2f)\n", mean(Y_ama), sd(Y_ama)))
cat(sprintf("  - Mean latent T: %.3f (SD: %.3f)\n", mean(T_ama), sd(T_ama)))

# Save
saveRDS(data_ama, "simulated_ama1_data.rds")
cat("  - Saved to: simulated_ama1_data.rds\n")

################################################################################
# 3. Simulate MSP1 Data
################################################################################

cat("\nStep 3: Simulating MSP1 antibody data\n")

# True parameter values (loosely based on Table 5 in manuscript)

# Conditional distribution parameters
mu0_msp <- -4.48
mu1_msp <- 1.26
sigma0_msp <- exp(-0.68)
sigma1_msp <- exp(-5.72)

# Beta distribution parameters (power law)
alpha0_msp <- 0.09
gamma_msp <- 0.28

beta0_msp <- 0.76
delta1_msp <- 0.11
zeta_msp <- 11.6
delta2_msp <- -0.06

cat("  MSP1 simulation parameters:\n")
cat(sprintf("    - Change point: %.1f years\n", zeta_msp))
cat(sprintf("    - Alpha exponent: γ = %.3f\n", gamma_msp))
cat(sprintf("    - Beta exponent shift: δ₂ = %.3f\n", delta2_msp))

# Simulate latent variables T
T_msp <- numeric(n_total)
for (i in seq_len(n_total)) {
  a <- ages[i]
  
  # Age-dependent Beta parameters
  alpha_a <- alpha0_msp * (a^gamma_msp)
  
  delta_a <- delta1_msp + (a > zeta_msp) * delta2_msp
  beta_a <- beta0_msp * (a^delta_a)
  
  T_msp[i] <- rbeta(1, alpha_a, beta_a)
}

# Simulate observed antibody levels Y | T
mu_msp <- T_msp * mu1_msp + (1 - T_msp) * mu0_msp
sigma_msp <- sqrt(T_msp * sigma1_msp^2 + (1 - T_msp) * sigma0_msp^2)
Y_msp <- rnorm(n_total, mean = mu_msp, sd = sigma_msp)

# Create data frame
data_msp <- tibble(
  age = ages,
  T = T_msp,
  msp_log = Y_msp
)

cat(sprintf("  - Simulated %d observations\n", nrow(data_msp)))
cat(sprintf("  - Mean log-antibody: %.2f (SD: %.2f)\n", mean(Y_msp), sd(Y_msp)))
cat(sprintf("  - Mean latent T: %.3f (SD: %.3f)\n", mean(T_msp), sd(T_msp)))

# Save
saveRDS(data_msp, "simulated_msp1_data.rds")
cat("  - Saved to: simulated_msp1_data.rds\n")

################################################################################
# 4. Summary Statistics
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("DATA GENERATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nAMA1 summary:\n")
cat(sprintf("  - N = %d\n", nrow(data_ama)))
cat(sprintf("  - Age: %.1f ± %.1f years (range: %.1f-%.1f)\n", 
            mean(data_ama$age), sd(data_ama$age), min(data_ama$age), max(data_ama$age)))
cat(sprintf("  - log(antibody): %.2f ± %.2f\n", mean(data_ama$ama_log), sd(data_ama$ama_log)))
cat(sprintf("  - Latent T: %.3f ± %.3f\n", mean(data_ama$T), sd(data_ama$T)))

cat("\nMSP1 summary:\n")
cat(sprintf("  - N = %d\n", nrow(data_msp)))
cat(sprintf("  - Age: %.1f ± %.1f years (range: %.1f-%.1f)\n", 
            mean(data_msp$age), sd(data_msp$age), min(data_msp$age), max(data_msp$age)))
cat(sprintf("  - log(antibody): %.2f ± %.2f\n", mean(data_msp$msp_log), sd(data_msp$msp_log)))
cat(sprintf("  - Latent T: %.3f ± %.3f\n", mean(data_msp$T), sd(data_msp$T)))

cat("\nGenerated files:\n")
cat("  - simulated_ama1_data.rds\n")
cat("  - simulated_msp1_data.rds\n\n")

cat("These datasets can now be used with the analysis scripts:\n")
cat("  - 01_fit_ama_by_age_group.R\n")
cat("  - 02_fit_msp_by_age_group.R\n")
cat("  - 03_joint_fit_ama.R\n")
cat("  - 04_joint_fit_msp.R\n\n")
