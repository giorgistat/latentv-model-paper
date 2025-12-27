################################################################################
# Helper Functions for Latent Variable Modeling of Antibody Response Data
# 
# This script contains shared functions used across the analysis scripts
# for the manuscript: "A flexible class of latent variable models for the 
# analysis of antibody response data"
#
# Functions included:
#   - Conditional normal density (Y|T)
#   - Log-sum-exp utility
#   - Age bracket label formatting
#   - Plotting utilities
#   - Data simulation utilities
################################################################################

# Required packages
required_packages <- c("dplyr", "tibble", "ggplot2", "scales", "purrr", "tidyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# Core Model Functions
################################################################################

#' Conditional Normal Density Y | T
#' 
#' @param t Vector of latent variable values in [0,1]
#' @param y Observed antibody concentration (log-scale)
#' @param mu1 Mean at T=0
#' @param mu2 Mean at T=1
#' @param sigma1 SD at T=0
#' @param sigma2 SD at T=1
#' @param log.p Logical, return log density?
#' @return Density values
dnorm_latent <- function(t, y, mu1, mu2, sigma1, sigma2, log.p = FALSE) {
  mu_t     <- t * mu2 + (1 - t) * mu1
  sigma2_t <- t * sigma2^2 + (1 - t) * sigma1^2
  
  dens <- dnorm(y, mean = mu_t, sd = sqrt(sigma2_t), log = log.p)
  
  if (log.p) {
    dens[sigma2_t <= 0] <- -Inf
  } else {
    dens[sigma2_t <= 0] <- 0
  }
  
  dens
}

#' Log-Sum-Exp for Numerical Stability
#' 
#' @param x Vector of log-values
#' @return log(sum(exp(x))) computed stably
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

################################################################################
# Data Organization Functions
################################################################################

#' Create Age Bracket Labels
#' 
#' @param age_breaks Vector of age break points
#' @return Character vector of bracket labels
make_bracket_labels <- function(age_breaks) {
  n <- length(age_breaks)
  labels <- character(n - 1)
  
  for (i in seq_len(n - 1)) {
    if (i < n - 1) {
      labels[i] <- sprintf("[%g, %g)", age_breaks[i], age_breaks[i + 1])
    } else {
      labels[i] <- sprintf("[%g, %g]", age_breaks[i], age_breaks[i + 1])
    }
  }
  
  labels
}

#' Assign Age Groups
#' 
#' @param age Numeric vector of ages
#' @param age_breaks Vector of age break points
#' @return Factor vector of age group assignments
assign_age_groups <- function(age, age_breaks) {
  labels <- make_bracket_labels(age_breaks)
  cut(age, breaks = age_breaks, labels = labels, 
      right = FALSE, include.lowest = TRUE)
}

################################################################################
# Plotting Utility Functions
################################################################################

#' Clean Y-Axis Theme
#' 
#' @param n_breaks Number of axis breaks
#' @param acc Accuracy for number formatting
#' @return List of ggplot2 theme elements
axis_y_clean <- function(n_breaks = 5, acc = NULL) {
  list(
    scale_y_continuous(
      breaks = pretty_breaks(n_breaks),
      labels = if (is.null(acc)) label_number() else label_number(accuracy = acc),
      expand = expansion(mult = c(0.02, 0.06))
    ),
    theme(
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 11),
      strip.text = element_text(size = 10),
      panel.spacing = unit(6, "pt")
    )
  )
}

################################################################################
# Simulation Functions
################################################################################

#' Simulate from Beta Distribution for Latent Variable T
#' 
#' @param n Number of observations
#' @param alpha Shape parameter 1
#' @param beta Shape parameter 2
#' @return Vector of T values in (0,1)
simulate_T_beta <- function(n, alpha, beta) {
  rbeta(n, alpha, beta)
}

#' Simulate Y | T under Linear Interpolation Model
#' 
#' @param t Vector of latent variable values
#' @param mu1 Mean at T=0
#' @param mu2 Mean at T=1
#' @param sigma1 SD at T=0
#' @param sigma2 SD at T=1
#' @return Vector of simulated Y values
simulate_Y_given_T <- function(t, mu1, mu2, sigma1, sigma2) {
  mu_t <- t * mu2 + (1 - t) * mu1
  sigma2_t <- t * sigma2^2 + (1 - t) * sigma1^2
  
  rnorm(length(t), mean = mu_t, sd = sqrt(sigma2_t))
}

#' Simulate Complete Dataset (T and Y)
#' 
#' @param n Number of observations
#' @param age Age vector
#' @param alpha_fun Function returning alpha(age)
#' @param beta_fun Function returning beta(age)
#' @param mu1 Mean at T=0
#' @param mu2 Mean at T=1
#' @param sigma1 SD at T=0
#' @param sigma2 SD at T=1
#' @return Data frame with age, T, and Y
simulate_serology_data <- function(n, age, alpha_fun, beta_fun, 
                                   mu1, mu2, sigma1, sigma2) {
  alpha_age <- alpha_fun(age)
  beta_age <- beta_fun(age)
  
  T_vals <- rbeta(n, alpha_age, beta_age)
  Y_vals <- simulate_Y_given_T(T_vals, mu1, mu2, sigma1, sigma2)
  
  tibble(
    age = age,
    T = T_vals,
    Y = Y_vals
  )
}

################################################################################
# Model Validation Functions
################################################################################

#' Create Histogram Envelope for Model Validation
#' 
#' @param y_obs Observed Y values
#' @param ages Corresponding ages
#' @param theta Parameter vector
#' @param simulate_fun Function to simulate Y|age
#' @param n_sims Number of simulation replicates
#' @param n_bins Number of histogram bins
#' @param probs Quantile probabilities for envelope
#' @return Tibble with histogram and envelope
histogram_envelope <- function(y_obs, ages, theta, simulate_fun, 
                               n_sims = 1000, n_bins = 20,
                               probs = c(0.025, 0.5, 0.975)) {
  n <- length(y_obs)
  
  # Simulate datasets
  y_sim_mat <- matrix(NA_real_, nrow = n, ncol = n_sims)
  for (s in seq_len(n_sims)) {
    y_sim_mat[, s] <- simulate_fun(ages, theta)
  }
  
  # Compute histogram range
  range_all <- range(c(y_obs, y_sim_mat))
  breaks <- seq(range_all[1], range_all[2], length.out = n_bins + 1)
  
  # Observed histogram
  h_obs <- hist(y_obs, breaks = breaks, plot = FALSE)
  mids <- h_obs$mids
  widths <- diff(breaks)
  dens_obs <- h_obs$counts / (n * widths)
  
  # Simulated histograms
  dens_mat <- matrix(NA_real_, nrow = length(mids), ncol = n_sims)
  for (s in seq_len(n_sims)) {
    h_s <- hist(y_sim_mat[, s], breaks = breaks, plot = FALSE)
    dens_mat[, s] <- h_s$counts / (n * widths)
  }
  
  # Compute quantiles
  dens_q <- apply(dens_mat, 1, quantile, probs = probs)
  
  tibble(
    mid = mids,
    density_obs = dens_obs,
    density_low = dens_q[1, ],
    density_med = dens_q[2, ],
    density_high = dens_q[3, ]
  )
}

################################################################################
# Model Comparison Functions
################################################################################

#' Compute BIC for Model Comparison
#' 
#' @param loglik Log-likelihood value
#' @param n_params Number of parameters
#' @param n_obs Number of observations
#' @return BIC value
compute_BIC <- function(loglik, n_params, n_obs) {
  -2 * loglik + n_params * log(n_obs)
}

#' Compute AIC for Model Comparison
#' 
#' @param loglik Log-likelihood value
#' @param n_params Number of parameters
#' @return AIC value
compute_AIC <- function(loglik, n_params) {
  -2 * loglik + 2 * n_params
}

cat("Helper functions loaded successfully.\n")
cat("Available functions:\n")
cat("  - dnorm_latent: Conditional density Y|T\n")
cat("  - log_sum_exp: Numerical stability utility\n")
cat("  - make_bracket_labels: Format age group labels\n")
cat("  - simulate_serology_data: Generate synthetic data\n")
cat("  - histogram_envelope: Model validation\n")
cat("  - compute_BIC, compute_AIC: Model comparison\n")
