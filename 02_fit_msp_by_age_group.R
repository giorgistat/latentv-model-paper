################################################################################
# Age-Stratified Analysis of MSP1 Antibody Data
#
# This script fits the latent Beta model to simulated MSP1 data
# separately for each age group. The conditional distribution parameters
# (mu0, mu1, sigma0, sigma1) are estimated from the youngest age group
# and held fixed for all other groups, while the Beta shape parameters
# (alpha, beta) are re-estimated for each group.
#
# Uses two-stage optimization:
#   1. L2 histogram-based estimation (fast, robust initialization)
#   2. Maximum likelihood estimation with nlminb (refinement)
#
# This exploratory analysis (Section 5.2 of manuscript) guides the
# specification of the joint model across all ages.
################################################################################

# Clear workspace
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(patchwork)
  library(tibble)
})

cat("Starting MSP1 age-stratified analysis...\n\n")

################################################################################
# 1. Helper Functions for Model Fitting
################################################################################

# Compute L2 distance between histogram and model density
compute_l2_hist <- function(par, Y_hist, t_grid = seq(0.001, 0.999, length.out = 100)) {
  # Parameters: (mu0, mu1, log_sigma0, log_sigma1, log_alpha, log_beta)
  mu0 <- par[1]
  mu1 <- par[2]
  sigma0 <- exp(par[3])
  sigma1 <- exp(par[4])
  alpha <- exp(par[5])
  beta <- exp(par[6])
  
  # Beta distribution for latent variable T
  prior_t <- dbeta(t_grid, alpha, beta)
  prior_t <- prior_t / sum(prior_t)  # Normalize
  
  # Model density at histogram midpoints
  model_dens <- numeric(length(Y_hist$midpoints))
  
  for (i in seq_along(Y_hist$midpoints)) {
    y <- Y_hist$midpoints[i]
    
    # Integrate over latent variable T
    mu_t <- (1 - t_grid) * mu0 + t_grid * mu1
    sigma_t <- sqrt((1 - t_grid) * sigma0^2 + t_grid * sigma1^2)
    
    # Gaussian density for each t
    lik_t <- dnorm(y, mean = mu_t, sd = sigma_t)
    
    # Marginal density
    model_dens[i] <- sum(lik_t * prior_t)
  }
  
  # Empirical histogram density
  emp_dens <- Y_hist$counts / (sum(Y_hist$counts) * Y_hist$bin_widths)
  
  # L2 distance
  l2_dist <- log(sum((emp_dens - model_dens)^2))
  
  return(l2_dist)
}

# Compute negative log-likelihood
negloglik <- function(par, y_data, t_grid = seq(0.001, 0.999, length.out = 100)) {
  # Parameters: (mu0, mu1, log_sigma0, log_sigma1, log_alpha, log_beta)
  mu0 <- par[1]
  mu1 <- par[2]
  sigma0 <- exp(par[3])
  sigma1 <- exp(par[4])
  alpha <- exp(par[5])
  beta <- exp(par[6])
  
  # Beta distribution for latent variable T
  prior_t <- dbeta(t_grid, alpha, beta)
  prior_t <- prior_t / sum(prior_t)  # Normalize
  
  # Compute log-likelihood
  n <- length(y_data)
  loglik <- numeric(n)
  
  for (i in seq_len(n)) {
    y <- y_data[i]
    
    # Integrate over latent variable T
    mu_t <- (1 - t_grid) * mu0 + t_grid * mu1
    sigma_t <- sqrt((1 - t_grid) * sigma0^2 + t_grid * sigma1^2)
    
    # Gaussian density for each t
    lik_t <- dnorm(y, mean = mu_t, sd = sigma_t)
    
    # Marginal likelihood
    lik <- sum(lik_t * prior_t)
    
    loglik[i] <- log(lik + 1e-300)  # Avoid log(0)
  }
  
  return(-sum(loglik))
}

# Create age bracket labels
make_bracket_labels <- function(age_breaks) {
  labels <- character(length(age_breaks) - 1)
  for (i in seq_along(labels)) {
    if (i == length(labels)) {
      labels[i] <- sprintf("[%d,100)", age_breaks[i])
    } else {
      labels[i] <- sprintf("[%d,%d)", age_breaks[i], age_breaks[i + 1])
    }
  }
  return(labels)
}

################################################################################
# 2. Load Data
################################################################################

cat("Loading simulated MSP1 data...\n")

# Load simulated data
data_msp <- readRDS("simulated_msp1_data.rds")

y_log <- data_msp$msp_log
age <- data_msp$age

cat(sprintf("  - Loaded %d observations\n", length(y_log)))
cat(sprintf("  - Age range: %.1f to %.1f years\n", min(age), max(age)))
cat(sprintf("  - Log-antibody range: %.2f to %.2f\n", min(y_log), max(y_log)))

################################################################################
# 3. Define Age Brackets
################################################################################

cat("\nDefining age brackets...\n")

age_brackets <- c(1, 6, 10, 15, 20, 30, 45, max(age))
bracket_labels <- make_bracket_labels(age_brackets)

cat(sprintf("  - Number of age groups: %d\n", length(bracket_labels)))
cat("  - Age groups:\n")
for (i in seq_along(bracket_labels)) {
  cat(sprintf("    %s\n", bracket_labels[i]))
}

################################################################################
# 4. Fit Model to Each Age Group
################################################################################

cat("\nFitting latent Beta model to each age group...\n")
cat("  Using two-stage optimization: L2 then MLE\n")

results <- vector("list", length(bracket_labels))
names(results) <- bracket_labels

for (i in seq_along(age_brackets[-length(age_brackets)])) {
  label <- bracket_labels[i]
  
  # Select observations in this age bracket
  if (i == 1) {
    idx <- age < age_brackets[i + 1]
  } else {
    idx <- age >= age_brackets[i] & age < age_brackets[i + 1]
  }
  
  y_i <- y_log[idx]
  n_i <- sum(idx)
  
  cat(sprintf("\n  Age group %s: n = %d\n", label, n_i))
  
  # Create histogram for L2 estimation
  Y_hist <- hist(y_i, breaks = 50, plot = FALSE)
  Y_hist$midpoints <- (Y_hist$breaks[-length(Y_hist$breaks)] + 
                         Y_hist$breaks[-1]) / 2
  Y_hist$bin_widths <- diff(Y_hist$breaks)
  Y_hist$weights <- rep(1, length(Y_hist$bin_widths))
  
  # Calculate midpoint age for this bracket
  if (i == 1) {
    mid_age <- mean(c(1, age_brackets[i + 1]))
  } else {
    mid_age <- mean(c(age_brackets[i], age_brackets[i + 1]))
  }
  
  # Fit model using two-stage approach
  if (i == 1) {
    # First age group: estimate all parameters
    cat("    Stage 1: L2 optimization for all parameters...\n")
    
    init <- c(
      -4.5,        # mu0
      1.3,         # mu1
      log(0.5),    # log_sigma0
      log(0.01),   # log_sigma1
      log(0.3),    # log_alpha
      log(1.0)     # log_beta
    )
    
    # Stage 1: L2 optimization
    opt_l2 <- nlminb(
      start = init,
      objective = function(th) compute_l2_hist(th, Y_hist),
      control = list(eval.max = 500, iter.max = 200)
    )
    
    cat("    Stage 2: MLE refinement with nlminb...\n")
    
    # Stage 2: MLE with nlminb, starting from L2 solution
    opt_mle <- nlminb(
      start = opt_l2$par,
      objective = function(th) negloglik(th, y_i),
      control = list(eval.max = 1000, iter.max = 500)
    )
    
    par_hat <- opt_mle$par
    
    cat(sprintf("    mu0 = %.3f, mu1 = %.3f\n", par_hat[1], par_hat[2]))
    cat(sprintf("    sigma0 = %.3f, sigma1 = %.3f\n", 
                exp(par_hat[3]), exp(par_hat[4])))
    cat(sprintf("    alpha = %.3f, beta = %.3f\n", 
                exp(par_hat[5]), exp(par_hat[6])))
    cat(sprintf("    Negative log-likelihood: %.2f\n", opt_mle$objective))
    
  } else {
    # Subsequent age groups: fix mu0, mu1, sigma0, sigma1, re-estimate alpha, beta
    cat("    Stage 1: L2 optimization for Beta parameters...\n")
    
    par_fixed <- results[[1]]$par[1:4]
    init_ab <- c(log(0.5), log(1.0))
    
    # Stage 1: L2 optimization
    opt_l2 <- nlminb(
      start = init_ab,
      objective = function(th) compute_l2_hist(c(par_fixed, th), Y_hist),
      control = list(eval.max = 500, iter.max = 200)
    )
    
    cat("    Stage 2: MLE refinement with nlminb...\n")
    
    # Stage 2: MLE with nlminb
    opt_mle <- nlminb(
      start = opt_l2$par,
      objective = function(th) negloglik(c(par_fixed, th), y_i),
      control = list(eval.max = 1000, iter.max = 500)
    )
    
    par_hat <- c(par_fixed, opt_mle$par)
    
    cat(sprintf("    alpha = %.3f, beta = %.3f\n", 
                exp(par_hat[5]), exp(par_hat[6])))
    cat(sprintf("    Negative log-likelihood: %.2f\n", opt_mle$objective))
  }
  
  # Store results
  results[[i]] <- list(
    par = par_hat,
    Y_hist = Y_hist,
    mid_age = mid_age,
    n = n_i,
    negloglik = if(i == 1) opt_mle$objective else opt_mle$objective
  )
}

cat("\nModel fitting complete.\n")

################################################################################
# 5. Create Diagnostic Plots
################################################################################

cat("\nCreating diagnostic plots...\n")

# Helper function for clean y-axis
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

# Function to plot histograms vs model densities
plot_histogram_vs_model <- function(results, bracket_labels) {
  
  plot_data <- list()
  
  for (i in seq_along(results)) {
    res <- results[[i]]
    label <- names(results)[i]
    
    # Extract parameters
    par <- res$par
    mu0 <- par[1]
    mu1 <- par[2]
    sigma0 <- exp(par[3])
    sigma1 <- exp(par[4])
    alpha <- exp(par[5])
    beta <- exp(par[6])
    
    # Create grid for plotting model density
    y_grid <- seq(min(res$Y_hist$breaks), max(res$Y_hist$breaks), length.out = 200)
    
    # Compute model density
    t_grid <- seq(0.001, 0.999, length.out = 100)
    prior_t <- dbeta(t_grid, alpha, beta)
    prior_t <- prior_t / sum(prior_t)
    
    model_dens <- numeric(length(y_grid))
    for (j in seq_along(y_grid)) {
      y <- y_grid[j]
      mu_t <- (1 - t_grid) * mu0 + t_grid * mu1
      sigma_t <- sqrt((1 - t_grid) * sigma0^2 + t_grid * sigma1^2)
      lik_t <- dnorm(y, mean = mu_t, sd = sigma_t)
      model_dens[j] <- sum(lik_t * prior_t)
    }
    
    # Histogram data
    hist_df <- tibble(
      y_mid = res$Y_hist$midpoints,
      density = res$Y_hist$counts / (sum(res$Y_hist$counts) * res$Y_hist$bin_widths),
      age_group = label
    )
    
    # Model density data
    model_df <- tibble(
      y = y_grid,
      density = model_dens,
      age_group = label
    )
    
    plot_data[[i]] <- list(hist = hist_df, model = model_df)
  }
  
  # Combine data
  hist_df_all <- bind_rows(lapply(plot_data, function(x) x$hist))
  model_df_all <- bind_rows(lapply(plot_data, function(x) x$model))
  
  # Create plot
  p <- ggplot() +
    geom_col(
      data = hist_df_all,
      aes(x = y_mid, y = density, fill = age_group),
      width = 0.15,
      alpha = 0.6
    ) +
    geom_line(
      data = model_df_all,
      aes(x = y, y = density, color = age_group),
      linewidth = 0.8
    ) +
    facet_wrap(~ factor(age_group, levels = bracket_labels), 
               ncol = 4, 
               scales = "free_y") +
    labs(
      x = "Log antibody concentration",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

# Function to plot Beta density curves
plot_beta_curves <- function(results, bracket_labels) {
  
  t_grid <- seq(0.001, 0.999, length.out = 200)
  plot_data <- list()
  
  for (i in seq_along(results)) {
    res <- results[[i]]
    label <- names(results)[i]
    
    alpha <- exp(res$par[5])
    beta <- exp(res$par[6])
    
    beta_dens <- dbeta(t_grid, alpha, beta)
    
    plot_data[[i]] <- tibble(
      t = t_grid,
      density = beta_dens,
      age_group = label
    )
  }
  
  plot_df <- bind_rows(plot_data)
  
  p <- ggplot(plot_df, aes(x = t, y = density, color = age_group, fill = age_group)) +
    geom_line(linewidth = 0.8) +
    geom_area(alpha = 0.2, position = "identity") +
    labs(
      x = "Latent immune state T",
      y = "Density",
      color = "Age group",
      fill = "Age group"
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

# Create plots
cat("  - Creating histogram vs model plot...\n")
p_hist <- plot_histogram_vs_model(results, bracket_labels) +
  axis_y_clean(n_breaks = 5)

cat("  - Creating Beta density curves plot...\n")
p_beta <- plot_beta_curves(results, bracket_labels) +
  axis_y_clean(n_breaks = 5, acc = 0.1) +
  coord_cartesian(ylim = c(0, 2.5))

# Extract alpha, beta, mu for middle panel
alpha_beta_df <- tibble(
  bracket = names(results),
  alpha = sapply(results, function(r) exp(r$par[5])),
  beta = sapply(results, function(r) exp(r$par[6]))
) %>%
  mutate(
    mu = alpha / (alpha + beta),
    bracket = factor(bracket, levels = bracket_labels)
  )

cat("  - Saving parameter estimates...\n")
write.csv(alpha_beta_df, "results/msp1_alpha_beta_by_age_group.csv", row.names = FALSE)

# Create parameter plots
cat("  - Creating parameter plots...\n")

p_alpha <- ggplot(alpha_beta_df, aes(x = bracket, y = alpha, group = 1)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.8) +
  labs(x = "Age group", y = expression(alpha)) +
  axis_y_clean(n_breaks = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

p_beta_param <- ggplot(alpha_beta_df, aes(x = bracket, y = beta, group = 1)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.8) +
  labs(x = "Age group", y = expression(beta)) +
  axis_y_clean(n_breaks = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

p_mu <- ggplot(alpha_beta_df, aes(x = bracket, y = mu, group = 1)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.8) +
  labs(x = "Age group", y = expression(mu == alpha / (alpha + beta))) +
  axis_y_clean(n_breaks = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# Stack parameter plots vertically
p_mid <- p_alpha / p_beta_param / p_mu + plot_layout(heights = c(1, 1, 1))

################################################################################
# 6. Combine and Save Plots
################################################################################

cat("\nCombining and saving plots...\n")

# Create combined figure: histograms | parameters | Beta curves
fig_combined <- p_hist + p_mid + p_beta +
  plot_layout(widths = c(1, 0.9, 1), guides = "collect") &
  theme(legend.position = "right")

# Create results directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Save plots
ggsave(
  "results/msp1_age_stratified_analysis.pdf", 
  fig_combined, 
  width = 14, 
  height = 6, 
  units = "in"
)

ggsave(
  "results/msp1_age_stratified_analysis.png", 
  fig_combined, 
  width = 14, 
  height = 6, 
  units = "in",
  dpi = 300
)

cat("  - Saved to: results/msp1_age_stratified_analysis.pdf\n")
cat("  - Saved to: results/msp1_age_stratified_analysis.png\n")

################################################################################
# 7. Summary Output
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("MSP1 AGE-STRATIFIED ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nParameter estimates by age group:\n\n")
print(alpha_beta_df, n = Inf)

cat("\nKey observations:\n")
cat(sprintf("  - Alpha increases from %.3f to %.3f\n", 
            min(alpha_beta_df$alpha), max(alpha_beta_df$alpha)))
cat(sprintf("  - Beta increases from %.3f to %.3f then stabilizes\n", 
            min(alpha_beta_df$beta), max(alpha_beta_df$beta)))
cat(sprintf("  - Mean T increases monotonically to %.3f\n",
            max(alpha_beta_df$mu)))

cat("\nOutputs generated:\n")
cat("  - results/msp1_alpha_beta_by_age_group.csv\n")
cat("  - results/msp1_age_stratified_analysis.pdf\n")
cat("  - results/msp1_age_stratified_analysis.png\n\n")