# ============================================================================
# Joint Model Across All Ages: AMA1
# ============================================================================
# Two-part model with change point τ:
#   - Age < τ: Mechanistic mixture (Dirac at 0 + Beta(α₂,1)), π(age) = 1-exp(-λ·age)
#   - Age ≥ τ: Regression model Beta(μ(age), φ), logit(μ) = η₀ + η₁·log(age)
#   - Continuity constraint at τ
#
# Outputs:
#   - Left panel: Histogram envelope (observed vs simulated) by age group
#   - Right panel: π(age) curve with change point τ marked
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

if (!dir.exists("results")) dir.create("results")

cat("\n=== AMA1 Joint Model Across All Ages ===\n\n")

# ============================================================================
# Load Data
# ============================================================================

cat("Loading data...\n")
data_ama <- readRDS("simulated_ama1_data.rds")

y <- data_ama$ama_log
age <- data_ama$age
age_min <- 1

cat(sprintf("  N = %d observations\n", length(y)))
cat(sprintf("  Age range: %.1f to %.1f years\n\n", min(age), max(age)))

# ============================================================================
# Optimization Setup
# ============================================================================

t_grid <- seq(0.001, 0.999, length.out = 200)

# Constants for first Beta component (Dirac at 0)
alpha1 <- 1
beta1_large <- 10000
prior_beta1 <- dbeta(t_grid, alpha1, beta1_large)
prior_beta1 <- prior_beta1 / sum(prior_beta1)

# L2 histogram objective
compute_l2_hist <- function(theta, Y_hist) {
  
  mu0 <- theta[1]
  mu1 <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha2 <- 1 + exp(theta[5])
  lambda <- exp(theta[6])
  phi <- exp(theta[7])
  eta1 <- theta[8]
  tau <- age_min + exp(theta[9])
  
  beta2 <- 1
  p0 <- 1
  
  prior_beta2 <- dbeta(t_grid, alpha2, beta2)
  prior_beta2 <- prior_beta2 / sum(prior_beta2)
  
  # Continuity at tau
  pi_tau <- p0 * exp(-lambda * tau)
  mu_tau_minus <- pi_tau * sum(t_grid * prior_beta1) + (1 - pi_tau) * sum(t_grid * prior_beta2)
  mu_tau_minus <- pmin(pmax(mu_tau_minus, 1e-6), 1 - 1e-6)
  eta0 <- qlogis(mu_tau_minus) - eta1 * log(tau)
  
  dens_model <- numeric(length(Y_hist$midpoints))
  
  for (j in seq_along(Y_hist$midpoints)) {
    y_j <- Y_hist$midpoints[j]
    
    # Find approximate age for this y value (use mean age as approximation)
    age_approx <- mean(age[abs(y - y_j) < 1])
    if (is.na(age_approx)) age_approx <- mean(age)
    
    if (age_approx < tau) {
      pi_i <- p0 * exp(-lambda * age_approx)
      prior_t <- pi_i * prior_beta1 + (1 - pi_i) * prior_beta2
    } else {
      eta_i <- eta0 + eta1 * log(age_approx)
      mu_i <- plogis(eta_i)
      alpha_i <- mu_i * phi
      beta_i <- (1 - mu_i) * phi
      prior_t <- dbeta(t_grid, alpha_i, beta_i)
      prior_t <- prior_t / sum(prior_t)
    }
    
    mu_t <- (1 - t_grid) * mu0 + t_grid * mu1
    sigma_t <- sqrt((1 - t_grid) * sigma0^2 + t_grid * sigma1^2)
    
    lik_t <- dnorm(y_j, mean = mu_t, sd = sigma_t)
    dens_model[j] <- sum(lik_t * prior_t)
  }
  
  # Empirical density (no weights)
  emp_dens <- Y_hist$counts / sum(Y_hist$counts)
  model_dens <- dens_model / sum(dens_model)
  
  # L2 on log scale
  sum((log(emp_dens + 1e-10) - log(model_dens + 1e-10))^2)
}

negloglik <- function(theta) {
  
  mu0 <- theta[1]
  mu1 <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha2 <- 1 + exp(theta[5])
  lambda <- exp(theta[6])
  phi <- exp(theta[7])
  eta1 <- theta[8]
  tau <- age_min + exp(theta[9])
  
  beta2 <- 1
  p0 <- 1
  
  prior_beta2 <- dbeta(t_grid, alpha2, beta2)
  prior_beta2 <- prior_beta2 / sum(prior_beta2)
  
  # Continuity at tau
  pi_tau <- p0 * exp(-lambda * tau)
  mu_tau_minus <- pi_tau * sum(t_grid * prior_beta1) + (1 - pi_tau) * sum(t_grid * prior_beta2)
  mu_tau_minus <- pmin(pmax(mu_tau_minus, 1e-6), 1 - 1e-6)
  eta0 <- qlogis(mu_tau_minus) - eta1 * log(tau)
  
  log_lik <- 0
  
  for (i in seq_along(y)) {
    age_i <- age[i]
    
    if (age_i < tau) {
      pi_i <- p0 * exp(-lambda * age_i)
      prior_t <- pi_i * prior_beta1 + (1 - pi_i) * prior_beta2
    } else {
      eta_i <- eta0 + eta1 * log(age_i)
      mu_i <- plogis(eta_i)
      alpha_i <- mu_i * phi
      beta_i <- (1 - mu_i) * phi
      prior_t <- dbeta(t_grid, alpha_i, beta_i)
      prior_t <- prior_t / sum(prior_t)
    }
    
    mu_t <- (1 - t_grid) * mu0 + t_grid * mu1
    sigma_t <- sqrt((1 - t_grid) * sigma0^2 + t_grid * sigma1^2)
    
    lik_t <- dnorm(y[i], mean = mu_t, sd = sigma_t)
    marg_lik <- sum(lik_t * prior_t)
    
    log_lik <- log_lik + log(max(marg_lik, 1e-300))
  }
  
  -log_lik
}

# ============================================================================
# Optimization
# ============================================================================

cat("Creating histogram for L2 optimization...\n")
Y_hist_all <- hist(y, breaks = 50, plot = FALSE)
Y_hist_all$midpoints <- (Y_hist_all$breaks[-length(Y_hist_all$breaks)] + Y_hist_all$breaks[-1]) / 2

cat("Stage 1: L2 histogram-based optimization...\n")
init <- c(-3.2, 0.75, log(0.75), log(0.09), log(0.5), log(0.15), log(4.5), -0.14, log(20))

opt_l2 <- nlminb(
  start = init,
  objective = function(th) compute_l2_hist(th, Y_hist_all),
  control = list(eval.max = 500, iter.max = 200)
)

cat(sprintf("  Converged: %d\n", opt_l2$convergence == 0))

cat("\nStage 2: MLE refinement...\n")
opt_mle <- nlminb(
  start = opt_l2$par,
  objective = negloglik,
  control = list(eval.max = 1000, iter.max = 500)
)

cat(sprintf("  Converged: %d\n", opt_mle$convergence == 0))

theta_hat <- opt_mle$par

# Extract estimates
tau_hat <- age_min + exp(theta_hat[9])
lambda_hat <- exp(theta_hat[6])
phi_hat <- exp(theta_hat[7])

cat(sprintf("\n  Change point τ = %.2f years\n", tau_hat))
cat(sprintf("  λ = %.3f\n", lambda_hat))
cat(sprintf("  φ = %.2f\n\n", phi_hat))

# ============================================================================
# Simulation Function
# ============================================================================

simulate_from_model <- function(age_vec, theta) {
  
  mu0 <- theta[1]
  mu1 <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha2 <- 1 + exp(theta[5])
  lambda <- exp(theta[6])
  phi <- exp(theta[7])
  eta1 <- theta[8]
  tau <- age_min + exp(theta[9])
  
  beta2 <- 1
  p0 <- 1
  
  prior_beta2 <- dbeta(t_grid, alpha2, beta2)
  prior_beta2 <- prior_beta2 / sum(prior_beta2)
  
  pi_tau <- p0 * exp(-lambda * tau)
  mu_tau_minus <- pi_tau * sum(t_grid * prior_beta1) + (1 - pi_tau) * sum(t_grid * prior_beta2)
  mu_tau_minus <- pmin(pmax(mu_tau_minus, 1e-6), 1 - 1e-6)
  eta0 <- qlogis(mu_tau_minus) - eta1 * log(tau)
  
  y_sim <- numeric(length(age_vec))
  
  for (i in seq_along(age_vec)) {
    age_i <- age_vec[i]
    
    if (age_i < tau) {
      pi_i <- p0 * exp(-lambda * age_i)
      prior_t <- pi_i * prior_beta1 + (1 - pi_i) * prior_beta2
    } else {
      eta_i <- eta0 + eta1 * log(age_i)
      mu_i <- plogis(eta_i)
      alpha_i <- mu_i * phi
      beta_i <- (1 - mu_i) * phi
      prior_t <- dbeta(t_grid, alpha_i, beta_i)
      prior_t <- prior_t / sum(prior_t)
    }
    
    t_val <- sample(t_grid, size = 1, prob = prior_t)
    mu_val <- (1 - t_val) * mu0 + t_val * mu1
    sigma_val <- sqrt((1 - t_val) * sigma0^2 + t_val * sigma1^2)
    
    y_sim[i] <- rnorm(1, mean = mu_val, sd = sigma_val)
  }
  
  y_sim
}

# ============================================================================
# Create Validation Envelope
# ============================================================================

cat("Creating histogram envelope validation plot...\n")

# Define age groups
age_breaks <- c(1, seq(6, floor(tau_hat), by = 4), ceiling(tau_hat), 20, 30, 45, 100)
age_breaks <- sort(unique(age_breaks))
age_labels <- paste0("[", head(age_breaks, -1), ", ", tail(age_breaks, -1), ")")

df_full <- data.frame(age = age, y = y)
df_full$age_group <- cut(df_full$age, breaks = age_breaks, labels = age_labels,
                         right = FALSE, include.lowest = TRUE)
df_full <- df_full[!is.na(df_full$age_group), ]

# Simulate envelopes
n_sim <- 100
cat(sprintf("  Simulating %d datasets...\n", n_sim))

# Process each age group
envelope_list <- list()

for (ag in unique(df_full$age_group)) {
  df_group <- df_full[df_full$age_group == ag, ]
  y_obs <- df_group$y
  ages <- df_group$age
  n <- length(y_obs)
  
  cat(sprintf("    Processing %s (n=%d)...\n", ag, n))
  
  # Simulate
  y_sim_mat <- matrix(NA, nrow = n, ncol = n_sim)
  for (s in 1:n_sim) {
    y_sim_mat[, s] <- simulate_from_model(ages, theta_hat)
  }
  
  # Common breaks with buffer
  range_obs <- range(y_obs, na.rm = TRUE)
  range_sim <- range(y_sim_mat, na.rm = TRUE)
  range_all <- c(min(range_obs[1], range_sim[1]), max(range_obs[2], range_sim[2]))
  range_buffer <- diff(range_all) * 0.1
  breaks <- seq(range_all[1] - range_buffer, range_all[2] + range_buffer, length.out = 21)
  
  # Observed
  h_obs <- hist(y_obs, breaks = breaks, plot = FALSE)
  mids <- h_obs$mids
  widths <- diff(breaks)
  dens_obs <- h_obs$counts / (n * widths)
  
  # Simulated
  dens_mat <- matrix(NA, nrow = length(mids), ncol = n_sim)
  for (s in 1:n_sim) {
    h_s <- hist(y_sim_mat[, s], breaks = breaks, plot = FALSE)
    dens_mat[, s] <- h_s$counts / (n * widths)
  }
  
  dens_q <- apply(dens_mat, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  envelope_list[[as.character(ag)]] <- data.frame(
    age_group = ag,
    mid = mids,
    density_obs = dens_obs,
    density_low = dens_q[1, ],
    density_med = dens_q[2, ],
    density_high = dens_q[3, ]
  )
}

envelope_data <- do.call(rbind, envelope_list)
rownames(envelope_data) <- NULL

p_envelope <- ggplot(envelope_data, aes(x = mid)) +
  geom_col(aes(y = density_obs), fill = "grey80", colour = "grey50", width = 0.15) +
  geom_ribbon(aes(ymin = density_low, ymax = density_high), alpha = 0.25, fill = "steelblue") +
  geom_line(aes(y = density_med), linewidth = 0.6) +
  facet_wrap(~ age_group, scales = "free_y") +
  labs(x = "log(AMA1)", y = "Density") +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    panel.spacing = unit(6, "pt")
  )

# ============================================================================
# Create π(age) Plot
# ============================================================================

cat("Creating π(age) plot...\n")

age_seq <- seq(1, tau_hat, length.out = 100)
pi_age <- 1 - exp(-lambda_hat * age_seq)

df_pi <- tibble(age = age_seq, pi = pi_age)

p_pi <- ggplot(df_pi, aes(x = age, y = pi)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = tau_hat, linetype = "dashed", colour = "red", linewidth = 0.7) +
  labs(x = "Age (years)", y = expression(pi(age))) +
  ylim(0, 1) +
  theme_bw(base_size = 10)

# ============================================================================
# Combine and Save
# ============================================================================

fig_combined <- p_envelope + p_pi + plot_layout(widths = c(2, 1))

print(fig_combined)

ggsave("results/ama1_validation_envelope.pdf", fig_combined, width = 14, height = 7)
ggsave("results/ama1_validation_envelope.png", fig_combined, width = 14, height = 7, dpi = 300)

# Save estimates
param_table <- tibble(
  parameter = c("mu0", "mu1", "log(sigma0)", "log(sigma1)", "log(alpha2-1)", 
                "log(lambda)", "log(phi)", "eta1", "log(tau-age_min)"),
  estimate = theta_hat
)

write.csv(param_table, "results/ama1_joint_model_estimates.csv", row.names = FALSE)

cat("\nFigures saved to: results/ama1_validation_envelope.{pdf,png}\n")
cat("Estimates saved to: results/ama1_joint_model_estimates.csv\n")
cat("\n=== Analysis Complete ===\n\n")