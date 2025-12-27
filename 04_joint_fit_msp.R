# ============================================================================
# Joint Model Across All Ages: MSP1
# ============================================================================
# Age-dependent Beta distribution with piecewise power law:
#   α(age) = α₀ · age^γ
#   β(age) = β₀ · age^δ(age)
#   where δ(age) = δ₁ (age ≤ ζ) or δ₁ + δ₂ (age > ζ)
#
# Outputs:
#   - Histogram envelope (observed vs simulated) by age group
#   - Age-dependent parameter plots: α(age), β(age), E[T]
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

if (!dir.exists("results")) dir.create("results")

cat("\n=== MSP1 Joint Model Across All Ages ===\n\n")

# ============================================================================
# Load Data
# ============================================================================

cat("Loading data...\n")
data_msp <- readRDS("simulated_msp1_data.rds")

y <- data_msp$msp_log
age <- data_msp$age

cat(sprintf("  N = %d observations\n", length(y)))
cat(sprintf("  Age range: %.1f to %.1f years\n\n", min(age), max(age)))

# ============================================================================
# Optimization Setup
# ============================================================================

t_grid <- seq(0.001, 0.999, length.out = 200)

# L2 histogram objective
compute_l2_hist <- function(theta, Y_hist) {
  
  mu0 <- theta[1]
  mu1 <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha0 <- exp(theta[5])
  gamma <- theta[6]
  beta0 <- exp(theta[7])
  delta1 <- theta[8]
  zeta <- theta[9]
  delta2 <- theta[10]
  
  dens_model <- numeric(length(Y_hist$midpoints))
  
  for (j in seq_along(Y_hist$midpoints)) {
    y_j <- Y_hist$midpoints[j]
    
    # Use mean age as approximation for histogram-based optimization
    age_approx <- mean(age[abs(y - y_j) < 1])
    if (is.na(age_approx)) age_approx <- mean(age)
    
    alpha_i <- alpha0 * age_approx^gamma
    
    if (age_approx <= zeta) {
      delta_i <- delta1
    } else {
      delta_i <- delta1 + delta2
    }
    
    beta_i <- beta0 * age_approx^delta_i
    
    prior_t <- dbeta(t_grid, alpha_i, beta_i)
    prior_t <- prior_t / sum(prior_t)
    
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
  alpha0 <- exp(theta[5])
  gamma <- theta[6]
  beta0 <- exp(theta[7])
  delta1 <- theta[8]
  zeta <- theta[9]
  delta2 <- theta[10]
  
  log_lik <- 0
  
  for (i in seq_along(y)) {
    age_i <- age[i]
    
    alpha_i <- alpha0 * age_i^gamma
    
    if (age_i <= zeta) {
      delta_i <- delta1
    } else {
      delta_i <- delta1 + delta2
    }
    
    beta_i <- beta0 * age_i^delta_i
    
    prior_t <- dbeta(t_grid, alpha_i, beta_i)
    prior_t <- prior_t / sum(prior_t)
    
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
init <- c(-4.5, 1.25, log(0.5), log(0.003), log(0.09), 0.28, log(0.75), 0.11, 11.6, -0.06)

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

# Extract key estimates
zeta_hat <- theta_hat[9]
alpha0_hat <- exp(theta_hat[5])
gamma_hat <- theta_hat[6]
beta0_hat <- exp(theta_hat[7])
delta1_hat <- theta_hat[8]
delta2_hat <- theta_hat[10]

cat(sprintf("\n  Change point ζ = %.2f years\n", zeta_hat))
cat(sprintf("  α₀ = %.3f, γ = %.3f\n", alpha0_hat, gamma_hat))
cat(sprintf("  β₀ = %.3f, δ₁ = %.3f, δ₂ = %.3f\n\n", beta0_hat, delta1_hat, delta2_hat))

# ============================================================================
# Simulation Function
# ============================================================================

simulate_from_model <- function(age_vec, theta) {
  
  mu0 <- theta[1]
  mu1 <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha0 <- exp(theta[5])
  gamma <- theta[6]
  beta0 <- exp(theta[7])
  delta1 <- theta[8]
  zeta <- theta[9]
  delta2 <- theta[10]
  
  y_sim <- numeric(length(age_vec))
  
  for (i in seq_along(age_vec)) {
    age_i <- age_vec[i]
    
    alpha_i <- alpha0 * age_i^gamma
    
    if (age_i <= zeta) {
      delta_i <- delta1
    } else {
      delta_i <- delta1 + delta2
    }
    
    beta_i <- beta0 * age_i^delta_i
    
    prior_t <- dbeta(t_grid, alpha_i, beta_i)
    prior_t <- prior_t / sum(prior_t)
    
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

age_breaks <- c(1, 6, 10, 15, 20, 30, 45, 100)
age_labels <- paste0("[", head(age_breaks, -1), ", ", tail(age_breaks, -1), ")")

df_full <- data.frame(age = age, y_log = y)
df_full$age_group <- cut(df_full$age, breaks = age_breaks, labels = age_labels,
                         right = FALSE, include.lowest = TRUE)
df_full <- df_full[!is.na(df_full$age_group), ]

n_sim <- 1000
cat(sprintf("  Simulating %d datasets...\n", n_sim))

# Process each age group
envelope_list <- list()

for (ag in unique(df_full$age_group)) {
  df_group <- df_full[df_full$age_group == ag, ]
  y_obs <- df_group$y_log
  ages <- df_group$age
  n <- length(y_obs)
  
  cat(sprintf("    Processing %s (n=%d)...\n", ag, n))
  
  y_sim_mat <- matrix(NA, nrow = n, ncol = n_sim)
  for (s in 1:n_sim) {
    y_sim_mat[, s] <- simulate_from_model(ages, theta_hat)
  }
  
  # Common breaks spanning observed and simulated data
  range_all <- range(c(y_obs, y_sim_mat))
  breaks <- seq(range_all[1], range_all[2], length.out = 21)
  
  h_obs <- hist(y_obs, breaks = breaks, plot = FALSE)
  mids <- h_obs$mids
  widths <- diff(breaks)
  dens_obs <- h_obs$counts / (n * widths)
  
  dens_mat <- matrix(NA, nrow = length(mids), ncol = n_sim)
  for (s in 1:n_sim) {
    h_s <- hist(y_sim_mat[, s], breaks = breaks, plot = FALSE)
    dens_mat[, s] <- h_s$counts / (n * widths)
  }
  
  dens_q <- apply(dens_mat, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
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
  labs(x = "log(MSP1)", y = "Density") +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    panel.spacing = unit(6, "pt")
  )

print(p_envelope)
ggsave("results/msp1_validation_envelope.pdf", p_envelope, width = 10, height = 7)
ggsave("results/msp1_validation_envelope.png", p_envelope, width = 10, height = 7, dpi = 300)

# ============================================================================
# Create Age-Dependent Parameter Plots
# ============================================================================

cat("Creating age-dependent parameter plots...\n")

age_seq <- seq(min(age), max(age), length.out = 200)

params_df <- tibble(
  age = age_seq,
  alpha = alpha0_hat * age_seq^gamma_hat,
  beta = beta0_hat * age_seq^(delta1_hat + ifelse(age_seq > zeta_hat, delta2_hat, 0)),
  mu_T = alpha / (alpha + beta)
)

p_alpha <- ggplot(params_df, aes(x = age, y = alpha)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  labs(x = "Age (years)", y = expression(alpha(age))) +
  theme_bw(base_size = 10)

p_beta <- ggplot(params_df, aes(x = age, y = beta)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = zeta_hat, linetype = "dashed", colour = "red", linewidth = 0.7) +
  labs(x = "Age (years)", y = expression(beta(age))) +
  theme_bw(base_size = 10)

p_mu_T <- ggplot(params_df, aes(x = age, y = mu_T)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = zeta_hat, linetype = "dashed", colour = "red", linewidth = 0.7) +
  labs(x = "Age (years)", y = expression(E[T])) +
  theme_bw(base_size = 10)

p_params <- p_alpha / p_beta / p_mu_T + plot_layout(heights = c(1, 1, 1))

print(p_params)
ggsave("results/msp1_age_dependent_params.pdf", p_params, width = 8, height = 8)
ggsave("results/msp1_age_dependent_params.png", p_params, width = 8, height = 8, dpi = 300)

# Save age-dependent parameters
write.csv(params_df, "results/msp1_age_dependent_parameters.csv", row.names = FALSE)

# Save estimates
param_table <- tibble(
  parameter = c("mu0", "mu1", "log(sigma0)", "log(sigma1)", "log(alpha0)", 
                "gamma", "log(beta0)", "delta1", "zeta", "delta2"),
  estimate = theta_hat
)

write.csv(param_table, "results/msp1_joint_model_estimates.csv", row.names = FALSE)

cat("\nFigures saved to:\n")
cat("  - results/msp1_validation_envelope.{pdf,png}\n")
cat("  - results/msp1_age_dependent_params.{pdf,png}\n")
cat("Parameter curves saved to: results/msp1_age_dependent_parameters.csv\n")
cat("Estimates saved to: results/msp1_joint_model_estimates.csv\n")
cat("\n=== Analysis Complete ===\n\n")