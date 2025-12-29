################################################################################
## Simulation Study: L2 vs MLE - Data Generation and Estimation Only
##
## This script:
## 1. Simulates data under different scenarios
## 2. Estimates parameters using both L2 and MLE
## 3. Records computation times
## 4. Saves all results to .rdata file
##
## No analysis, plotting, or summary statistics - just raw simulation output
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})


################################################################################
## INTEGRATION GRID
################################################################################

# t_grid is a fine grid over [0,1] used for numerical integration ONLY
# - NOT used for data simulation (we use rbeta() directly)
# - Primary method: adaptive integrate() function (more accurate)
# - Fallback method: trapezoidal rule on t_grid (when integrate() fails)
# 
# Grid of 200 points provides good balance between accuracy and speed
t_grid <- seq(0.001, 0.999, length.out = 200)

################################################################################
## SIMULATION PARAMETERS
################################################################################

# Sample sizes to investigate
sample_sizes <- c(100, 500, 1000, 5000)

# Number of simulation replicates per scenario
n_sims <- 1000

################################################################################
## SIMULATION SCENARIOS
################################################################################

scenarios <- list(
  # Scenario 1: Low transmission
  list(
    name = "low_transmission",
    mu0 = -3.0,
    mu1 = 1.0,
    sigma0 = 0.8,
    sigma1 = 0.3,
    alpha = 0.5,
    beta = 3.0,
    description = "Low transmission: U-shaped Beta concentrated near T=0"
  ),
  
  # Scenario 2: High transmission
  list(
    name = "high_transmission",
    mu0 = -3.0,
    mu1 = 1.0,
    sigma0 = 0.8,
    sigma1 = 0.3,
    alpha = 3.0,
    beta = 0.5,
    description = "High transmission: U-shaped Beta concentrated near T=1"
  ),
  
  # Scenario 3: Intermediate
  list(
    name = "intermediate",
    mu0 = -3.0,
    mu1 = 1.0,
    sigma0 = 0.8,
    sigma1 = 0.3,
    alpha = 2.0,
    beta = 2.0,
    description = "Intermediate transmission: symmetric unimodal Beta"
  ),
  
  # Scenario 4: Strongly bimodal
  list(
    name = "bimodal",
    mu0 = -3.5,
    mu1 = 1.5,
    sigma0 = 0.7,
    sigma1 = 0.2,
    alpha = 0.5,
    beta = 0.5,
    description = "Strongly bimodal: U-shaped with equal weights at extremes"
  )
)

################################################################################
## HELPER FUNCTIONS
################################################################################

# Compute model density at a point y for given Beta parameters
# Uses adaptive numerical integration (integrate() function)
# t_grid is used as fallback if integrate() fails
compute_lbm_density <- function(y_val, theta, t_grid) {
  mu0    <- theta[1]
  mu1    <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha  <- exp(theta[5])
  beta   <- exp(theta[6])
  
  # Integrand: f(y|t) * g(t)
  integrand <- function(t) {
    mu_t     <- (1 - t) * mu0 + t * mu1
    sigma2_t <- (1 - t) * sigma0^2 + t * sigma1^2
    if (any(sigma2_t <= 0)) return(rep(0, length(t)))
    dnorm(y_val, mean = mu_t, sd = sqrt(sigma2_t)) * dbeta(t, alpha, beta)
  }
  
  # Try adaptive integration first (more accurate)
  result <- tryCatch({
    integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-6)$value
  }, error = function(e) {
    # Fallback: use trapezoidal rule on t_grid
    K  <- length(t_grid)
    dt <- t_grid[2] - t_grid[1]
    vals <- integrand(t_grid)
    # Trapezoidal rule
    dt * (0.5 * vals[1] + sum(vals[2:(K-1)]) + 0.5 * vals[K])
  })
  
  result
}

# L2 criterion: sum of squared differences between empirical and model densities
# Uses histogram approximation for computational efficiency
L2_criterion_lbm <- function(theta, mids, dens_obs, t_grid) {
  
  # Compute model density at each bin midpoint
  dens_model <- sapply(mids, function(m) {
    compute_lbm_density(m, theta, t_grid)
  })
  
  # L2 contribution: sum of squared differences
  log(sum((dens_obs - dens_model)^2))
}

# Log-likelihood for LBM
# Uses adaptive integration with t_grid as fallback for numerical stability
loglik_lbm <- function(theta, y, t_grid) {
  
  mu0    <- theta[1]
  mu1    <- theta[2]
  sigma0 <- exp(theta[3])
  sigma1 <- exp(theta[4])
  alpha  <- exp(theta[5])
  beta   <- exp(theta[6])
  
  # Check for valid parameters
  if (alpha <= 0 || beta <= 0 || sigma0 <= 0 || sigma1 <= 0) return(-1e10)
  
  loglik <- 0
  n      <- length(y)
  
  for (i in seq_len(n)) {
    y_i <- y[i]
    
    # Integrand: f(y|t) * g(t)
    integrand <- function(t) {
      mu_t     <- (1 - t) * mu0 + t * mu1
      sigma2_t <- (1 - t) * sigma0^2 + t * sigma1^2
      if (any(sigma2_t <= 0)) return(rep(0, length(t)))
      dnorm(y_i, mean = mu_t, sd = sqrt(sigma2_t)) * dbeta(t, alpha, beta)
    }
    
    # Primary: Use adaptive integration (more accurate but can fail)
    int_result <- tryCatch({
      integrate(integrand, lower = 0, upper = 1, 
                rel.tol = 1e-8, abs.tol = 1e-10,
                subdivisions = 500)
    }, error = function(e) {
      # Fallback: Grid-based trapezoidal integration on t_grid
      K  <- length(t_grid)
      dt <- t_grid[2] - t_grid[1]
      vals <- integrand(t_grid)
      # Trapezoidal rule: (h/2) * (f(0) + 2*sum(f(interior)) + f(1))
      list(value = dt * (0.5 * vals[1] + sum(vals[2:(K-1)]) + 0.5 * vals[K]))
    })
    
    dens_y_i <- int_result$value
    
    if (is.finite(dens_y_i) && dens_y_i > 0) {
      loglik <- loglik + log(dens_y_i)
    } else {
      return(-1e10)
    }
  }
  
  loglik
}

################################################################################
## DATA SIMULATION FUNCTION
################################################################################

simulate_lbm_data <- function(n, mu0, mu1, sigma0, sigma1, alpha, beta) {
  
  # Sample latent states directly from Beta distribution (no discretization)
  T_vals <- rbeta(n, alpha, beta)
  
  # For each individual, compute conditional mean and variance
  mu_T     <- (1 - T_vals) * mu0 + T_vals * mu1
  sigma2_T <- (1 - T_vals) * sigma0^2 + T_vals * sigma1^2
  
  # Sample observed outcomes
  Y_vals <- rnorm(n, mean = mu_T, sd = sqrt(sigma2_T))
  
  return(Y_vals)
}

################################################################################
## ESTIMATION FUNCTIONS
################################################################################

# L2 estimation
estimate_L2 <- function(y, theta_init, t_grid, n_bins = NULL) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Use Sturges rule for number of bins if not provided
  if (is.null(n_bins)) {
    n_bins <- ceiling(log2(length(y)) + 1)
  }
  
  # Create histogram
  y_range <- range(y)
  breaks  <- seq(y_range[1], y_range[2], length.out = n_bins + 1)
  widths  <- diff(breaks)
  mids    <- (breaks[-1] + breaks[-length(breaks)]) / 2
  
  h <- hist(y, breaks = breaks, plot = FALSE)
  dens_obs <- h$counts / (length(y) * widths)
  
  # Optimize L2 criterion
  fit_L2 <- optim(
    par     = theta_init,
    fn      = L2_criterion_lbm,
    mids    = mids,
    dens_obs = dens_obs,
    t_grid  = t_grid,
    method  = "BFGS",
    control = list(maxit = 300, trace = 0)
  )
  
  # End timing
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    par = fit_L2$par,
    value = fit_L2$value,
    convergence = fit_L2$convergence,
    n_bins = n_bins,
    time_seconds = elapsed_time
  )
}

# MLE estimation (starting from L2 estimates)
estimate_MLE <- function(y, theta_init, t_grid) {
  
  # Start timing
  start_time <- Sys.time()
  
  neg_loglik <- function(theta) {
    -loglik_lbm(theta, y = y, t_grid = t_grid)
  }
  
  fit_MLE <- optim(
    par     = theta_init,
    fn      = neg_loglik,
    method  = "BFGS",
    control = list(maxit = 500, reltol = 1e-8, trace = 0)
  )
  
  # End timing
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    par = fit_MLE$par,
    value = -fit_MLE$value,  # Return log-likelihood (not negative)
    convergence = fit_MLE$convergence,
    time_seconds = elapsed_time
  )
}

################################################################################
## SINGLE SIMULATION REPLICATE
################################################################################

run_single_sim <- function(scenario, n, sim_id, t_grid) {
  
  # Extract true parameters
  mu0_true    <- scenario$mu0
  mu1_true    <- scenario$mu1
  sigma0_true <- scenario$sigma0
  sigma1_true <- scenario$sigma1
  alpha_true  <- scenario$alpha
  beta_true   <- scenario$beta
  
  # True parameter vector (on log scale for sigma, alpha, beta)
  theta_true <- c(mu0_true, mu1_true, 
                  log(sigma0_true), log(sigma1_true),
                  log(alpha_true), log(beta_true))
  
  # Simulate data
  y <- simulate_lbm_data(n, mu0_true, mu1_true, sigma0_true, sigma1_true,
                         alpha_true, beta_true)
  
  # Initial values for estimation (add some noise to true values)
  theta_init <- theta_true + rnorm(6, 0, 0.1)
  
  # L2 estimation
  fit_L2 <- tryCatch({
    estimate_L2(y, theta_init, t_grid)
  }, error = function(e) {
    list(par = rep(NA, 6), value = NA, convergence = 1, n_bins = NA, time_seconds = NA)
  })
  
  # MLE estimation (starting from L2)
  fit_MLE <- tryCatch({
    estimate_MLE(y, fit_L2$par, t_grid)
  }, error = function(e) {
    list(par = rep(NA, 6), value = NA, convergence = 1, time_seconds = NA)
  })
  
  # Transform to original scale
  theta_L2 <- c(fit_L2$par[1:2], exp(fit_L2$par[3:6]))
  theta_MLE <- c(fit_MLE$par[1:2], exp(fit_MLE$par[3:6]))
  theta_true_original <- c(mu0_true, mu1_true, sigma0_true, sigma1_true, 
                           alpha_true, beta_true)
  
  # Return results
  tibble(
    scenario_name = scenario$name,
    n = n,
    sim_id = sim_id,
    
    # L2 estimates (original scale)
    mu0_L2 = theta_L2[1],
    mu1_L2 = theta_L2[2],
    sigma0_L2 = theta_L2[3],
    sigma1_L2 = theta_L2[4],
    alpha_L2 = theta_L2[5],
    beta_L2 = theta_L2[6],
    
    # MLE estimates (original scale)
    mu0_MLE = theta_MLE[1],
    mu1_MLE = theta_MLE[2],
    sigma0_MLE = theta_MLE[3],
    sigma1_MLE = theta_MLE[4],
    alpha_MLE = theta_MLE[5],
    beta_MLE = theta_MLE[6],
    
    # True values
    mu0_true = mu0_true,
    mu1_true = mu1_true,
    sigma0_true = sigma0_true,
    sigma1_true = sigma1_true,
    alpha_true = alpha_true,
    beta_true = beta_true,
    
    # Convergence flags
    converged_L2 = fit_L2$convergence == 0,
    converged_MLE = fit_MLE$convergence == 0,
    
    # Timing (seconds)
    time_L2 = fit_L2$time_seconds,
    time_MLE = fit_MLE$time_seconds,
    
    # Additional info
    n_bins_L2 = fit_L2$n_bins,
    L2_criterion = fit_L2$value,
    loglik_MLE = fit_MLE$value
  )
}

################################################################################
## MAIN SIMULATION LOOP
################################################################################

cat("================================================================================\n")
cat("SIMULATION STUDY: L2 vs MLE\n")
cat("================================================================================\n\n")

cat("Configuration:\n")
cat("  Number of scenarios:", length(scenarios), "\n")
cat("  Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("  Replicates per scenario-sample size:", n_sims, "\n")
cat("  Total simulations:", length(scenarios) * length(sample_sizes) * n_sims, "\n\n")

# Initialize results storage
all_results <- list()
counter <- 1

# Loop over scenarios
for (s in seq_along(scenarios)) {
  scenario <- scenarios[[s]]
  
  cat("--------------------------------------------------------------------------------\n")
  cat("Scenario", s, "of", length(scenarios), ":", scenario$name, "\n")
  cat(scenario$description, "\n")
  cat("Parameters: mu0=", scenario$mu0, ", mu1=", scenario$mu1, 
      ", sigma0=", scenario$sigma0, ", sigma1=", scenario$sigma1,
      ", alpha=", scenario$alpha, ", beta=", scenario$beta, "\n")
  
  # Loop over sample sizes
  for (n in sample_sizes) {
    
    cat("  Sample size n =", n, "... ")
    start_time <- Sys.time()
    
    # Run simulations
    results_n <- lapply(1:n_sims, function(sim_id) {
      run_single_sim(scenario, n, sim_id, t_grid)
    })
    
    # Combine results for this n
    results_n <- bind_rows(results_n)
    
    elapsed <- difftime(Sys.time(), start_time, units = "secs")
    cat("Done in", round(elapsed, 1), "seconds\n")
    
    all_results[[counter]] <- results_n
    counter <- counter + 1
  }
}

# Combine all results
results_df <- bind_rows(all_results)

################################################################################
## SAVE RESULTS
################################################################################

# Generate random number for filename
random_id <- sample(10000:99999, 1)
filename <- sprintf("simulation_results_%d.rdata", random_id)

# Save results and metadata
save(results_df, scenarios, sample_sizes, n_sims, file = filename)

cat("\n================================================================================\n")
cat("SIMULATION COMPLETE\n")
cat("================================================================================\n")
cat("Results saved to:", filename, "\n")
cat("Total rows in results:", nrow(results_df), "\n\n")

# Quick convergence check
cat("--- Convergence Summary ---\n")
convergence_summary <- results_df %>%
  group_by(scenario_name, n) %>%
  summarise(
    prop_L2_converged = mean(converged_L2, na.rm = TRUE),
    prop_MLE_converged = mean(converged_MLE, na.rm = TRUE),
    .groups = "drop"
  )
print(convergence_summary)

cat("\n")
cat("================================================================================\n")
cat("Simulation complete. Run analysis script to compute summaries.\n")
cat("================================================================================\n")