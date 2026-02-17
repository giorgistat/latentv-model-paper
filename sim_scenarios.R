################################################################################
# Vaccination Contamination Simulation Study - Single Iteration
# ADAPTED L2 AND LOGLIK FUNCTIONS WITH T-GRID INTEGRATION
################################################################################

rm(list=ls())

################################################################################
# TRUE PARAMETERS
################################################################################

TRUE_PARAMS <- list(
  mu_0 = -3.5,
  mu_1 = 1.0,
  sigma_0 = 0.8,
  sigma_1 = 0.3,
  eta_0 = -2.0,
  eta_1 = 0.4,
  phi = 3.0
)

SAMPLE_SIZES <- c(1000, 2000, 5000)

# Global t-grid for integration
t_grid <- seq(0.001, 0.999, length.out = 50)

################################################################################
# HELPER FUNCTIONS
################################################################################

mu_age <- function(age, eta_0, eta_1) {
  plogis(eta_0 + eta_1 * log(age))
}

alpha_age <- function(age, eta_0, eta_1, phi) {
  mu_age(age, eta_0, eta_1) * phi
}

beta_age <- function(age, eta_0, eta_1, phi) {
  (1 - mu_age(age, eta_0, eta_1)) * phi
}

################################################################################
# DATA GENERATION FUNCTIONS (unchanged)
################################################################################

generate_S1 <- function(n, rho, params = TRUE_PARAMS) {
  ages <- runif(n, 1, 50)
  vaccinated <- runif(n) < rho
  T_vals <- numeric(n)
  
  for(i in 1:n) {
    alpha <- alpha_age(ages[i], params$eta_0, params$eta_1, params$phi)
    beta <- beta_age(ages[i], params$eta_0, params$eta_1, params$phi)
    T_nat <- rbeta(1, alpha, beta)
    
    if(vaccinated[i]) {
      delta <- rbeta(1, 1, 9)
      T_vals[i] <- min(1, T_nat + delta)
    } else {
      T_vals[i] <- T_nat
    }
  }
  
  mu_cond <- (1 - T_vals) * params$mu_0 + T_vals * params$mu_1
  sigma_cond <- sqrt((1 - T_vals) * params$sigma_0^2 + T_vals * params$sigma_1^2)
  Y <- rnorm(n, mu_cond, sigma_cond)
  
  data.frame(age = ages, T = T_vals, Y = Y)
}

generate_S2 <- function(n, delta, params = TRUE_PARAMS) {
  ages <- runif(n, 1, 50)
  in_target_age <- (ages >= 5) & (ages <= 15)
  T_vals <- numeric(n)
  
  for(i in 1:n) {
    alpha <- alpha_age(ages[i], params$eta_0, params$eta_1, params$phi)
    beta <- beta_age(ages[i], params$eta_0, params$eta_1, params$phi)
    
    if(in_target_age[i]) {
      mu_i <- mu_age(ages[i], params$eta_0, params$eta_1)
      mu_reduced <- (mu_i * exp(-delta)) / (mu_i * exp(-delta) + 1 - mu_i)
      alpha_red <- mu_reduced * params$phi
      beta_red <- (1 - mu_reduced) * params$phi
      T_vals[i] <- rbeta(1, alpha_red, beta_red)
    } else {
      T_vals[i] <- rbeta(1, alpha, beta)
    }
  }
  
  mu_cond <- (1 - T_vals) * params$mu_0 + T_vals * params$mu_1
  sigma_cond <- sqrt((1 - T_vals) * params$sigma_0^2 + T_vals * params$sigma_1^2)
  Y <- rnorm(n, mu_cond, sigma_cond)
  
  data.frame(age = ages, T = T_vals, Y = Y)
}

generate_S3 <- function(n, pi, params = TRUE_PARAMS) {
  ages <- runif(n, 1, 50)
  naive <- runif(n) < pi
  T_vals <- numeric(n)
  
  for(i in 1:n) {
    if(naive[i]) {
      T_vals[i] <- 0
    } else {
      alpha <- alpha_age(ages[i], params$eta_0, params$eta_1, params$phi)
      beta <- beta_age(ages[i], params$eta_0, params$eta_1, params$phi)
      T_vals[i] <- rbeta(1, alpha, beta)
    }
  }
  
  mu_cond <- (1 - T_vals) * params$mu_0 + T_vals * params$mu_1
  sigma_cond <- sqrt((1 - T_vals) * params$sigma_0^2 + T_vals * params$sigma_1^2)
  Y <- rnorm(n, mu_cond, sigma_cond)
  
  data.frame(age = ages, T = T_vals, Y = Y)
}

################################################################################
# L2 HISTOGRAM OBJECTIVE - ADAPTED
################################################################################

fit_L2 <- function(data, n_bins = 30) {
  
  # Create histogram
  breaks <- seq(min(data$Y) - 0.1, max(data$Y) + 0.1, length.out = n_bins + 1)
  hist_data <- hist(data$Y, breaks = breaks, plot = FALSE)
  
  Y_hist <- list(
    midpoints = hist_data$mids,
    counts = hist_data$counts,
    breaks = breaks
  )
  
  objective <- function(theta) {
    mu_0 <- theta[1]
    mu_1 <- mu_0 + exp(theta[2])
    sigma_0 <- exp(theta[3])
    sigma_1 <- exp(theta[4])
    eta_0 <- theta[5]
    eta_1 <- theta[6]
    phi <- exp(theta[7])
    
    dens_model <- numeric(length(Y_hist$midpoints))
    
    for (j in seq_along(Y_hist$midpoints)) {
      y_j <- Y_hist$midpoints[j]
      
      # Use mean age as approximation for this y-value
      age_approx <- mean(data$age)
      
      # Compute prior for T at this age
      mu_t_age <- mu_age(age_approx, eta_0, eta_1)
      alpha_t <- mu_t_age * phi
      beta_t <- (1 - mu_t_age) * phi
      
      prior_t <- dbeta(t_grid, alpha_t, beta_t)
      prior_t <- prior_t / sum(prior_t)
      
      # Conditional distribution Y|T
      mu_y_t <- (1 - t_grid) * mu_0 + t_grid * mu_1
      sigma_y_t <- sqrt((1 - t_grid) * sigma_0^2 + t_grid * sigma_1^2)
      
      lik_t <- dnorm(y_j, mean = mu_y_t, sd = sigma_y_t)
      dens_model[j] <- sum(lik_t * prior_t)
    }
    
    # Empirical density
    emp_dens <- Y_hist$counts / sum(Y_hist$counts)
    model_dens <- dens_model / sum(dens_model)
    
    # L2 on log scale (more stable)
    sum((log(emp_dens + 1e-10) - log(model_dens + 1e-10))^2)
  }
  
  init <- c(TRUE_PARAMS$mu_0, 
            log(TRUE_PARAMS$mu_1 - TRUE_PARAMS$mu_0), 
            log(TRUE_PARAMS$sigma_0), 
            log(TRUE_PARAMS$sigma_1),
            TRUE_PARAMS$eta_0, 
            TRUE_PARAMS$eta_1, 
            log(TRUE_PARAMS$phi))
  
  start_time <- Sys.time()
  res <- tryCatch({
    nlminb(init, objective, control = list(eval.max = 5000, iter.max = 5000))
  }, error = function(e) {
    list(par = rep(NA, 7), convergence = 1)
  })
  time_L2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  list(
    mu_0 = res$par[1],
    mu_1 = res$par[1] + exp(res$par[2]),
    sigma_0 = exp(res$par[3]),
    sigma_1 = exp(res$par[4]),
    eta_0 = res$par[5],
    eta_1 = res$par[6],
    phi = exp(res$par[7]),
    converged = (res$convergence == 0),
    time = time_L2
  )
}

################################################################################
# NEGATIVE LOG-LIKELIHOOD - ADAPTED
################################################################################

fit_MLE <- function(data, init) {
  
  negloglik <- function(theta) {
    mu_0 <- theta[1]
    mu_1 <- mu_0 + exp(theta[2])
    sigma_0 <- exp(theta[3])
    sigma_1 <- exp(theta[4])
    eta_0 <- theta[5]
    eta_1 <- theta[6]
    phi <- exp(theta[7])
    
    log_lik <- 0
    
    for (i in 1:nrow(data)) {
      age_i <- data$age[i]
      y_i <- data$Y[i]
      
      # Compute prior for T at this age
      mu_t_age <- mu_age(age_i, eta_0, eta_1)
      alpha_t <- mu_t_age * phi
      beta_t <- (1 - mu_t_age) * phi
      
      prior_t <- dbeta(t_grid, alpha_t, beta_t)
      prior_t <- prior_t / sum(prior_t)
      
      # Conditional distribution Y|T
      mu_y_t <- (1 - t_grid) * mu_0 + t_grid * mu_1
      sigma_y_t <- sqrt((1 - t_grid) * sigma_0^2 + t_grid * sigma_1^2)
      
      lik_t <- dnorm(y_i, mean = mu_y_t, sd = sigma_y_t)
      marg_lik <- sum(lik_t * prior_t)
      
      log_lik <- log_lik + log(max(marg_lik, 1e-300))
    }
    
    -log_lik
  }
  
  # Ensure valid starting values
  if(init$mu_1 <= init$mu_0) {
    init$mu_1 <- init$mu_0 + 0.5
  }
  if(init$sigma_0 <= 0) init$sigma_0 <- 0.1
  if(init$sigma_1 <= 0) init$sigma_1 <- 0.1
  if(init$phi <= 0) init$phi <- 0.5
  
  start_params <- c(init$mu_0, 
                    log(init$mu_1 - init$mu_0), 
                    log(init$sigma_0), 
                    log(init$sigma_1),
                    init$eta_0, 
                    init$eta_1, 
                    log(init$phi))
  
  start_time <- Sys.time()
  res <- tryCatch({
    nlminb(start_params, negloglik, control = list(eval.max = 2000, iter.max = 2000))
  }, error = function(e) {
    list(par = rep(NA, 7), convergence = 1)
  })
  time_MLE <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  list(
    mu_0 = res$par[1],
    mu_1 = res$par[1] + exp(res$par[2]),
    sigma_0 = exp(res$par[3]),
    sigma_1 = exp(res$par[4]),
    eta_0 = res$par[5],
    eta_1 = res$par[6],
    phi = exp(res$par[7]),
    converged = (res$convergence == 0),
    time = time_MLE
  )
}

################################################################################
# RUN ALL SCENARIOS - MATRIX STORAGE
################################################################################

all_results <- list()

for(n in SAMPLE_SIZES) {
  
  cat(sprintf("\n\n===== SAMPLE SIZE n = %d =====\n", n))
  
  # Scenario 1
  for(rho in c(0.2, 0.5, 0.8)) {
    scenario_name <- sprintf("S1_rho%.1f_n%d", rho, n)
    cat(sprintf("\n%s: Generating data...", scenario_name))
    
    data <- generate_S1(n, rho)
    
    cat(" Fitting L2...")
    fit_l2 <- fit_L2(data)
    
    cat(" Fitting MLE...")
    fit_mle <- fit_MLE(data, fit_l2)
    
    all_results[[scenario_name]] <- c(
      mu_0 = fit_mle$mu_0,
      mu_1 = fit_mle$mu_1,
      sigma_0 = fit_mle$sigma_0,
      sigma_1 = fit_mle$sigma_1,
      eta_0 = fit_mle$eta_0,
      eta_1 = fit_mle$eta_1,
      phi = fit_mle$phi,
      n = n,
      rho = rho,
      delta = NA,
      pi = NA,
      converged = as.numeric(fit_mle$converged),
      time = fit_mle$time
    )
    
    cat(sprintf(" Done.\n"))
  }
  
  # Scenario 2
  for(delta in c(0.4, 0.7, 1.2)) {
    scenario_name <- sprintf("S2_delta%.1f_n%d", delta, n)
    cat(sprintf("\n%s: Generating data...", scenario_name))
    
    data <- generate_S2(n, delta)
    
    cat(" Fitting L2...")
    fit_l2 <- fit_L2(data)
    
    cat(" Fitting MLE...")
    fit_mle <- fit_MLE(data, fit_l2)
    
    all_results[[scenario_name]] <- c(
      mu_0 = fit_mle$mu_0,
      mu_1 = fit_mle$mu_1,
      sigma_0 = fit_mle$sigma_0,
      sigma_1 = fit_mle$sigma_1,
      eta_0 = fit_mle$eta_0,
      eta_1 = fit_mle$eta_1,
      phi = fit_mle$phi,
      n = n,
      rho = NA,
      delta = delta,
      pi = NA,
      converged = as.numeric(fit_mle$converged),
      time = fit_mle$time
    )
    
    cat(sprintf(" Done.\n"))
  }
  
  # Scenario 3
  for(pi in c(0.10, 0.25, 0.40)) {
    scenario_name <- sprintf("S3_pi%.2f_n%d", pi, n)
    cat(sprintf("\n%s: Generating data...", scenario_name))
    
    data <- generate_S3(n, pi)
    
    cat(" Fitting L2...")
    fit_l2 <- fit_L2(data)
    
    cat(" Fitting MLE...")
    fit_mle <- fit_MLE(data, fit_l2)
    
    all_results[[scenario_name]] <- c(
      mu_0 = fit_mle$mu_0,
      mu_1 = fit_mle$mu_1,
      sigma_0 = fit_mle$sigma_0,
      sigma_1 = fit_mle$sigma_1,
      eta_0 = fit_mle$eta_0,
      eta_1 = fit_mle$eta_1,
      phi = fit_mle$phi,
      n = n,
      rho = NA,
      delta = NA,
      pi = pi,
      converged = as.numeric(fit_mle$converged),
      time = fit_mle$time
    )
    
    cat(sprintf(" Done.\n"))
  }
}

# Convert to matrix
results_matrix <- do.call(rbind, all_results)
rownames(results_matrix) <- names(all_results)

cat("\n\nResults Matrix:\n")
print(results_matrix)

# Save
output <- list(
  true_params = TRUE_PARAMS,
  sample_sizes = SAMPLE_SIZES,
  results_matrix = results_matrix
)

random_id <- paste0(sample(c(letters, 0:9), 8, replace = TRUE), collapse = "")
filename <- sprintf("sim_iter_id%s.rds", random_id)

saveRDS(output, filename)
cat(sprintf("\n\nSaved: %s\n", filename))