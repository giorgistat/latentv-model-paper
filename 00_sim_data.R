############################################################
## Synthetic dataset simulation for AMA1 and MSP
## using MLEs from fitted models (no real data used)
##
## Requires:
##  - experiment/article/ams/joint_fit_full.rds   (AMA1 full-lik fit)
##  - experiment/article/msp/joint_fit_msp.rds    (MSP full-lik fit)
##
## Output:
##  - experiment/article/synthetic/synthetic_AMA1_MSP_dataset.csv
##
## Variables simulated:
##  - id
##  - age
##  - T_ama (latent seroreactivity for AMA1)
##  - log_AMA1, ama_norm
##  - T_msp (latent seroreactivity for MSP)
##  - log_MSP, msp_norm
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(tibble)
})

set.seed(123)

############################################################
## 1. Load fitted models (MLEs)
############################################################

# AMA1 full-likelihood fit (theta length 12)
fit_ama <- readRDS("ama_fit.rds")
theta_hat_ama <- fit_ama$par

# MSP full-likelihood fit (theta length 10)
fit_msp <- readRDS("msp_fit.rds")
theta_hat_msp <- fit_msp$par

############################################################
## 2. Common helpers
############################################################

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

t_grid  <- seq(0.01, 1 - 0.01, length.out = 200)
age_cut <- 15

############################################################
## 3. Reconstruct AMA1 model from theta_hat_ama
##
## Model:
##  - T | age < 15 ~ pi(a) * Beta(1, 10000) + (1 - pi(a)) * Beta(alpha2, 1)
##        pi(a) = exp(-r1 * a), p0 = 1
##  - T | age >= 15 ~ Beta(mean = logit^{-1}(eta0 + eta1 * log(a)),
##                         precision = phi)
##  - Y | T = t ~ N(mu_t, sigma_t^2)
##      mu_t     = t * mu2 + (1 - t) * mu1
##      sigma_t^2 = t * sigma2^2 + (1 - t) * sigma1^2
############################################################

## Unpack AMA parameters
mu1_ama    <- theta_hat_ama[1]
mu2_ama    <- theta_hat_ama[2]
sigma1_ama <- exp(theta_hat_ama[3])
sigma2_ama <- exp(theta_hat_ama[4])

alpha1_ama <- 1
beta1_ama  <- 10000
alpha2_ama <- 1 + exp(theta_hat_ama[7])
beta2_ama  <- 1

p0_ama   <- 1
r1_ama   <- exp(theta_hat_ama[10])
eta1_ama <- theta_hat_ama[11]
phi_ama  <- exp(theta_hat_ama[12])

## Precompute truncated Beta(1,10000) on t_grid
log_dbeta1_const  <- dbeta(t_grid, alpha1_ama, beta1_ama, log = TRUE)
lognorm1_const    <- log_sum_exp(log_dbeta1_const)
norm_dbeta1_const <- exp(log_dbeta1_const - lognorm1_const)
ET1_trunc_const   <- sum(t_grid * norm_dbeta1_const)

## Truncated Beta(alpha2_ama,1) on t_grid
log_dbeta2_hat  <- dbeta(t_grid, alpha2_ama, beta2_ama, log = TRUE)
lognorm2_hat    <- log_sum_exp(log_dbeta2_hat)
norm_dbeta2_hat <- exp(log_dbeta2_hat - lognorm2_hat)
ET2_trunc_hat   <- sum(t_grid * norm_dbeta2_hat)

## Continuity at 15 to get eta0_ama
pi15_hat <- p0_ama * exp(-r1_ama * age_cut)
m15_hat  <- pi15_hat * ET1_trunc_const + (1 - pi15_hat) * ET2_trunc_hat
m15_hat  <- pmin(pmax(m15_hat, 1e-6), 1 - 1e-6)

eta0_ama <- qlogis(m15_hat) - eta1_ama * log(age_cut)

## AMA1: simulate T | age
simulate_T_ama_given_age <- function(a_vec) {
  n <- length(a_vec)
  T_out <- numeric(n)

  # Below 15: mixture of two Betas
  idx_low  <- which(a_vec < age_cut)
  idx_high <- which(a_vec >= age_cut)

  if (length(idx_low) > 0) {
    a_low <- a_vec[idx_low]
    pi_a  <- p0_ama * exp(-r1_ama * a_low)
    pi_a  <- pmin(pmax(pi_a, 0), 1)

    u_mix <- runif(length(a_low))
    use_comp1 <- (u_mix < pi_a)

    # Component 1: Beta(1, 10000)
    T_out[idx_low[use_comp1]] <-
      rbeta(sum(use_comp1), alpha1_ama, beta1_ama)
    # Component 2: Beta(alpha2_ama, 1)
    T_out[idx_low[!use_comp1]] <-
      rbeta(sum(!use_comp1), alpha2_ama, beta2_ama)
  }

  # Above / equal 15: single Beta with mean logit^{-1}(eta0 + eta1 log a), precision phi_ama
  if (length(idx_high) > 0) {
    a_high <- a_vec[idx_high]
    eta_a  <- eta0_ama + eta1_ama * log(a_high)
    m_a    <- plogis(eta_a)
    m_a    <- pmin(pmax(m_a, 1e-6), 1 - 1e-6)

    alpha_a <- m_a * phi_ama
    beta_a  <- (1 - m_a) * phi_ama

    T_out[idx_high] <- rbeta(length(idx_high), alpha_a, beta_a)
  }

  T_out
}

## AMA1: simulate Y | T
simulate_Y_ama_given_T <- function(T_vec) {
  mu_t     <- T_vec * mu2_ama + (1 - T_vec) * mu1_ama
  sigma2_t <- T_vec * sigma2_ama^2 + (1 - T_vec) * sigma1_ama^2
  sigma_t  <- sqrt(pmax(sigma2_t, 1e-12))
  rnorm(length(T_vec), mean = mu_t, sd = sigma_t)
}

############################################################
## 4. Reconstruct MSP model from theta_hat_msp
##
## Model:
##  - T | age ~ Beta(alpha_age, beta_age)
##    alpha_age = exp(alpha_intercept + alpha_slope * log(age))
##    beta_age  = exp(beta_intercept +
##                    (beta_slope + 1(age > change_point)*beta_diff) * log(age))
##
##  - Y | T = t ~ N(mu_t, sigma_t^2) with same mixture structure
############################################################

mu1_msp    <- theta_hat_msp[1]
mu2_msp    <- theta_hat_msp[2]
sigma1_msp <- exp(theta_hat_msp[3])
sigma2_msp <- exp(theta_hat_msp[4])

alpha_intercept_msp <- theta_hat_msp[5]
alpha_slope_msp     <- theta_hat_msp[6]
beta_intercept_msp  <- theta_hat_msp[7]
beta_slope_msp      <- theta_hat_msp[8]
change_point_msp    <- theta_hat_msp[9]
beta_diff_msp       <- theta_hat_msp[10]

## MSP: simulate T | age
simulate_T_msp_given_age <- function(a_vec) {
  n <- length(a_vec)

  alpha_age <- exp(alpha_intercept_msp + alpha_slope_msp * log(a_vec))
  beta_age  <- exp(
    beta_intercept_msp +
      (beta_slope_msp + (a_vec > change_point_msp) * beta_diff_msp) * log(a_vec)
  )

  rbeta(n, alpha_age, beta_age)
}

## MSP: simulate Y | T
simulate_Y_msp_given_T <- function(T_vec) {
  mu_t     <- T_vec * mu2_msp + (1 - T_vec) * mu1_msp
  sigma2_t <- T_vec * sigma2_msp^2 + (1 - T_vec) * sigma1_msp^2
  sigma_t  <- sqrt(pmax(sigma2_t, 1e-12))
  rnorm(length(T_vec), mean = mu_t, sd = sigma_t)
}

############################################################
## 5. Simulate full synthetic dataset
############################################################

simulate_AMA1_MSP_dataset <- function(n = 2000,
                                      age_min = 1,
                                      age_max = 60,
                                      seed = 2025,
                                      id_prefix = "syn") {
  set.seed(seed)

  # 1) Simulate ages (no real data used)
  age <- runif(n, min = age_min, max = age_max)

  # 2) AMA1 latent T and outcome
  T_ama   <- simulate_T_ama_given_age(age)
  log_AMA <- simulate_Y_ama_given_T(T_ama)
  ama_val <- exp(log_AMA)

  # 3) MSP latent T and outcome
  T_msp   <- simulate_T_msp_given_age(age)
  log_MSP <- simulate_Y_msp_given_T(T_msp)
  msp_val <- exp(log_MSP)

  tibble(
    id        = paste0(id_prefix, seq_len(n)),
    age       = age,
    T_ama     = T_ama,
    log_AMA1  = log_AMA,
    ama_norm  = ama_val,
    T_msp     = T_msp,
    log_MSP   = log_MSP,
    msp_norm  = msp_val
  )
}

## Run the simulation
synthetic_dat <- simulate_AMA1_MSP_dataset(
  n       = 3000,   # choose sample size
  age_min = 1,
  age_max = 60,
  seed    = 2025
)

## Save synthetic dataset for GitHub
write.csv(
  synthetic_dat,
  file = "sim_data.csv",
  row.names = FALSE
)
