############################################################
## Joint AMA1 model + full-likelihood fit in Y
## + histogram–envelope diagnostics by age group
##
##
## Model:
##  - Y = log(AMA1)
##  - T | age < 15 ~ pi(a) * Beta(1, 10000) + (1 - pi(a)) * Beta(alpha2, 1)
##        with  pi(a) = exp(-r1 * a)  (p0 = 1)
##  - T | age >= 15 ~ Beta(mean = logit^{-1}(eta0 + eta1 * log(a)),
##                         precision = phi)
##    where eta0 is fixed by continuity at age 15 using truncated Betas
##
##  - Y | T = t ~ N(mu_t, sigma_t^2),
##        mu_t     = t * mu2 + (1 - t) * mu1
##        sigma_t^2 = t * sigma2^2 + (1 - t) * sigma1^2
##
## Fitting:
##  - Full log-likelihood with numerical integration over T on t_grid.
##
## Outputs:
##  1. joint_fit_full.rds : optim() result for full log-likelihood
##  2. Param table printed to console
##  3. Histogram–envelope by age group
##     [1, 6), [6, 10), [10, 15), [15, 20), [20, 30), [30, 45), [45, 100)
##  4. Final plot: E[log(AMA1)|age] (curve) + empirical means (points)
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

set.seed(123)

############################################################
## 1. Load and prepare data
############################################################

redhot <- read.csv("sim_data.csv")

# Same filter as before: > 0 and age > 1
sel <- redhot$ama_norm > 0 & redhot$age > 1
redhot <- redhot[sel, ]

y_raw   <- redhot$ama_norm
age_raw <- redhot$age

# log-transform AMA1
y1 <- log(y_raw)
y  <- y1

age_cut <- 15
t_grid  <- seq(0.01, 1 - 0.01, length.out = 100)

## Collapse ages to unique values and counts (for speed)
age_unique <- sort(unique(age_raw))
age_counts <- as.numeric(table(age_raw)[as.character(age_unique)])
n_age      <- length(age_unique)
n_total    <- length(age_raw)

## Map each individual to an index in age_unique
age_index <- match(age_raw, age_unique)

############################################################
## 2. Conditional Normal mixture Y | T
############################################################

dnorm_mixture <- function(t, y, mu1, mu2, sigma1, sigma2, log.p = FALSE) {
  mu_t     <- t * mu2 + (1 - t) * mu1
  sigma2_t <- t * sigma2^2 + (1 - t) * sigma1^2
  dens     <- dnorm(y, mean = mu_t, sd = sqrt(sigma2_t), log = log.p)
  if (log.p) {
    dens[sigma2_t <= 0] <- -Inf
  } else {
    dens[sigma2_t <= 0] <- 0
  }
  dens
}

############################################################
## 3. Helper: log-sum-exp
############################################################

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

############################################################
## 4. Precompute truncated Beta(1, 10000) over t_grid
############################################################

alpha1 <- 1
beta1  <- 10000

log_dbeta1_const  <- dbeta(t_grid, alpha1, beta1, log = TRUE)
lognorm1_const    <- log_sum_exp(log_dbeta1_const)
norm_dbeta1_const <- exp(log_dbeta1_const - lognorm1_const)

# Truncated mean of Beta(1, 10000) on t_grid
ET1_trunc_const <- sum(t_grid * norm_dbeta1_const)

############################################################
## 5. Full log-likelihood for the AMA1 T-model
##
## theta (length 12):
##  [1]  mu1
##  [2]  mu2
##  [3]  log(sigma1)
##  [4]  log(sigma2)
##  [5]  (unused; alpha1 fixed)
##  [6]  (unused; beta1 fixed)
##  [7]  alpha2 via alpha2 = 1 + exp(theta[7])
##  [8]  (unused; beta2 fixed = 1)
##  [9]  (unused; p0 fixed = 1)
##  [10] r1 via r1 = exp(theta[10])
##  [11] eta1 (slope on log(age) for age >= 15)
##  [12] log(phi) via phi = exp(theta[12])
############################################################

loglik_ama <- function(theta,
                       y,
                       age_unique,
                       age_index,
                       t_grid,
                       age_cut = 15) {

  mu1    <- theta[1]
  mu2    <- theta[2]
  sigma1 <- exp(theta[3])
  sigma2 <- exp(theta[4])

  alpha2 <- 1 + exp(theta[7])
  beta2  <- 1

  p0  <- 1
  r1  <- exp(theta[10])

  eta1 <- theta[11]
  phi  <- exp(theta[12])

  # Truncated Beta(alpha2,1) on t_grid
  log_dbeta2  <- dbeta(t_grid, alpha2, beta2, log = TRUE)
  lognorm2    <- log_sum_exp(log_dbeta2)
  norm_dbeta2 <- exp(log_dbeta2 - lognorm2)
  ET2_trunc   <- sum(t_grid * norm_dbeta2)

  # Continuity at age 15 via truncated means
  pi15      <- p0 * exp(-r1 * age_cut)
  m15_minus <- pi15 * ET1_trunc_const + (1 - pi15) * ET2_trunc
  m15_minus <- pmin(pmax(m15_minus, 1e-6), 1 - 1e-6)

  # For age >= 15: mean m(a) = logit^{-1}(eta0 + eta1 * log(a))
  eta0 <- qlogis(m15_minus) - eta1 * log(age_cut)

  # Prior over T | unique ages (truncated on t_grid, normalised)
  K  <- length(t_grid)
  prior_age_mat <- matrix(NA_real_, nrow = length(age_unique), ncol = K)

  for (u in seq_along(age_unique)) {
    a_u <- age_unique[u]

    if (a_u < age_cut) {
      pi_u <- p0 * exp(-r1 * a_u)
      prior_age_mat[u, ] <- pi_u * norm_dbeta1_const + (1 - pi_u) * norm_dbeta2
    } else {
      eta_u <- eta0 + eta1 * log(a_u)
      m_u   <- plogis(eta_u)

      alpha_u <- m_u * phi
      beta_u  <- (1 - m_u) * phi

      log_dbeta_u <- dbeta(t_grid, alpha_u, beta_u, log = TRUE)
      lognorm_u   <- log_sum_exp(log_dbeta_u)
      prior_age_mat[u, ] <- exp(log_dbeta_u - lognorm_u)
    }
  }

  # Now compute log-likelihood over all individuals
  loglik <- 0
  n      <- length(y)

  for (i in seq_len(n)) {
    u_i <- age_index[i]

    lik_y_given_t <- dnorm_mixture(
      t_grid, y[i],
      mu1, mu2,
      sigma1, sigma2,
      log.p = FALSE
    )

    dens_y_i <- sum(lik_y_given_t * prior_age_mat[u_i, ])

    if (is.finite(dens_y_i) && dens_y_i > 0) {
      loglik <- loglik + log(dens_y_i)
    } else {
      return(-1e10)
    }
  }

  loglik
}

neg_loglik_ama <- function(theta) {
  -loglik_ama(theta,
              y          = y,
              age_unique = age_unique,
              age_index  = age_index,
              t_grid     = t_grid,
              age_cut    = age_cut)
}

############################################################
## 6. Initial values for theta and optimisation (full likelihood)
############################################################

theta_init <- c(
  -3.2824216,  # mu1
  0.7757835,  # mu2
  -0.2794324,  # log(sigma1)
  -2.6948222,  # log(sigma2)
  0.0,        # unused
  0.0,        # unused
  -1.0373180,  # alpha2 = 1 + exp(...)
  0.0,        # unused (beta2 fixed)
  0.0,        # unused (p0 fixed)
  -1.8952039,  # r1 = exp(...)
  0.05,       # eta1
  log(10)     # log(phi)
)

fit <- optim(
  par     = theta_init,
  fn      = neg_loglik_ama,
  method  = "BFGS",
  control = list(
    maxit  = 500,
    reltol = 1e-6,
    trace  = 1,
    REPORT = 10
  )
)

cat("Converged (AMA1, full lik):", fit$convergence == 0, "\n")
cat("Maximum log-lik (AMA1):", -fit$value, "\n")

theta_hat <- fit$par
print(theta_hat)

############################################################
## 7. Parameters on original scale
############################################################

mu1_hat    <- theta_hat[1]
mu2_hat    <- theta_hat[2]
sigma1_hat <- exp(theta_hat[3])
sigma2_hat <- exp(theta_hat[4])

alpha1_hat <- 1
beta1_hat  <- 10000
alpha2_hat <- 1 + exp(theta_hat[7])
beta2_hat  <- 1

p0_hat     <- 1
r1_hat     <- exp(theta_hat[10])

eta1_hat   <- theta_hat[11]
phi_hat    <- exp(theta_hat[12])

# Truncated Beta(alpha2_hat,1) on t_grid
log_dbeta2_hat  <- dbeta(t_grid, alpha2_hat, beta2_hat, log = TRUE)
lognorm2_hat    <- log_sum_exp(log_dbeta2_hat)
norm_dbeta2_hat <- exp(log_dbeta2_hat - lognorm2_hat)
ET2_trunc_hat   <- sum(t_grid * norm_dbeta2_hat)

# Continuity at 15
pi15_hat <- p0_hat * exp(-r1_hat * age_cut)
m15_hat  <- pi15_hat * ET1_trunc_const + (1 - pi15_hat) * ET2_trunc_hat
m15_hat  <- pmin(pmax(m15_hat, 1e-6), 1 - 1e-6)
eta0_hat <- qlogis(m15_hat) - eta1_hat * log(age_cut)

param_table <- tibble(
  parameter = c("mu1", "mu2",
                "sigma1", "sigma2",
                "alpha1", "beta1",
                "alpha2", "beta2",
                "p0", "r1",
                "eta0", "eta1",
                "phi"),
  estimate  = c(mu1_hat, mu2_hat,
                sigma1_hat, sigma2_hat,
                alpha1_hat, beta1_hat,
                alpha2_hat, beta2_hat,
                p0_hat, r1_hat,
                eta0_hat, eta1_hat,
                phi_hat)
)

print(param_table)

############################################################
## 8. E[Y | age] function under fitted model
############################################################

expected_Y_given_age <- function(age_vec, theta) {

  mu1    <- theta[1]
  mu2    <- theta[2]
  sigma1 <- exp(theta[3])
  sigma2 <- exp(theta[4])

  alpha2 <- 1 + exp(theta[7])
  beta2  <- 1

  p0   <- 1
  r1   <- exp(theta[10])
  eta1 <- theta[11]
  phi  <- exp(theta[12])

  # Truncated Beta(alpha2,1) on t_grid
  log_dbeta2  <- dbeta(t_grid, alpha2, beta2, log = TRUE)
  lognorm2    <- log_sum_exp(log_dbeta2)
  norm_dbeta2 <- exp(log_dbeta2 - lognorm2)
  ET2_trunc   <- sum(t_grid * norm_dbeta2)

  # Continuity at 15 via truncated means
  pi15      <- p0 * exp(-r1 * age_cut)
  m15_minus <- pi15 * ET1_trunc_const + (1 - pi15) * ET2_trunc
  m15_minus <- pmin(pmax(m15_minus, 1e-6), 1 - 1e-6)
  eta0      <- qlogis(m15_minus) - eta1 * log(age_cut)

  # Conditional mean of Y given T
  mu_t_grid <- t_grid * mu2 + (1 - t_grid) * mu1

  EY <- numeric(length(age_vec))

  for (i in seq_along(age_vec)) {
    a_i <- age_vec[i]

    if (a_i < age_cut) {
      # Below 15: mixture of Beta(1, 10000) and Beta(alpha2, 1)
      pi_i    <- p0 * exp(-r1 * a_i)
      prior_t <- pi_i * norm_dbeta1_const + (1 - pi_i) * norm_dbeta2
    } else {
      # Above 15: single Beta with logit-mean and precision phi
      eta_i <- eta0 + eta1 * log(a_i)
      m_i   <- plogis(eta_i)

      alpha_i <- m_i * phi
      beta_i  <- (1 - m_i) * phi

      log_dbeta_i <- dbeta(t_grid, alpha_i, beta_i, log = TRUE)
      lognorm_i   <- log_sum_exp(log_dbeta_i)
      prior_t     <- exp(log_dbeta_i - lognorm_i)
    }

    EY[i] <- sum(mu_t_grid * prior_t)
  }

  EY
}

############################################################
## 9. Histogram–envelope by age group (same classes as before)
##
## Age groups: [1, 6), [6, 10), [10, 15),
##              [15, 20), [20, 30), [30, 45), [45, 100)
############################################################

make_bracket_labels <- function(breaks) {
  paste0("[", head(breaks, -1), ", ", tail(breaks, -1), ")")
}

age_group_breaks <- c(1, 6, 10, 15, 20, 30, 45, 100)
age_group_labels <- make_bracket_labels(age_group_breaks)

df_full <- tibble(
  age   = age_raw,
  y_log = y1
) %>%
  mutate(
    age_group = cut(
      age,
      breaks = age_group_breaks,
      labels = age_group_labels,
      right = FALSE,          # [a, b)
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(age_group))

# Simulator from the fitted model: Y | age
simulate_Y_given_age_vec <- function(age_vec, n_sims = 1000) {

  ages_unique  <- sort(unique(age_vec))
  mu_t_grid    <- t_grid * mu2_hat + (1 - t_grid) * mu1_hat
  sigma_t_grid <- sqrt(t_grid * sigma2_hat^2 + (1 - t_grid) * sigma1_hat^2)

  # Precompute prior over T | age for each unique age
  prior_age_list <- vector("list", length(ages_unique))
  names(prior_age_list) <- as.character(ages_unique)

  for (k in seq_along(ages_unique)) {
    a <- ages_unique[k]
    if (a < age_cut) {
      pi_a    <- p0_hat * exp(-r1_hat * a)
      prior_t <- pi_a * norm_dbeta1_const + (1 - pi_a) * norm_dbeta2_hat
    } else {
      eta_a <- eta0_hat + eta1_hat * log(a)
      m_a   <- plogis(eta_a)
      alpha_a <- m_a * phi_hat
      beta_a  <- (1 - m_a) * phi_hat

      log_dbeta_a <- dbeta(t_grid, alpha_a, beta_a, log = TRUE)
      lognorm_a   <- log_sum_exp(log_dbeta_a)
      prior_t     <- exp(log_dbeta_a - lognorm_a)
    }
    prior_age_list[[as.character(a)]] <- prior_t
  }

  n      <- length(age_vec)
  y_mat  <- matrix(NA_real_, nrow = n, ncol = n_sims)

  for (s in seq_len(n_sims)) {
    t_idx <- integer(n)
    for (i in seq_len(n)) {
      prior_t <- prior_age_list[[as.character(age_vec[i])]]
      t_idx[i] <- sample(seq_along(t_grid), size = 1, prob = prior_t)
    }
    mu_vec    <- mu_t_grid[t_idx]
    sigma_vec <- sigma_t_grid[t_idx]
    y_mat[, s] <- rnorm(n, mean = mu_vec, sd = sigma_vec)
  }

  y_mat
}

# Function to build histogram envelope for one age group
hist_envelope_one_group <- function(df_group,
                                    n_sims = 1000,
                                    n_bins = 20,
                                    probs = c(0.025, 0.5, 0.975)) {
  y_obs <- df_group$y_log
  ages  <- df_group$age
  n     <- length(y_obs)

  # 1) simulate from the fitted model for this age group
  y_sim_mat <- simulate_Y_given_age_vec(ages, n_sims = n_sims)

  # 2) choose common breaks that span both observed and simulated values
  range_all <- range(c(y_obs, y_sim_mat))
  breaks    <- seq(range_all[1], range_all[2], length.out = n_bins + 1)

  # 3) observed histogram
  h_obs   <- hist(y_obs, breaks = breaks, plot = FALSE)
  mids    <- h_obs$mids
  widths  <- diff(breaks)
  dens_obs <- h_obs$counts / (n * widths)

  # 4) simulated histograms using the SAME breaks
  dens_mat <- matrix(NA_real_, nrow = length(mids), ncol = n_sims)
  for (s in seq_len(n_sims)) {
    h_s <- hist(y_sim_mat[, s], breaks = breaks, plot = FALSE)
    dens_mat[, s] <- h_s$counts / (n * widths)
  }

  # 5) envelope
  dens_q <- apply(dens_mat, 1, quantile, probs = probs)

  tibble(
    mid          = mids,
    density_obs  = dens_obs,
    density_low  = dens_q[1, ],
    density_med  = dens_q[2, ],
    density_high = dens_q[3, ]
  )
}

# Build envelope for all groups
n_sims_env <- 1000
n_bins_env <- 20

df_env_all <- df_full %>%
  group_by(age_group) %>%
  group_modify(~{
    env_df <- hist_envelope_one_group(.x,
                                      n_sims = n_sims_env,
                                      n_bins = n_bins_env)
    env_df
  }) %>%
  ungroup()

fig_env <- ggplot(df_env_all, aes(x = mid)) +
  geom_col(aes(y = density_obs),
           fill = "grey80", colour = "grey50") +
  geom_ribbon(aes(ymin = density_low, ymax = density_high),
              alpha = 0.25) +
  geom_line(aes(y = density_med)) +
  facet_wrap(~ age_group, scales = "free_y") +
  labs(
    x = "log(AMA1)",
    y = "Density",
    title = "Histogram–envelope diagnostic by age group (AMA1, full lik fit)"
  ) +
  theme_bw()

fig_env

############################################################
## 10. Empirical mean of Y by age class and plot (base R)
############################################################

# Empirical means in 1-year age bins
age_breaks <- seq(floor(min(age_raw)), ceiling(max(age_raw)), by = 1)

age_df <- tibble(
  age   = age_raw,
  y_log = y1
) %>%
  mutate(
    age_class = cut(age,
                    breaks = age_breaks,
                    include.lowest = TRUE,
                    right = FALSE)
  ) %>%
  group_by(age_class) %>%
  summarise(
    age_mid = mean(age),
    y_mean  = mean(y_log),
    .groups = "drop"
  )

# Model-based E[Y | age]
age_seq <- seq(min(age_raw), max(age_raw), length.out = 200)
EY_seq  <- expected_Y_given_age(age_seq, theta_hat)

# y-limits that include both model curve and empirical means
ylim_range <- range(c(EY_seq, age_df$y_mean), na.rm = TRUE)

plot(
  age_seq, EY_seq,
  type = "l",
  xlab = "Age (years)",
  ylab = "E[log(AMA1) | age]",
  main = "Expected log-AMA1 vs age (full-likelihood fit)",
  ylim = ylim_range
)

points(
  age_df$age_mid, age_df$y_mean,
  pch = 16,
  col = "red"
)

############################################################
## End of AMA1 full-likelihood script
############################################################
