############################################################
## MSP model + full-likelihood fit in Y
## + histogram–envelope diagnostics by age group
##
##
## Model:
##  - Y = log(MSP)
##  - T | age ~ Beta(alpha_age, beta_age)
##    with
##      alpha_age = exp(alpha_intercept + alpha_slope * log(age))
##      beta_age  = exp(beta_intercept +
##                      (beta_slope + 1(age > change_point)*beta_diff) * log(age))
##
##  - Y | T = t ~ N(mu_t, sigma_t^2),
##        mu_t      = t * mu2 + (1 - t) * mu1
##        sigma_t^2 = t * sigma2^2 + (1 - t) * sigma1^2
##
## Likelihood:
##  - For each i:
##      f_Y(y_i | age_i) = ∫ N(y_i | mu_t(t), sigma_t^2(t))
##                               Beta(t | alpha_age_i, beta_age_i) dt
##    approximated on a fixed grid t_grid.
##
## Outputs:
##  1. joint_fit_msp.rds : optim() result for full likelihood
##  2. Param table printed to console
##  3. Histogram–envelope by age group
##     [1, 6), [6, 10), [10, 15), [15, 20), [20, 30), [45, 100)
##  4. Final plot: E[log(MSP)|age] (curve) + empirical means (points)
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

# Filter and log-transform MSP (consistent with your bootstrap code)
sel <- redhot$msp_norm > 0 & redhot$age > 0
redhot <- redhot[sel, ]

y_raw   <- redhot$msp_norm
age_raw <- redhot$age

# log-transform MSP
y1 <- log(y_raw)
y  <- y1  # alias

t_grid <- seq(0.01, 1 - 0.01, length.out = 100)
n_total <- length(y)

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
## 3. MSP log-likelihood (as in your bootstrap code)
############################################################

loglik_norm_mixture_alpha_age <- function(theta,
                                          y,
                                          age,
                                          t.grid = t_grid) {

  mu1 <- theta[1]
  mu2 <- theta[2]
  sigma1 <- exp(theta[3])
  sigma2 <- exp(theta[4])

  alpha_intercept <- theta[5]
  alpha_slope     <- theta[6]
  beta_intercept  <- theta[7]
  beta_slope      <- theta[8]
  change_point    <- theta[9]
  beta_diff       <- theta[10]

  loglik <- 0

  for (i in seq_along(y)) {
    a_i <- age[i]

    alpha_age <- exp(alpha_intercept + alpha_slope * log(a_i))
    beta_age  <- exp(beta_intercept +
                       (beta_slope + (a_i > change_point) * beta_diff) * log(a_i))

    prior_t <- dbeta(t.grid, alpha_age, beta_age, log = FALSE)

    lik_y_given_t <- dnorm_mixture(
      t        = t.grid,
      y        = y[i],
      mu1      = mu1,
      mu2      = mu2,
      sigma1   = sigma1,
      sigma2   = sigma2,
      log.p    = FALSE
    )

    dens_y <- sum(lik_y_given_t * prior_t)

    if (is.finite(dens_y) && dens_y > 0) {
      loglik <- loglik + log(dens_y)
    } else {
      return(-1e10)
    }
  }

  loglik
}

neg_loglik <- function(theta, y, age, t.grid = t_grid) {
  -loglik_norm_mixture_alpha_age(theta, y, age, t.grid)
}

############################################################
## 4. Optimisation
############################################################

# Starting values inferred from your theta_true used in bootstrap
theta_init <- c(
  -4.1019356,  1.1644738,  # mu1, mu2
  -0.4873009, -5.9488632,  # log(sigma1), log(sigma2)
  -0.3054093,  0.3390226,  # alpha_intercept, alpha_slope
  0.5554429,  0.1481352,  # beta_intercept, beta_slope
  11.3299291, -0.0662268   # change_point, beta_diff
)

fit <- optim(
  par     = theta_init,
  fn      = neg_loglik,
  y       = y,
  age     = age_raw,
  t.grid  = t_grid,
  method  = "BFGS",
  control = list(
    maxit  = 1000,
    reltol = 1e-6,
    trace  = 1,
    REPORT = 10
  )
)


cat("Converged (MSP):", fit$convergence == 0, "\n")
cat("Maximum log-lik (MSP):", -fit$value, "\n")

theta_hat <- fit$par
print(theta_hat)

############################################################
## 5. Parameters on original scale
############################################################

mu1_hat  <- theta_hat[1]
mu2_hat  <- theta_hat[2]
sigma1_hat <- exp(theta_hat[3])
sigma2_hat <- exp(theta_hat[4])

alpha_intercept_hat <- theta_hat[5]
alpha_slope_hat     <- theta_hat[6]
beta_intercept_hat  <- theta_hat[7]
beta_slope_hat      <- theta_hat[8]
change_point_hat    <- theta_hat[9]
beta_diff_hat       <- theta_hat[10]

param_table <- tibble(
  parameter = c("mu1", "mu2",
                "sigma1", "sigma2",
                "alpha_intercept", "alpha_slope",
                "beta_intercept", "beta_slope",
                "change_point", "beta_diff"),
  estimate  = c(mu1_hat, mu2_hat,
                sigma1_hat, sigma2_hat,
                alpha_intercept_hat, alpha_slope_hat,
                beta_intercept_hat, beta_slope_hat,
                change_point_hat, beta_diff_hat)
)

print(param_table)

############################################################
## 6. Simulator Y | age under fitted MSP model
############################################################

simulate_y_msp <- function(theta, age, t.grid = t_grid) {
  mu1 <- theta[1]
  mu2 <- theta[2]
  sigma1 <- exp(theta[3])
  sigma2 <- exp(theta[4])
  alpha_intercept <- theta[5]
  alpha_slope     <- theta[6]
  beta_intercept  <- theta[7]
  beta_slope      <- theta[8]
  change_point    <- theta[9]
  beta_diff       <- theta[10]

  # Age-specific beta parameters with logarithmic growth
  alpha_age <- exp(alpha_intercept + alpha_slope * log(age))
  beta_age  <- exp(beta_intercept +
                     (beta_slope + (age > change_point) * beta_diff) * log(age))

  # Sample T from age-specific Beta distribution
  t_val <- rbeta(length(age), alpha_age, beta_age)

  # Sample Y from normal mixture given T
  mu_t     <- t_val * mu2 + (1 - t_val) * mu1
  sigma2_t <- t_val * sigma2^2 + (1 - t_val) * sigma1^2

  rnorm(length(age), mean = mu_t, sd = sqrt(sigma2_t))
}

simulate_Y_given_age_vec <- function(age_vec, theta, n_sims = 1000) {
  n <- length(age_vec)
  y_mat <- matrix(NA_real_, nrow = n, ncol = n_sims)
  for (s in seq_len(n_sims)) {
    y_mat[, s] <- simulate_y_msp(theta, age_vec)
  }
  y_mat
}

############################################################
## 7. Histogram–envelope by age group (as for AMA1)
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

hist_envelope_one_group <- function(df_group,
                                    theta,
                                    n_sims = 1000,
                                    n_bins = 20,
                                    probs = c(0.025, 0.5, 0.975)) {
  y_obs <- df_group$y_log
  ages  <- df_group$age
  n     <- length(y_obs)

  # 1) simulate from the fitted model for this age group
  y_sim_mat <- simulate_Y_given_age_vec(ages, theta, n_sims = n_sims)

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

n_sims_env <- 1000
n_bins_env <- 20

df_env_all <- df_full %>%
  group_by(age_group) %>%
  group_modify(~{
    env_df <- hist_envelope_one_group(.x,
                                      theta = theta_hat,
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
    x = "log(MSP)",
    y = "Density",
    title = "Histogram–envelope diagnostic by age group (MSP)"
  ) +
  theme_bw()

fig_env

############################################################
## 8. E[log(MSP) | age] and empirical means (base R)
############################################################

expected_Y_given_age <- function(age_vec,
                                 theta,
                                 t.grid = t_grid) {

  mu1 <- theta[1]
  mu2 <- theta[2]

  alpha_intercept <- theta[5]
  alpha_slope     <- theta[6]
  beta_intercept  <- theta[7]
  beta_slope      <- theta[8]
  change_point    <- theta[9]
  beta_diff       <- theta[10]

  mu_t_grid <- t.grid * mu2 + (1 - t.grid) * mu1
  EY <- numeric(length(age_vec))

  for (i in seq_along(age_vec)) {
    a <- age_vec[i]

    alpha_age <- exp(alpha_intercept + alpha_slope * log(a))
    beta_age  <- exp(beta_intercept +
                       (beta_slope + (a > change_point) * beta_diff) * log(a))

    prior_t_raw <- dbeta(t.grid, alpha_age, beta_age)
    prior_t     <- prior_t_raw / sum(prior_t_raw)

    EY[i] <- sum(mu_t_grid * prior_t)
  }

  EY
}

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

# Model-based E[log(MSP) | age]
age_seq <- seq(min(age_raw), max(age_raw), length.out = 200)
EY_seq  <- expected_Y_given_age(age_seq, theta_hat)

# y-limits that include both model curve and empirical means
ylim_range <- range(c(EY_seq, age_df$y_mean), na.rm = TRUE)

plot(
  age_seq, EY_seq,
  type = "l",
  xlab = "Age (years)",
  ylab = "E[log(MSP) | age]",
  main = "Expected log-MSP vs age (age-dependent Beta model)",
  ylim = ylim_range
)

points(
  age_df$age_mid, age_df$y_mean,
  pch = 16,
  col = "red"
)

############################################################
## End of MSP script
############################################################
