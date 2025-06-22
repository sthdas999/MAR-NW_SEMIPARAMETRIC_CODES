
# semiparametric_simulation.R
# Simulation study for testing independence in semiparametric regression under MAR
set.seed(123)
library(np)        # for kernel regression
library(stats)     # for logistic and basic distributions
library(dplyr)     # for data wrangling

# Parameters
n <- 500                     # Sample size
replications <- 1000        # Number of Monte Carlo repetitions
beta0 <- 1.0                 # True intercept
beta1 <- 2.0                 # Coefficient for X
beta2 <- 1.5                 # Coefficient for Z
sigma2 <- 1.0                # Variance of epsilon
gamma1 <- 1.0                # Slope in MAR model
bandwidths <- seq(0.1, 1, 0.1)
missing_rates <- c(0.1, 0.3, 0.5)

# Generate covariate Z: choose "normal" or "t5"
generate_Z <- function(type = "normal") {
  if (type == "normal") {
    rnorm(n)
  } else if (type == "t5") {
    rt(n, df = 5)
  } else {
    stop("Invalid distribution type")
  }
}

# MAR missingness mechanism
generate_missing <- function(X, target_rate, gamma1 = 1) {
  logit_fn <- function(g0) {
    mean(1 / (1 + exp(-(g0 + gamma1 * X)))) - target_rate
  }
  gamma0 <- uniroot(logit_fn, c(-10, 10))$root
  prob <- 1 / (1 + exp(-(gamma0 + gamma1 * X)))
  rbinom(length(X), 1, prob)
}

# Data generator
simulate_data <- function(dist_Z = "normal", mar_rate = 0.1) {
  X <- runif(n)
  Z <- generate_Z(dist_Z)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  Y <- beta0 + beta1 * X + beta2 * Z + epsilon
  M <- generate_missing(X, mar_rate, gamma1)
  data.frame(X, Z, Y, M)
}

# Nadaraya-Watson estimator
nw_estimate <- function(X, Y, eval_points, h) {
  n <- length(X)
  m_hat <- numeric(length(eval_points))
  for (j in seq_along(eval_points)) {
    weights <- (1 - ((X - eval_points[j]) / h)^2) * (abs((X - eval_points[j]) / h) <= 1)
    denom <- sum(weights)
    if (denom > 0) {
      m_hat[j] <- sum(weights * Y) / denom
    } else {
      m_hat[j] <- NA
    }
  }
  m_hat
}

# Main simulation loop (example for one rep)
run_one_simulation <- function(dist_Z = "normal", mar_rate = 0.3, h = 0.2) {
  dat <- simulate_data(dist_Z, mar_rate)
  obs <- dat %>% filter(M == 1)

  # Step 1: Estimate mY(X) and mZ(X) using NW
  mY <- nw_estimate(obs$X, obs$Y, obs$X, h)
  mZ <- nw_estimate(obs$X, obs$Z, obs$X, h)

  # Step 2: Residuals
  epsY <- obs$Y - mY
  epsZ <- obs$Z - mZ

  # Step 3: Estimate beta
  beta_hat <- sum(epsY * epsZ) / sum(epsZ^2)

  # Step 4: Transform Y and estimate m(X)
  Y_tilde <- obs$Y - beta_hat * obs$Z
  m_X <- nw_estimate(obs$X, Y_tilde, obs$X, h)

  # Step 5: Impute missing Y
  missing <- dat %>% filter(M == 0)
  m_miss <- nw_estimate(obs$X, Y_tilde, missing$X, h)
  Y_imp <- c(obs$Y, m_miss)
  Z_all <- c(obs$Z, missing$Z)
  X_all <- c(obs$X, missing$X)

  # Final dataset
  data.frame(X = X_all, Z = Z_all, Y = Y_imp)
}

# Example execution for 1 bandwidth and missing rate
res <- run_one_simulation(dist_Z = "normal", mar_rate = 0.3, h = 0.2)
head(res)
