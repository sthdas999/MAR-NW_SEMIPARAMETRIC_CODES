# Required libraries
library(MASS)  # For mvrnorm if needed
library(stats)

# Set parameters
beta0 <- 1.0
beta1 <- 0.5
beta2 <- 0.5
sigma2 <- 1.0
n <- 500
n_sim <- 10000
gamma1 <- 0.7  # Moderate MAR dependence

# Function to generate data
generate_data <- function(n, beta0, beta1, beta2, sigma2, gamma0, gamma1, z_dist = "normal") {
  X <- runif(n, 0, 1)
  if (z_dist == "normal") {
    Z <- rnorm(n, 0, 1)
  } else if (z_dist == "t5") {
    Z <- rt(n, df = 5)
  }
  epsilon <- rnorm(n, 0, sqrt(sigma2))
  Y <- beta0 + beta1 * X + beta2 * Z + epsilon
  
  # MAR missingness
  logit_p <- gamma0 + gamma1 * X
  p_miss <- exp(logit_p) / (1 + exp(logit_p))
  M <- rbinom(n, 1, p_miss)
  
  list(X = X, Z = Z, Y = Y, M = M)
}

# Calibrate gamma0 for desired missing rate
find_gamma0 <- function(target_missing, gamma1) {
  f <- function(gamma0) {
    X <- runif(1e5, 0, 1)
    p <- exp(gamma0 + gamma1 * X) / (1 + exp(gamma0 + gamma1 * X))
    mean(p) - target_missing
  }
  uniroot(f, c(-10, 10))$root
}

# Simulate test statistic under null or alternative
simulate_test <- function(data) {
  # Placeholder: replace with actual test statistic computation
  # Example: simple t-test or custom semiparametric test
  # Return a numeric test statistic
  abs(mean(data$Y[!data$M]))  # Example: mean of observed Y
}

# Run simulations
run_simulations <- function(n_sim, n, beta0, beta1, beta2, sigma2, gamma0, gamma1, z_dist, null = TRUE) {
  test_stats <- numeric(n_sim)
  for (i in 1:n_sim) {
    data <- generate_data(n, beta0, beta1, beta2, sigma2, gamma0, gamma1, z_dist)
    if (null) {
      data$Y <- beta0 + beta1 * data$X + beta2 * data$Z + rnorm(n, 0, sqrt(sigma2))
    }
    test_stats[i] <- simulate_test(data)
  }
  test_stats
}

# Example workflow
target_missing <- 0.1
gamma0 <- find_gamma0(target_missing, gamma1)

# Simulate under null for Type-I error calibration
null_stats_Tn2 <- run_simulations(n_sim, n, beta0, beta1, beta2, sigma2, gamma0, gamma1, "normal", null = TRUE)
null_cutoff_Tn2 <- quantile(null_stats_Tn2, 0.95)

# Simulate under alternative for power
alt_stats_Tn2 <- run_simulations(n_sim, n, beta0, beta1, beta2, sigma2, gamma0, gamma1, "normal", null = FALSE)
empirical_power_Tn2 <- mean(alt_stats_Tn2 > null_cutoff_Tn2)

# Output
cat("Size-adjusted empirical power for Tn2:", empirical_power_Tn2 * 100, "%\n")
