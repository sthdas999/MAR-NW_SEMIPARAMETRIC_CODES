
# =====================================
# Combined R Code for Power Curves
# Examples 1 & 2 (Tn1, Tn2, Tn3 Statistics)
# MAR Mechanism with NW and ILLS Estimation
# =====================================

# Load Libraries
library(KernSmooth)
library(MASS)
library(np)
library(dplyr)
library(ggplot2)

# Function: U-Statistics Tn1, Tn2, Tn3
Tn1 <- function(x, y) {
  n <- length(x)
  idx <- combn(n, 2)
  sgns <- apply(idx, 2, function(i) sign((x[i[1]] - x[i[2]]) * (y[i[1]] - y[i[2]])))
  return(mean(sgns))
}
Tn2 <- function(x, y) {
  n <- length(x)
  idx <- combn(n, 4)
  kernel <- apply(idx, 2, function(i) {
    ax <- sign(abs(x[i[1]] - x[i[2]]) + abs(x[i[3]] - x[i[4]]) - abs(x[i[1]] - x[i[3]]) - abs(x[i[2]] - x[i[4]]))
    ay <- sign(abs(y[i[1]] - y[i[2]]) + abs(y[i[3]] - y[i[4]]) - abs(y[i[1]] - y[i[3]]) - abs(y[i[2]] - y[i[4]]))
    ax * ay
  })
  return(mean(kernel))
}
Tn3 <- function(x, y) {
  n <- length(x)
  y_ord <- y[order(x)]
  y_diff <- y_ord
  y_diff[1] <- y_ord[2] - 3*y_ord[2] + 3*y_ord[1] - y_ord[1]
  y_diff[2] <- y_ord[3] - 3*y_ord[3] + 3*y_ord[2] - y_ord[1]
  for (i in 3:n) {
    y_diff[i] <- y_ord[i] - 3*y_ord[i-1] + 3*y_ord[i-2] - y_ord[i-3]
  }
  return(cor(x, y_diff, method = "spearman"))
}

# Estimation functions
NW_estimation <- function(x, y) {
  h <- dpill(x, y)
  yhat <- ksmooth(x, y, "normal", bandwidth = h)$y
  return(yhat)
}

ILLS_estimation <- function(x, y) {
  fit <- npreg(y ~ x, regtype = "ll", ckertype = "epanechnikov")
  return(fitted(fit))
}

# Simulation setup
simulate_example <- function(example_id, missing_prop, method = c("NW", "ILLS"), nrep = 100, n = 100) {
  method <- match.arg(method)
  Tn1_power <- Tn2_power <- Tn3_power <- numeric(nrep)
  beta0 <- 2
  for (i in 1:nrep) {
    x <- rnorm(n)
    z <- rnorm(n)
    if (example_id == 1) {
      eps <- rnorm(n, 0, sqrt((1 + 5*x)/100))
    } else {
      l_x <- 1 / (3 * x)
      l_x[is.infinite(l_x)] <- 1000
      eps <- (rchisq(n, df = l_x) - l_x) / (10 * sqrt(2 * l_x))
    }
    y <- beta0 * z + eps
    m_x <- ifelse(method == "NW", NW_estimation(x, y), ILLS_estimation(x, y))
    prob_missing <- exp(m_x) / (1 + exp(m_x))
    m_idx <- sample(1:n, size = floor(missing_prop * n), prob = prob_missing)
    observed <- setdiff(1:n, m_idx)
    x_obs <- x[observed]
    y_obs <- y[observed]
    z_obs <- z[observed]
    beta_hat <- lm(y_obs ~ z_obs)$coefficients[2]
    yhat <- beta_hat * z + eps
    m_hat <- ifelse(method == "NW", NW_estimation(x[observed], yhat[observed]), ILLS_estimation(x[observed], yhat[observed]))
    residuals <- yhat - m_hat
    Tn1_power[i] <- Tn1(x[observed], residuals[observed])
    Tn2_power[i] <- Tn2(x[observed], residuals[observed])
    Tn3_power[i] <- Tn3(x[observed], residuals[observed])
  }
  return(list(Tn1 = Tn1_power, Tn2 = Tn2_power, Tn3 = Tn3_power))
}

# Run for all scenarios (set manually or automate)
# [Code snippet only; simulate_example() is now available for batch runs.]
