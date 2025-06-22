
# Load necessary packages
library(MASS)       # for mvrnorm
library(KernSmooth) # for kernel regression
library(energy)     # for distance correlation (if needed)
set.seed(123)

# Simulation Parameters
n <- 200                      # Sample size
beta <- 1.5                   # True value of beta
rho <- 0.5                    # Correlation between X and Z
gamma <- 1                   # Control parameter for dependence in alternatives

# Generate correlated covariates (X, Z)
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
covariates <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
X <- covariates[, 1]
Z <- covariates[, 2]

# Nonlinear regression function m(x)
m <- function(x) { 1.5 * sin(pi * x) }

# Error generation
generate_error <- function(X, alt = FALSE, gamma = 1) {
  if (alt) {
    # Under alternative: dependent error
    e <- rnorm(n, mean = gamma * cos(2 * pi * X), sd = 1)
  } else {
    # Under null: independent error
    e <- rnorm(n, mean = 0, sd = 1)
  }
  return(e)
}

# Generate Y under null or alternative
alt <- TRUE  # set to FALSE for null hypothesis
epsilon <- generate_error(X, alt = alt, gamma = gamma)
Y <- Z * beta + m(X) + epsilon

# MAR missingness: Missing depends on X only
p_miss <- 0.8 - 0.5 * X
delta <- rbinom(n, 1, prob = pmin(pmax(p_miss, 0.1), 0.9))  # observed indicator

# Observed data
X_obs <- X[delta == 1]
Z_obs <- Z[delta == 1]
Y_obs <- Y[delta == 1]

# Estimate beta using Robinson’s method
library(np)
h <- dpill(X_obs, Y_obs)  # bandwidth selection
mY_hat <- ksmooth(X_obs, Y_obs, kernel = "normal", bandwidth = h, x.points = X_obs)$y
mZ_hat <- ksmooth(X_obs, Z_obs, kernel = "normal", bandwidth = h, x.points = X_obs)$y

resY <- Y_obs - mY_hat
resZ <- Z_obs - mZ_hat
beta_hat <- sum(resY * resZ) / sum(resZ^2)

# Impute missing Y using NW estimator on residuals
Y_tilde <- Y_obs - Z_obs * beta_hat
mX_hat <- ksmooth(X_obs, Y_tilde, kernel = "normal", bandwidth = h, x.points = X)$y
Y_imp <- Z * beta_hat + mX_hat  # imputed full Y

# Compute third-order differences on ordered pairs
df <- data.frame(X = X, Y_imp = Y_imp)
df_sorted <- df[order(df$X), ]
Y_star <- with(df_sorted, c(
  df_sorted$Y_imp[1] - 3 * df_sorted$Y_imp[1] + 3 * df_sorted$Y_imp[1] - df_sorted$Y_imp[1],  # padding
  sapply(3:(n-1), function(i) {
    df_sorted$Y_imp[i+1] - 3 * df_sorted$Y_imp[i] + 3 * df_sorted$Y_imp[i-1] - df_sorted$Y_imp[i-2]
  }),
  0, 0  # padding for symmetry
))

# Compute Kendall's tau between X and Y_star
test_result <- cor.test(df_sorted$X, Y_star, method = "kendall")
cat("Kendall's tau:", test_result$estimate, "\n")
cat("p-value:", test_result$p.value, "\n")
