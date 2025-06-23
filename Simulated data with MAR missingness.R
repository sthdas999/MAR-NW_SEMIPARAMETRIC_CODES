# -------------------------------
# Simulated data with MAR missingness
# Nadaraya-Watson estimate at x = 0.5
# -------------------------------

# No external libraries needed for this basic task

# Set seed for reproducibility
set.seed(123)

# Simulate data
X <- seq(0.1, 1.0, by = 0.1)
epsilon <- rnorm(length(X), mean = 0, sd = 0.1)
Y <- sin(2 * pi * X) + epsilon

# Induce MAR missingness
prob_obs <- 0.7 + 0.2 * X
delta <- rbinom(length(X), size = 1, prob = prob_obs)

# Combine into data frame and print
sim_data <- data.frame(X = round(X, 2), Y = round(Y, 3), delta)
print(sim_data)

# Epanechnikov kernel function
K <- function(u) {
  w <- 0.75 * (1 - u^2)
  w[abs(u) > 1] <- 0
  return(w)
}

# Set bandwidth (from CV result in your example)
h <- 0.2

# Target x value
x0 <- 0.5

# Compute kernel weights (only for observed data)
u <- (X - x0) / h
weights <- K(u) * delta  # zero weight if Y is missing

# Compute NW estimate
numerator <- sum(weights * Y)
denominator <- sum(weights)

m_hat <- numerator / denominator

# Display result
cat(sprintf("Estimated m(%.1f): %.3f\n", x0, m_hat))
