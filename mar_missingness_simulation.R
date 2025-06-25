
# Step 1: Simulate X
set.seed(123)
n <- 1000
X <- rnorm(n)

# Step 2: Set desired missing rate and fix gamma1
target_missing_rate <- 0.3
gamma1 <- 1

# Step 3: Function to solve for gamma0
missing_rate_diff <- function(gamma0) {
  p_miss <- exp(gamma0 + gamma1 * X) / (1 + exp(gamma0 + gamma1 * X))
  mean(p_miss) - target_missing_rate
}

# Step 4: Solve for gamma0
gamma0 <- uniroot(missing_rate_diff, lower = -10, upper = 10)$root

# Step 5: Compute missingness probabilities
p_miss <- exp(gamma0 + gamma1 * X) / (1 + exp(gamma0 + gamma1 * X))

# Step 6: Generate delta (1 = observed, 0 = missing)
delta <- rbinom(n, 1, prob = 1 - p_miss)

# Step 7: Confirm actual missing rate
actual_missing_rate <- mean(delta == 0)

# Output results
cat("Estimated gamma0:", round(gamma0, 4), "\n")
cat("Fixed gamma1:", gamma1, "\n")
cat("Target missing rate:", target_missing_rate, "\n")
cat("Achieved missing rate:", round(actual_missing_rate, 4), "\n")
