# Tn3_simulation.R
# Simulation study for degenerate U-statistic T_{n,3}

# Load necessary libraries
library(stats)
library(graphics)
library(utils)

set.seed(123)
n_rep <- 100
n <- 100
B <- 2000

Tn3 <- function(x, y) {
  n <- length(x)
  combs <- combn(n, 4)
  stat <- 0
  for (i in 1:ncol(combs)) {
    idx <- combs[, i]
    xi <- x[idx[1]]
    xj <- x[idx[2]]
    yk <- y[idx[3]]
    yl <- y[idx[4]]
    stat <- stat + sign((xi - xj) * (yk - yl))
  }
  return(stat / ncol(combs))
}

pvals <- numeric(n_rep)
Tn3_vals <- numeric(n_rep)

for (r in 1:n_rep) {
  x <- rnorm(n)
  y <- rnorm(n)
  stat.obs <- Tn3(x, y)
  stat.boot <- replicate(B, Tn3(x, sample(y)))
  pvals[r] <- mean(stat.boot >= stat.obs)
  Tn3_vals[r] <- n * stat.obs
}

avg_pval <- mean(pvals)
cat("Average P-value under H0:", avg_pval, "\\n")

# Save plots
png("hist_Tn3.png", width=600, height=400)
hist(Tn3_vals, breaks=20, col='lightgreen', main="Histogram of n*T_{n,3} under H0",
     xlab="n * T_{n,3}")
dev.off()

png("qq_Tn3.png", width=600, height=400)
qqnorm(Tn3_vals)
qqline(Tn3_vals)
title("QQ Plot of n*T_{n,3} under H0")
dev.off()
