# ---------------------------------------------
# Simulation Study for T_{n,2} U-Statistic (Degeneracy Order 2)
# ----------------------------
# Load Required Base Libraries
# ----------------------------
library(utils)     # for combn()
library(stats)     # for rnorm(), density(), qqnorm()
library(graphics)  # for hist(), lines(), qqline()

# ----------------------------
# Define the Test Statistic T_{n,2}
# ----------------------------
Tn2 <- function(x, y) {
  n <- length(x)
  stopifnot(length(y) == n)
  
  # Generate all 4-element combinations
  idx <- utils::combn(n, 4)
  
  # Apply the kernel function
  a_xy <- apply(idx, 2, function(i) {
    xv <- x[i]; yv <- y[i]
    
    ax <- sign(abs(xv[1] - xv[2]) + abs(xv[3] - xv[4]) -
                 abs(xv[1] - xv[3]) - abs(xv[2] - xv[4]))
    ay <- sign(abs(yv[1] - yv[2]) + abs(yv[3] - yv[4]) -
                 abs(yv[1] - yv[3]) - abs(yv[2] - yv[4]))
    
    ax * ay
  })
  
  # Return mean of the kernel evaluations
  mean(a_xy, na.rm = TRUE)
}

# ----------------------------
# Simulation Settings
# ----------------------------
set.seed(2025)
n <- 100              # Sample size
B <- 2000             # Number of bootstrap replicates
nrep <- 100           # Number of simulation runs

Tn2_values <- numeric(nrep)
pvals <- numeric(nrep)

# ----------------------------
# Simulation Loop (under H0)
# ----------------------------
for (rep in 1:nrep) {
  # Generate data under null: x and y independent
  x <- rnorm(n)
  y <- rnorm(n)
  
  # Compute observed statistic
  stat.obs <- Tn2(x, y)
  
  # Bootstrap null distribution by permuting y
  stat.boot <- replicate(B, {
    y_perm <- sample(y)
    Tn2(x, y_perm)
  })
  
  # Bootstrap p-value
  pval <- mean(stat.boot >= stat.obs)
  
  # Store results
  Tn2_values[rep] <- stat.obs
  pvals[rep] <- pval
}

# ----------------------------
# Results Summary
# ----------------------------
cat("Average bootstrap p-value across", nrep, "runs:", round(mean(pvals), 4), "\n")

# ----------------------------
# Histogram of n * T_{n,2}
# ----------------------------
hist(n * Tn2_values,
     breaks = 30,
     col = "skyblue",
     main = expression(paste("Histogram of ", n, "·T"[n,2], " under H0")),
     xlab = expression(n * T[n,2]),
     probability = TRUE)
lines(density(n * Tn2_values), col = "darkblue", lwd = 2)

# ----------------------------
# QQ Plot to Assess Normality (Expected to Deviate)
# ----------------------------
qqnorm(n * Tn2_values,
       main = expression("QQ Plot of " ~ n * T[n,2]))
qqline(n * Tn2_values, col = "red", lwd = 2)

# ----------------------------
# End of Script
# ----------------------------
