# Required libraries
library(combinat)   # for combn() permutations
library(magrittr)   # for pipe `%>%` if preferred

# Define Tn2: test statistic for T_{n,2}
Tn2 <- function(x, y) {
  n <- length(x)
  stopifnot(length(y) == n)
  
  idx <- utils::combn(n, 4)
  a_xy <- apply(idx, 2, function(i) {
    xv <- x[i]; yv <- y[i]
    ax <- sign(abs(xv[1] - xv[2]) + abs(xv[3] - xv[4]) -
                 abs(xv[1] - xv[3]) - abs(xv[2] - xv[4]))
    ay <- sign(abs(yv[1] - yv[2]) + abs(yv[3] - yv[4]) -
                 abs(yv[1] - yv[3]) - abs(yv[2] - yv[4]))
    ax * ay
  })
  
  mean(a_xy, na.rm = TRUE)
}

# Example data - replace with your observed sample
set.seed(123)
n <- 100
x <- rnorm(n)
y <- rnorm(n)

# Compute observed test statistic
stat.obs <- Tn2(x, y)

# Bootstrap (permutation) calibration
B <- 2000
stat.boot <- replicate(B, {
  y_perm <- sample(y)
  Tn2(x, y_perm)
})

# 95th percentile critical value
crit.95 <- quantile(stat.boot, 0.95)

# Empirical p-value
p.value <- mean(stat.boot >= stat.obs)

# Display results
cat("Observed Tn2:         ", round(stat.obs, 4), "\n",
    "Bootstrap 95% cutoff: ", round(crit.95, 4), "\n",
    "Bootstrap p-value:    ", round(p.value, 4), "\n")
