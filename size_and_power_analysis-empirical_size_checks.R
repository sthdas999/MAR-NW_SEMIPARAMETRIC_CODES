# empirical_size.R -- compute empirical size for six tests via simulation

library(MonteCarlo)  # for simulation grid and table

# 1. Define simulation functions for each test under H0
#    Replace `test_stat_*` with your actual test implementations.
t1 <- function(n) { test_stat_1 <- rnorm(n); decision <- my_test1(test_stat_1); list(decision = decision) }
t2 <- function(n) { test_stat_2 <- rnorm(n); decision <- my_test2(test_stat_2); list(decision = decision) }
t3 <- function(n) { test_stat_3 <- rnorm(n); decision <- my_test3(test_stat_3); list(decision = decision) }
t4 <- function(n) { test_stat_4 <- rnorm(n); decision <- my_test4(test_stat_4); list(decision = decision) }
t5 <- function(n) { test_stat_5 <- rnorm(n); decision <- my_test5(test_stat_5); list(decision = decision) }
t6 <- function(n) { test_stat_6 <- rnorm(n); decision <- my_test6(test_stat_6); list(decision = decision) }

# 2. Combine into a single simulation function
sim_null <- function(n) {
  c1 <- t1(n)$decision
  c2 <- t2(n)$decision
  c3 <- t3(n)$decision
  c4 <- t4(n)$decision
  c5 <- t5(n)$decision
  c6 <- t6(n)$decision
  list(Test1 = c1, Test2 = c2, Test3 = c3,
       Test4 = c4, Test5 = c5, Test6 = c6)
}

# 3. Define grid (just sample size here)
param_list <- list(n = 20)  # or multiple n's if you want

# 4. Run Monte Carlo simulation
set.seed(123)
MC0 <- MonteCarlo(func = sim_null, nrep = 5000, param_list = param_list)

# 5. Summarize: proportion of rejections per test
emp_sizes <- sapply(summary(MC0)$decision, function(x) x["mean"]) * 100
names(emp_sizes) <- paste0("Test", 1:6)

# 6. Optionally, recalibrate critical values if outside [4.5, 5.5]
recalibrate <- emp_sizes < 4.5 | emp_sizes > 5.5
recalibrated <- emp_sizes  # placeholder

if (any(recalibrate)) {
  # e.g., find new critical value via simulated null quantile
  # stats_null <- replicate(5000, test_stat_under_H0())
  # critical_new <- quantile(stats_null, 0.95)
  # Use critical_new in rerunning power or size
}

# 7. Print LaTeX table
library(xtable)
tab <- data.frame(
  Test = names(emp_sizes),
  Nominal = rep(5.00, 6),
  Observed = round(emp_sizes, 2)
)
print(xtable(tab, caption = "Empirical Type?I error (%) under the null (5,000 sims)"),
      include.rownames = FALSE)
