source('R/ttest.R')

ns <- seq(10, 90, 10)
x <- HybridPowerTtest$new(
  parallel=T,
  ns=ns,
  n_prior = 1000,
  n_MC = 1000,
  alt='two.sided',
  alpha=.05,
  prior = 'uniform',
  prior_lower = 0.1,
  prior_upper = 0.2
)

x$assurance
x$hybrid_powers
x$plot_power(x$hybrid_powers)
