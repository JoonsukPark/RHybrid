source('R/ttest.R')

ns <- seq(10, 90, 10)
x <- HybridPowerTtest$new(
  parallel=T,
  ns=ns,
  n_prior = 10,
  n_MC = 10,
  alt='two.sided',
  alpha=.05,
  prior = 'uniform',
  prior_lower = 0.1,
  prior_upper = 0.2
)

x$assurance()
x$generate_hybrid_power()
x$generate_hybrid_power()$plot_power()
