library(R6)
library(ggplot2)

HybridPower <- R6Class(
  'HybridPower',
  public = list(
    parallel = FALSE,
    ns = c(),
    n_prior = NULL,
    n_MC = NULL,
    prior = NULL,
    alpha = NULL,
    alt = NULL,
    prior_mu=NULL,
    prior_sd=NULL,
    prior_lower=NULL,
    prior_upper=NULL,
    powers = c(),
    
    initialize = function(
      parallel = FALSE,
      ns = c(),
      n_prior = 1,
      n_MC = 1,
      prior = 'normal',
      alpha = 0.05,
      alt = 'two.sided',
      prior_mu=NULL,
      prior_sd=NULL,
      prior_lower=NULL,
      prior_upper=NULL
    ) {
      # Validate inputs
      if (length(ns) == 0) stop('Input sample sizes!')
      else {
        for (i in 1:length(ns)) {
          if (!is.numeric(ns[i]) | ns[i] %% 1 != 0 | ns[i] <= 0)
            stop('Invalid sample size!')
        }
      }
      if (n_prior %% 1 != 0 | n_prior <= 0)
        stop('Invalid # draws from prior!')
      if (n_MC %% 1 != 0 | n_MC <= 0)
        stop('Invalid # Monte Carlo simulations!')
      if (
        !is.character(prior) |
        length(prior) != 1 |
        !(prior == 'normal' | prior == 'uniform')
      )
        stop('Invalid prior!')
      if (alpha <= 0 | alpha >= 1)
        stop('Alpha should be between 0 and 1!')
      if (alt != 'one.sided' & alt != 'two.sided')
        stop('Alternative hypothesis should be either \'one.sided\' or \'two.sided\'!')
      
      self$parallel <- parallel
      self$ns <- ns
      self$n_prior <- n_prior
      self$n_MC <- n_MC
      self$prior <- prior
      self$alpha <- alpha
      self$alt <- alt
      if (prior == 'normal') {
        self$prior_mu <- prior_mu
        self$prior_sd <- prior_sd
      }
      else if (prior == 'uniform') {
        self$prior_lower <- prior_lower
        self$prior_upper <- prior_upper
      }
    },
    
    print = function(){
      cat('HybridPower Instance Description: \n\n')
      cat('Parallelize: ', self$parallel, '\n\n')
      cat('Sample sizes: ', self$ns, '\n')
      cat('Draws from prior: ', self$n_prior, '\n')
      cat('# Monte Carlo simulations: ', self$n_MC, '\n\n')
      cat('Type of prior: ', self$prior, '\n')
      if (self$prior == 'normal') {
        cat('Prior mean: ', self$prior_mu, '\n')
        cat('Prior sd: ', self$prior_sd, '\n\n')
      }
      else if (self$prior == 'uniform') {
        cat('Prior lower bound: ', self$prior_lower, '\n')
        cat('Prior upper bound: ', self$prior_upper, '\n\n')
      }
      cat('Alternative Hypothesis: ', self$alt, '\n')
      cat('Level of significance: ', self$alpha, '\n')
    }
  )
)

plot_powers <- function(powers, n_bins=30) {
  hist(
    powers,
    prob=T,
    breaks = n_bins,
    main='Distribution of power',
    xlab='Power',
    ylab='Proportions'
  )
}

# HybridPower$new(
#   prior_mu=0.5,
#   prior_sd=0.1,
#   ns=c(10,20,30),
#   parallel=T
# )
