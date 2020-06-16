library(R6)
library(ggplot2)

is_numeric <- function(x) return(is.numeric(x) & length(x) == 1)

HybridPower <- R6Class(
  'HybridPower',
  public = list(
    es = NULL,
    design = NULL,
    parallel = FALSE,
    ns = c(),
    n_prior = NULL,
    n_MC = NULL,
    prior = NULL,
    alpha = NULL,
    alt = NULL,
    powers = c(),
    hybrid_powers = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
    prior_lower = NULL,
    prior_upper = NULL,

    initialize = function(
      parallel = FALSE,
      ns = c(),
      n_prior = 1,
      n_MC = 1,
      prior = 'normal',
      alpha = 0.05,
      alt = 'two.sided'
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
        !(prior %in% c('normal', 'uniform', 'beta', 'truncnorm'))
      )
        stop('Invalid prior!')
      if (alpha <= 0 | alpha >= 1)
        stop('Alpha should be between 0 and 1.')
      if (alt != 'one.sided' & alt != 'two.sided')
        stop('Alternative hypothesis should be either \'one.sided\' or \'two.sided\'!')

      self$parallel <- parallel
      self$ns <- sort(ns)
      self$n_prior <- n_prior
      self$n_MC <- n_MC
      self$prior <- prior
      self$alpha <- alpha
      self$alt <- alt
    },

    print = function(){
      cat('HybridPower Instance Description: \n\n')
      cat('Parallelize: ', self$parallel, '\n\n')
      cat('Sample sizes: ', self$ns, '\n')
      cat('Draws from prior: ', self$n_prior, '\n')
      cat('# Monte Carlo simulations: ', self$n_MC, '\n\n')
      cat('Type of prior: ', self$prior, '\n')
      cat('Alternative Hypothesis: ', self$alt, '\n')
      cat('Level of significance: ', self$alpha, '\n')
    }
  )
)
