source('R/HybridPower.R')

HybridPowerSLR <- R6Class(
  'HybridPowerSLR',
  inherit = HybridPower,
  public = list(
    r2 = NULL,
    hybrid_powers = NULL,
    prior_alpha = NULL,
    prior_beta = NULL,
    prior_upper = NULL,
    prior_lower = NULL,
    prior_mu = NULL,
    prior_sd = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      alt = 'two.sided',
      r2 = NULL,
      prior_alpha = NULL,
      prior_beta = NULL,
      prior_upper = NULL,
      prior_lower = NULL,
      prior_mu = NULL,
      prior_sd = NULL,
      assurance_props = NULL
    ) {
      super$initialize(
        parallel = FALSE,
        ns=ns,
        n_prior=n_prior,
        n_MC=n_MC,
        prior=prior,
        alpha=alpha,
        alt=alt,
        assurance_props=assurance_props
      )
      if (!is.null(r2)) {
        if (!(is.numeric(r2)) | r2 > 1 | r2 < 0)
          stop('Invalid r2')
        else
          self$r2 <- r2
      }
      if (!(is.null(prior))) {
        if (prior == 'beta') {
          if (!(is_numeric(prior_pi_1_alpha) &
                is_numeric(prior_pi_1_beta)))
            stop('Invalid input type for the priors')
          if (prior_alpha <= 0)
            stop('prior_alpha should be positive')
          if (prior_beta <= 0)
            stop('prior_beta should be positive')
          self$prior_alpha <- prior_alpha
          self$prior_beta <- prior_beta
        }
        else if (prior == 'uniform') {
          if ((prior_pi_1_lower > prior_pi_1_upper) |
               prior_pi_1_lower < 0 | prior_pi_1_upper > 1
              )
            stop('Invalid limits for the uniform prior(s)')
          self$prior_lower <- prior_lower
          self$prior_upper <- prior_upper
        }
        else if (prior == 'truncnorm') {
          if (!(is_numeric(prior_mu) &
                is_numeric(prior_sd))
          )
            stop('Invalid input type for the priors')
          if (prior_mu <= 0 | prior_mu >= 1)
            stop('Prior means must be between 0 and 1')
          if (prior_sd <= 0)
            stop('Prior standard deviations must be positive')
          self$prior_mu <- prior_mu
          self$prior_sd <- prior_sd
        }
      }
    },

    print = function() {
      super$print()
      cat('R-squared under H_1: ', self$r2, '\n')
      if (self$prior == 'beta') {
        cat('Alpha: ', self$prior_alpha, '\n')
        cat('Beta: ', self$prior_beta, '\n')
      }
      else if (self$prior == 'uniform') {
        cat('Lower (r2): ', self$prior_lower, '\n')
        cat('Upper (r2): ', self$prior_upper, '\n')
      }
      else if (self$prior == 'truncnorm') {
        cat('Prior mean (r2): ', self$prior_mu, '\n')
        cat('Prior sd (r2): ', self$prior_sd, '\n')
      }
      cat('Study design: Simple Regression\n')
    },

    classical_power = function(n=self$ns, r2=self$r2) {
      return(1-pf(qf(1-self$alpha, df1=1, df2=n-2), df1=1, df2=n-2, ncp=r2/(1-r2)*n))
    }
  ),

  private = list(
    draw_prior_es = function() {
      if (self$prior == 'beta') {
        return(rbeta(self$n_prior, self$prior_alpha, self$prior_beta))
      }
      else if (self$prior == 'uniform') {
        return(runif(self$n_prior, self$prior_lower, self$prior_upper))
      }
      else if (self$prior == 'truncnorm') {
        return(
          truncnorm::rtruncnorm(
            n=self$n_prior,
            a=0,
            b=1,
            mean=self$prior_mu,
            sd=self$prior_sd
          )
        )
      }
    },

    generate_hybrid_power = function(n) {
      return(self$classical_power(n, r2=private$draw_prior_es()))
    }
  )
)
