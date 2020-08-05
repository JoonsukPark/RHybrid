source('R/HybridPower.R')

hp_cor <- R6Class(
  'hp_cor',
  inherit = HybridPower,
  public = list(
    rho = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      alpha = 0.05,
      alt = 'one.sided',
      prior = NULL,
      prior_a = NULL,
      prior_b = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      prior_mu = NULL,
      prior_sigma = NULL,
      rho = NULL,
      quantiles=NULL,
      assurance_level_props=NULL
    ) {
      if (!(is.null(prior))) {
        if (prior == 'truncnorm' | prior == 'uniform') {
          if (!(is.null(prior_lower))) {
            if (prior_lower < 0) {
              stop('Lower bound of rho cannot be less than 0')
            }
          }
          else {
            prior_lower <- 0
          }
          if (!(is.null(prior_upper))) {
            if (prior_upper > 1) {
              stop('Upper bound of rho cannot be greater than 1')
            }
          }
          else {
            prior_upper <- 1
          }
        }
        if (prior == 'truncnorm') {
          if (abs(prior_mu) > 1) {
            stop('The absolute value of prior_mu cannot be greater than 1')
          }
          if (abs(prior_sigma) <= 0) {
            stop('prior_sigma cannot be negative')
          }
        }
        else if (prior == 'beta') {}
        else {
          stop('Invalid prior')
        }
      }
      super$initialize(
        parallel = parallel,
        ns = ns,
        n_prior = n_prior,
        n_MC = n_MC,
        prior = prior,
        prior_a = prior_a,
        prior_b = prior_b,
        prior_lower = prior_lower,
        prior_upper = prior_upper,
        prior_mu = prior_mu,
        prior_sigma = prior_sigma,
        alpha = alpha,
        alt = alt,
        quantiles = quantiles,
        assurance_level_props=assurance_level_props
      )
      if (!(is.null(rho))) {
        if (!(is.numeric(rho)))
          stop('rho must be numeric')
        if (rho > 1 | rho < -1)
          stop('rho must be between -1 and 1')
        self$rho <- abs(rho)
      }
    },

    print = function() {
      super$print()
      cat('Test type: Bivariate correlation\n')
      cat('Point effect size: ', self$rho, '\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'beta') {
          cat('Prior a: ', self$prior_a, '\n')
          cat('Prior b: ', self$prior_b, '\n')
        }
        else if (self$prior == 'truncnorm') {
          cat('Prior mu: ', self$prior_mu, '\n')
          cat('Prior sigma: ', self$prior_sigma, '\n')
        }
        else {
          cat('Prior lower: ', self$prior_lower, '\n')
          cat('Prior upper: ', self$prior_upper, '\n')
        }
      }
    },

    classical_power = function(rho=self$rho, n=self$ns) {
      if (self$alt == 'one.sided')
        rho = abs(rho)
      crit_t <- ifelse(
        self$alt == 'two.sided',
        qt(1-self$alpha/2, df=n-2),
        qt(1-self$alpha, df=n-2)
      )
      crit_r <- crit_t / sqrt(crit_t^2 + n - 2)
      r_z <- atanh(rho) + rho/(2*n-2)
      crit_r_z <- atanh(crit_r) + crit_r/(2*n-2)
      if (self$alt == 'two.sided') {
        return(
          pnorm(-crit_r_z, r_z, 1/sqrt(n-3)) +
            1 - pnorm(crit_r_z, r_z, 1/sqrt(n-3))
        )
      }
      else {
        return(
          1 - pnorm(crit_r_z, r_z, 1/sqrt(n-3))
        )
      }
    }
  ),

  private = list(
    draw_prior_es = function() {
      if (self$prior == 'beta') {
        return(
          rbeta(self$n_prior, self$prior_a, self$prior_b)
        )
      }
      else if (self$prior == 'truncnorm') {
        return(
          truncnorm::rtruncnorm(
            self$n_prior,
            mean = self$prior_mu,
            sd = self$prior_sigma,
            a = self$prior_lower,
            b = self$prior_upper
          )
        )
      }
      else {
        return(runif(self$n_prior, self$lower, self$upper))
      }
    },

    generate_hybrid_power = function(n, es) {
      return(
        sapply(es, FUN=self$classical_power, n=n)
      )
    }
  )
)
