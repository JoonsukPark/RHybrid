source('R/HybridPower.R')

HybridPowerCorrelation <- R6Class(
  'HybridPowerCorrelation',
  inherit = HybridPower,
  public = list(
    hybrid_powers = NULL,
    prior = NULL,
    prior_a = NULL,
    prior_b = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
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
      prior_mu = NULL,
      prior_sd = NULL,
      rho = NULL,
      assurance_props=NULL
    ) {
      super$initialize(
        parallel = parallel,
        ns = ns,
        n_prior = n_prior,
        n_MC = n_MC,
        prior = prior,
        alpha = alpha,
        alt = alt,
        assurance_props = assurance_props
      )
      if (!(is.null(prior))) {
        if (prior == 'beta') {
          if (is.null(prior_a) | is.null(prior_a))
            stop('Provide prior_a and prior_b first')
          if (prior_a <= 0 | prior_b <= 0)
            stop('Parameters for the beta prior must be positive')
          self$prior_a <- prior_a
          self$prior_b <- prior_b
        }
        else if (prior == 'normal' | prior == 'truncnorm') {
          if (is.null(prior_mu) | is.null(prior_sd))
            stop('Provide prior_mu and prior_sd first')
          if (prior_sd <= 0)
            stop('Prior standard deviation must be positive')
          if (prior_mu <= -1 | prior_mu >= 1)
            stop('Prior mean of correlation coefficients must be between -1 and 1')
          self$prior_mu <- prior_mu
          self$prior_sd <- prior_sd
        }
      }
      if (!(is.null(rho))) {
        if (!(is.numeric(rho)))
          stop('rho must be numeric')
        if (rho > 1 | rho < 0)
          stop('rho must be between 0 and 1')
        self$rho <- rho
      }
    },

    print = function() {
      super$print()
      cat('Test type: Bivariate correlation\n')
      cat('Point effect size: ', self$rho, '\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'beta') {
          cat('Prior_a: ', self$prior_a, '\n')
          cat('Prior_b: ', self$prior_b, '\n')
        }
        else if (self$prior == 'normal') {
          cat('Prior_mu: ', self$prior_mu, '\n')
          cat('Prior_sd: ', self$prior_sd, '\n')
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
      else {
        return(
          truncnorm::rtruncnorm(
            self$n_prior,
            mean = self$prior_mu,
            sd = self$prior_sd,
            a = -1,
            b = 1
          )
        )
      }
    },

    generate_hybrid_power = function(n) {
      rhos <- private$draw_prior_es()
      return(
        sapply(rhos, FUN=self$classical_power, n=n)
      )
    }
  )
)
