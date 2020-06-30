source('R/HybridPower.R')

hp_slr <- R6Class(
  'hp_slr',
  inherit = HybridPower,
  public = list(
    r2 = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      alt = 'two.sided',
      r2 = NULL,
      prior_a = NULL,
      prior_b = NULL,
      prior_upper = NULL,
      prior_lower = NULL,
      prior_mu = NULL,
      prior_sigma = NULL,
      quantiles = NULL,
      assurance_level_props=NULL
    ) {
      if (!(is.null(prior))) {
        if (!(prior %in% c('truncnorm','beta','uniform')))
          stop('Prior must be one of truncnorm, beta or uniform')
      }
      super$initialize(
        parallel = FALSE,
        ns=ns,
        n_prior=n_prior,
        n_MC=n_MC,
        prior=prior,
        prior_a = prior_a,
        prior_b = prior_b,
        prior_upper = prior_upper,
        prior_lower = prior_lower,
        prior_mu = prior_mu,
        prior_sigma = prior_sigma,
        alpha=alpha,
        alt=alt,
        quantiles=quantiles,
        assurance_level_props=assurance_level_props
      )
      if (!is.null(r2)) {
        if (!(is.numeric(r2)) | r2 > 1 | r2 < 0)
          stop('Invalid r2')
        else
          self$r2 <- r2
      }
    },

    print = function() {
      super$print()
      cat('R-squared under H_1: ', self$r2, '\n')
      if (self$prior == 'beta') {
        cat('Alpha: ', self$prior_a, '\n')
        cat('Beta: ', self$prior_b, '\n')
      }
      else if (self$prior == 'uniform') {
        cat('Lower (r2): ', self$prior_lower, '\n')
        cat('Upper (r2): ', self$prior_upper, '\n')
      }
      else if (self$prior == 'truncnorm') {
        cat('Prior mean (r2): ', self$prior_mu, '\n')
        cat('Prior sd (r2): ', self$prior_sigma, '\n')
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
        return(rbeta(self$n_prior, self$prior_a, self$prior_b))
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
            sd=self$prior_sigma
          )
        )
      }
    },

    generate_hybrid_power = function(n) {
      return(self$classical_power(n, r2=private$draw_prior_es()))
    }
  )
)
