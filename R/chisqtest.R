source('R/HybridPower.R')

hp_chisq <- R6Class(
  'hp_chisq',
  inherit = HybridPower,
  public = list(
    p_0 = NULL,
    p_1 = NULL,
    hybrid_powers = NULL,
    prior_mu = NULL,
    prior_sigma = NULL,
    prior_lower = NULL,
    prior_upper = NULL,
    prior_a = NULL,
    prior_b = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      p_0 = NULL,
      p_1 = NULL,
      alt = 'two.sided',
      prior_mu = NULL,
      prior_sigma = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      prior_a = NULL,
      prior_b = NULL,
      quantiles = NULL,
      assurance_level_props = NULL
    ) {
      if (!(is.null(prior))) {
        if (!(prior %in% c('dirichlet', 'beta', 'uniform', 'truncnorm')))
          stop('Invalid prior')
      }
      super$initialize(
        parallel = FALSE,
        ns=ns,
        n_prior=n_prior,
        n_MC=n_MC,
        prior=prior,
        prior_mu = prior_mu,
        prior_sigma = prior_sigma,
        prior_lower = prior_lower,
        prior_upper = prior_upper,
        prior_a = prior_a,
        prior_b = prior_b,
        alpha=alpha,
        alt=alt,
        quantiles=quantiles,
        assurance_level_props=assurance_level_props
      )
      if (sum(p_0) != 1)
        stop('p_0 must sum to 1')
      for (i in 1:length(p_0)) {
        if (!(is.numeric(p_0[i])))
          stop('Elements of p_0 must be numeric')
      }
      if (!(is.null(p_1))) {
        if (sum(p_1) != 1)
          stop('p_1 must sum to 1')
        if (length(p_0) != length(p_1))
          stop('p_0 and p_1 should have same lengths')
        for (i in 1:length(p_1)) {
          if (!(is.numeric(p_1[i])))
            stop('Elements of p_1 must be numeric')
        }
      }
      self$p_0 <- p_0
      self$p_1 <- p_1
    },

    print = function() {
      super$print()
      cat('Proportions under H_0: ', self$p_0, '\n')
      cat('Proportions under H_1: ', self$p_1, '\n\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'dirichlet') {
          cat('Prior alpha: ', self$prior_a, '\n')
        }
        else if (self$prior == 'beta') {
          cat('Prior alpha: ', self$prior_a, '\n')
          cat('Prior beta: ', self$prior_b, '\n\n')
        }
        else if (self$prior == 'normal') {
          cat('Prior mean: ', self$prior_mu, '\n')
          cat('Prior sd: ', self$prior_sigma, '\n\n')
        }
        else if (self$prior == 'uniform') {
          cat('Prior lower bound: ', self$prior_lower, '\n')
          cat('Prior upper bound: ', self$prior_upper, '\n\n')
        }
      }
      cat('Test type: Chi^2\n')
    },

    classical_power = function(n=self$ns, p_1=self$p_1) {
      return(1 - pchisq(qchisq(1-self$alpha, df=length(p_1)-1), df=length(p_1)-1, ncp=sum((p_1 - self$p_0)^2/self$p_0)*n))
    }
  ),

  private = list(
    draw_prior_es = function() {
      if (self$prior == 'beta') {
        es <- rbeta(self$n_prior, self$prior_a, self$prior_b)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'truncnorm') {
        es <- truncnorm::rtruncnorm(self$n_prior, mean=self$prior_mu, sd=self$prior_sigma, a=0, b=1)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'uniform') {
        es <- runif(self$n_prior, self$prior_lower, self$prior_upper)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'dirichlet') {
        es <- gtools::rdirichlet(self$n_prior, self$prior_a)
        return(es)
      }
    },

    generate_hybrid_power = function(n) {
      return(apply(private$draw_prior_es(), 1, FUN=self$classical_power, n=n))
    }
  )
)
