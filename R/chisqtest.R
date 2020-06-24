source('R/HybridPower.R')

HybridPowerChisqTest <- R6Class(
  'HybridPowerChisqTest',
  inherit = HybridPower,
  public = list(
    p_0 = NULL,
    p_1 = NULL,
    hybrid_powers = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
    prior_lower = NULL,
    prior_upper = NULL,
    prior_alpha = NULL,
    prior_beta = NULL,

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
      prior_sd = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      prior_alpha = NULL,
      prior_beta = NULL
    ) {
      super$initialize(
        parallel = FALSE,
        ns,
        n_prior,
        n_MC,
        prior,
        alpha,
        alt
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

      if (!(is.null(prior))) {
        if (length(prior_alpha) > 1) {
          if (prior != 'dirichlet')
            stop('The prior type must be dirichlet for (# cells) > 2')
          if (length(prior_alpha) != length(p_0))
            stop('Lengths of prior_alpha and p_0 should match')
          for (i in 1:length(prior_alpha)) {
            if (!(is_numeric(prior_alpha[i])))
              stop('Invalid input type for the priors')
            if (length(prior_alpha[i]) != 1)
              stop('prior_alpha should be a scalar')
            if (prior_alpha[i] <= 0)
              stop('Elements of prior_alpha must be positive')
          }
          self$prior_alpha <- prior_alpha
        }
        else {
          if (prior == 'dirichlet')
            stop('Cannot use a dirichlet prior in this case')
          if (prior == 'beta') {
            if (is.null(prior_alpha) | is.null(prior_beta))
              stop('Please specify prior parameters')
            if (!(is_numeric(prior_alpha) & is_numeric(prior_beta)))
              stop('Invalid input type for the priors')
            if (length(prior_alpha) != 1)
              stop('prior_alpha should be a scalar')
            if (length(prior_beta) != 1)
              stop('prior_beta should be a scalar')
            if (prior_alpha <= 0)
              stop('prior_alpha should be positive')
            if (prior_beta <= 0)
              stop('prior_beta should be positive')
            self$prior_alpha <- prior_alpha
            self$prior_beta <- prior_beta
          }
          if (prior == 'truncnorm') {
            if (is.null(prior_mu) | is.null(prior_sd))
              stop('Please specify prior parameters')
            if (!(is_numeric(prior_mu) & is_numeric(prior_sd)))
              stop('Invalid input type for the priors')
            if (length(prior_mu) != 1)
              stop('prior_mu should be a scalar')
            if (length(prior_sd) != 1)
              stop('prior_sd should be a scalar')
            if (prior_mu <= 0 | prior_mu >= 1)
              stop('prior_mu should be between 0 and 1')
            if (prior_sd <= 0)
              stop('prior sd should be positive')
            self$prior_mu <- prior_mu
            self$prior_sd <- prior_sd
          }
          else if (prior == 'uniform') {
            if (is.null(prior_lower) | is.null(prior_upper))
              stop('Please specify prior parameters')
            if (prior_lower < 0 | prior_lower > 1)
              stop('The lower bound should be between 0 and 1')
            if (prior_upper < 0 | prior_upper > 1)
              stop('The upper bound should be between 0 and 1')
            if (prior_lower >= prior_upper)
              stop('The lower bound should be smaller than the upper bound')
            self$prior_lower <- prior_lower
            self$prior_upper <- prior_upper
          }
        }
      }
    },

    print = function() {
      super$print()
      cat('Proportions under H_0: ', self$p_0, '\n')
      cat('Proportions under H_1: ', self$p_1, '\n\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'dirichlet') {
          cat('Prior alpha: ', self$prior_alpha, '\n')
        }
        else if (self$prior == 'beta') {
          cat('Prior alpha: ', self$prior_alpha, '\n')
          cat('Prior beta: ', self$prior_beta, '\n\n')
        }
        else if (self$prior == 'normal') {
          cat('Prior mean: ', self$prior_mu, '\n')
          cat('Prior sd: ', self$prior_sd, '\n\n')
        }
        else if (self$prior == 'uniform') {
          cat('Prior lower bound: ', self$prior_lower, '\n')
          cat('Prior upper bound: ', self$prior_upper, '\n\n')
        }
      }
      cat('Test type: Chi^2\n')
    },

    classical_power = function(n=self$ns, p=self$p_1) {
      return(1 - pchisq(qchisq(1-self$alpha, df=length(p)-1), df=length(p)-1, ncp=sum((p - self$p_0)^2/self$p_0)*n))
    },

    hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        self$output <- mclapply(self$ns, private$generate_hybrid_power)
        private$melt_output()
        return(self$output)
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$generate_hybrid_power(self$ns[i])
        }
        self$output <- res
        private$melt_output()
        return(self$output)
      }
    }
  ),

  private = list(
    draw_prior_es = function() {
      if (self$prior == 'beta') {
        es <- rbeta(self$n_prior, self$prior_alpha, self$prior_beta)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'truncnorm') {
        es <- truncnorm::rtruncnorm(self$n_prior, mean=self$prior_mu, sd=self$prior_sd, a=0, b=1)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'uniform') {
        es <- runif(self$n_prior, self$prior_lower, self$prior_upper)
        return(cbind(es, 1-es))
      }
      else if (self$prior == 'dirichlet') {
        es <- gtools::rdirichlet(self$n_prior, self$prior_alpha)
        return(es)
      }
    },

    generate_hybrid_power = function(n) {
      ps <- private$draw_prior_es()
      return(apply(ps, 1, FUN=self$classical_power, n=n))
    }
  )
)

# Classical power analysis
x_classical <- HybridPowerChisqTest$new(
  ns = seq(10, 90, 10),
  p_0 = c(1/3, 1/3, 1/3),
  p_1 = c(1/4, 1/4, 1/2)
)
x_classical$classical_power()

# Hybrid power analysis
x_hybrid <- HybridPowerChisqTest$new(
  prior='beta',
  parallel = T,
  ns = seq(10, 90, 10),
  n_prior=1000,
  prior_alpha = 1,
  prior_beta = 2,
  p_0 = c(1/2, 1/2),
  p_1 = c(1/3, 2/3)
)

x_hybrid$classical_power()
x_hybrid$hybrid_power()
x_hybrid$assurance()
x_hybrid$boxplot()
