source('R/HybridPower.R')

hp_sign <- R6Class(
  'hp_sign',
  inherit = HybridPower,
  public = list(
    p_0 = NULL,
    p_1 = NULL,
    MC = F,
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
      p_0 = 0.5,
      p_1 = NULL,
      alt = 'two.sided',
      MC=F,
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
        if (!(prior %in% c('truncnorm','beta','uniform')))
          stop('Prior must be one of truncnorm, beta or uniform')
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
      if (!(is.numeric(p_0)) | p_0 > 1 | p_0 < 0)
        stop('Invalid p_0')
      if (!(is.null(p_1))) {
        if (!(is.numeric(p_1)) | p_1 > 1 | p_1 < 0)
          stop('Invalid p_1')
      }

      self$p_0 <- p_0
      self$p_1 <- p_1
      self$MC <- MC
    },

    print = function() {
      super$print()
      cat('p under H_0: ', self$p_0, '\n')
      if (!(is.null(self$p_1)))
        cat('p under H_1: ', self$p_1, '\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'beta') {
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
      cat('Test type: Sign Test\n')
    },

    classical_power = function(n=self$ns, p_0=self$p_0, p_1=self$p_1) {
      if (is.null(p_1))
        stop('Provide a value of p_1')
      if (self$MC) {
        return(
          mean(
            sapply(rbinom(self$n_MC, n, p_1), FUN=private$binom_test, n=n, p_0=p_0)
          )
        )
      }
      else {
        mu_0 <- n*p_0
        mu_1 <- n*p_1
        sigma_0 <- sqrt(mu_0*(1-p_0))
        sigma_1 <- sqrt(mu_1*(1-p_1))
        rho <- sigma_1 / sigma_0
        if (self$alt == 'two.sided') {
          z_lower <- ((mu_0 - mu_1) / sigma_0 + qnorm(self$alpha/2))/rho
          z_upper <- ((mu_0 - mu_1) / sigma_0 + qnorm(1-self$alpha/2))/rho
          return(pnorm(z_lower) + 1 - pnorm(z_upper))
        }
        else if (self$alt == 'one.sided') {
          return(1-pnorm(qnorm(1-self$alpha), abs(mu_0 - mu_1) / sigma_0, 1/rho))
        }
      }
    },

    hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        if (!(cores)) cores <- detectCores()
        self$output <- parallel::mclapply(self$ns, private$generate_hybrid_power)
        private$melt_output()
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$generate_hybrid_power(self$ns[i])
        }
        self$output <- res
        private$melt_output()
      }
      return(self$output)
    }
  ),

  private = list(
    binom_test = function(x, n, p_0) {
      if (self$alt == 'two.sided') {
        crit_lower <- qbinom(self$alpha/2, n, p_0)-1
        crit_upper <- qbinom(1-self$alpha/2, n, p_0)
        return((x <= crit_lower) | (x >= crit_upper))
      }
      else if (self$alt == 'one.sided') {
        return(x >= qbinom(1-self$alpha/2, n, p_0))
      }
    },

    draw_prior_es = function() {
      if (self$prior == 'beta') {
        return(
          rbeta(self$n_prior, self$prior_a, self$prior_b)
        )
      }
      else if (self$prior == 'truncnorm') {
        return(
          truncnorm::rtruncnorm(self$n_prior, mean=self$prior_mu, sd=self$prior_sigma, a=0, b=1)
        )
      }
      else if (self$prior == 'uniform') {
        return(
          runif(self$n_prior, self$prior_lower, self$prior_upper)
        )
      }
    },

    generate_hybrid_power = function(n) {
      es <- private$draw_prior_es()
      return(
        sapply(es, FUN=self$classical_power, n=n, p_0=self$p_0)
      )
    }
  )
)
