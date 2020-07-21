source('R/HybridPower.R')

hp_oneway_anova <- R6Class(
  'hp_oneway_anova',
  inherit = HybridPower,
  public = list(
    mu = NULL,
    es = NULL,
    sd = 1,
    k = NULL,
    design = NULL,
    rho = 0,
    epsilon = 1,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=10,
      n_MC=10,
      prior=NULL,
      alpha = 0.05,
      mu = NULL,
      prior_mu = NULL,
      prior_sigma = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      sd = 1,
      design = 'fe',
      rho = 0,
      epsilon = 1,
      alt = 'two.sided',
      quantiles = NULL,
      assurance_level_props=NULL
    ) {
      if (!(is.null(prior))) {
        if (!(prior %in% c('normal', 'uniform')))
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
        alpha=alpha,
        alt=alt,
        quantiles=quantiles,
        assurance_level_props=assurance_level_props
      )
      if (!is.null(prior)) {
        if (!(prior %in% c('normal','uniform')))
          stop('Invalid type of prior')
      }
      if (!(design %in% c('fe', 'rm')))
        stop('Design should be one of \'fe\' or \'rm\'!')
      if (rho >= 1 | rho < 0)
        stop('rho should be between 0 and 1!')
      if (epsilon > 1 | epsilon < 0)
        stop('epsilon should be between 0 and 1!')

      self$mu <- mu
      self$design <- design
      self$rho <- rho
      self$sd <- sd
      self$epsilon <- epsilon

      if (!(is.null(prior))) {
        if (prior == 'normal')
          self$k <- length(prior_mu)
        else if (prior == 'uniform') {
          self$k <- length(prior_lower)
        }
      }
      else
        self$k <- length(mu)
      if (design == 'fe')
        self$ns <- self$ns * self$k
    },

    print = function() {
      super$print()
      cat('Fixed means for classical power analysis: ', self$mu, '\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'normal') {
          cat('Prior means: ', self$prior_mu, '\n')
          cat('Prior sds: ', self$prior_sigma, '\n\n')
        }
        else if (self$prior == 'uniform') {
          cat('Prior lower bounds: ', self$prior_lower, '\n')
          cat('Prior upper bounds: ', self$prior_upper, '\n\n')
        }
      }
      cat('Test type: One-way ANOVA\n')
      cat('Study design: ', self$design, '\n')
      if (self$design == 'fe')
        cat('Fixed effects One-way ANOVA')
      else
        cat('Repeated measure One-way ANOVA')
    },

    classical_power = function(mu = self$mu, n=self$ns, f2=NULL) {
      if (is.null(mu))
        stop('Input effect size is null')
      else {
        if (is.null(f2))
          f2 <- private$compute_f(mu)^2
        if (self$design == 'fe') {
          ncp <- f2*n
          df1 <- self$k-1
          df2 <- n - self$k
        }
        else {
          u <- self$k / (1-self$rho)
          ncp <- f2*n*u*self$epsilon
          df1 <- (self$k-1)*self$epsilon
          df2 <- (n-1)*(self$k-1)*self$epsilon
          print(c(f2, n, ncp, df1, df2))
        }
        return(private$compute_f_prob(f2, ncp, df1, df2))
      }
    }
  ),

  private = list(
    compute_f = function(means) {
      return(sqrt(var(means)*(self$k-1)/self$k)/self$sd)
    },

    compute_f_prob = function(f, ncp, df1, df2) {
      crit <- qf(
        1-self$alpha,
        df1 = df1,
        df2 = df2
      )
      return(
        1-pf(
          crit,
          ncp = ncp,
          df1 = df1,
          df2 = df2
        )
      )
    },

    draw_prior_es = function() {
      means <- vector()
      if (self$prior == 'normal') {
        for (i in 1:self$k) {
          means <- cbind(
            means,
            rnorm(self$n_prior, self$prior_mu[i], self$prior_sigma[i])
          )
        }
      }
      else if (self$prior == 'uniform') {
        for (i in 1:length(self$k)) {
          means <- cbind(
            means,
            rnorm(self$n_prior, self$prior_lower[i], self$prior_upper[i])
          )
        }
      }
      return(means)
    },

    generate_hybrid_power = function(n) {
      return(
        apply(
          private$draw_prior_es(),
          1,
          FUN=self$classical_power,
          n=n
        )
      )
    }
  )
)
