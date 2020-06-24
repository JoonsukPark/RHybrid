source('R/HybridPower.R')

HybridPowerOnewayANOVA <- R6Class(
  'HybridPowerOnewayANOVA',
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
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      mu = NULL,
      prior_mu = c(),
      prior_sigma = c(),
      prior_lower = c(),
      prior_upper = c(),
      sd = 1,
      design = 'fe',
      rho = 0,
      epsilon = 1,
      alt = 'two.sided',
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
        if (prior == 'normal') {
          if (length(prior_mu) == 1)
            stop('Specify more than 2 groups\' prior means')
          if (sum(prior_sigma < 0) > 0)
            stop('Prior standard deviations must be nonnegative')
          if (length(prior_mu) == 1 | (length(prior_mu) > 1 & length(prior_sigma) > 1)) {
            if (length(prior_mu) != length(prior_sigma))
              stop('Lengths of prior means and sds should be identical')
          }
          else
            prior_sigma <- rep(prior_sigma, length(prior_mu))
          self$prior_mu <- prior_mu
          self$prior_sigma <- prior_sigma
          self$k <- length(prior_mu)
        }
        else if (prior == 'uniform') {
          if (length(prior_lower) == 1)
            stop('Specify more than 2 groups\' prior means')
          if (length(prior_upper) == 1)
            stop('Specify more than 2 groups\' prior means')
          if (length(prior_lower) != length(prior_upper))
            stop('Lengths of prior lower and upper bounds should be identical')
          self$prior_lower <- prior_lower
          self$prior_upper <- prior_upper
          self$k <- length(prior_lower)
        }
      }
      else
        self$k <- length(mu)
    },

    print = function() {
      super$print()
      cat('Fixed means for classical power analysis: ', self$mu, '\n')
      if (self$prior == 'normal') {
        cat('Prior means: ', self$prior_mu, '\n')
        cat('Prior sds: ', self$prior_sigma, '\n\n')
      }
      else if (self$prior == 'uniform') {
        cat('Prior lower bounds: ', self$prior_lower, '\n')
        cat('Prior upper bounds: ', self$prior_upper, '\n\n')
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
          ncp <- f2*n*u
          df1 <- (self$k-1)*self$epsilon
          df2 <- (n-1)*(self$k-1)*self$epsilon
        }
        return(private$compute_f_prob(f2, ncp, df1, df2))
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
