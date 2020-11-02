source('R/HybridPower.R')

hp_oneway_anova <- R6Class(
  'hp_oneway_anova',
  inherit = hp,
  public = list(
    mu = NULL,
    es = NULL,
    sigma = 1,
    k = NULL,
    design = NULL,
    rho = 0,
    epsilon = 1,

    initialize = function(
      parallel = FALSE,
      cores=NULL,
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
      sigma = 1,
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
      if (length(prior_sigma) == 1) prior_sigma <- rep(prior_sigma, length(prior_mu))
      if (length(sigma) == 1) {
        if(is.null(mu)) sigma <- rep(sigma, length(prior_mu))
        else sigma <- rep(sigma, length(mu))
      }
      super$initialize(
        cores=cores,
        parallel = parallel,
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
      if (!(is.numeric(sigma))) {
        stop('sigma must be numeric')
      }
      if (length(sigma) == 1) {
        if (sigma <= 0) {
          stop('sigma must be greater than 0')
        }
      }
      else if (length(unique(sigma)) == 1) {
        sigma <- sigma[1]
        if (sigma <= 0) {
          stop('sigma must be greater than 0')
        }
      }

      else {
        if (design != 'fe') stop('Welch ANOVA power analysis is only supported for fixed effects design')
        if (!(is.null(mu))) {
          if (length(sigma) != length(mu)) stop('Lengths of mu and sigma do not match')
        }
        if (!(is.null(prior_mu))) {
          if (length(sigma) != length(prior_mu)) stop('Lengths of prior_mu and sigma do not match')
        }
        for (i in 1:length(sigma)) {
          if (sigma[i] <= 0) {
            stop('Every sigma must be greater than 0')
          }
        }
      }

      self$mu <- mu
      self$design <- design
      self$rho <- rho
      self$sigma <- sigma
      self$epsilon <- epsilon

      if (!(is.null(prior))) {
        if (prior == 'normal') {
          if (!(is.null(prior_mu))) {
            if (length(prior_mu) != length(prior_sigma)) {
              stop('Lengths of prior_mu and prior_sigma do not match')
            }
          }
          self$k <- length(prior_mu)
        }
        else if (prior == 'uniform') {
          if ((!(is.null(prior_lower))) | (!(is.null(prior_upper))))
            stop('Please provide both prior_lower and prior_upper')
          else {
            if (length(prior_lower) != length(prior_upper)) {
              stop('Lengths of prior_lower and prior_upper do not match')
            }
          }
          self$k <- length(prior_lower)
        }
      }
      else {
        self$k <- length(mu)
      }
    },

    print = function() {
      super$print()
      if (!(is.null(self$mu))) {
        cat('Fixed means for classical power analysis: ', self$mu, '\n')
      }
      cat('Standard deviation(s): ', self$sigma, '\n')
      cat('Sample sizes: ', self$ns/2, '\n')
      if (!(is.null(self$prior))) {
        if (self$prior == 'normal') {
          cat('Prior means: ', self$prior_mu, '\n')
          cat('Prior sigmas: ', self$prior_sigma, '\n\n')
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
        if (length(self$sigma) == 1) {
          if (is.null(f2))
            f2 <- private$compute_f(mu)^2
          if (self$design == 'fe') {
            ncp <- f2*n*self$k
            df1 <- self$k-1
            df2 <- n - self$k
          }
          else {
            u <- self$k / (1-self$rho)
            ncp <- f2*n*u*self$epsilon
            df1 <- (self$k-1)*self$epsilon
            df2 <- (n-1)*(self$k-1)*self$epsilon
          }
          return(private$compute_f_prob(f2, ncp, df1, df2))
        }
        else {
          if (self$parallel) {
            if (self$parallel) {
              if (is.null(self$cores)) cl <- parallel::makeCluster(parallel::detectCores()-1)
              else cl <- parallel::makeCluster(self$cores)
              doParallel::registerDoParallel(cl)
              res <- unlist(parallel::parLapply(cl, n, fun=private$simulate_welch, mu=mu))
              parallel::stopCluster(cl)
              return(res)
            }
          }
          else {
            return(sapply(n, private$simulate_welch, mu=mu))
          }
        }
      }
    }
  ),

  private = list(
    simulate_anova_uneq_var = function(i, n, mu) {
      len_mu <- length(mu)
      data <- vector()
      group <- vector()
      for (i in 1:len_mu) {
        data <- c(data, rnorm(n, mu[i], self$sigma[i]))
        group <- c(group, rep(i, n))
      }
      df <- data.frame(
        data = data,
        group = factor(group)
      )
      res <- oneway.test(data ~ group, data=df, var.equal=F)$p.value < self$alpha
      return(res)
    },

    simulate_welch = function(n, mu) {
      return(mean(sapply(1:self$n_MC, private$simulate_anova_uneq_var, n=n, mu=mu)))
    },

    compute_f = function(means) {
      return(sqrt(var(means)*(self$k-1)/self$k/self$sigma^2))
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

    generate_hybrid_power = function(n, es) {
      return(
        apply(
          es,
          1,
          FUN=self$classical_power,
          n=n
        )
      )
    }
  )
)
