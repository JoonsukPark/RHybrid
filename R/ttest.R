source('R/HybridPower.R')

hp_ttest <- R6Class(
  'hp_ttest',
  inherit = hp,
  public = list(
    es = NULL,
    design = NULL,
    d = NULL,
    sigma = 1,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      alt = 'two.sided',
      d = NULL,
      prior_mu = NULL,
      prior_sigma = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      design = 'one.sample',
      sigma = 1,
      quantiles = NULL,
      assurance_level_props=NULL
    ) {
      if (!(is.null(prior))) {
        if (!(prior %in% c('normal', 'uniform')))
          stop('Invalid prior')
      }
      super$initialize(
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
      self$design <- design
      if (!(is.null(d))) {
        if (!(is.numeric(d)))
          stop('Input effect size should be numeric')
        if (length(d) != 1)
          stop('Cohen\'s d must be a single number!')
      }
      if (!(is.null(d)))
        self$d <- abs(d)

      if (is.numeric(sigma) & length(sigma) <= 2) {
        for (i in 1:length(sigma)) {
          if (sigma[i] <= 0)
            stop('Invalid sigma')
        }
        self$sigma <- sigma
      }
      else
        stop('Invalid sigma')
      if (!(is.null(prior)) & is.null(sigma))
        stop('sigma cannot be null for hybrid power calculation')
    },

    print = function() {
      super$print()
      if (!(is.null(self$d)))
        cat('Cohen\'s d: ', self$d, '\n')
      if (!(is.null(self$sigma)))
        cat('Standard deviation for data: ', self$sigma, '\n')
      if (!(is.null(prior))) {
        if (self$prior == 'normal') {
          cat('Prior mean: ', self$prior_mu, '\n')
          cat('Prior sigma: ', self$prior_sigma, '\n\n')
        }
        else if (self$prior == 'uniform') {
          cat('Prior lower bound: ', self$prior_lower, '\n')
          cat('Prior upper bound: ', self$prior_upper, '\n\n')
        }
      }
      cat('Test type: t-test\n')
      cat('Study design: ', self$design, '\n')
    },

    classical_power = function(d = self$d, n = self$ns, cores=NULL) {
      if (is.null(d))
        stop('Effect size is not provided!')
      else {
        if (length(self$sigma) == 1 | (length(self$sigma) == 2 & self$sigma[1] == self$sigma[2])) {
          if (length(self$sigma) > 1) sigma <- self$sigma[1]
          else sigma <- self$sigma
          return(
            power.t.test(
              n = n,
              delta = d,
              sd = 1,
              sig.level = self$alpha,
              alt = self$alt,
              type = self$design,
              power=NULL
            )$power
          )
        }
        else {
          if (self$parallel) {
            if (is.null(cores)) cl <- parallel::makeCluster(parallel::detectCores()-1)
            else cl <- parallel::makeCluster(cores)
            doParallel::registerDoParallel(cl)
            res <- unlist(parallel::parLapply(cl, n, fun=private$mc_ttest))
            parallel::stopCluster(cl)
            return(res)
          }
          else {
            return(unlist(lapply(n, private$mc_ttest)))
          }
        }
      }
    }
  ),

  private = list(
    draw_prior_es = function() {
      if (self$prior == 'normal') {
        return(
          rnorm(self$n_prior, self$prior_mu, self$prior_sigma)
        )
      }
      else if (self$prior == 'uniform') {
        return(
          runif(self$n_prior, self$prior_lower, self$prior_upper)
        )
      }
      else
        stop('Invalid prior type')
    },

    sim_ttest = function(i, d, n) {
      x <- rnorm(n, 0, self$sigma[1])
      y <- rnorm(n, d, self$sigma[2])
      return(t.test(x, y, var.equal=F, paired=F)$p.value < self$alpha)
    },

    mc_ttest = function(n, d=self$d) {
      return(mean(unlist(lapply(1:self$n_MC, private$sim_ttest, d=d, n=n))))
    },

    generate_hybrid_power = function(es, n) {
      if (length(self$sigma) == 1 | (length(self$sigma) == 2 & self$sigma[1] == self$sigma[2])) {
        return(sapply(es, FUN=self$classical_power, n=n))
      }
      else {
          return(sapply(es, FUN=private$mc_ttest, n=n))
      }
    }
  )
)
