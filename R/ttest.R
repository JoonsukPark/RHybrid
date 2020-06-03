setwd('~/HybridPower')
source('HybridPower.R')

HybridPowerTtest <- R6Class(
  'HybridPowerTtest',
  inherit = HybridPower,
  public = list(
    es = NULL,
    design = NULL,
    hybrid_powers = NULL,
    
    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior='normal',
      alpha = 0.05,
      alt = 'two.sided',
      prior_mu,
      prior_sd,
      prior_lower,
      prior_upper,
      design = 'one.sample'
    ) {
      super$initialize(
        parallel = FALSE,
        ns,
        n_prior,
        n_MC,
        prior,
        alpha,
        alt,
        prior_mu,
        prior_sd,
        prior_lower,
        prior_upper
      )
      self$design <- design
    },
    
    print = function() {
      super$print()
      cat('Test type: t-test\n')
      cat('Study design: ', self$design, '\n')
    },
    
    classical_power = function(n) {
      delta <- ifelse(
        self$prior == 'normal',
        self$prior_mu,
        (self$prior_lower + self$prior_upper)/2
      )
      return(
        power.t.test(
          n = n,
          delta = delta,
          sig.level = self$alpha,
          alt = self$alt,
          type = self$design,
          power=NULL
        )$power
      )
    },
    
    draw_prior_es = function() {
      if (self$prior == 'normal') {
        return(
          rnorm(self$n_prior, self$prior_mu, self$prior_sd)
        )
      }
      else if (self$prior == 'uniform') {
        return(
          runif(self$n_prior, self$prior_lower, self$prior_upper)
        )
      }
    },

    hybrid_power = function(n) {
      es <- self$draw_prior_es()
      return(
        power.t.test(
          n = n,
          delta = es,
          sig.level = self$alpha,
          alt = self$alt,
          type = self$design,
          power = NULL
        )$power
      )
    },
    
    assurance = function(n) {
      return(
        mean(self$hybrid_power(n))
      )
    },
    
    generate_hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        return(mclapply(self$ns, self$hybrid_power))
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- self$hybrid_power(self$ns[i])
        }
        return(res)
      }
    }
  )
)
