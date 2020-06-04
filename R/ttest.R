source('R/HybridPower.R')
library(reshape2)
library(dplyr)

HybridPowerTtest <- R6Class(
  'HybridPowerTtest',
  inherit = HybridPower,
  public = list(
    es = NULL,
    design = NULL,

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

    classical_power = function() {
      delta <- ifelse(
        self$prior == 'normal',
        self$prior_mu,
        (self$prior_lower + self$prior_upper)/2
      )
      powers <- power.t.test(
        n = self$ns,
        delta = delta,
        sig.level = self$alpha,
        alt = self$alt,
        type = self$design,
        power=NULL
      )$power
      power_df <- data.frame(self$ns)
      colnames(power_df) <- 'n'
      power_df$power <- powers
      return(power_df)
    }
  ),

  private = list(
    hybrid_power = function(n) {
      es <- private$draw_prior()
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
    }
  )
)

