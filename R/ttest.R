setwd('~/HybridPower')
source('HybridPower.R')

HybridPowerTtest <- R6Class(
  'HybridPowerTtest',
  inherit = HybridPower,
  public = list(
    es = NULL,
    design = NULL,
    hybrid_powers = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
    prior_lower = NULL,
    prior_upper = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior='normal',
      alpha = 0.05,
      alt = 'two.sided',
      prior_mu = NULL,
      prior_sd = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      design = 'one.sample'
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
      self$design <- design
      if (prior == 'normal') {
        self$prior_mu <- prior_mu
        self$prior_sd <- prior_sd
      }
      else if (prior == 'uniform') {
        self$prior_lower <- prior_lower
        self$prior_upper <- prior_upper
      }
    },

    print = function() {
      super$print()
      if (self$prior == 'normal') {
        cat('Prior mean: ', self$prior_mu, '\n')
        cat('Prior sd: ', self$prior_sd, '\n\n')
      }
      else if (self$prior == 'uniform') {
        cat('Prior lower bound: ', self$prior_lower, '\n')
        cat('Prior upper bound: ', self$prior_upper, '\n\n')
      }
      cat('Test type: t-test\n')
      cat('Study design: ', self$design, '\n')
    },

    classical_power = function() {
      delta <- ifelse(
        self$prior == 'normal',
        self$prior_mu,
        (self$prior_lower + self$prior_upper)/2
      )
      return(
        power.t.test(
          n = self$ns,
          delta = delta,
          sig.level = self$alpha,
          alt = self$alt,
          type = self$design,
          power=NULL
        )$power
      )
    },

    generate_hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        return(
          private$melt_powers(
            mclapply(self$ns, private$hybrid_power)
          )
        )
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$hybrid_power(self$ns[i])
        }
        return(private$melt_powers(res))
      }
    },

    assurances = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        return(mclapply(self$ns, private$assurance))
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$assurance(self$ns[i])
        }
        return(res)
      }
    },

    plot_power = function(power_df) {
      p <- ggplot(power_df, aes(x=factor(n), y=power)) + geom_boxplot()
      p <- p + xlab('Sample Size') + ylab('Power') + ggtitle('Distributions of Power')
      p <- p + stat_summary(fun=mean, geom="point", shape=5, size=4)
      p
    }
  ),

  private = list(

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
      es <- private$draw_prior_es()
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

    melt_powers = function(power_list) {
      powers <- data.frame(power_list)
      colnames(powers) = self$ns
      return(
        suppressMessages(
          melt(powers, variable.name='n', value.name = 'power')
        )
      )
    },

    assurance = function(n) {
      return(
        mean(private$hybrid_power(n))
      )
    }
  )
)

# x <- HybridPowerTtest$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior_mu = 0.3,
#   prior_sd = 0.1
# )
#
# x$classical_power()
# x$generate_hybrid_power()
# x$assurances()
# x$plot_power(x$generate_hybrid_power())
#
