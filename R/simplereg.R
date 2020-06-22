source('R/HybridPower.R')

HybridPowerSLR <- R6Class(
  'HybridPowerSLR',
  inherit = HybridPower,
  public = list(
    hybrid_powers = NULL,
    prior_alpha = NULL,
    prior_beta = NULL,
    prior_upper = NULL,
    prior_lower = NULL,
    prior_mu = NULL,
    prior_sd = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior='beta',
      alpha = 0.05,
      alt = 'two.sided',

      prior_alpha = NULL,
      prior_beta = NULL,

      prior_upper = NULL,
      prior_lower = NULL,

      prior_mu = NULL,
      prior_sd = NULL
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
      if (prior == 'beta') {
        if (!(is_numeric(prior_pi_1_alpha) &
              is_numeric(prior_pi_1_beta)))
          stop('Invalid input type for the priors')
        if (prior_alpha <= 0)
          stop('prior_alpha should be positive')
        if (prior_beta <= 0)
          stop('prior_beta should be positive')
        self$prior_alpha <- prior_alpha
        self$prior_beta <- prior_beta
      }
      else if (prior == 'uniform') {
        if ((prior_pi_1_lower > prior_pi_1_upper) |
             prior_pi_1_lower < 0 | prior_pi_1_upper > 1
            )
          stop('Invalid limits for the uniform prior(s)')
        self$prior_lower <- prior_lower
        self$prior_upper <- prior_upper
      }
      else if (prior == 'truncnorm') {
        if (!(is_numeric(prior_mu) &
              is_numeric(prior_sd))
        )
          stop('Invalid input type for the priors')
        if (prior_mu <= 0 | prior_mu >= 1)
          stop('Prior means must be between 0 and 1')
        if (prior_sd <= 0)
          stop('Prior standard deviations must be positive')
        self$prior_mu <- prior_mu
        self$prior_sd <- prior_sd
      }
    },

    print = function() {
      super$print()
      if (self$prior == 'beta') {
        cat('Alpha: ', self$prior_alpha, '\n')
        cat('Beta: ', self$prior_beta, '\n')
      }
      else if (self$prior == 'uniform') {
        cat('Lower (r2): ', self$prior_lower, '\n')
        cat('Upper (r2): ', self$prior_upper, '\n')
      }
      else if (self$prior == 'truncnorm') {
        cat('Prior mean (r2): ', self$prior_mu, '\n')
        cat('Prior sd (r2): ', self$prior_sd, '\n')
      }
      cat('Study design: Simple Regression\n')
    },

    classical_power = function(n, r2) {
      f2 <- r2 / (1-r2)
      return(
        1-pf(qf(1-self$alpha, df1=1, df2=n-2), df1=1, df2=n-2, ncp=f2*n)
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
      if (self$prior == 'beta') {
        return(rbeta(self$n_prior, self$prior_alpha, self$prior_beta))
      }
      else if (self$prior == 'uniform') {
        return(runif(self$n_prior, self$prior_lower, self$prior_upper))
      }
      else if (self$prior == 'truncnorm') {
        return(
          truncnorm::rtruncnorm(
            n=self$n_prior,
            a=0,
            b=1,
            mean=self$prior_mu,
            sd=self$prior_sd
          )
        )
      }
    },

    hybrid_power = function(n) {
      return(self$classical_power(n, r2=private$draw_prior_es()))
    },

    melt_powers = function(power_list) {
      powers <- data.frame(power_list)
      colnames(powers) = self$ns
      return(
        suppressMessages(
          reshape2::melt(powers, variable.name='n', value.name = 'power')
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

# x2 <- HybridPowerSLR$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior = 'truncnorm',
#   prior_mu = .2,
#   prior_sd = .2,
#   alt = 'two.sided'
# )
#
# x2$classical_power(n=30, r2=0.3)
# x2$generate_hybrid_power()
# x2$assurances()
# x2$plot_power(x2$generate_hybrid_power())
