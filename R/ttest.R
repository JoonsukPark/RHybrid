setwd('~/RHybrid/R')
source('HybridPower.R')

HybridPowerTtest <- R6Class(
  'HybridPowerTtest',
  inherit = HybridPower,
  public = list(
    es = NULL,
    design = NULL,
    hybrid_powers = NULL,
    d = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
    prior_lower = NULL,
    prior_upper = NULL,

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
      if (!(is.null(d))) {
        if (!(is.numeric(d)))
          stop('Input effect size should be numeric')
        if (length(d) != 1)
          stop('Cohen\'s d must be a single number!')
      }
      self$d <- d

      if (!(is.null(prior))) {
        if (prior == 'normal') {
          if (!(is.numeric(prior_mu)) | !(is.numeric(prior_sd)))
            stop('prior parameters must be numeric')
          if (prior_sd <= 0)
            stop('prior_sd must be positive')
          self$prior_mu <- prior_mu
          self$prior_sd <- prior_sd
        }
        else if (prior == 'uniform') {
          if (!(is.numeric(prior_lower)) | !(is.numeric(prior_upper)))
            stop('prior parameters must be numeric')
          if (prior_lower > prior_upper)
            stop('The lower bound cannot be greater than the upper bound')
          self$prior_lower <- prior_lower
          self$prior_upper <- prior_upper
        }
      }
    },

    print = function() {
      super$print()
      if (!(is.null(self$d)))
        cat('Cohen\'s d: ', self$d, '\n')
      if (!(is.null(prior))) {
        if (self$prior == 'normal') {
          cat('Prior mean: ', self$prior_mu, '\n')
          cat('Prior sd: ', self$prior_sd, '\n\n')
        }
        else if (self$prior == 'uniform') {
          cat('Prior lower bound: ', self$prior_lower, '\n')
          cat('Prior upper bound: ', self$prior_upper, '\n\n')
        }
      }
      cat('Test type: t-test\n')
      cat('Effect size type: Cohen\'s d')
      cat('Study design: ', self$design, '\n')
    },

    classical_power = function(d = self$d, n = self$ns) {
      if (is.null(d))
        stop('Input effect size not provided!')
      else {
        return(
          power.t.test(
            n = n,
            delta = d,
            sig.level = self$alpha,
            alt = self$alt,
            type = self$design,
            power=NULL
          )$power
        )
      }
    },

    generate_hybrid_power = function(cores=NULL) {
      if (is.null(self$prior))
        stop('Specify a prior first')
      else {
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
      else
        stop('Invalid prior type')
    },

    hybrid_power = function(n) {
      return(
        sapply(
          private$draw_prior_es(),
          FUN=self$classical_power,
          n=n
        )
      )
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
      return(mean(private$hybrid_power(n)))
    }
  )
)

power_classical <- HybridPowerTtest$new(
  ns = seq(10, 90, 10),
  d = 0.5
)
power_classical$classical_power()

power_hybrid <- HybridPowerTtest$new(
  ns = seq(10, 90, 10),
  n_prior=1000,
  prior = 'normal',
  prior_mu = 0.3,
  prior_sd = 0.1
)

# This should generate an error
power_hybrid$classical_power()

power_hybrid$generate_hybrid_power()
power_hybrid$assurances()
powers <- power_hybrid$generate_hybrid_power()
power_hybrid$plot_power(powers)

power_both <- HybridPowerTtest$new(
  d = 0.5,
  ns = seq(10, 90, 10),
  n_prior=1000,
  prior = 'normal',
  prior_mu = 0.3,
  prior_sd = 0.1
)

power_both$classical_power()
power_both$generate_hybrid_power()
power_both$assurances()
powers <- power_both$generate_hybrid_power()
power_both$plot_power(powers)
