setwd('~/RHybrid/R')
source('HybridPower.R')

HybridPowerSignTest <- R6Class(
  'HybridPowerSigntest',
  inherit = HybridPower,
  public = list(
    p_0 = NULL,
    MC = F,
    hybrid_powers = NULL,
    prior_mu = NULL,
    prior_sd = NULL,
    prior_lower = NULL,
    prior_upper = NULL,
    prior_alpha = NULL,
    prior_beta = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior='beta',
      alpha = 0.05,
      p_0 = 0.5,
      alt = 'two.sided',
      MC=F,
      prior_mu = 0.5,
      prior_sd = 0.1,
      prior_lower = 0,
      prior_upper = 1,
      prior_alpha = 1,
      prior_beta = 1
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
      self$p_0 <- p_0
      self$MC <- MC
      if (prior == 'beta') {
        if (!(is_numeric(prior_alpha) & is_numeric(prior_beta)))
          stop('Invalid input type for the priors')
        if (prior_alpha <= 0)
          stop('prior_alpha should be positive')
        if (prior_beta <= 0)
          stop('prior_beta should be positive')
        self$prior_alpha <- prior_alpha
        self$prior_beta <- prior_beta
      }
      if (prior == 'truncnorm') {
        if (prior_mu <= 0 | prior_mu >= 1)
          stop('prior_mu should be between 0 and 1')
        if (prior_sd <= 0)
          stop('prior sd should be positive')
        self$prior_mu <- prior_mu
        self$prior_sd <- prior_sd
      }
      else if (prior == 'uniform') {
        if (prior_lower < 0 | prior_lower > 1)
          stop('The lower bound should be between 0 and 1')
        if (prior_upper < 0 | prior_upper > 1)
          stop('The upper bound should be between 0 and 1')
        if (prior_lower >= prior_upper)
          stop('The lower bound should be smaller than the upper bound')
        self$prior_lower <- prior_lower
        self$prior_upper <- prior_upper
      }
    },

    print = function() {
      super$print()
      if (self$prior == 'beta') {
        cat('Prior alpha: ', self$prior_alpha, '\n')
        cat('Prior beta: ', self$prior_beta, '\n\n')
      }
      else if (self$prior == 'normal') {
        cat('Prior mean: ', self$prior_mu, '\n')
        cat('Prior sd: ', self$prior_sd, '\n\n')
      }
      else if (self$prior == 'uniform') {
        cat('Prior lower bound: ', self$prior_lower, '\n')
        cat('Prior upper bound: ', self$prior_upper, '\n\n')
      }
      cat('Test type: Sign Test\n')
    },

    classical_power = function(n, p, p_0) {
      if (self$MC) {
        return(
          mean(
            sapply(rbinom(self$n_MC, n, p), FUN=private$binom_test, n=n, p_0=p_0)
          )
        )
      }

      else {
        mu_0 <- n*p_0
        mu_1 <- n*p
        sigma_0 <- sqrt(mu_0*(1-p_0))
        sigma_1 <- sqrt(mu_1*(1-p))
        rho <- sigma_1 / sigma_0
        if (self$alt == 'two.sided') {
          z_lower <- ((mu_0 - mu_1) / sigma_0 + qnorm(self$alpha/2))/rho
          z_upper <- ((mu_0 - mu_1) / sigma_0 + qnorm(1-self$alpha/2))/rho
          return(pnorm(z_lower) + 1 - pnorm(z_upper))
        }
        else if (self$alt == 'greater') {
          z_upper <- ((mu_0 - mu_1) / sigma_0 + qnorm(1-self$alpha/2))/rho
          return(1 - pnorm(z_upper))
        }
        else {
          z_lower <- ((mu_0 - mu_1) / sigma_0 + qnorm(self$alpha/2))/rho
          return(pnorm(z_lower))
        }
      }
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
    binom_test = function(x, n, p_0) {
      if (self$alt == 'two.sided') {
        crit_lower <- qbinom(self$alpha/2, n, p_0)-1
        crit_upper <- qbinom(1-self$alpha/2, n, p_0)
        return((x <= crit_lower) | (x >= crit_upper))
      }
      else if (self$alt == 'greater') {
        return(x >= qbinom(1-self$alpha/2, n, p_0))
      }
      else {
        return(x <= qbinom(self$alpha/2, n, p_0)-1)
      }
    },

    draw_prior_es = function() {
      if (self$prior == 'beta') {
        return(
          rbeta(self$n_prior, self$prior_alpha, self$prior_beta)
        )
      }
      else if (self$prior == 'truncnorm') {
        return(
          truncnorm::rtruncnorm(self$n_prior, mean=self$prior_mu, sd=self$prior_sd, a=0, b=1)
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
        sapply(es, FUN=self$classical_power, n=n, p_0=self$p_0)
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
      return(
        mean(private$hybrid_power(n))
      )
    }
  )
)

x <- HybridPowerSignTest$new(
  prior='uniform',
  parallel = T,
  ns = seq(10, 90, 10),
  n_prior=1000,
  n_MC = 100,
  p_0 = 0.5,
  MC=F
)

x$classical_power(n=50, p=0.65, p_0=0.5)
x$generate_hybrid_power()
x$assurances()
x$plot_power(x$generate_hybrid_power())

