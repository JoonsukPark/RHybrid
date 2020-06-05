setwd('~/HybridPower')
source('HybridPower.R')

library(reshape2)

HybridPowerOnewayANOVA <- R6Class(
  'HybridPowerOnewayANOVA',
  inherit = HybridPower,
  public = list(
    es = NULL,
    hybrid_powers = NULL,
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
      prior='normal',
      alpha = 0.05,
      prior_mu = c(),
      prior_sd = c(),
      prior_lower = c(),
      prior_upper = c(),
      sd = 1,
      design = 'fe',
      rho = 0,
      epsilon = 1,
      alt = 'one.sided'
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

      if (!(design %in% c('fe', 'rm')))
        stop('Design should be one of \'fe\' or \'rm\'!')
      if (self$prior == 'normal') {
          if (length(prior_mu) == 1)
            stop('Specify more than 2 groups\' prior means')
          if (length(prior_sd) == 1)
            stop('Specify more than 2 groups\' prior means')
          if (length(prior_mu) != length(prior_sd))
            stop('Lengths of prior means and sds should be identical')
      }
      else if (prior == 'uniform') {
        if (length(prior_lower) == 1)
          stop('Specify more than 2 groups\' prior means')
        if (length(prior_upper) == 1)
          stop('Specify more than 2 groups\' prior means')
        if (length(prior_lower) != length(prior_upper))
          stop('Lengths of prior lower and upper bounds should be identical')
      }
      if (rho > 1 | rho < 0)
        stop('rho should be between 0 and 1!')
      if (epsilon > 1 | epsilon < 0)
        stop('epsilon should be between 0 and 1!')

      self$design <- design
      self$rho <- rho
      self$sd <- sd
      self$epsilon <- epsilon
      if (prior == 'normal') {
        self$prior_mu <- prior_mu
        self$prior_sd <- prior_sd
        self$k <- length(prior_mu)
      }
      else if (prior == 'uniform') {
        self$prior_lower <- prior_lower
        self$prior_upper <- prior_upper
        self$k <- length(prior_lower)
      }
    },

    print = function() {
      super$print()
      if (self$prior == 'normal') {
        cat('Prior means: ', self$prior_mu, '\n')
        cat('Prior sds: ', self$prior_sd, '\n\n')
      }
      else if (self$prior == 'uniform') {
        cat('Prior lower bounds: ', self$prior_lower, '\n')
        cat('Prior upper bounds: ', self$prior_upper, '\n\n')
      }
      cat('Test type: ANOVA\n')
      cat('Study design: ', self$design, '\n')
      if (self$design == 'fe')
        cat('Fixed effects One-way ANOVA')
      else
        cat('Repeated measure One-way ANOVA')
    },

    classical_power = function() {

      if (self$prior == 'normal')
        f <- private$compute_f(self$prior_mu)
      else
        f <- private$compute_f((self$prior_lower + self$prior_upper)/2)
      if (self$design == 'fe') {
        ncp <- f^2*self$ns
        df1 <- self$k-1
        df2 <- self$ns - self$k
      }
      else {
        u <- self$k / (1-self$rho)
        ncp <- f^2*self$ns*u
        df1 <- (self$k-1)*self$epsilon
        df2 <- (self$ns-1)*(self$k-1)*self$epsilon
      }
      return(private$compute_f_prob(f, ncp, df1, df2))
    },

    generate_hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        return(
          private$melt_powers(mclapply(self$ns, private$hybrid_power)))
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
            rnorm(self$n_prior, self$prior_mu[i], self$prior_sd[i])
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
      return(apply(means, 1, private$compute_f))
    },

    hybrid_power = function(n) {
      f <- private$draw_prior_es()
      if (self$design == 'fe') {
        ncp <- f^2*n
        df1 <- self$k-1
        df2 <- n - self$k
      }
      else {
        u <- self$k / (1-self$rho)
        ncp <- f^2*n*u
        df1 <- (self$k-1)*self$epsilon
        df2 <- (n-1)*(self$k-1)*self$epsilon
      }
      return(private$compute_f_prob(f, ncp, df1, df2))
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

# x <- HybridPowerOnewayANOVA$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior_mu = c(2, 2.2),
#   prior_sd = c(0.5, 0.5),
#   sd = 1,
#   rho = 0.1,
#   design='rm'
# )
#
# x$classical_power()
# x$generate_hybrid_power()
# x$assurances()
# x$plot_power(x$generate_hybrid_power())
