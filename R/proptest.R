setwd('~/RHybrid/R')
source('HybridPower.R')

HybridPowerProp <- R6Class(
  'HybridPowerProp',
  inherit = HybridPower,
  public = list(
    hybrid_powers = NULL,
    design = NULL,

    pi_1 = NULL,
    pi_2 = NULL,
    c = NULL,

    prior_pi_1_alpha = NULL,
    prior_pi_1_beta = NULL,
    prior_pi_2_alpha = NULL,
    prior_pi_2_beta = NULL,

    prior_pi_1_upper = NULL,
    prior_pi_1_lower = NULL,
    prior_pi_2_upper = NULL,
    prior_pi_2_lower = NULL,

    prior_pi_1_mu = NULL,
    prior_pi_1_sd = NULL,
    prior_pi_2_mu = NULL,
    prior_pi_2_sd = NULL,

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior='beta',
      alpha = 0.05,
      alt = 'two.sided',
      design = 'one.sample',
      pi_1 = NULL,
      pi_2 = NULL,
      c = NULL,

      prior_pi_1_alpha = NULL,
      prior_pi_1_beta = NULL,
      prior_pi_2_alpha = NULL,
      prior_pi_2_beta = NULL,

      prior_pi_1_upper = NULL,
      prior_pi_1_lower = NULL,
      prior_pi_2_upper = NULL,
      prior_pi_2_lower = NULL,

      prior_pi_1_mu = NULL,
      prior_pi_1_sd = NULL,
      prior_pi_2_mu = NULL,
      prior_pi_2_sd = NULL
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
      if (prior == 'beta') {
        if (!(is_numeric(prior_pi_1_alpha) &
             is_numeric(prior_pi_1_beta) &
             (is.null(prior_pi_2_mu) | is_numeric(prior_pi_2_mu)) &
             (is.null(prior_pi_2_sd) | is_numeric(prior_pi_2_sd))
          )
        )
          stop('Invalid input type for the priors')
        if (prior_pi_1_alpha <= 0)
          stop('prior_pi_1_alpha should be positive')
        if (prior_pi_1_beta <= 0)
          stop('prior_pi_1_beta should be positive')
        if (prior_pi_2_alpha <= 0)
          stop('prior_pi_2_alpha should be positive')
        if (prior_pi_2_beta <= 0)
          stop('prior_pi_2_beta should be positive')
        self$prior_pi_1_alpha <- prior_pi_1_alpha
        self$prior_pi_1_beta <- prior_pi_1_beta
        self$prior_pi_2_alpha <- prior_pi_2_alpha
        self$prior_pi_2_beta <- prior_pi_2_beta
      }
      else if (prior == 'uniform') {
        if ((prior_pi_1_lower > prior_pi_1_upper) |
            prior_pi_1_lower < 0 | prior_pi_1_upper > 1 |
            (prior_pi_2_lower > prior_pi_2_upper) |
            prior_pi_2_lower < 0 | prior_pi_2_upper > 1)
          stop('Invalid limits for the uniform prior(s)')
        self$prior_pi_1_lower <- prior_pi_1_lower
        self$prior_pi_1_upper <- prior_pi_1_upper
        self$prior_pi_2_lower <- prior_pi_2_lower
        self$prior_pi_2_upper <- prior_pi_2_upper
      }
      else if (prior == 'truncnorm') {
        if (!(is_numeric(prior_pi_1_mu) &
              is_numeric(prior_pi_1_sd) &
              (is.null(prior_pi_2_mu) | is_numeric(prior_pi_2_mu)) &
              (is.null(prior_pi_2_sd) | is_numeric(prior_pi_2_sd))
          )
        )
          stop('Invalid input type for the priors')
        if (prior_pi_1_mu <= 0 | prior_pi_1_mu >= 1)
          stop('Prior means must be between 0 and 1')
        if (!(is.null(prior_pi_2_mu))) {
          if (prior_pi_2_mu <= 0 | prior_pi_2_mu >= 1 )
            stop('Prior means must be between 0 and 1')
        }
        if (prior_pi_1_sd <= 0)
          stop('Prior standard deviations must be positive')
        if  (!(is.null(prior_pi_2_sd))) {
          if (prior_pi_2_sd <= 0)
            stop('Prior standard deviations must be positive')
        }
        self$prior_pi_1_mu <- prior_pi_1_mu
        self$prior_pi_1_sd <- prior_pi_1_sd
        self$prior_pi_2_mu <- prior_pi_2_mu
        self$prior_pi_2_sd <- prior_pi_2_sd
      }
      if (!(is.null(pi_1))) {
        if (length(pi_1) == 1) {
          if (!(is.numeric(pi_1) | pi_1 < 0 | pi_1 > 1))
            stop('pi_1 is not a valid input')
        }
        else {
          stop('pi_1 should be a constant')
        }
        if (!(is.null(pi_2))) {
          if (length(pi_1) == 1) {
            if (!(is.numeric(pi_1) | pi_1 < 0 | pi_1 > 1))
              stop('pi_2 is not a valid input')
          }
          else {
            stop('pi_2 should be a constant')
          }
        }
      }
      else if (is.null(pi_1) & !(is.null(pi_2))) {
        stop('Single inputs should go into pi_1.')
      }
      if (!(is.null(c))) {
        if (c >= 1 | c <= 0)
          stop('c should be between 0 and 1')
      }

      self$pi_1 <- pi_1
      self$pi_2 <- pi_2
      self$c <- c
    },

    print = function() {
      super$print()
      if (self$prior == 'beta') {
        cat('Alpha (p_1): ', self$prior_pi_1_alpha, '\n')
        cat('Beta (p_1): ', self$prior_pi_1_beta, '\n')
        cat('Alpha (p_2): ', self$prior_pi_2_alpha, '\n')
        cat('Beta (p_2): ', self$prior_pi_2_beta, '\n')
      }
      else if (self$prior == 'uniform') {
        cat('Lower (p_1): ', self$prior_pi_1_lower, '\n')
        cat('Upper (p_1): ', self$prior_pi_1_upper, '\n')
        cat('Lower (p_2): ', self$prior_pi_2_lower, '\n')
        cat('Upper (p_2): ', self$prior_pi_2_upper, '\n')
      }
      else if (self$prior == 'truncnorm') {
        cat('Prior mean (p_1): ', self$prior_pi_1_mu, '\n')
        cat('Prior sd (p_1): ', self$prior_pi_1_sd, '\n')
        cat('Prior mean (p_2): ', self$prior_pi_2_mu, '\n')
        cat('Prior sd (p_2): ', self$prior_pi_2_sd, '\n')
      }
      cat('p_1: ',self$p_1, '\n')
      cat('p_2: ',self$p_2, '\n')
      cat('c: ',self$c, '\n')
      cat('Study design: ', self$design, '\n')
    },

    classical_power = function(n, pi_1, pi_2 = NULL) {
      if (is.null(pi_2)) {
        sd_0 <- sqrt(self$c*(1-self$c)/n)
        mu_1 <- (pi_1 - self$c)
        sd_1 <- sqrt(pi_1*(1-pi_1)/n)
        if (self$alt == 'two.sided') {
          crit_lower <- qnorm(self$alpha/2, 0, sd_0)
          crit_upper <- qnorm(1-self$alpha/2, 0, sd_0)
          return(pnorm(crit_lower, mu_1, sd_1) + 1 - pnorm(crit_upper, mu_1, sd_1))
        }
        else if (self$alt == 'less') {
          crit <- qnorm(self$alpha/2)
          return(pnorm(crit, mu_es, sd_es))
        }
        else {
          crit <- qnorm(1-self$alpha/2)
          return(1 - pnorm(crit, mu_es, sd_es))
        }
      }
      else {
        return(
          power.prop.test(
            n=n,
            p1=pi_1,
            p2=pi_2,
            sig.level=self$alpha,
            alt = self$alt
          )$power
        )
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
    draw_prior_es = function() {
      if (self$prior == 'beta') {
        if (is.null(self$prior_pi_2_alpha) | is.null(self$prior_pi_2_beta)) {
          return(rbeta(self$n_prior, self$prior_pi_1_alpha, self$prior_pi_1_beta))
        }
        else {
          return(
            cbind(
              rbeta(self$n_prior, self$prior_pi_1_alpha, self$prior_pi_1_beta),
              rbeta(self$n_prior, self$prior_pi_2_alpha, self$prior_pi_2_beta)
            )
          )
        }
      }
      else if (self$prior == 'uniform') {
        if (is.null(self$prior_pi_2_lower) | is.null(self$prior_pi_2_upper)) {
          return(runif(self$n_prior, self$prior_pi_1_lower, self$prior_pi_1_upper))
        }
        else {
          return(
            cbind(
              runif(self$n_prior, self$prior_pi_1_lower, self$prior_pi_1_upper),
              runif(self$n_prior, self$prior_pi_2_lower, self$prior_pi_2_upper)
            )
          )
        }
      }
      else if (self$prior == 'truncnorm') {
        if (is.null(self$prior_pi_2_mu) | is.null(self$prior_pi_2_sd)) {
          return(
            truncnorm::rtruncnorm(
              n=self$n_prior,
              a=0,
              b=1,
              mean=self$prior_pi_1_mu,
              sd=self$prior_pi_1_sd
            )
          )
        }
        else {
          return(
            cbind(
              truncnorm::rtruncnorm(
                n=self$n_prior,
                a=0,
                b=1,
                mean=self$prior_pi_1_mu,
                sd=self$prior_pi_1_sd
              ),
              truncnorm::rtruncnorm(
                n=self$n_prior,
                a=0,
                b=1,
                mean=self$prior_pi_2_mu,
                sd=self$prior_pi_2_sd
              )
            )
          )
        }
      }
    },

    hybrid_power = function(n) {
      es <- private$draw_prior_es()
      if (is.null(dim(es))) return(self$classical_power(n, pi_1=es, pi_2=NULL))
      else return(self$classical_power(n, pi_1=es[,1], pi_2=es[,2]))
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

x <- HybridPowerProp$new(
  ns = seq(10, 90, 10),
  n_prior=10,
  prior = 'truncnorm',
  prior_pi_1_mu = .6,
  prior_pi_1_sd = .1,
  c = 0.5,
  alt = 'two.sided'
)

# x$classical_power(n=200, pi_1 = 0.6)
# x$generate_hybrid_power()
# x$assurances()
# x$plot_power(x2$generate_hybrid_power())
