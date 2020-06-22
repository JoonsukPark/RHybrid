source('R/HybridPower.R')

HybridPowerTtest <- R6Class(
  'HybridPowerTtest',
  inherit = HybridPower,
  public = list(
    es = NULL,
    design = NULL,
    d = NULL,
    sd = 1,

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
      sd = 1,
      assurance_props = c(0.5)
    ) {
      super$initialize(
        parallel = FALSE,
        ns=ns,
        n_prior=n_prior,
        n_MC=n_MC,
        prior=prior,
        alpha=alpha,
        alt=alt,
        assurance_props=assurance_props
      )
      self$design <- design
      if (!(is.null(d))) {
        if (!(is.numeric(d)))
          stop('Input effect size should be numeric')
        if (length(d) != 1)
          stop('Cohen\'s d must be a single number!')
      }
      self$d <- d

      if (is.numeric(sd)) {
        if (sd > 0)
          self$sd <- sd
      }
      else
        stop('Invalid sd')
      if (!(is.null(prior)) & is.null(sd))
        stop('sd cannot be null for hybrid power calculation')

      if (!(is.null(prior))) {
        if (prior == 'normal') {
          if (!(is.numeric(prior_mu)) | !(is.numeric(prior_sigma)))
            stop('prior parameters must be numeric')
          if (prior_sigma <= 0)
            stop('prior_sigma must be positive')
          self$prior_mu <- prior_mu
          self$prior_sigma <- prior_sigma
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
      if (!(is.null(self$sd)))
        cat('Standard deviation for data: ', self$sd, '\n')
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

    hybrid_power = function(cores=NULL) {
      if (is.null(self$prior))
        stop('Specify a prior first')
      else {
        if (self$parallel) {
          library(parallel)
          if (!(cores)) cores <- detectCores()
          self$output <- private$melt_output(
            mclapply(self$ns, private$generate_hybrid_power)
          )
        }
        else {
          res <- list()
          for (i in 1:length(self$ns)) {
            res[[i]] <- private$generate_hybrid_power(self$ns[i])
          }
          self$output <- private$melt_output(res)
        }
        return(self$output)
      }
    },

    assurance = function() {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      return(summarise(group_by(self$output, n), assurance = mean(power), .groups='keep'))
    },

    assurance_level = function(props=self$assurance_props) {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      if (is.null(props))
        stop('Provide target proportions')
      for (i in 1:length(props))
        if (!(is.numeric(props[i])) | props[i] > 1 | props[i] < 0)
          stop('Invalid proportion(s)')
      props <- sort(props)
      res <- summarise(group_by(self$output, n), quantile(power, probs=props[1]), .groups='keep')
      if (length(props) > 1) {
        for (i in 2:length(props)) {
          res <- left_join(res, summarise(group_by(self$output, n), quantile(power, probs=props[i]), .groups='keep'), by='n')
        }
      }
      col_names <- c('n', props)
      colnames(res) <- col_names
      return(res)
    },

    boxplot = function() {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      p <- ggplot(self$output, aes(x=factor(n), y=power)) + geom_boxplot()
      p <- p + xlab('Sample Size') + ylab('Power') + ggtitle('Distributions of Power')
      p <- p + stat_summary(fun=mean, geom="point", shape=5, size=4)
      p
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

    generate_hybrid_power = function(n) {
      return(
        sapply(
          private$draw_prior_es(),
          FUN=self$classical_power,
          n=n
        )
      )
    },

    melt_output = function(power_list) {
      output <- data.frame(power_list)
      colnames(output) = self$ns
      return(
        suppressMessages(
          reshape2::melt(output, variable.name='n', value.name = 'power')
        )
      )
    }
  )
)

