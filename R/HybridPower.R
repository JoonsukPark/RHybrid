library(R6)
library(ggplot2)
library(dplyr)
library(reshape2)
library(parallel)

is_numeric <- function(x) return(is.numeric(x) & length(x) == 1)

hp <- R6Class(
  'HybridPower',
  public = list(
    es = NULL,
    design = NULL,
    parallel = FALSE,
    cores = NULL,
    ns = c(),
    n_prior = NULL,
    n_MC = NULL,
    prior = NULL,
    alpha = NULL,
    alt = NULL,
    output = NULL,
    hybrid_output = NULL,
    prior_mu = NULL,
    prior_sigma = NULL,
    prior_lower = NULL,
    prior_upper = NULL,
    prior_a = NULL,
    prior_b = NULL,
    quantiles = c(0, .25, .5, .75, 1),
    assurance_level_props = NULL,

    initialize = function(
      parallel = FALSE,
      cores = NULL,
      ns = c(),
      n_MC = 100,
      alpha = 0.05,
      alt = 'two.sided',
      quantiles = c(0, .25, .5, .75, 1),
      assurance_level_props = NULL,
      prior = NULL,
      n_prior = 100,
      prior_mu = NULL,
      prior_sigma = NULL,
      prior_lower = NULL,
      prior_upper = NULL,
      prior_a = NULL,
      prior_b = NULL,
      validate = T
    ) {
      # Validate inputs

      if (!(is.null(prior))) {
        if (
          !is.character(prior) |
          length(prior) != 1 |
          !(prior %in% c('normal', 'uniform', 'beta', 'truncnorm', 'dirichlet'))
        )
          stop('Invalid prior!')
        if (validate) {
          if (prior == 'normal' | prior == 'truncnorm') {
            if (is.null(prior_mu))
              stop('prior_mu is missing')
            if (is.null(prior_sigma))
              stop('prior_sigma is missing')
            if (!(is.numeric(prior_mu)) | !(is.numeric(prior_sigma)))
              stop('prior parameters must be numeric')
            if (length(prior_sigma) == 1) {
              if (prior_sigma < 0)
                stop('prior_sigma must be nonnegative')
            }
            else {
              for (i in 1:length(prior_sigma)) {
                if (prior_sigma[i] < 0)
                  stop('All prior_sigmas must be nonnegative')
              }
            }
            self$prior_mu <- prior_mu
            self$prior_sigma <- prior_sigma
            if (prior == 'truncnorm') {
              if (is.null(prior_lower) != is.null(prior_upper)) {
                stop('prior_lower and prior_upper must be both NULL or specified!')
              }
              if (!(is.null(prior_lower))) {
                if (prior_lower >= prior_upper) {
                  stop('prior_lower must be smaller than prior_upper')
                }
              }
              self$prior_lower = prior_lower
              self$prior_upper = prior_upper
            }
          }
          else if (prior == 'uniform') {
            if (is.null(prior_lower))
              stop('prior_lower is missing')
            if (is.null(prior_upper))
              stop('prior_upper is missing')
            if (!(is.numeric(prior_lower)) | !(is.numeric(prior_upper)))
              stop('prior parameters must be numeric')
            if (prior_lower > prior_upper)
              stop('The lower bound cannot be greater than the upper bound')
            self$prior_lower <- prior_lower
            self$prior_upper <- prior_upper
          }
          else if (prior == 'beta') {
            if (is.null(prior_a))
              stop('prior_a is missing')
            if (is.null(prior_b))
              stop('prior_b is missing')
            if (is.null(prior_a) | is.null(prior_b))
              stop('Provide prior_a and prior_b first')
            if (prior_a <= 0 | prior_b <= 0)
              stop('Parameters for the beta prior must be positive')
            self$prior_a <- prior_a
            self$prior_b <- prior_b
          }
          else if (prior == 'dirichlet') {
            if (is.null(prior_a))
              stop('prior_a is missing')
            for (i in 1:length(prior_a)) {
              if (!(is_numeric(prior_a[i])))
                stop('Invalid input type for the priors')
              if (length(prior_a[i]) != 1)
                stop('prior_a should be a scalar')
              if (prior_a[i] <= 0)
                stop('Elements of prior_a must be positive')
            }
            self$prior_a <- prior_a
          }
        }
      }

      if (length(ns) == 0) stop('Input sample sizes!')
      else {
        for (i in 1:length(ns)) {
          if (!is.numeric(ns[i]) | ns[i] %% 1 != 0 | ns[i] <= 0)
            stop('Invalid sample size!')
        }
      }
      if (n_prior %% 1 != 0 | n_prior <= 0)
        stop('Invalid # draws from prior!')
      if (n_MC %% 1 != 0 | n_MC <= 0)
        stop('Invalid # Monte Carlo simulations!')
      if (alpha <= 0 | alpha >= 1)
        stop('Alpha should be between 0 and 1.')
      if (!(alt %in% c('one.sided', 'two.sided', 'greater', 'less')))
        stop('Alternative hypothesis should be either \'one.sided\' or \'two.sided\'!')
      if (!(is.null(quantiles))) {
        if (length(quantiles) == 1) {
          if (!(is.numeric(quantiles) & quantiles <= 1 & quantiles >= 0))
            stop('Invalid quantile values')
        }
        else {
          for (i in 1:length(quantiles)) {
            if (!(is.numeric(quantiles[i]) & (quantiles[i] <= 1) & (quantiles[i] >= 0)))
              stop('Invalid quantile values')
          }
        }
      }
      else
        quantiles <- c(0, .25, .5, .75, 1)

      if (!(is.null(assurance_level_props))) {
        if (length(assurance_level_props) == 1) {
          if (!(is.numeric(assurance_level_props) & assurance_level_props <= 1 & assurance_level_props >= 0))
            stop('Invalid proportions for assurance levels')
        }
        else {
          for (i in 1:length(assurance_level_props)) {
            if (!(is.numeric(assurance_level_props[i]) & (assurance_level_props[i] <= 1) & (assurance_level_props[i] >= 0)))
              stop('Invalid proportions for assurance levels')
          }
        }
      }
      else
        assurance_level_props <- NULL

      self$parallel <- parallel
      self$cores <- cores
      self$ns <- sort(ns)
      self$n_prior <- n_prior
      self$n_MC <- n_MC
      self$prior <- prior
      self$alpha <- alpha
      self$alt <- alt
      self$quantiles <- quantiles
      self$assurance_level_props <- assurance_level_props
    },

    print = function(){
      cat('HybridPower Instance Description: \n\n')
      cat('Parallelize: ', self$parallel, '\n\n')
      cat('Sample sizes: ', self$ns, '\n')
      cat('Draws from prior: ', self$n_prior, '\n')
      cat('# Monte Carlo simulations: ', self$n_MC, '\n\n')
      cat('Type of prior: ', self$prior, '\n')
      cat('Alternative Hypothesis: ', self$alt, '\n')
      cat('Level of significance: ', self$alpha, '\n')
      cat('Proportions for assurance levels: ', self$assurance_level_props, '\n')
    },

    assurance = function() {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      return(summarise(group_by(self$output, n), assurance = mean(power), .groups='keep'))
    },

    power_quantiles = function(props=self$quantiles) {
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

    assurance_level = function(props=self$assurance_level_props) {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      if (is.null(props))
        stop('Provide target proportions')
      for (i in 1:length(props))
        if (!(is.numeric(props[i])) | props[i] > 1 | props[i] < 0)
          stop('Invalid proportion(s)')
      props <- sort(props)
      res <- summarise(group_by(self$output, n), private$compute_assurance_level(power, prop=props[1]), .groups='keep')
      if (length(props) > 1) {
        for (i in 2:length(props)) {
          res <- left_join(res, summarise(group_by(self$output, n), private$compute_assurance_level(power, prop=props[i]), .groups='keep'), by='n')
        }
      }
      col_names <- c('n', props)
      colnames(res) <- col_names
      return(res)
    },

    hybrid_power = function(cores=NULL) {
      if (is.null(self$prior)) {
        stop('Specify a prior first')
      }
      else {
        es <- private$draw_prior_es()
        if (self$parallel) {
          if (is.null(cores)) cl <- parallel::makeCluster(detectCores()-1)
          else cl <- parallel::makeCluster(cores)
          doParallel::registerDoParallel(cl)
          self$output <- parallel::parLapply(cl, self$ns, fun=private$generate_hybrid_power, es=es)
          private$melt_output()
          stopCluster(cl)
        }
        else {
          res <- list()
          for (i in 1:length(self$ns)) {
            res[[i]] <- private$generate_hybrid_power(self$ns[i], es=es)
          }
          self$output <- res
          private$melt_output()
        }
        cat('\nExample output:\n\n')
        print(head(self$output))
        cat('\n...\n')
        cat('\nFor the complete list of power values, access $output!\n')
      }
    },

    boxplot = function() {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      p <- ggplot(self$output, aes(x=factor(n), y=power)) + geom_boxplot()
      p <- p + xlab('Sample Size') + ylab('Power') + ggtitle('Distributions of Power')
      p <- p + stat_summary(fun=mean, geom='point', shape=5, size=4)
      p
    }
  ),

  private = list(
    melt_output = function() {
      output <- data.frame(self$output)
      colnames(output) = self$ns
      self$output <- suppressMessages(
        reshape2::melt(output, variable.name='n', value.name = 'power')
      )
    },

    compute_assurance_level = function(vec, prop) {
      return(sum(vec >= prop) / length(vec))
    }
  )
)
