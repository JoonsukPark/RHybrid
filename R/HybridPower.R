library(R6)
library(ggplot2)
library(dplyr)
library(reshape2)
library(parallel)

is_numeric <- function(x) return(is.numeric(x) & length(x) == 1)

HybridPower <- R6Class(
  'HybridPower',
  public = list(
    es = NULL,
    design = NULL,
    parallel = FALSE,
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
    assurance_props = NULL,

    initialize = function(
      parallel = FALSE,
      ns = c(),
      n_prior = 1,
      n_MC = 1,
      prior = NULL,
      alpha = 0.05,
      alt = 'two.sided',
      assurance_props = NULL
    ) {
      # Validate inputs
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
      if (!(is.null(prior))) {
        if (
          !is.character(prior) |
          length(prior) != 1 |
          !(prior %in% c('normal', 'uniform', 'beta', 'truncnorm', 'dirichlet'))
        )
          stop('Invalid prior!')
      }
      if (alpha <= 0 | alpha >= 1)
        stop('Alpha should be between 0 and 1.')
      if (alt != 'one.sided' & alt != 'two.sided')
        stop('Alternative hypothesis should be either \'one.sided\' or \'two.sided\'!')
      if (!(is.null(assurance_props))) {
        if (length(assurance_props) == 1) {
          if (!(is.numeric(assurance_props) & assurance_props <= 1 & assurance_props >= 0))
            stop('Invalid assurance level')
        }
        else {
          for (i in 1:length(assurance_props)) {
            if (!(is.numeric(assurance_props[i]) & (assurance_props[i] <= 1) & (assurance_props[i] >= 0)))
              stop('Invalid assurance proportions')
          }
        }
      }
      else
        assurance_props <- c(.25, .5, .75)

      self$parallel <- parallel
      self$ns <- sort(ns)
      self$n_prior <- n_prior
      self$n_MC <- n_MC
      self$prior <- prior
      self$alpha <- alpha
      self$alt <- alt
      self$assurance_props <- assurance_props
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

<<<<<<< HEAD
    hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        if (is.null(cores)) cores <- detectCores()
        self$output <- mclapply(self$ns, private$generate_hybrid_power)
        private$melt_output()
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$generate_hybrid_power(self$ns[i])
        }
        self$output <- res
        private$melt_output()
      }
      return(self$output)
    },

=======
>>>>>>> 642f19fa7a415fd4fafeeed7834d571453bd004b
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
    }
  )
)
