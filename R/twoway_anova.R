source('R/HybridPower.R')

hp_twoway_anova <- R6Class(
  'hp_twoway_anova',
  inherit = HybridPower,
  public = list(
    es = NULL,
    cellmeans = NULL,
    sd = 1,
    k = NULL,
    m = NULL,
    design = NULL,
    rho = NULL,
    epsilon = NULL,
    fe_factor = c(),
    rm_factor = c(),
    rm_dim = c(),

    initialize = function(
      parallel = FALSE,
      ns=c(),
      n_prior=1,
      n_MC=1,
      prior=NULL,
      alpha = 0.05,
      cellmeans=NULL,
      prior_mu = c(),
      prior_sigma = c(),
      prior_lower = c(),
      prior_upper = c(),
      sd = 1,
      design = 'fe',
      rho = NULL,
      epsilon = NULL,
      alt = 'one.sided',
      quantiles = NULL,
      assurance_level_props=NULL
    ) {
      super$initialize(
        parallel = FALSE,
        ns=ns,
        n_prior=n_prior,
        n_MC=n_MC,
        alpha=alpha,
        alt=alt,
        quantiles=quantiles,
        assurance_level_props=assurance_level_props
      )
      self$prior <- prior
      self$fe_factor <- which(design == 'fe')
      self$rm_factor <- which(design == 'rm')
      if (!(is.null(cellmeans))) {
        if (length(dim(cellmeans)) != 2)
          stop('Input a 2-dimensional cell means matrix')
        dims <- dim(cellmeans)
        self$k <- dims[self$fe_factor]
        self$m <- dims[self$rm_factor]
        for(i in 1:dims[1]) {
          for(j in 1:dims[2]) {
            if (!(is.numeric(cellmeans[i][j])))
              stop('Invalid cell means')
          }
        }
        self$cellmeans <- cellmeans
      }

      if (length(design) != 2)
        stop('The vector of designs should have a length of 2')
      if (sum(design == 'rm') == 2)
        stop('At least one of the factors should be fixed effect')
      self$rm_dim <- sum(design == 'rm')
      for (i in 1:2) {
        if (!(design[i] %in% c('fe', 'rm')))
          stop('Design should be one of \'fe\' or \'rm\'!')
      }
      if (!(is.null(prior))) {
        self$prior <- prior
        if (prior == 'normal') {
          if (!(is.matrix(prior_mu)) | !(length(dim(prior_mu)) == 2))
            stop('Input a 2-dimensional matrix for prior means!')
          if (sum(dim(prior_mu) == 1) > 0)
            stop('At least 1 of the IVs has a single level')
          if (!(is.matrix(prior_sigma)) | !(length(dim(prior_sigma)) == 2))
            stop('Input a 2-dimensional matrix for prior standard deviations!')
          if (sum(dim(prior_sigma) == 1) > 0)
            stop('At least 1 of the IVs has a single level')
          if (sum(prior_sigma <= 0) > 0)
            stop('Prior standard deviations must be strictly positive')
          if (sum((dim(prior_mu) == dim(prior_sigma)) == FALSE) > 0)
            stop('Dimensions of prior means and sds should be identical')
        }
        else if (prior == 'uniform') {
          if (!(is.matrix(prior_lower)) | !(length(dim(prior_lower)) == 2))
            stop('Input a 2-dimensional matrix for prior lower bounds!')
          if (sum(dim(prior_lower) == 1) > 0)
            stop('At least 1 of the IVs has a single level')
          if (!(is.matrix(prior_upper)) | !(length(dim(prior_upper)) == 2))
            stop('Input a 2-dimensional matrix for prior upper bounds!')
          if (sum(dim(prior_upper) == 1) > 0)
            stop('At least 1 of the IVs has a single level')
          if (sum(prior_lower > prior_upper) > 0)
            stop('At least one lower bound is greater than the upper bound')
          if (sum((dim(prior_lower) == dim(prior_upper)) == FALSE) > 0)
            stop('Dimensions of prior lower/upper bounds must be identical')
        }
      }

      if (!(design[1] == 'fe' & design[2] == 'fe')) {
        if (length(rho) != self$rm_dim)
          stop('length of rho should match that of repeated measure factors')
        if (length(epsilon) != self$rm_dim)
          stop('length of epsilon should match that of repeated measure factors')
      }
      if (self$rm_dim > 0) {
        for (i in 1:self$rm_dim)
          if (rho[i] >= 1 | rho[i] < 0)
            stop('rho should be between 0 and 1!')
        if (epsilon[i] > 1 | epsilon[i] < 0)
          stop('epsilon should be between 0 and 1!')
      }
      if (sd <=0) stop('Standard deviation must be strictly positive')

      self$design <- design
      self$rho <- rho
      self$sd <- sd
      self$epsilon <- epsilon
      if (!(is.null(prior))) {
        if (prior == 'normal') {
          self$prior_mu <- prior_mu
          self$prior_sigma <- prior_sigma
          dims <- dim(prior_mu)
          self$k <- dims[self$fe_factor]
          self$m <- dims[self$rm_factor]
        }
        else if (prior == 'uniform') {
          self$prior_lower <- prior_lower
          self$prior_upper <- prior_upper
          dims <- dim(prior_lower)
          self$k <- dims[self$fe_factor]
          self$m <- dims[self$rm_factor]
        }
      }
    },

    print = function() {
      super$print()
      if (!(is.null(self$cellmeans))) {
        cat('Cell means: ', self$cellmeans, '\n\n')
      }
      if (!(is.null(self$prior))) {
        if (self$prior == 'normal') {
          cat('\nPrior means: \n\n')
          cat(self$prior_mu[1,], '\n')
          cat(self$prior_mu[2,], '\n\n')
          cat('Prior standard deviations: \n\n')
          cat(self$prior_sigma[1,], '\n')
          cat(self$prior_sigma[2,], '\n')
          cat('Implied f values: ',
              private$compute_f(self$prior_mu, 1), ' ',
              private$compute_f(self$prior_mu, 2), '\n\n'
          )
        }
        else if (self$prior == 'uniform') {
          cat('\nPrior lower bounds: \n\n')
          cat(self$prior_lower[1,], '\n')
          cat(self$prior_lower[2,], '\n\n')
          cat('Prior upper bounds: \n\n')
          cat(self$prior_upper[1,], '\n')
          cat(self$prior_upper[2,], '\n\n')
          cat('Implied f values: ',
              private$compute_f((self$prior_lower + self$prior_upper)/2, 1), ' ',
              private$compute_f((self$prior_lower + self$prior_upper)/2, 2), '\n\n'
          )
        }
      }
      cat('Test type: Two-way ANOVA\n')
      cat('Study design: ', self$design, '\n')
      if (!(is.null(self$rho)))
        cat('rho: ', self$rho, '\n')
      if (!(is.null(self$epsilon)))
        cat('epsilon: ', self$epsilon, '\n')
    },

    classical_power = function(means=self$cellmeans, n=self$ns) {
      powers <- list()
      if (is.null(means)) {
        if (self$prior == 'normal')
          means <- self$prior_mu
        else
          means <- (self$prior_lower + self$prior_upper) / 2
      }
      if (self$design[1] != self$design[2]) {
        for (dim in 1:2) {
          if (self$design[dim] == 'fe') {
            f <- private$compute_f(means, dim)
            u <- self$m / (1 + (self$m-1)*self$rho)
            ncp <- f^2*n*u
            df1 <- self$k-1
            df2 <- n - self$k
            powers[[dim]] <- private$compute_f_prob(f, ncp, df1, df2)
          }
          else {
            f <- private$compute_f(means, dim)
            u <- self$m / (1 - self$rho)
            ncp <- f^2*n*u*self$epsilon
            df1 <- (self$m-1)*self$epsilon
            df2 <- (n - self$k)*(self$m-1)*self$epsilon
            powers[[dim]] <- private$compute_f_prob(f, ncp, df1, df2)
          }
        }
        f <- private$compute_f_mat(means)
        u <- self$m / (1 - self$rho)
        ncp <- f^2*n*u*self$epsilon
        df1 <- (self$k-1)*(self$m-1)*self$epsilon
        df2 <- (n - self$k)*(self$m-1)*self$epsilon
        powers[[3]] <- private$compute_f_prob(f, ncp, df1, df2)
      }
      else {
        for (dim in 1:2) {
          f <- private$compute_f(means, dim)
          ncp <- f^2*n
          df1 <- self$k[dim]-1
          df2 <- n - self$k[1]*self$k[2]
          powers[[dim]] <- private$compute_f_prob(f, ncp, df1, df2)
        }
        f <- private$compute_f_mat(means)
        ncp <- f^2*n
        df1 <- (self$k[1]-1)*(self$k[2]-1)
        df2 <- n - self$k[1]*self$k[2]
        powers[[3]] <- private$compute_f_prob(f, ncp, df1, df2)
      }
      return(powers)
    },

    hybrid_power = function(cores=NULL) {
      if (self$parallel) {
        library(parallel)
        if (!(cores)) cores <- detectCores()
        self$output <- melt(
          private$gather_powers(
            mclapply(
              self$ns, private$generate_hybrid_power)
          ),
          id.vars = 'n',
          variable.name = 'type',
          value.name = 'power'
        )
        return(self$output)
      }
      else {
        res <- list()
        for (i in 1:length(self$ns)) {
          res[[i]] <- private$generate_hybrid_power(self$ns[i])
        }
        self$output <- melt(private$gather_powers(res),
          id.vars = 'n',
          variable.name = 'type',
          value.name = 'power'
        )
        return(self$output)
      }
    },

    boxplot = function(power_df=self$output, factors='all') {
      if (sum(factors == 'all') > 0)
        p <- ggplot(power_df, aes(x=factor(n), y=power, color=type))
      else {
        tryCatch(
          {
            power_df <- power_df %>% filter(type %in% factors)
            p <- ggplot(power_df, aes(x=factor(n), y=power, color=type))
          }, error = function(e) stop('Invalid IV type!')
        )
      }
      p <- p + geom_boxplot()
      p <- p + xlab('Sample Size') + ylab('Power') + ggtitle('Distributions of Power')
      p
    },

    assurance = function() {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      return(summarise(group_by(self$output, n, type), assurance = mean(power), .groups='keep'))
    },

    assurance_level = function(props=self$quantiles) {
      if (is.null(self$output))
        stop('Run hybrid_power() first')
      if (is.null(props))
        stop('Provide target proportions')
      for (i in 1:length(props))
        if (!(is.numeric(props[i])) | props[i] > 1 | props[i] < 0)
          stop('Invalid proportion(s)')
      props <- sort(props)
      res <- summarise(group_by(self$output, n, type), quantile(power, probs=props[1]), .groups='keep')
      if (length(props) > 1) {
        for (i in 2:length(props)) {
          res <- left_join(res, summarise(group_by(self$output, n, type), quantile(power, probs=props[i]), .groups='keep'), by=c('n', 'type'))
        }
      }
      col_names <- c('n', 'type', props)
      colnames(res) <- col_names
      return(res)
    }
  ),

  private = list(

    compute_f = function(means, dim) {
      means <- apply(means, dim, mean)
      n_means <- length(means)
      return(sqrt(var(means)*(n_means-1)/n_means)/self$sd)
    },

    compute_f_mat = function(means) {
      means <- as.vector(means)
      n_means <- length(means)
      return(sqrt(var(means)*(n_means-1)/n_means)/self$sd)
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
      means <- list()
      if (self$prior == 'normal') {
        prior_mu <- as.vector(self$prior_mu)
        prior_sigma <- as.vector(self$prior_sigma)
        if (length(self$k) == 1)
          len <- self$k*self$m
        else len <- prod(self$k)
        means <- list()
        for (i in 1:self$n_prior) {
          means[[i]] <- matrix(rnorm(len, prior_mu, prior_sigma), ncol = ncol(self$prior_mu))
        }
      }
      else if (self$prior == 'uniform') {
        prior_lower <- as.vector(self$prior_lower)
        prior_upper <- as.vector(self$prior_upper)
        len <- self$k*self$m
        means <- list()
        for (i in 1:self$n_prior) {
          means[[i]] <- matrix(runif(len, prior_lower, prior_upper), ncol = ncol(self$prior_mu))
        }
      }
      return(means)
    },

    generate_hybrid_power = function(n) {
      means <- private$draw_prior_es()
      return(
        matrix(
          unlist(
            t(sapply(means, self$classical_power, n=n))
          ),
          nrow=self$n_prior
        )
      )
    },

    gather_powers = function(power_list) {
      powers <- data.frame()
      for (i in 1:length(self$ns)) {
        power_df <- cbind(self$ns[i], power_list[[i]])
        powers <- rbind(powers, power_df)
      }
      if (self$design[1] != self$design[2])
        colnames(powers) = c('n', self$design, 'interaction')
      else
        colnames(powers) = c('n', 'fe1', 'fe2', 'interaction')
      return(powers)
    },

    gather_assurances = function(assurance_list) {
      res <- data.frame()
      for (i in 1:length(assurance_list)) {
        res <- rbind(res, c(self$ns[i], assurance_list[[i]]))
      }
      colnames(res) = c('n', self$design, 'interaction')
      return(res)
    }
  )
)
