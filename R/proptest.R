source('R/HybridPower.R')

hp_prop <- R6Class(
  'hp_prop',
  inherit = HybridPower,
  public = list(
    hybrid_powers = NULL,
    design = NULL,
    n_MC = NULL,
    exact=FALSE,
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
      c = NULL,
      exact=FALSE,

      pi_1 = NULL,
      pi_2 = NULL,

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
      assurance_level_props = NULL,
      quantiles = NULL
    ) {
      if (!(is.null(prior))) {
        if (!(prior %in% c('beta', 'uniform', 'truncnorm')))
          stop('Invalid prior')
      }
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
      self$design <- design
      self$n_MC <- n_MC

      if (!(is.null(pi_1))) {
        if (!(is.numeric(pi_1)))
          stop('pi_1 is not numeric')
        if (pi_1 > 1 | pi_1 < 0)
          stop('Invalid value of pi_1')
        self$pi_1 <- pi_1
      }
      if (!(is.null(pi_2))) {
        if (!(is.numeric(pi_2)))
          stop('pi_2 is not numeric')
        if (pi_2 > 1 | pi_2 < 0)
          stop('Invalid value of pi_2')
        self$pi_2 <- pi_2
      }
      if (!(is.null(prior))) {
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
      }

      if (!(is.logical(exact)))
        stop('\'Exact\' must be logical')
      if (!(is.null(c))) {
        if (c >= 1 | c <= 0)
          stop('c should be between 0 and 1')
      }
      self$c <- c
      self$exact <- exact
    },

    print = function() {
      super$print()
      if (!(is.null(self$pi_1)))
        cat('Population proportion 1: ', self$pi_1, '\n')
      if (!(is.null(self$pi_2)))
        cat('Population proportion 2: ', self$pi_2, '\n')
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
      if (!(is.null(self$c)))
        cat('Population proportion under the null: ',self$c, '\n')
      cat('Study design: ', self$design, '\n')
    },

    classical_power = function(n=self$ns, pi_1=self$pi_1, pi_2 = self$pi_2, exact=T) {
      if (is.null(pi_2)) {
        mu_1 <- abs(pi_1 - self$c)
        sd_0 <- sqrt(self$c*(1-self$c)/n)
        sd_1 <- sqrt(pi_1*(1-pi_1)/n)
        if (self$alt == 'two.sided') {
          crit <- qnorm(1-self$alpha/2, 0, sd_0)
          return(pnorm(-crit, mu_1, sd_1) + 1 - pnorm(crit, mu_1, sd_1))
        }
        else {
          crit <- qnorm(1-self$alpha, 0, sd_0)
          return(1 - pnorm(crit, mu_1, sd_1))
        }
      }
      else {
        if (!(self$exact) & (!(is.null(pi_2)))) {
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
        else {
          if ((length(pi_1) != 1) | (length(pi_2) != 1))
            stop('Invalid pi_1 or pi_2 for a Fisher\'s exact test!')
          if (pi_1 > pi_2)
            alt <- 'less'
          else
            alt <- 'greater'
          return(sapply(n, private$sim_exact_test, n_MC = self$n_MC, pi_1=pi_1, pi_2=pi_2))
        }
      }
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

    classical_power2 = function(pi, n=self$ns, exact=FALSE) {
      pi_1 = pi[1]
      pi_2 = pi[2]
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
          return(pnorm(crit, mu_1, sd_1))
        }
        else {
          crit <- qnorm(1-self$alpha/2)
          return(1 - pnorm(crit, mu_1, sd_1))
        }
      }
      else {
        if (!(self$exact) & pi_2) {
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
        else if (self$exact & pi_2){
          if ((length(pi_1) != 1) | (length(pi_1) != 1))
            stop('Invalid pi_1 or pi_2 for a Fisher\'s exact test!')
          return(sapply(n, private$sim_exact_test, n_MC = self$n_MC, pi_1=pi_1, pi_2=pi_2))
        }
      }
    },

    is_significant = function(n, x) {
      return(
        fisher.test(
          x=matrix(c(x[2], n-x[2], x[1], n-x[1]), ncol=2),
          alt=self$alt,
          simulate.p.value = FALSE
        )$p.value < self$alpha
      )
    },

    sim_exact_test = function(n, n_MC, pi_1, pi_2) {
      X <- cbind(rbinom(n_MC, n, pi_1), rbinom(n_MC, n, pi_2))
      return(
        mean(apply(X, 1, private$is_significant, n=n))
      )
    },

    generate_hybrid_power = function(n) {
      es <- private$draw_prior_es()
      if (is.null(dim(es))) return(self$classical_power(n, pi_1=es, pi_2=NULL, exact=self$exact))
      else {
        return(
          apply(es, 1, private$classical_power2, n=n, exact=self$exact)
        )
      }
    }
  )
)
