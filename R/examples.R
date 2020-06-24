# library(hybridpower)
#
# #######################
# ## Example 1: t-test ##
# #######################
#
# # Classical power analysis
# power_classical <- HybridPowerTtest$new(
#   ns = seq(10, 90, 10),
#   d = 0.5,
#   alpha=.01
# )
# power_classical$classical_power()
#
# # B-C hybrid power analysis
#
# # Note that assurance_props is now another input to the initializer (default is 0.5, which is the median)
# # Also note that you can change the significance level (default is still alpha=0.05)
# power_hybrid <- HybridPowerTtest$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior = 'normal',
#   prior_mu = 0.3,
#   prior_sigma = 0.1,
#   alpha=0.1,
#   assurance_props = c(.5, .75)
# )
#
# # This should generate an error because you haven't specified the 'd' here.
# power_hybrid$classical_power()
#
# # Now it will be OK
# power_hybrid <- HybridPowerTtest$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior = 'normal',
#   prior_mu = 0.3,
#   prior_sigma = 0.1,
#   alpha=0.1,
#   assurance_props = c(.5, .75),
#   d = 0.1
# )
# power_hybrid$classical_power()
#
# # This saves the generated hybrid power values at self$output (in this case, power_hybrid$output)
# power_hybrid$hybrid_power()
#
# # You can also do this because hybrid_power() saves the power values as well as returning them as output values
# powers <- power_hybrid$hybrid_power()
#
# # You can retrieve the saved power values like this
# power_hybrid$output
#
# # These must be run after running hybrid_power()
# # because they assume that power values are already saved in self$output
# # Otherwise, they would raise an error.
#
# power_hybrid$assurance()
# power_hybrid$boxplot()
# power_hybrid$assurance()
# power_hybrid$assurance_level()
#
# # Also you can change the internal variables whenever you want
# power_hybrid$alpha <- 0.05
# power_hybrid$classical_power()
# power_hybrid$hybrid_power()
# power_hybrid$output
# power_hybrid$assurance()
# power_hybrid$boxplot()
# power_hybrid$assurance()
# power_hybrid$assurance_level()
#
#
# ############################
# ## Example 2: 1-way ANOVA ##
# ############################
#
# # Classical power analysis
# power_classical <- HybridPowerOnewayANOVA$new(
#   ns = seq(10, 90, 10),
#   mu = c(2, 2.2),
#   sd = 1,
#   design='fe'
# )
# power_classical$classical_power()
#
# # If you provide the value of f^2, it will compute power based on that
# # Otherwise, it will compute power based on the provided cell means
# power_classical$classical_power(f2=0.5)
#
# # B-C hybrid power analysis
# power_hybrid <- HybridPowerOnewayANOVA$new(
#   ns = seq(10, 90, 10),
#   prior_mu = c(2, 2.5),
#   prior_sigma = c(0, .2),
#   sd = 1,
#   design='fe',
#   prior = 'normal',
#   n_prior = 1000
# )
#
# # This should generate an error
# power_hybrid$classical_power()
# power_hybrid$hybrid_power()
# power_hybrid$assurance()
# power_hybrid$assurance_level()
# power_hybrid$boxplot()
#
# # Now it should be OK
# power_both <- HybridPowerOnewayANOVA$new(
#   ns = seq(10, 90, 10),
#   mu = c(2, 2.2),
#   prior_mu = c(2, 2.5),
#   prior_sigma = c(0, .2),
#   sd = 1,
#   design='fe',
#   prior = 'normal',
#   n_prior = 1000,
#   assurance_props = c(.2, .5, .8)
# )
#
# power_both$classical_power()
# power_both$hybrid_power()
# power_both$assurance_level()
# power_both$boxplot()
#
# ############################
# ## Example 3: 2-way ANOVA ##
# ############################
#
# x_classical <- HybridPowerTwowayANOVA$new(
#   ns = c(10, 20, 30, 40),
#   design = c('fe', 'rm'),
#   cellmeans = matrix(c(2, 2.2, 2, 1.8, 2, 2.4), nrow=3),
#   sd = 1,
#   rho = 0.2,
#   epsilon=1
# )
# x_classical$classical_power()
#
# x_hybrid <- HybridPowerTwowayANOVA$new(
#   parallel = T,
#   ns = c(10, 20, 30, 40),
#   n_prior=10,
#   design = c('fe', 'rm'),
#   prior_mu = matrix(c(2, 2.2, 2, 1.8, 2, 2.4), nrow=3),
#   prior_sigma = matrix(rep(.5, 6), nrow=3),
#   sd = 2,
#   rho = 0.2,
#   epsilon=1,
#   assurance_props = c(.2, .5, .8)
# )
# x_hybrid$hybrid_power()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$assurance_level()
# Classical power analysis
#
# ######################################
# ## Example 4: bivariate correlation ##
# ######################################
#
# Classical power analysis (minimal example)
#
# x_classical <- HybridPowerCorrelation$new(
#   ns = seq(10, 90, 10),
#   rho = .7,
#   alt = 'two.sided'
# )
#
# x_classical$classical_power()
#
# # Hybrid power analysis
#
# x_hybrid <- HybridPowerCorrelation$new(
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   rho = .5,
#   prior_mu = .3,
#   prior_sd = .1,
#   prior = 'truncnorm',
#   alt = 'two.sided'
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$assurance_level()
#
# # Classical power analysis (minimal example)
#
# x_classical <- HybridPowerCorrelation$new(
#   ns = seq(10, 90, 10),
#   rho = .7,
#   alt = 'two.sided'
# )
#
# x_classical$classical_power()
#
# # Hybrid power analysis
#
# x_hybrid <- HybridPowerCorrelation$new(
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   rho = .5,
#   prior_mu = .3,
#   prior_sd = .1,
#   prior = 'truncnorm',
#   alt = 'two.sided',
#   assurance_props = seq(.2, .5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$assurance_level()
#
# #################################
# ## Example 5: Proportions test ##
# #################################
#
# x <- HybridPowerProp$new(
#   parallel=T,
#   ns = seq(30, 90, 10),
#   n_prior=100,
#   prior = 'truncnorm',
#   prior_pi_1_mu = .6,
#   prior_pi_1_sd = .1,
#   prior_pi_2_mu = .5,
#   prior_pi_2_sd = .1,
#   c = 0.5,
#   n_MC = 100,
#   alt = 'two.sided',
#   exact=T,
#   pi_1 = 0.5,
#   pi_2 = 0.7
# )
#
# x$classical_power()
# x$hybrid_power()
# x$assurance()
# x$boxplot()
#
# ##########################
# ## Example 6: Sign test ##
# ##########################
#
# x_classical <- HybridPowerSignTest$new(
#   ns = seq(10, 90, 10),
#   p_0 = 0.5,
#   p_1 = 0.6
# )
# x_classical$classical_power()
#
# x_hybrid <- HybridPowerSignTest$new(
#   prior='uniform',
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   n_MC = 100,
#   p_0 = 0.5,
#   p_1 = 0.6,
#   MC=F,
#   assurance_props = c(.2, .5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
#
# #########################################
# ## Example 7: Simple linear regression ##
# #########################################
#
# x_classical <- HybridPowerSLR$new(
#   ns = seq(10, 90, 10),
#   r2 = 0.1
# )
# x_classical$classical_power()
#
#
# x_hybrid <- HybridPowerSLR$new(
#   ns = seq(10, 90, 10),
#   r2 = 0.1,
#   n_prior=1000,
#   prior = 'truncnorm',
#   prior_mu = .2,
#   prior_sd = .2,
#   alt = 'two.sided'
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
#
# #######################################
# ## Example 8: Chi-square test of GOF ##
# #######################################
#
# x_classical <- HybridPowerChisqTest$new(
#   ns = seq(10, 90, 10),
#   p_0 = c(1/3, 1/3, 1/3),
#   p_1 = c(1/4, 1/4, 1/2)
# )
# x_classical$classical_power()
#
# # Hybrid power analysis
# x_hybrid <- HybridPowerChisqTest$new(
#   prior='beta',
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior_alpha = 1,
#   prior_beta = 2,
#   p_0 = c(1/2, 1/2),
#   p_1 = c(1/3, 2/3)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
