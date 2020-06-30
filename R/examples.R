# library(hybridpower)
#
# #######################
# ## Example 1: t-test ##
# #######################
#
# # Classical power analysis
# x_classical <- hp_ttest$new(
#   ns = seq(10, 90, 10),
#   d = 0.5,
#   alpha=.01
# )
# x_classical$classical_power()
#
# # This should fail
# x_classical$assurance_level()
#
# # B-C hybrid power analysis
#
# # Note that quantiles is now another input to the initializer (default is 0.5, which is the median)
# # Also note that you can change the significance level (default is still alpha=0.05)
# # Now the default quantiles are those from Tukey (0, .25, .5, .75, 1)
# x_hybrid <- hp_ttest$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior = 'normal',
#   prior_mu = 0.3,
#   prior_sigma = 0.1,
#   alpha=0.1,
#   assurance_level_props = c(.5, .8)
# )
#
# # This should generate an error because you haven't specified the 'd' here.
# x_hybrid$classical_power()
#
# # Now it will be OK
# x_hybrid <- hp_ttest$new(
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior = 'normal',
#   prior_mu = 0.3,
#   prior_sigma = 0.1,
#   alpha=0.1,
#   d = 0.1,
#   assurance_level_props = c(.5, .8)
# )
# x_hybrid$classical_power()
#
# # This saves the generated hybrid power values at self$output (in this case, x_hybrid$output)
# x_hybrid$hybrid_power()
#
# # You can also do this because hybrid_power() saves the power values as well as returning them as output values
# powers <- x_hybrid$hybrid_power()
#
# # You can retrieve the saved power values like this
# x_hybrid$output
#
# # These must be run after running hybrid_power()
# # because they assume that power values are already saved in self$output
# # Otherwise, they would raise an error.
#
# x_hybrid$assurance()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$assurance_level()
#
# # Also you can change the internal variables whenever you want
# x_hybrid$alpha <- 0.05
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$output
# x_hybrid$assurance()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$assurance_level()
#
#
# ############################
# ## Example 2: 1-way ANOVA ##
# ############################
#
# # Classical power analysis
# x_classical <- hp_oneway_anova$new(
#   ns = seq(10, 90, 10),
#   mu = c(2, 2.2),
#   sd = 1,
#   design='fe'
# )
# x_classical$classical_power()
#
# # If you provide the value of f^2, it will compute power based on that
# # Otherwise, it will compute power based on the provided cell means
# x_classical$classical_power(f2=0.5)
#
# # B-C hybrid power analysis
# x_hybrid <- hp_oneway_anova$new(
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
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$boxplot()
#
# # Now it should be OK
# x_hybrid <- hp_oneway_anova$new(
#   ns = seq(10, 90, 10),
#   mu = c(2, 2.2),
#   prior_mu = c(2, 2.5),
#   prior_sigma = c(0, .2),
#   sd = 1,
#   design='fe',
#   prior = 'normal',
#   n_prior = 1000,
#   quantiles = c(.2, .5, .8),
#   assurance_level_props = c(.5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$power_quantiles()
# x_hybrid$boxplot()
# x_hybrid$assurance_level()
#
# ############################
# ## Example 3: 2-way ANOVA ##
# ############################
#
# x_classical <- hp_twoway_anova$new(
#   ns = c(10, 20, 30, 40),
#   design = c('fe', 'rm'),
#   cellmeans = matrix(c(2, 2.2, 2, 1.8, 2, 2.4), nrow=3),
#   sd = 1,
#   rho = 0.2,
#   epsilon=1
# )
# x_classical$classical_power()
#
# x_hybrid <- hp_twoway_anova$new(
#   parallel = T,
#   ns = c(10, 20, 30, 40),
#   n_prior=10,
#   design = c('fe', 'rm'),
#   prior = 'normal',
#   prior_mu = matrix(c(2, 2.2, 2, 1.8, 2, 2.4), nrow=3),
#   prior_sigma = matrix(rep(.5, 6), nrow=3),
#   sd = 2,
#   rho = 0.2,
#   epsilon=1,
#   quantiles = c(.2, .5, .8)
# )
# x_hybrid$hybrid_power()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$assurance_level()
#
# #####################################
# # Example 4: bivariate correlation ##
# #####################################
#
# # Classical power analysis (minimal example)
#
# x_classical <- hp_cor$new(
#   ns = seq(10, 90, 10),
#   rho = .7,
#   alt = 'two.sided'
# )
# x_classical$classical_power()
#
# # Hybrid power analysis
#
# x_hybrid <- hp_cor$new(
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   rho = .5,
#   prior_mu = .3,
#   prior_sigma = .1,
#   prior = 'truncnorm',
#   alt = 'two.sided',
#   assurance_level_props = c(.5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$boxplot()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$assurance_level()
#
# #################################
# ## Example 5: Proportions test ##
# #################################
#
# x <- hp_prop$new(
#   parallel=T,
#   ns = seq(30, 90, 10),
#   n_prior=10,
#   prior = 'truncnorm',
#   prior_pi_1_mu = .6,
#   prior_pi_1_sd = .1,
#   prior_pi_2_mu = .5,
#   prior_pi_2_sd = .1,
#   c = 0.5,
#   n_MC = 10,
#   alt = 'two.sided',
#   exact=T,
#   pi_1 = 0.5,
#   pi_2 = 0.7,
#   assurance_level_props = c(.5, .8)
# )
#
# x$classical_power()
# x$hybrid_power()
# x$assurance()
# x$boxplot()
# x$assurance_level()
#
# ##########################
# ## Example 6: Sign test ##
# ##########################
#
# x_classical <- hp_sign$new(
#   ns = seq(10, 90, 10),
#   p_0 = 0.5,
#   p_1 = 0.6
# )
# x_classical$classical_power()
#
# x_hybrid <- hp_sign$new(
#   prior='uniform',
#   prior_lower = 0.1,
#   prior_upper = 0.3,
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   n_MC = 100,
#   p_0 = 0.5,
#   p_1 = 0.6,
#   MC=F,
#   quantiles = c(.2, .5, .8),
#   assurance_level_props = c(.5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
# x_hybrid$assurance_level()
#
# #########################################
# ## Example 7: Simple linear regression ##
# #########################################
#
# x_classical <- hp_slr$new(
#   ns = seq(10, 90, 10),
#   r2 = 0.1
# )
# x_classical$classical_power()
#
#
# x_hybrid <- hp_slr$new(
#   ns = seq(10, 90, 10),
#   r2 = 0.1,
#   n_prior=1000,
#   prior = 'truncnorm',
#   prior_mu = .2,
#   prior_sigma = .2,
#   alt = 'two.sided',
#   assurance_level_props = c(.5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
# x_hybrid$assurance_level()
#
# #######################################
# ## Example 8: Chi-square test of GOF ##
# #######################################
#
# x_classical <- hp_chisq$new(
#   ns = seq(10, 90, 10),
#   p_0 = c(1/3, 1/3, 1/3),
#   p_1 = c(1/4, 1/4, 1/2)
# )
# x_classical$classical_power()
#
# # Hybrid power analysis
# x_hybrid <- hp_chisq$new(
#   prior='beta',
#   parallel = T,
#   ns = seq(10, 90, 10),
#   n_prior=1000,
#   prior_a = 1,
#   prior_b = 2,
#   p_0 = c(1/2, 1/2),
#   p_1 = c(1/3, 2/3),
#   assurance_level_props = c(.5, .8)
# )
#
# x_hybrid$classical_power()
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$boxplot()
# x_hybrid$power_quantiles()
# x_hybrid$assurance_level()
