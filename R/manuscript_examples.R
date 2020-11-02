# library(hybridpower)
# SEED <- 1234
#
# # Case 1: independent samples t-test example
# set.seed(SEED)
# x_classical <- hp_ttest$new(
#   ns = seq(10, 90, 10),
#   design = 'two.sample',
#   d = 0.5
# )
# x_classical$classical_power()
#
# # Bayesian-classical hybrid power
# set.seed(SEED)
# x_hybrid <- hp_ttest$new(
#   parallel=T,
#   ns = seq(10, 90, 10),
#   design = 'two.sample',
#   n_prior = 1000,
#   prior = 'normal',
#   prior_mu = 0.5,
#   prior_sigma = 0.1
# )
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
# x_hybrid$power_quantiles()
# x_hybrid$boxplot()
#
# # Both can be calculated using a single instance
# set.seed(SEED)
# x_both <- hp_ttest$new(
#   parallel=T,
#   ns = seq(10, 90, 10),
#   d = 0.5,
#   design = 'two.sample',
#   n_prior = 1000,
#   prior = 'normal',
#   prior_mu = 0.5,
#   prior_sigma = 0.1,
#   assurance_level_props = c(.5, .8),
#   quantiles = c(0, .2, .5, .7, 1) #override defaults
# )
#
# x_both$hybrid_power()
# x_both$assurance()
# x_both$power_quantiles()
# x_both$boxplot()
# x_both$assurance_level()
#
# ## Case 2: sign test example
# set.seed(SEED)
# x_classical <- hp_sign$new(
#   ns = seq(50, 150, 10),
#   p_0 = 0.5,
#   p_1 = 0.3
# )
# x_classical$classical_power()
#
# set.seed(SEED)
# x_hybrid <- hp_sign$new(
#   parallel=T,
#   ns = seq(50, 150, 10),
#   p_0 = 0.5,
#   prior = 'uniform',
#   prior_lower = 0.2,
#   prior_upper = 0.4,
#   n_prior = 5000
# )
# x_hybrid$hybrid_power()
# x_hybrid$assurance()
#
# # Case 3: Welch ANOVA example
# set.seed(SEED)
# x_classical_welch <- hp_oneway_anova$new(
#   parallel=T,
#   ns = seq(40, 120, 20),
#   mu = c(2, 2.5, 2.2),
#   sigma = c(1, 1.1, 1.2),
#   design = 'fe',
#   n_MC=1000
# )
# x_classical_welch$classical_power()
#
# set.seed(SEED)
# x_hybrid_welch <- hp_oneway_anova$new(
#   parallel=T,
#   ns = seq(40, 120, 20),
#   prior_mu = c(2, 2.5, 2.2),
#   prior_sigma = c(.1, .2, .15),
#   sigma = c(1.0, 1.1, 1.2),
#   design = 'fe',
#   prior = 'normal',
#   n_prior = 1000,
#   n_MC = 1000
# )
# begin <- Sys.time()
# x_hybrid_welch$hybrid_power()
# end <- Sys.time()
# print(end-begin)
# x_hybrid_welch$assurance()
# pdf('fig3.pdf', width=6, height=4)
# x_hybrid_welch$boxplot()
# dev.off()
