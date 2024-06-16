## evaluate relationship between borealization and abundance
library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(mgcv)
library(ggpubr)

source("./scripts/stan_utils.R")

theme_set(theme_bw())


# reload data

# load borealization DFA trend

trend <- read.csv("./output/dfa_trend.csv")

# add rolling means and lags to figure borealization effects over differ time scales
trend <- trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, trend = estimate) %>%
  mutate(trend2 = zoo::rollmean(trend, 2, fill = NA, align = "right"),
         trend3 = zoo::rollmean(trend, 3, fill = NA, align= "right"),
         trend3_lag1 = lag(trend3, 1),
         trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1),
         trend2_lag2 = lag(trend2, 2))

abundance_male <- read.csv("./output/male_imputed_weighted_log_cpue.csv")

acf(abundance_male$log_mean[abundance_male$year %in% 1975:2019])

ggplot(abundance_male, aes(year, log_mean)) +
  geom_point() + 
  geom_line()

male_dat <- left_join(trend, abundance_male)
# dat <- left_join(trend, check)

# fill in log_mean_lag1 for 2021
male_dat$log_mean_lag1[male_dat$year == 2021] <- 
  mean(c(male_dat$log_mean_lag1[male_dat$year == 2020], male_dat$log_mean_lag1[male_dat$year == 2022]))

# and save
write.csv(male_dat, "male_dat.csv", row.names = F)

# exploratory models in mgcv
# using previous-year abundance and borealization index at different lags

male_boreal_mod1 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod1)

male_boreal_mod2 <- gam(log_mean ~ log_mean_lag1 + s(trend2) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod2)

male_boreal_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend3) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod3)

male_boreal_mod4 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod4)

male_boreal_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod5)

male_boreal_mod6 <- gam(log_mean ~  log_mean_lag1 + s(trend_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod6)

male_boreal_mod7 <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod7)
# plot(male_boreal_mod7, resid = T, pch = 19, se = T)

male.aicc <- MuMIn::AICc(male_boreal_mod1, male_boreal_mod2, male_boreal_mod3, male_boreal_mod4, 
                           male_boreal_mod5, male_boreal_mod6, male_boreal_mod7)

male.aicc

# fit best male model in brms---------

male_form <- bf(log_mean ~ log_mean_lag1 + s(trend3_lag1) + s(year, k = 4))

fit_male <- brm(male_form, 
                  data = male_dat,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.999, max_treedepth = 14),
                seed = 99)

saveRDS(fit_male, file = "output/fit_male.rds")

# diagnostics
male_brm <- readRDS("./output/fit_male.rds")
check_hmc_diagnostics(male_brm$fit)
neff_lowest(male_brm$fit)
rhat_highest(male_brm$fit)
summary(male_brm)
bayes_R2(male_brm)


## analyze female abundance -------------------

abundance_female <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year) %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

# save to plot 

write.csv(abundance_female, "./output/female_imputed_weighted_log_cpue.csv")


female_dat <- left_join(trend, abundance_female)

# fill in log_mean_lag1 for 2021
female_dat$log_mean_lag1[female_dat$year == 2021] <- 
  mean(c(female_dat$log_mean_lag1[female_dat$year == 2020], female_dat$log_mean_lag1[female_dat$year == 2022]))

female_boreal_mod1 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod1)

female_boreal_mod2 <- gam(log_mean ~ log_mean_lag1 + s(trend2) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod2)

female_boreal_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend3) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod3)

female_boreal_mod4 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod4)

female_boreal_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod5)

female_boreal_mod6 <- gam(log_mean ~  log_mean_lag1 + s(trend_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod6)

female_boreal_mod7 <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod7)
# plot(female_boreal_mod7, resid = T, pch = 19, se = T)

female.aicc <- MuMIn::AICc(female_boreal_mod1, female_boreal_mod2, female_boreal_mod3, female_boreal_mod4, 
                           female_boreal_mod5, female_boreal_mod6, female_boreal_mod7)

female.aicc

# fit best female model in brms---------

female_form <- bf(log_mean ~ log_mean_lag1 + s(trend3_lag1) + s(year, k = 4))

fit_female <- brm(female_form, 
                  data = female_dat,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.9999, max_treedepth = 14),
                  seed = 99)

saveRDS(fit_female, file = "output/fit_female.rds")

# diagnostics
female_brm <- readRDS("./output/fit_female.rds")
check_hmc_diagnostics(female_brm$fit)
neff_lowest(female_brm$fit)
rhat_highest(female_brm$fit)
summary(female_brm)
bayes_R2(female_brm)
