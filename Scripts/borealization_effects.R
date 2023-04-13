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
         trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1),
         trend2_lag2 = lag(trend2, 2))

# load male snow crab abundance and process
abundance_male <- read.csv("./output/male_30-95_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year) %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

acf(abundance_male$log_mean[abundance_male$year %in% 1975:2019])

ggplot(abundance_male, aes(year, log_mean)) +
  geom_point() + 
  geom_line()

dat <- left_join(trend, abundance_male)

# fill in log_mean_lag1 for 2021
dat$log_mean_lag1[dat$year == 2021] <- mean(dat$log_mean_lag1[dat$year == 2020], dat$log_mean_lag1[dat$year == 2022])

# exploratory models in mgcv
# using previous-year abundance and borealization index at different lags

male_mod1 <- gam(log_mean ~ s(log_mean_lag1, k = 4), data = dat)
summary(male_mod1)
plot(male_mod1, se = F, resid = T, pch = 19)

# previous-year abundance can be a linear term for simplicity

male_mod2 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1), data = dat)
summary(male_mod2)
plot(male_mod2, se = F, resid = T, pch = 19)

male_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend2), data = dat)
summary(male_mod3)

male_mod4 <- gam(log_mean ~ log_mean_lag1 + s(trend3), data = dat)
summary(male_mod4)

male_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2), data = dat)
summary(male_mod5)

MuMIn::AICc(male_mod1, male_mod2, male_mod3, male_mod4, male_mod5) # model 2 is best

# add smoothed year term to best male model
male_mod2a <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year), data = dat)
summary(male_mod2a)
plot(male_mod2a, se = T, resid = T, pch = 19)

# fit best male model in brms---------

male_form <- bf(log_mean ~ log_mean_lag1 + s(trend2_lag1) + s(year, k = 4))

fit_male <- brm(male_form, 
                  data = dat,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(fit_male, file = "output/fit_male.rds")

# diagnostics
male_brm <- readRDS("./output/fit_male.rds")
check_hmc_diagnostics(male_brm$fit)
neff_lowest(male_brm$fit)
rhat_highest(male_brm$fit)
summary(male_brm)
bayes_R2(male_brm)

check <- plot(conditional_effects(male_brm), ask = FALSE, points = T)

# get vector of points conditioned on year and lagged abundance effects
new.dat <- data.frame(year = c(1976:2019, 2021, 2022),
                      log_mean_lag1 = dat$log_mean_lag1[dat$year %in% c(1976:2019, 2021, 2022)],
                      trend2_lag1 = mean(dat$trend2_lag1[dat$year %in% c(1976:2019, 2021, 2022)]))

pred <- posterior_epred(male_brm, newdata = new.dat)

points_dat <- data.frame(year = c(1976:2019, 2021, 2022),
                         trend2_lag1 = dat$trend2_lag1[dat$year %in% c(1976:2019, 2021, 2022)],
                         log_mean = colMeans(pred))


# plot male_brm

## 95% CI
ce1s_1 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]

male_brm_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (lag 1-2)", y = "Log CPUE") 

print(male_brm_plot)

ggsave("./figs/male_brm_borealization_plot.png", width = 6, height = 4, units = 'in')

## analyze female abundance -------------------

abundance_female <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year) %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

dat <- left_join(trend, abundance_female)

# fill in log_mean_lag1 for 2021
dat$log_mean_lag1[dat$year == 2021] <- mean(dat$log_mean_lag1[dat$year == 2020], dat$log_mean_lag1[dat$year == 2022])

female_mod1 <- gam(log_mean ~ s(log_mean_lag1, k = 4), data = dat)
summary(female_mod1)
plot(female_mod1, se = F, resid = T, pch = 19)

female_mod2 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1), data = dat)
summary(female_mod2)
plot(female_mod2, se = F, resid = T, pch = 19)

female_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend2), data = dat)
summary(female_mod3)

female_mod4 <- gam(log_mean ~ log_mean_lag1 + s(trend3), data = dat)
summary(female_mod4)
plot(female_mod4, se = F, resid = T, pch = 19)

female_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2), data = dat)
summary(female_mod5)

MuMIn::AICc(female_mod1, female_mod2, female_mod3, female_mod4, female_mod5)

# add smoothed year effect to best female model
female_mod4a <- gam(log_mean ~ log_mean_lag1 + s(trend3) + s(year), data = dat)
summary(female_mod4a)
plot(female_mod4a, se = F, resid = T, pch = 19)

# fit best female model in brms---------

female_form <- bf(log_mean ~ log_mean_lag1 + s(trend3) + s(year, k = 4))

fit_female <- brm(female_form, 
                  data = dat,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(fit_female, file = "output/fit_female.rds")

# diagnostics
female_brm <- readRDS("./output/fit_female.rds")
check_hmc_diagnostics(female_brm$fit)
neff_lowest(female_brm$fit)
rhat_highest(female_brm$fit)
summary(female_brm)
bayes_R2(female_brm)

plot(conditional_effects(female_brm), ask = FALSE, points = T)

# plot female_brm

## 95% CI
ce1s_1 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend3
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend3[["lower__"]]

female_brm_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (lag 0-2)", y = "Log CPUE") 

print(female_brm_plot)

ggsave("./figs/female_brm_borealization_plot.png", width = 6, height = 4, units = 'in')

## second-best model for females (using trend2_lag1, as for males)-------
female_form <- bf(log_mean ~ log_mean_lag1 + s(trend2_lag1) + s(year, k = 4))

fit_female <- brm(female_form, 
                  data = dat,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(fit_female, file = "output/fit_female_2.rds")

# diagnostics
female_brm <- readRDS("./output/fit_female_2.rds")
check_hmc_diagnostics(female_brm$fit)
neff_lowest(female_brm$fit)
rhat_highest(female_brm$fit)
summary(female_brm)
bayes_R2(female_brm)

plot(conditional_effects(female_brm), ask = FALSE, points = T)

# plot female_brm

## 95% CI
ce1s_1 <- conditional_effects(female_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(female_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(female_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]

female_brm_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Log CPUE") 

print(female_brm_plot)

## combine into multivariate brms model ----

# reload and combine abundance data
abundance_male <- read.csv("./output/male_30-95_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean_male = mean) %>%
  arrange(year) %>%
  mutate(log_mean_male_lag1 = lag(log_mean_male, 1)) %>%
  select(-sd)

abundance_female <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean_female = mean) %>%
  arrange(year) %>%
  mutate(log_mean_female_lag1 = lag(log_mean_female, 1)) %>%
  select(-sd)

# fill in missing 2021 lag abundance
abundance_female$log_mean_female_lag1[abundance_female$year==2021] <- 
  mean(abundance_female$log_mean_female_lag1[abundance_female$year==2020],
       abundance_female$log_mean_female_lag1[abundance_female$year==2022])

abundance_male$log_mean_male_lag1[abundance_male$year==2021] <- 
  mean(abundance_male$log_mean_male_lag1[abundance_male$year==2020],
       abundance_male$log_mean_male_lag1[abundance_male$year==2022])

# combine with trend
dat <- left_join(trend, abundance_female) %>%
  left_join(., abundance_male) %>%
  mutate(trend2 = zoo::rollmean(trend, 2, fill = NA, align = "right"),
         trend2_lag1 = lag(trend2, 1))
  

head(dat)

# examine correlation between male and female abundance
ggplot(dat, aes(log_mean_male, log_mean_female)) +
  geom_point()

ggplot(dat, aes(year, log_mean_female)) +
  geom_point() +
  geom_line()

ggplot(dat, aes(year, log_mean_male)) +
  geom_point() +
  geom_line()

ggplot(dat, aes(log_mean_female)) +
  geom_histogram(bins = 12, fill = "grey", color = "black")

# set up brms model

female_form <- bf(log_mean_female ~ log_mean_female_lag1 + s(trend3) + s(year, k = 4))
male_form <- bf(log_mean_male ~ log_mean_male_lag1 + s(trend2_lag1) + s(year, k = 4))

# scale all time series 
dat$log_mean_female <- as.vector(scale(dat$log_mean_female))
dat$log_mean_female_lag1 <- as.vector(scale(dat$log_mean_female_lag1))

dat$log_mean_male <- as.vector(scale(dat$log_mean_male))
dat$log_mean_male_lag1 <- as.vector(scale(dat$log_mean_male_lag1))

dat$trend2_lag1 <- as.vector(scale(dat$trend2_lag1))



fit_both <- brm(female_form + male_form, 
           data = dat,
           cores = 4, chains = 4, iter = 5000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(fit_both, file = "output/multivariate_brm.rds")


# diagnostics
fit_both <- readRDS("./output/multivariate_brm.rds")
check_hmc_diagnostics(fit_both$fit)
neff_lowest(fit_both$fit)
rhat_highest(fit_both$fit)

bayes_R2(fit_both)

plot(conditional_effects(fit_both), ask = FALSE, points = T)
