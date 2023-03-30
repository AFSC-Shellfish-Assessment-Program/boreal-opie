## modify model to estimate missing abundance in 2020 ----------------------------
library(tidyverse)
library(mice)
library(rstan)
library(brms)
library(bayesplot)
library(mgcv)

source("./scripts/stan_utils.R")

theme_set(theme_bw())


# reload data

# load borealization DFA trend

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

# abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

abun1 <- read.csv("./output/male3059_drop5_df.csv", row.names = 1) %>%
  mutate(size = "30-59") %>%
  dplyr::select(-wtd_log_mean, -unwtd_log_mean, -n_stations)

abun2 <- read.csv("./output/male6095_drop5_df.csv", row.names = 1) %>%
  mutate(size = "60-95") %>%
  dplyr::select(-wtd_log_mean, -unwtd_log_mean, -n_stations)

# add in NAs
xtra <- data.frame(year = 2020,
                   imp_log_mean = NA,
                   imp_sd = NA,
                   size = "30-59")

abun1 <- rbind(abun1, xtra) %>%
  arrange(year)

xtra <- xtra %>%
  mutate(size = "60-95")

abun2 <- rbind(abun2, xtra) %>%
  arrange(year)

# clean up and combine

abundance <- rbind(abun1, abun2) %>%
  rename(log_mean = imp_log_mean) %>%
  dplyr::select(year, log_mean, size) %>%
  pivot_longer(cols = c(-year, -size)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = size, values_from = value)

dat <- left_join(trend, abundance)

# cross-correlation for age classes: strongest at lag 1
ccf(as.vector(dat[dat$year %in% 1975:2019,3]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf

# cross-correlation for borealization index v. small size class - strongest at lag 1!
ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,3]))$acf

# cross-correlation for borealization index v. large size class - strongest at lag 0!
ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf

dat$trend2 <- zoo::rollmean(dat$trend, 2, fill = NA, align = "right")
dat$trend3 <- zoo::rollmean(dat$trend, 3, fill = NA, align= "right")

str(dat)

# cross-correlation for smoothed borealization index v. small size class - strongest at lag 0
ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,3]))

# cross-correlation for smo-thed borealization index v. small size class - strongest at lag 0
ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,4]))

# set up lag trend for plotting
dat <- dat %>%
  mutate(trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1))
  
ggplot(dat, aes(trend, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend_lag1, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend_lag1, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2_lag1, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend3, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

plot.dat <- dat %>%
  pivot_longer(cols = -year)

# exploratory GAMs - 60-95 size class as the lag of 30-59 + borealization effect
dat <- dat %>%
  rename(small_male = `30-59`,
         large_male = `60-95`) %>%
  mutate(small_male_lag1 = lag(small_male, 1))
 
# model small males 
sm_mod1 <- gam(small_male ~ s(trend_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod1)
plot(sm_mod1, resid = T, se = F, pch = 19)

sm_mod2 <- gam(small_male ~ s(trend2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod2)
plot(sm_mod2, resid = T, se = F, pch = 19)

sm_mod3 <- gam(small_male ~ s(trend, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod3)
plot(sm_mod3, resid = T, se = F, pch = 19)

sm_mod4 <- gam(small_male ~ s(trend2_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod4)
plot(sm_mod4, resid = T, se = F, pch = 19)

ggplot(dat, aes(trend2_lag1, small_male)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

sm_mod5 <- gam(small_male ~ s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod5)
plot(sm_mod5, resid = T, se = F, pch = 19)

dat <- dat %>%
  mutate(trend_lag2 = lag(trend, 2))

sm_mod6 <- gam(small_male ~ s(trend_lag2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod6)
plot(sm_mod6, resid = T, se = F, pch = 19)

sm_mod7 <- gam(small_male ~ s(trend_lag1, k = 4) + s(trend_lag2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod7)
plot(sm_mod7, resid = T, se = F, pch = 19)


MuMIn::AICc(sm_mod1, sm_mod2, sm_mod3, sm_mod4, sm_mod5, sm_mod6, sm_mod7) # best model is sm_mod4

## brms versions of small male models --------------------

# model 1
form <- bf(small_male ~ s(trend_lag1, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm1.rds")

# model 2
form <- bf(small_male ~ s(trend2, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm2.rds")

# model 3
form <- bf(small_male ~ s(trend, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm3.rds")

# model 4
form <- bf(small_male ~ s(trend2_lag1, k = 4))

## fit
brm <- brm(form,
                         data = dat,
                         cores = 4, chains = 4, iter = 3000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm4.rds")

# model 5
form <- bf(small_male ~ s(trend3, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm5.rds")

# model 6
form <- bf(small_male ~ s(trend_lag2, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm6.rds")

# model 7
form <- bf(small_male ~ s(trend_lag1, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(brm, file = "output/sm_male_brm7.rds")

## brms small male model selection ---------------------------

# load model objects
sm_brm1 <- readRDS("./output/sm_male_brm1.rds")
sm_brm2 <- readRDS("./output/sm_male_brm2.rds")
sm_brm3 <- readRDS("./output/sm_male_brm3.rds")
sm_brm4 <- readRDS("./output/sm_male_brm4.rds")
sm_brm5 <- readRDS("./output/sm_male_brm5.rds")
sm_brm6 <- readRDS("./output/sm_male_brm6.rds")
sm_brm7 <- readRDS("./output/sm_male_brm7.rds")

# calculate loo values for each
sm_brm1 <- add_criterion(sm_brm1, "loo")
sm_brm2 <- add_criterion(sm_brm2, "loo")
sm_brm3 <- add_criterion(sm_brm3, "loo")
sm_brm4 <- add_criterion(sm_brm4, "loo")
sm_brm5 <- add_criterion(sm_brm5, "loo")
sm_brm6 <- add_criterion(sm_brm6, "loo")
sm_brm7 <- add_criterion(sm_brm7, "loo")

# and compare the models
compare <- loo_compare(sm_brm1, sm_brm2, sm_brm3, sm_brm4, sm_brm5, sm_brm6, sm_brm7, criterion = "loo")

compare # same result as for mgcv - model 4 is best

# diagnostics for best model
brm <- readRDS("./output/sm_male_brm4.rds")
check_hmc_diagnostics(brm$fit)
neff_lowest(brm$fit)
rhat_highest(brm$fit)

bayes_R2(brm)

plot(conditional_effects(brm), ask = FALSE)


# plot

## 95% CI
ce1s_1 <- conditional_effects(brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]


sm_brm <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (mean lag 1-2)", y = "Mean log CPUE") +
  geom_text(data = filter(dat, year != 2020), aes(x = trend2_lag1, y = small_male, label = year), size = 3)

print(sm_brm)

ggsave("./Figs/borealization_abundance_regression_small_male_no_imputation__no_ar.png", width = 6, height = 4, units = 'in')


## large male models in gam -----------------------------------------------

large_mod1 <- gam(large_male ~ s(small_male_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod1)
plot(large_mod1, resid = T, se = F, pch = 19)

large_mod2 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod2)
plot(large_mod2, resid = T, se = F, pch = 19)

large_mod3 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend2, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod3)
plot(large_mod3, resid = T, se = F, pch = 19)

large_mod4 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod4)
plot(large_mod4, resid = T, se = F, pch = 19)

large_mod5 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod5)
plot(large_mod5, resid = T, se = F, pch = 19)

large_mod6 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend2_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod6)
plot(large_mod6, resid = T, se = F, pch = 19)

MuMIn::AICc(large_mod1, large_mod2, large_mod3, large_mod4, large_mod5, large_mod6) # mod 2 is best

## brms versions of large male models --------------------

# model 1
form <- bf(large_male ~ s(small_male_lag1, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm1.rds")

# model 2
form <- bf(large_male ~ s(small_male_lag1, k = 4) + s(trend, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm2.rds")

# model 3
form <- bf(large_male ~ s(small_male_lag1, k = 4) + s(trend2, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm3.rds")

# model 4
form <- bf(large_male ~ s(small_male_lag1, k = 4) + s(trend3, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm4.rds")

# model 5
form <- bf(large_male ~ s(small_male_lag1, k = 4) + s(trend_lag1, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm5.rds")

# model 6
form <- bf(large_male ~ s(small_male_lag1, k = 4) + s(trend2_lag1, k = 4))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/lrg_male_brm6.rds")

## brms large male model selection ---------------------------

# load model objects
lrg_brm1 <- readRDS("./output/lrg_male_brm1.rds")
lrg_brm2 <- readRDS("./output/lrg_male_brm2.rds")
lrg_brm3 <- readRDS("./output/lrg_male_brm3.rds")
lrg_brm4 <- readRDS("./output/lrg_male_brm4.rds")
lrg_brm5 <- readRDS("./output/lrg_male_brm5.rds")
lrg_brm6 <- readRDS("./output/lrg_male_brm6.rds")

# calculate loo values for each
lrg_brm1 <- add_criterion(lrg_brm1, "loo")
lrg_brm2 <- add_criterion(lrg_brm2, "loo", moment_match = T)
lrg_brm3 <- add_criterion(lrg_brm3, "loo", moment_match = T)
lrg_brm4 <- add_criterion(lrg_brm4, "loo")
lrg_brm5 <- add_criterion(lrg_brm5, "loo")
lrg_brm6 <- add_criterion(lrg_brm6, "loo")

# and compare the models
compare <- loo_compare(lrg_brm1, lrg_brm2, lrg_brm3, lrg_brm4, lrg_brm5, lrg_brm6, criterion = "loo")

compare # same result as for mgcv - model 2 is best

# diagnostics
lrg_brm <- readRDS("./output/lrg_male_brm2.rds")
check_hmc_diagnostics(lrg_brm$fit)
neff_lowest(lrg_brm$fit)
rhat_highest(lrg_brm$fit)

bayes_R2(lrg_brm)

plot(conditional_effects(lrg_brm), ask = FALSE)

# plot lrg_brm

## 95% CI
ce1s_1 <- conditional_effects(lrg_brm, effect = "trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(lrg_brm, effect = "trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(lrg_brm, effect = "trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend[["lower__"]]


lrg_brm2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Mean log CPUE") 

print(lrg_brm2)


ce1s_1 <- conditional_effects(lrg_brm, effect = "small_male_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(lrg_brm, effect = "small_male_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(lrg_brm, effect = "small_male_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))

dat_ce <- ce1s_1$small_male_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$small_male_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$small_male_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$small_male_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$small_male_lag1[["lower__"]]


lrg_brm1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Small male log CPUE (previous year)", y = "Mean log CPUE") 

print(lrg_brm1)

## combine and plot -----------------------------------

sm_abun <- ggplot(abun1, aes(year, imp_log_mean)) +
  geom_l