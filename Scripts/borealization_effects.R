## modify model to estimate missing abundance in 2020 ----------------------------
library(tidyverse)
library(mice)
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

trend <- trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

# abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

# abun1 <- read.csv("./output/male3059_drop5_df.csv", row.names = 1) %>%
#   mutate(size = "30-59") %>%
#   dplyr::select(-wtd_log_mean, -unwtd_log_mean, -n_stations)
# 
# abun2 <- read.csv("./output/male6095_drop5_df.csv", row.names = 1) %>%
#   mutate(size = "60-95") %>%
#   dplyr::select(-wtd_log_mean, -unwtd_log_mean, -n_stations)
# 
# # add in NAs
# xtra <- data.frame(year = 2020,
#                    imp_log_mean = NA,
#                    imp_sd = NA,
#                    size = "30-59")
# 
# abun1 <- rbind(abun1, xtra) %>%
#   arrange(year)
# 
# xtra <- xtra %>%
#   mutate(size = "60-95")
# 
# abun2 <- rbind(abun2, xtra) %>%
#   arrange(year)
# 
# abun1 <- read.csv("./output/male3059_drop5_df_simple.csv", row.names = 1) %>%
#   mutate(size = "30-59")
# 
# abun2 <- read.csv("./output/male_60-95_drop5_df_simple.csv", row.names = 1) %>%
#   mutate(size = "60-95")


abundance <- read.csv("./output/male_30-95_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year)

# # clean up and combine
# 
# abundance <- rbind(abun1, abun2) %>%
#   rename(log_mean = mean) %>%
#   dplyr::select(year, log_mean, size) %>%
#   pivot_longer(cols = c(-year, -size)) %>%
#   dplyr::select(-name) %>%
#   pivot_wider(names_from = size, values_from = value)

dat <- left_join(trend, abundance)

# # cross-correlation for age classes: strongest at lag 1 & 2
# ccf(as.vector(dat[dat$year %in% 1975:2019,3]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf
# 
# # autocorrelation for small size class
# acf(as.vector(dat[dat$year %in% 1975:2019,3]))$acf
# 
# # cross-correlation for borealization index v. small size class - nada
# ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,3]))$acf

# cross-correlation for borealization index v. large size class - nada
# ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf

dat$trend2 <- zoo::rollmean(dat$trend, 2, fill = NA, align = "right")
dat$trend3 <- zoo::rollmean(dat$trend, 3, fill = NA, align= "right")

# str(dat)
# 
# # cross-correlation for smoothed borealization index v. small size class - nada
# ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,3]))
# 
# # cross-correlation for smo0thed borealization index v. large size class - nada
# ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,4]))
# 
# set up lag trend for plotting
dat <- dat %>%
  mutate(trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1),
         trend2_lag2 = lag(trend2, 2))
#   
# ggplot(dat, aes(trend, `30-59`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend_lag1, `30-59`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend2, `30-59`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend, `60-95`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend_lag1, `60-95`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend2, `60-95`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend2_lag1, `60-95`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# ggplot(dat, aes(trend3, `60-95`)) +
#   geom_text(aes(label = year)) +
#   geom_smooth(method = "gam", se = F)
# 
# plot.dat <- dat %>%
#   pivot_longer(cols = -year)
# 
# # exploratory GAMs - 60-95 size class as the lag of 30-59 + borealization effect
# dat <- dat %>%
#   rename(small_male = `30-59`,
#          large_male = `60-95`) %>%
#   mutate(small_male_lag1 = lag(small_male, 1),
#          small_male_lag2 = lag(small_male, 2))

dat <- dat %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

ggplot(dat, aes(year, trend)) +
  geom_point() +
  geom_line()

# fill in log_mean_lag1 for 2021
dat$log_mean_lag1[dat$year == 2021] <- mean(dat$log_mean_lag1[dat$year == 2020], dat$log_mean_lag1[dat$year == 2022])

male_mod1 <- gam(log_mean ~ s(log_mean_lag1, k = 4), data = dat)
summary(male_mod1)
plot(male_mod1, se = F, resid = T, pch = 19)

male_mod2 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1), data = dat)
summary(male_mod2)
plot(male_mod2, se = F, resid = T, pch = 19)

male_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend2), data = dat)
summary(male_mod3)

male_mod4 <- gam(log_mean ~ log_mean_lag1 + s(trend3), data = dat)
summary(male_mod4)

male_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2), data = dat)
summary(male_mod5)

MuMIn::AICc(male_mod1, male_mod2, male_mod3, male_mod4, male_mod5)

# add smoothed year term
male_mod2a <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year), data = dat)
summary(male_mod2a)
plot(male_mod2a, se = T, resid = T, pch = 19)


## analyze female abundance -------------------

abundance <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year)

# # clean up and combine
# 
# abundance <- rbind(abun1, abun2) %>%
#   rename(log_mean = mean) %>%
#   dplyr::select(year, log_mean, size) %>%
#   pivot_longer(cols = c(-year, -size)) %>%
#   dplyr::select(-name) %>%
#   pivot_wider(names_from = size, values_from = value)

dat <- left_join(trend, abundance)


dat$trend2 <- zoo::rollmean(dat$trend, 2, fill = NA, align = "right")
dat$trend3 <- zoo::rollmean(dat$trend, 3, fill = NA, align= "right")


dat <- dat %>%
  mutate(trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1),
         trend2_lag2 = lag(trend2, 2))


dat <- dat %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

ggplot(dat, aes(year, trend)) +
  geom_point() +
  geom_line()

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

female_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2), data = dat)
summary(female_mod5)

MuMIn::AICc(female_mod1, female_mod2, female_mod3, female_mod4, female_mod5) 

## combine into multivariate brms model ----

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

ggplot(dat, aes(log_mean_female)) +
  geom_histogram(bins = 12, fill = "grey", color = "black")

# set up brms model

female_form <- bf(log_mean_female ~ log_mean_female_lag1 + s(trend2_lag1) + s(year, k = 4))
male_form <- bf(log_mean_male ~ log_mean_male_lag1 + s(trend2_lag1) + s(year, k = 4))

# scale all time series 
dat$log_mean_female <- as.vector(scale(dat$log_mean_female))
dat$log_mean_female_lag1 <- as.vector(scale(dat$log_mean_female_lag1))

dat$log_mean_male <- as.vector(scale(dat$log_mean_male))
dat$log_mean_male_lag1 <- as.vector(scale(dat$log_mean_male_lag1))

dat$trend2_lag1 <- as.vector(scale(dat$trend2_lag1))


# begin with sex-specific models
fit_female <- brm(female_form, 
           data = dat,
           cores = 4, chains = 4, iter = 2000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99, max_treedepth = 14))

saveRDS(fit_female, file = "output/fit_female.rds")


fit_both <- brm(female_form + male_form, 
           data = dat,
           cores = 4, chains = 4, iter = 5000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(fit_both, file = "output/multivariate_brm.rds")

# model small males 
sm_mod1 <- gam(small_male ~ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod1)
plot(sm_mod1, resid = T, se = F, pch = 19)

sm_mod2 <- gam(small_male ~ s(trend_lag1) + s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod2)
plot(sm_mod2, resid = T, se = F, pch = 19)

sm_mod3 <- gam(small_male ~ s(trend2) + s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod3)
plot(sm_mod3, resid = T, se = F, pch = 19)

sm_mod4 <- gam(small_male ~ s(trend)+ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod4)
plot(sm_mod4, resid = T, se = F, pch = 19)

sm_mod5 <- gam(small_male ~ s(trend2_lag1)+ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod5)
plot(sm_mod5, resid = T, se = F, pch = 19)

sm_mod6 <- gam(small_male ~ s(trend3)+ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod6)
plot(sm_mod6, resid = T, se = F, pch = 19)


ggplot(dat, aes(trend2_lag1, small_male)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)


dat <- dat %>%
  mutate(trend_lag2 = lag(trend, 2))

sm_mod7 <- gam(small_male ~ s(trend_lag2, k = 4)+ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod7)
plot(sm_mod7, resid = T, se = F, pch = 19)

sm_mod8 <- gam(small_male ~ s(trend_lag1, k = 4) + s(trend_lag2, k = 4)+ s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(sm_mod8)
plot(sm_mod8, resid = T, se = F, pch = 19)


MuMIn::AICc(sm_mod1, sm_mod2, sm_mod3, sm_mod4, sm_mod5, sm_mod6, sm_mod7, sm_mod8) # best model is sm_mod3

sm_mod9 <- gam(small_male ~ s(trend), dat = dat, na.action = "na.omit")

sm_mod10 <- gam(small_male ~ s(trend2), dat = dat, na.action = "na.omit")

sm_mod11 <- gam(small_male ~ s(trend3), dat = dat, na.action = "na.omit")

sm_mod12 <- gam(small_male ~ s(trend_lag1), dat = dat, na.action = "na.omit")

sm_mod13 <- gam(small_male ~ s(trend2_lag1), dat = dat, na.action = "na.omit")

sm_mod14 <- gam(small_male ~ s(trend_lag2), dat = dat, na.action = "na.omit")

MuMIn::AICc(sm_mod9, sm_mod10, sm_mod11, sm_mod12, sm_mod13, sm_mod14) # sm_mod13 (trend2_lag1) is best

summary(sm_mod13)

plot(sm_mod13, pch = 19, se = F, resid = T)


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


sm_brm_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (mean lag 1-2)", y = "Mean log CPUE") +
  geom_text(data = filter(dat, year != 2020), aes(x = trend2_lag1, y = small_male, label = year), size = 3)

print(sm_brm_plot)

ggsave("./Figs/borealization_abundance_regression_small_male_no_imputation__no_ar.png", width = 6, height = 4, units = 'in')


## large male models in gam -----------------------------------------------
dat <- dat %>%
  mutate(large_male_lag1 = lag(large_male))

large_mod1 <- gam(large_male ~ s(trend), dat = dat, na.action = "na.omit")
summary(large_mod1)

large_mod2 <- gam(large_male ~ s(trend2), dat = dat, na.action = "na.omit")
summary(large_mod2)

large_mod3 <- gam(large_male ~ s(trend3), dat = dat, na.action = "na.omit")
summary(large_mod3)

large_mod4 <- gam(large_male ~ s(trend_lag1), dat = dat, na.action = "na.omit")

large_mod5 <- gam(large_male ~ s(trend2_lag1), dat = dat, na.action = "na.omit")

large_mod6 <- gam(large_male ~ s(trend_lag2), dat = dat, na.action = "na.omit")

MuMIn::AICc(large_mod1, large_mod2, large_mod3, large_mod4, large_mod5, large_mod6)

summary(large_mod6)
summary(large_mod5)


summary(large_mod1)
plot(large_mod1, resid = T, se = F, pch = 19)

large_mod2 <- gam(large_male ~ s(large_male_lag1) + s(small_male_lag1), dat = dat, na.action = "na.omit")
summary(large_mod2)
plot(large_mod2, resid = T, se = F, pch = 19)

large_mod3 <- gam(large_male ~ s(large_male_lag1) + s(small_male_lag1) + s(small_male_lag2), dat = dat, na.action = "na.omit")
summary(large_mod3)
plot(large_mod3, resid = T, se = F, pch = 19)

large_mod4 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2), dat = dat, na.action = "na.omit")
summary(large_mod4)
plot(large_mod4, resid = T, se = F, pch = 19)

# find best way to model autocorrelation
MuMIn::AICc(large_mod1, large_mod2, large_mod3, large_mod4) # mod4 is best


large_mod5 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend), dat = dat, na.action = "na.omit")
summary(large_mod5)
plot(large_mod5, resid = T, se = F, pch = 19)

large_mod6 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend2), dat = dat, na.action = "na.omit")
summary(large_mod6)
plot(large_mod5, resid = T, se = F, pch = 19)

large_mod7 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend3), dat = dat, na.action = "na.omit")
summary(large_mod7)
plot(large_mod7, resid = T, se = F, pch = 19)

large_mod8 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend_lag1), dat = dat, na.action = "na.omit")
summary(large_mod8)
plot(large_mod8, resid = T, se = F, pch = 19)

large_mod9 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend2_lag1), dat = dat, na.action = "na.omit")
summary(large_mod9)
plot(large_mod9, resid = T, se = F, pch = 19)

large_mod10 <- gam(large_male ~ s(small_male_lag1) + s(small_male_lag2) +
                    s(trend_lag2), dat = dat, na.action = "na.omit")
summary(large_mod10)
plot(large_mod10, resid = T, se = F, pch = 19)

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

lrg_brm_plot2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Mean log CPUE") 

print(lrg_brm_plot2)


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


lrg_brm_plot1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Small male log CPUE (previous year)", y = "Mean log CPUE") 

print(lrg_brm_plot1)

## combine and plot -----------------------------------

# make abundance plots for each time series
sm_abun_plot <- ggplot(abun1, aes(year, imp_log_mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = imp_log_mean - 2*imp_sd,
                    ymax = imp_log_mean + 2*imp_sd)) +
  labs(y = "Log CPUE") +
  theme(axis.title.x = element_blank())

sm_abun_plot

lrg_abun_plot <- ggplot(abun2, aes(year, imp_log_mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = imp_log_mean - 2*imp_sd,
                    ymax = imp_log_mean + 2*imp_sd)) +
  labs(y = "Log CPUE") +
  theme(axis.title.x = element_blank())

lrg_abun_plot

empty_plot <- ggplot + theme_void()

png("./Figs/combined_borealization_effects_plot.png", width = 9, height = 5, units = 'in', res = 300)

ggarrange(ggarrange(sm_abun_plot, sm_brm_plot, empty_plot, ncol = 3, labels = c("a", "b", ""), widths = c(3,3,2)),
                  ggarrange(lrg_abun_plot, lrg_brm_plot1, lrg_brm_plot2, ncol = 3, labels = c("c", "d", "e")),
          nrow = 2)

dev.off()


##

dat <- dat %>%
  mutate(trend2_lag2 = lag(trend2_lag1))


combined_mod1 <- gam(large_male ~ small_male_lag1 + s(trend2_lag2) + s(trend), dat = dat, na.action = "na.omit")
summary(combined_mod1)
plot(combined_mod1, resid = T, se = F, pch = 19, all.terms = T)

## brms version

form <- bf(large_male ~ small_male_lag1 + trend2_lag2 + trend + ar(time = year, p = 1, cov = TRUE))

## fit
brm <- brm(form,
           data = dat,
           cores = 4, chains = 4, iter = 3000,
           save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(brm, file = "output/combined_male_brm1.rds")

plot(conditional_effects(brm), ask = F)
