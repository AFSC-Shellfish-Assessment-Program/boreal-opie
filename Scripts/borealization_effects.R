# Effects of Bering Sea borealization on snow crab abundance

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(mgcv)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

# load borealization DFA trend 

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

# clean up and combine

abundance <- abundance %>%
  rename(year = AKFIN_SURVEY_YEAR) %>%
  mutate(log_abundance = log(ABUNDANCE_MIL), .keep = "unused")

dat <- left_join(trend, abundance)

# estimate 2020 abundance as mid value

dat$log_abundance[dat$year == 2020] <-
  mean(dat$log_abundance[dat$year %in% c(2019, 2021)])

plot.dat <- dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)


dat <- dat %>%
  mutate(trend2 = zoo::rollmean(trend, 2, fill = "NA", align = "right"),
         log_abundance_lag1 = lead(log_abundance, n = 1),
         abundance_change = log_abundance_lag1 - log_abundance)


ggplot(dat, aes(abundance_change)) +
  geom_histogram(fill = "grey", color = "black")

# intitial models in mgcv

mod1 <- gam(abundance_change ~ log_abundance + s(trend), data = dat)

summary(mod1)
plot(mod1)


mod2 <- gam(abundance_change ~ log_abundance + s(trend2), data = dat)

summary(mod2)
plot(mod2)


mod3 <- gam(abundance_change ~ s(log_abundance) + s(trend2), data = dat)

summary(mod3)
plot(mod3)
MuMIn::AICc(mod1, mod2, mod3) # mod1 is slightly better

## fit brms models ---------------------

trend_formula <-  bf(abundance_change ~ log_abundance + s(trend))

trend2_formula <-  bf(abundance_change ~ log_abundance + s(trend2))

trend3_formula <-  bf(abundance_change ~ s(log_abundance) + s(trend2))

# drop NAs
dat <- na.omit(dat)

## Show default priors
get_prior(trend_formula, dat)

## fit
trend_brm <- brm(trend_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 2000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(trend_brm, file = "output/trend_brm.rds")

trend_brm  <- add_criterion(trend_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(trend_brm, file = "output/trend_brm.rds")

# diagnostics
trend_brm <- readRDS("./output/trend_brm.rds")
check_hmc_diagnostics(trend_brm$fit)
neff_lowest(trend_brm$fit)
rhat_highest(trend_brm$fit)
summary(trend_brm)
bayes_R2(trend_brm)

plot(trend_brm$criteria$loo, "k")

plot(conditional_effects(trend_brm), ask = FALSE)

y <- dat$abundance_change
yrep_trend_brm  <- fitted(trend_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_trend_brm[sample(nrow(yrep_trend_brm), 25), ]) +
  ggtitle("trend_brm")

trace_plot(trend_brm$fit)

# and the model with two-year trend as covariate

## fit
trend2_brm <- brm(trend2_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 2000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(trend2_brm, file = "output/trend2_brm.rds")

trend2_brm  <- add_criterion(trend2_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(trend2_brm, file = "output/trend2_brm.rds")

# diagnostics
trend2_brm <- readRDS("./output/trend2_brm.rds")
check_hmc_diagnostics(trend2_brm$fit)
neff_lowest(trend2_brm$fit)
rhat_highest(trend2_brm$fit)
summary(trend2_brm)
bayes_R2(trend2_brm)

plot(trend2_brm$criteria$loo, "k")

plot(conditional_effects(trend2_brm), ask = FALSE)

y <- dat$abundance_change
yrep_trend2_brm  <- fitted(trend2_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_trend2_brm[sample(nrow(yrep_trend2_brm), 25), ]) +
  ggtitle("trend2_brm")

trace_plot(trend2_brm$fit)

# model comparison
loo(trend_brm, trend2_brm)
# trend 2 is only slightly better

# plot both to compare

# plot trend_brm

## 95% CI
ce1s_1 <- conditional_effects(trend_brm, effect = "trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(trend_brm, effect = "trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(trend_brm, effect = "trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$trend), amount = 0.01),
                          rep(NA, 100-length(unique(dat$trend))))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization trend", y = "Abundance change") +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g1)

ggsave("./Figs/borealization_abundance_regression.png", width = 6, height = 4, units = 'in')

# plot trend2 for comparison

## 95% CI
ce1s_1 <- conditional_effects(trend2_brm, effect = "trend2", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(trend2_brm, effect = "trend2", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(trend2_brm, effect = "trend2", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$trend2), amount = 0.01),
                          rep(NA, 100-length(unique(dat$trend2))))


g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization trend", y = "Abundance change") +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g2)


# compare with a model that looks at abundance change in the current year, from the previous year

dat <- left_join(trend, abundance)

# estimate 2020 abundance as mid value

dat$log_abundance[dat$year == 2020] <-
  mean(dat$log_abundance[dat$year %in% c(2019, 2021)])


dat <- dat %>%
  mutate(log_abundance_lead1 = lag(log_abundance, n = 1),
         abundance_change2 = log_abundance - log_abundance_lead1)

trend_lead_abundance_formula <-  bf(abundance_change2 ~ log_abundance_lead1 + s(trend))


# drop NAs
dat <- na.omit(dat)

## Show default priors
get_prior(trend_lead_abundance_formula, dat)

## fit
trend_lead_abundance_brm <- brm(trend_lead_abundance_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 2000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(trend_lead_abundance_brm, file = "output/trend_lead_abundance_brm.rds")

trend_lead_abundance_brm  <- add_criterion(trend_lead_abundance_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(trend_lead_abundance_brm, file = "output/trend_lead_abundance_brm.rds")

# diagnostics
trend_lead_abundance_brm <- readRDS("./output/trend_lead_abundance_brm.rds")
check_hmc_diagnostics(trend_lead_abundance_brm$fit)
neff_lowest(trend_lead_abundance_brm$fit)
rhat_highest(trend_lead_abundance_brm$fit)
summary(trend_lead_abundance_brm)
bayes_R2(trend_lead_abundance_brm) # not as good a model as the previous two

plot(trend_lead_abundance_brm$criteria$loo, "k")

plot(conditional_effects(trend_lead_abundance_brm), ask = FALSE) # not as good a model

y <- dat$abundance_change2
yrep_trend_lead_abundance_brm  <- fitted(trend_lead_abundance_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_trend_lead_abundance_brm[sample(nrow(yrep_trend_lead_abundance_brm), 25), ]) +
  ggtitle("trend_lead_abundance_brm")

trace_plot(trend_lead_abundance_brm$fit)

