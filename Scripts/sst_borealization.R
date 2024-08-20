# create operating sst-borealization model
# for attribution of recent borealization events

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

# load borealization DFA trend and sst anomaly time series

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate)

sst <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv")

sst <- sst %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(sst.anomaly = annual.anomaly.unsmoothed)

dat <- left_join(trend, sst)

ggplot(dat, aes(sst.anomaly, trend)) +
  geom_point() +
  geom_smooth(method = "gam", color = "blue", se = F) +
  geom_smooth(method = "lm", color = "red", se = F)

# fit brms model ------------------------------

sst_boreal_formula <-  bf(trend ~ s(sst.anomaly))

## Show default priors
get_prior(sst_boreal_formula, dat)

## fit 
sst_boreal_brm <- brm(sst_boreal_formula,
                      data = dat,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.9999, max_treedepth = 10))

saveRDS(sst_boreal_brm, file = "output/sst_boreal_brm.rds")

# run diagnostics
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")
check_hmc_diagnostics(sst_boreal_brm$fit)
neff_lowest(sst_boreal_brm$fit)
rhat_highest(sst_boreal_brm$fit)
summary(sst_boreal_brm)
bayes_R2(sst_boreal_brm)