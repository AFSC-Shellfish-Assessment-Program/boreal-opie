# Attribution of Bering Sea borealization - two-year running mean SST and borealization index

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

mc.cores = parallel::detectCores()

# load borealization DFA trend and sst anomaly time series

trend <- read.csv("./output/dfa_trend.csv") %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate) %>%
  mutate(two.yr.index = zoo::rollmean(trend, 2, align = "right", fill = NA))

sst <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(year, annual.anomaly.unsmoothed, annual.anomaly.two.yr.running.mean) %>%
  rename(annual.sst = annual.anomaly.unsmoothed,
         two.yr.sst = annual.anomaly.two.yr.running.mean) %>%
  mutate(test.2yr = zoo::rollmean(annual.sst, 2, align = "right", fill = NA)) 
# so we want right-aligned running means, attribution study uses left-aligned!

dat <- left_join(trend, sst) %>%
  select(year, test.2yr, two.yr.index) %>% # using right-alinged sst!
  rename(two.yr.sst = test.2yr)

ggplot(dat, aes(two.yr.sst, two.yr.index)) +
  geom_point() +
  geom_smooth(method = "gam", color = "blue", se = F) +
  geom_smooth(method = "lm", color = "red", se = F)

# fit brms model ------------------------------

sst_boreal_2yr_formula <-  bf(two.yr.index ~ s(two.yr.sst))

## Show default priors
get_prior(sst_boreal_2yr_formula, dat)

## fit 
sst_boreal_2yr_brm <- brm(sst_boreal_2yr_formula,
                      data = dat,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.9999, max_treedepth = 10))

saveRDS(sst_boreal_2yr_brm, file = "output/sst_boreal_2yr_running_mean_brm.rds")

# run diagnostics
sst_boreal_2yr_brm <- readRDS("./output/sst_boreal_2yr_running_mean_brm.rds")
check_hmc_diagnostics(sst_boreal_2yr_brm$fit)
neff_lowest(sst_boreal_2yr_brm$fit)
rhat_highest(sst_boreal_2yr_brm$fit)
summary(sst_boreal_2yr_brm)
bayes_R2(sst_boreal_2yr_brm)

# plot
## 95% CI
ce1s_1 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$two.yr.sst
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$two.yr.sst[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$two.yr.sst[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$two.yr.sst[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$two.yr.sst[["lower__"]]
# dat_ce[["rug.anom"]] <- c(jitter(unique(dat$two.yr.sst), amount = 0.01),
# rep(NA, 100-length(unique(dat$two.yr.sst))))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  labs(x = "Two-year running mean SST anomaly", y = "Two-year running mean borealization index") +
  geom_text(data = dat, aes(two.yr.sst, two.yr.index, label = year), size = 3)
# geom_rug](aes(x=rug.anom, y=NULL)) 

print(g1)
