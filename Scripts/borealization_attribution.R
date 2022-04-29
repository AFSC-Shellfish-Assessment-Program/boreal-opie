# Attribution of Bering Sea borealization

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

sst <- read.csv("./Data/ersst_anomaly_time_series.csv")

sst <- sst %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(sst.anomaly = annual.anomaly.unsmoothed)

dat <- left_join(trend, sst)

ggplot(dat, aes(sst.anomaly, trend)) +
  geom_point()

# fit brms model ------------------------------

sst_boreal_formula <-  bf(trend ~ s(sst.anomaly))

## Show default priors
get_prior(sst_boreal_formula, dat)

## fit 
sst_boreal_brm <- brm(sst_boreal_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(sst_boreal_brm, file = "output/sst_boreal_brm.rds")

# run diagnostics
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")
check_hmc_diagnostics(sst_boreal_brm$fit)
neff_lowest(sst_boreal_brm$fit)
rhat_highest(sst_boreal_brm$fit)
summary(sst_boreal_brm)
bayes_R2(sst_boreal_brm)

# save high-resolution predicted values for expected return time
ce1s_1 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.025, 0.975), resolution = 5000)

write.csv(ce1s_1$sst.anomaly, "./output/sst_borealization_predicted_relationship.csv", row.names = F)

# plot
## 95% CI
ce1s_1 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$sst.anomaly
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$sst.anomaly[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$sst.anomaly[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$sst.anomaly[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$sst.anomaly[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$sst.anomaly), amount = 0.01),
                          rep(NA, 100-length(unique(dat$sst.anomaly))))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  labs(x = "SST anomaly wrt 1950-1999", y = "Borealization trend") +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g1)

## estimate FAR -------------------------

# load model object for predicting preindustrial and historical probabilities
mod <- readRDS("./output/Eastern_Bering_sea_rolling_window_binomial2.rds")

far_pred <- data.frame()

  
# set up sst anomalies
ersst.temp <- sst %>%
  rename(annual.anomaly.1yr = sst.anomaly)
  
## setup new data
nd <- data.frame(period = c("historical", "preindustrial"),
                   ersst.year = rep(ersst.temp$year, each = 2),
                   annual.anomaly.1yr = rep(ersst.temp$annual.anomaly.1yr, each = 2),
                   N = 1000,
                   model_fac = NA)
  
nd_pre <- nd[nd$period == "preindustrial", ]
nd_his <- nd[nd$period == "historical", ]
  
## make predictions
## exclude random effects for model_fac
pre_pp <- posterior_epred(mod, newdata = nd_pre, re_formula = NA)
his_pp <- posterior_epred(mod, newdata = nd_his, re_formula = NA)
  
## Calc probabilities
## These are our posterior probabilities to use for FAR calculation
pre_prob <- pre_pp / unique(nd$N)
his_prob <- his_pp / unique(nd$N)
  
  
## Calc FAR
far <- 1 - (pre_prob / his_prob)
range(far, na.rm = TRUE)
  
far_pred <- rbind(far_pred,
                         data.frame(annual.anomaly.1yr = nd_pre$annual.anomaly.1yr,
                                    year = nd_pre$ersst.year,
                                    far = apply(far, 2, mean),
                                    lower = apply(far, 2, quantile, probs = 0.025),
                                    upper = apply(far, 2, quantile, probs = 0.975)))
  

g2 <- ggplot(far_pred) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_line(aes(x = year, y = far), size = 0.2, color = "red3") +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
  ylab("Fraction of attributable risk") +
  xlab("Year")

print(g2)


# aside - calculate risk ratio

RR <- far_pred %>%
  filter(year >= 2014) %>%
  mutate(RR = 1/(1-far))



## fit attribution model in brms----------------------

# NB - could add uncertainty in FAR estimates into this model

far_add <- far_pred %>%
  select(year, far)

dat <- left_join(dat, far_add)

attribution_boreal_formula <-  bf(trend ~ s(far))

## Show default priors
get_prior(attribution_boreal_formula, dat)

## fit 
attribution_boreal_brm <- brm(attribution_boreal_formula,
                      data = filter(dat, year >= 1994), # limit to positive far
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(attribution_boreal_brm, file = "output/attribution_boreal_brm.rds")

# run diagnostics
attribution_boreal_brm <- readRDS("./output/attribution_boreal_brm.rds")
check_hmc_diagnostics(attribution_boreal_brm$fit)
neff_lowest(attribution_boreal_brm$fit)
rhat_highest(attribution_boreal_brm$fit)
summary(attribution_boreal_brm)
bayes_R2(attribution_boreal_brm)

# plot
## 95% CI
ce1s_1 <- conditional_effects(attribution_boreal_brm, effect = "far", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(attribution_boreal_brm, effect = "far", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(attribution_boreal_brm, effect = "far", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$far
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$far[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$far[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$far[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$far[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$far[dat$year >= 1994]), amount = 0.01),
                          rep(NA, 100-length(unique(dat$far[dat$year >= 1994]))))


g3 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  labs(x = "Fraction of attributable risk", y = "Borealization trend") +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g3)


# save
png("./Figs/borealization_attribution.png", width = 4, height = 8, units = 'in', res = 300)

ggpubr::ggarrange(g1,
                  g2,
                  g3,
                  ncol = 1,
                  labels = "auto")

dev.off()

