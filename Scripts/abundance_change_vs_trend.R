## evaluate abundance change approach

library(tidyverse)
library(mice)
library(rstan)
library(brms)
library(bayesplot)

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
  rename(year = YEAR) %>%
  mutate(log_abundance = log(ABUNDANCE), .keep = "unused")

dat <- left_join(trend, abundance)

plot.dat <- dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)


# estimate 2020 with mice ---------------------------------------
library(mice)

# load prediction time series
pred_dat <- read.csv("./Data/imputation_data.csv")

# clean up and remove survey data - can consider these in a lagged application if needed
pred_dat <- pred_dat[1:42,1:8] %>%
  select(-survey_female_mat_biomass, -survey_male_mat_biomass)

# and log transform

pred_dat[,2:6] <- log(pred_dat[,2:6])

dat <- abundance %>%
  rbind(.,
        data.frame(year = 2020,
                   log_abundance = NA)) %>%
  arrange(year) %>%
  left_join(., pred_dat)

# examine correlations!
cor(dat[,-1], use = "p")

imp <- mice(data = dat, method = "norm", m = 100, seed = 957)

# pull out 2020 estimates
estimated_2020 <- NA


for(i in 1:100){
  
  estimated_2020[i] <- complete(imp, i)$log_abundance[41]
  
}


# plot estimate value compared with observed time series
estimated <- data.frame(year = 2020,
                        log_abundance = mean(estimated_2020),
                        LCI = quantile(estimated_2020, 0.025),
                        UCI = quantile(estimated_2020, 0.975))


abundance.plot <- abundance %>%
  mutate(LCI = NA,
         UCI = NA)

abundance.plot <- rbind(abundance.plot, estimated) %>%
  arrange(year)

ggplot(abundance.plot, aes(year, log_abundance)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  ylab("log(immature snow crab abundance)") +
  theme(axis.title.x = element_blank())


abundance.plot2 <- abundance.plot %>%
  mutate(lead_abundance = lead(log_abundance),
         abundance_change = (lead_abundance - log_abundance)) %>%
  left_join(., trend)

ggplot(abundance.plot2, aes(trend, abundance_change)) +
  geom_text(aes(label = year))

# fit exploratory gam
mod <- mgcv::gam(abundance_change ~ log_abundance + s(trend), data = abundance.plot2)
summary(mod)

plot(mod, resid=T, pch = 19, se = F)
