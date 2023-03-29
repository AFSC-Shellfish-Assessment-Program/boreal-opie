# sst-borealization relationship

library(tidyverse)
library(mgcv)

# load dfa trend and immature abundance

trend <- read.csv("./Output/dfa_trend.csv")

trend <- trend %>%
  rename(year = t,
         boreal.trend = estimate) %>%
  select(year, boreal.trend)

sst <- read.csv("./Data/ersst_time_series.csv")

# use annual and winter unsmoothed

sst <- sst %>%
  select(year, annual.unsmoothed, winter.unsmoothed)

dat <- left_join(trend, sst)

plot.dat <- dat %>%
  pivot_longer(cols = c(-year, -boreal.trend))

ggplot(plot.dat, aes(value, boreal.trend)) +
  geom_point() +
  facet_wrap(~name, scales = "free_x")

mod1 <- gamm(boreal.trend ~ s(annual.unsmoothed), correlation = corAR1(), data = dat)

mod2 <- gamm(boreal.trend ~ s(winter.unsmoothed), correlation = corAR1(), data = dat)

MuMIn::AICc(mod1, mod2)

summary(mod1$gam)
