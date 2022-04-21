# Effects of Bering Sea borealization on snow crab abundance

library(tidyverse)
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

