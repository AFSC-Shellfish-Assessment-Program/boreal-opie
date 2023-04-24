# plot harvest time series for talk / paper

library(tidyverse)

dat <- read.csv("./data/harvest.csv")

# convert to mt
dat <- dat %>%
  mutate(tons = Harvest / 2205)

ggplot(dat, aes(Season, tons)) +
  geom_point() +
  geom_line()

dat <- read.csv("./data/value.csv") %>%
  select(Season, Ex.vessel)

ggplot(dat, aes(Season, Ex.vessel)) +
  geom_point() +
  geom_line()

# control for inflation
infl <- read.csv("./data/inflation_data.csv") 

# get adjuster for 2022 value

infl <- infl %>%
  mutate(value_2022 = infl$cpi[infl$year == 2022]/ cpi)  %>%
  select(-cpi) %>%
  rename(Season = year)

dat <- left_join(dat, infl) %>%
  mutate(adj_value = Ex.vessel * value_2022)

dat$adj_value[dat$Season == 2023] <- 0
  
ggplot(dat, aes(Season, adj_value/1e6)) +
  geom_point() +
  geom_line()

# add number of vessels each year
vessels <- read.csv("./data/vessel_count.csv")

dat <- left_join(dat, vessels) %>%
  mutate(value_per_boat = adj_value / Vessels)

dat$value_per_boat[dat$Season == 2023] <- 0

ggplot(dat, aes(Season, value_per_boat/1e6)) +
  geom_point() +
  geom_line() +
  labs(y = "Million $USD per boat") +
  theme(axis.title.x = element_blank())

ggsave("./figs/value_per_vessel.png", width = 3.5, height = 2.5, units = 'in')
