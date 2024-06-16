# plot catch, bycatch, and survey biomass for extended data

library(tidyverse)

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set theme
theme_set(theme_bw())

# load data
biomass <- read.csv("./data/survey_biomass.csv")

# change to kt
biomass_kt <- biomass
biomass_kt[,2:7] <- biomass[,2:7]/1000


catch <- read.csv("./data/catch_bycatch_data.csv")

# vector of years (including 2020, which is missing from survey biomass)
data <- data.frame(year = 1988:2022) %>%
  left_join(., biomass_kt) %>%
  left_join(., catch)


plot <- data %>%
  select(year, male_lt95, immature_female, trawl_bycatch) %>%
  rename(`Male biomass` = male_lt95,
         `Female biomass` = immature_female,
         `Trawl bycatch` = trawl_bycatch) %>%
  pivot_longer(cols = -year)


ggplot(plot, aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cb[c(2,4,6)]) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab(expression(Biomass ~(10^3~t)))

ggsave("./figs/Extended_Data_Fig_6.png", width = 6, height = 4, units = 'in')
