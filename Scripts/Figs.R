# all the script required to make Figs 1 & 2 for proposed ms.

library(tidyverse)

theme_set(theme_bw())

# abundance (Fig. 1b)

abundance <- read.csv("./data/Abundance_Snow_table.csv", row.names = 1) %>%
  select(Year, total.abundance) %>%
  filter(Year >=1988)

xtra <- data.frame(Year = 2020,
                   total.abundance = NA)

abundance <- rbind(abundance, xtra)

ggplot(abundance, aes(Year, total.abundance)) +
  geom_line() +
  geom_point() +
  theme(axis.title.x = element_blank()) +
  labs(y = expression(Snow~crab~abundance~(10^6)))


biomass <- read.csv("./data/Biomass_Snow_table.csv", row.names = 1) %>%
  select(Year, total.biomass) %>%
  filter(Year >=1988)

xtra <- data.frame(Year = 2020,
                   total.biomass = NA)

biomass <- rbind(biomass, xtra)

ggplot(biomass, aes(Year, total.biomass)) +
  geom_line() +
  geom_point() 

       