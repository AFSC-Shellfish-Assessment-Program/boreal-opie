# plot sst and crab abundance

library(tidyverse)

theme_set(theme_bw())

# load sst

sst <- read.csv("./Data/ersst_time_series.csv")

# clean up
sst <- sst %>%
  filter(year >= 1972) %>%
  select(year, annual.unsmoothed) %>%
  rename(annual_sst = annual.unsmoothed)

# load snow crab abundance

abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

# clean up and combine

abundance <- abundance %>%
  rename(year = AKFIN_SURVEY_YEAR) %>%
  mutate(log_abundance = log(ABUNDANCE_MIL), .keep = "unused")

xtra <- data.frame(year = 2020, 
                   log_abundance = mean(abundance$log_abundance[abundance$year %in% c(2019, 2021)]))

abundance <- rbind(abundance, xtra)

# combine

data <- left_join(sst, abundance) %>%
  pivot_longer(cols = -year) %>%
  filter(year < 2022)


ggplot(data, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank())

ggsave("./Figs/sst_abundance.png", width = 4, height = 5, units = 'in')
