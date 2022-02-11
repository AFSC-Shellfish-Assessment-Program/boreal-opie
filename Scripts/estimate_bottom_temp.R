library(tidyverse)
library(lubridate)

dat <- read.csv("./data/haul.csv")

head(dat)

# extract gear temp, date, year, gis_station

dat <- dat %>%
  mutate(julian=yday(parse_date_time(START_TIME, "d-m-y", "US/Alaska"))) %>%
  select(julian, GEAR_TEMPERATURE, SURVEY_YEAR, GIS_STATION) %>%
  rename(bottom.temp = GEAR_TEMPERATURE, year = SURVEY_YEAR, station = GIS_STATION)

# check for  which stations are sampled each year
ff <- function(x) sum(!is.na(x))

check <- tapply(dat$julian, list(dat$year, dat$station), ff)
View(check)

# quite a few are sampled only once!

check <- dat %>%
  group_by(station) %>%
  summarise(count = ff(julian)) %>% 
  arrange(count)

check            

ggplot(check, aes(count)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

# so we can cut off above 20...
# but it also looks like we have some that were sampled > once in a year (> 46 samples in 46-year time series)
# might make the multiple imputation more complicated (i.e., if we want to use a year x station table for imputation!)

# first, cut off all the n < 20 sets

keep <- check$station[check$count > 20]

dat <- dat %>%
  filter(station %in% keep)

# look at these stations sampled more than once a year 

check <- dat %>%
  group_by(year, station) %>%
  summarize(count = n()) %>%
  filter(count > 1)

View(check)

# aha - I think these are the BBRKC retows!

# easiest/best not to use these for temperature
# bad to throw away information, but these are out of the norm re julian:station combinations,
# so might be too complicated to try to analyze (in addition to the difficulties in mice)

# need to go through the situations where 2 tows are made and use the first