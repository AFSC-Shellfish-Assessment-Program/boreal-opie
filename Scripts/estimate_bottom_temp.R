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

# make sure there are a max of 2 tows per station/year combo
max(check$count) # yes

# aha - I think these are the BBRKC retows!

# easiest/best not to use these for temperature
# bad to throw away information, but these are out of the norm re julian:station combinations,
# so might be too complicated to try to analyze (in addition to the difficulties in mice)

# need to go through the situations where 2 tows are made and use the first

# set up an index for the ones we want to drop
dat$station_year = paste(dat$station, dat$year, sep = "_")

check$station_year = paste(check$station, check$year, sep = "_")

# now pick these out of dat

doubles <- dat %>%
  filter(station_year %in% check$station_year)

# confirm this has the double sample events

check.double <- doubles %>%
  group_by(station_year) %>%
  summarise(count = n())

unique(check.double$count) # good...

# now select the maximum date for each station-year
# these are the ones to drop from dat

drop.check <- doubles %>%
  group_by(station_year) %>%
  summarise(julian = max(julian)) %>% # picks the second day in each station-year combo
  mutate(station_year_julian = paste(station_year, julian, sep = "_")) 

# and drop these from dat
dat$station_year_julian <- paste(dat$station_year, dat$julian, sep = "_")
  
dat <- dat %>%
  filter(!station_year_julian %in% drop.check$station_year_julian)

# check
c.dat <- dat %>%
  group_by(year, station) %>%
  summarise(count = n())

unique(c.dat$count) # looks good

# now set up for multiple imputation