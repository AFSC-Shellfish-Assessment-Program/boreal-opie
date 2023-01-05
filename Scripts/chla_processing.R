# Goal: create an index of annual chla size fraction,
# controlling for spatial-temporal differences in sampling among years

#Data source includes summer/fall depth integrated (top 50m) chla biomass  
  #from the Bering Sea fall BASIS survey 2003-2018
#NOTES: chl-a size fraction ratio is avg of >10um/avg total chla

library(tidyverse)
library(mgcv)
library(lubridate)

#Author: Erin Fedewa
#Last update: 1/5/2023

# load data
dat <- read.csv("./Data/chla_raw.csv")

###############################################
#data wrangling
dat %>%
  separate(gear_time, c("date", "time"), sep= "\\s", extra = "drop", fill = "right") %>%
  mutate(date = mdy(date),
        julian = yday(date),
        year = as.factor(year)) -> dat.new
           
#Spatiotemporal data exploration:

#Stations sampled each yr (EBS & NBS)
dat.new %>%
  group_by(year) %>%
  summarise(num_stations = length(unique(station_ID)))

#Stations sampled each yr (EBS ONLY)
dat.new %>%
  filter(gear_lat < 60) %>%
  group_by(year) %>%
  summarise(num_stations = length(unique(station_ID)))
#No data for 2013/2017, need to throw out 2019 

#Plot stations w/ EBS/NBS border in red 
dat.new %>%
  group_by(year, gear_lat, gear_long) %>%
  distinct(station_ID) %>%
  ggplot() +
  geom_point(aes(x = gear_long, y = gear_lat), size=.5) +
  geom_hline(yintercept=60, color = "red") +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~year)

#EBS stations sampled 50% of the years
dat.new %>%
  filter(gear_lat < 60,
         year != 2019) %>%
  group_by(gear_lat, gear_long) %>%
#Need to come back to this......


#EBS Mean date sampled 
dat.new %>%
  filter(gear_lat < 60,
         year != 2019) %>%
  group_by(year) %>%
  summarise(mean_date = mean(julian, na.rm=T),
            min_date = min(julian, na.rm=T),
            max_date = max(julian, na.rm=T)) 
#plot
dat.new %>%
  filter(gear_lat < 60,
         year != 2019) %>%
  ggplot(aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year) # big differences 

#To do: spatially subset by bsierp regions? How to deliniate EBS vrs NBS
#Need to filter for EBS, pull 2019, filter for stations sampled in >50% of years?

###############################################
# fit model to chl-a ratio to control for lat,long and julian day and 
# estimate annual chl-a

#ALL data- not spatially subset yet 
chla.mod <- gam(avg_ratio ~ te(gear_long, gear_lat) + s(julian, k = 5) + year, 
                  data = dat.new)

summary(chla.mod)
plot(chla.mod)

#Create new data frame with average lat, long, and day and predict ratio for each year
new.chla.dat <- data.frame(julian = mean(dat.new$julian, na.rm=T),
                             lat = mean(dat.new$gear_lat),
                             long = mean(dat.new$gear_long),
                             year = unique(dat.new$year))

chla.pred <- predict(chla.mod, newdata = new.chla.dat, se = T)


pred.dat <- data.frame(year = new.chla.dat$year,
                             chla_ratio = chla.pred$fit,
                             LCI = chla.pred$fit - 1.96 * chla.pred$se.fit,
                             UCI = chla.pred$fit + 1.96 * chla.pred$se.fit,
                             group = "chla") %>%
  mutate(year = as.numeric(as.character(year)))

ggplot(pred.dat, aes(year, log_abundance)) +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  facet_wrap(~group, scales = "free_y", ncol = 1)


# Save product
write.csv(pred.dat, "./Data/summarized_chla.csv")