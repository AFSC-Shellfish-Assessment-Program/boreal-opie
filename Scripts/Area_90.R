# notes ----

# Determine core habitat of immature snow crab in EBS
  #1) Determine stations that compose average core habitat across the long-term timeseries 
  #2) Determine average bottom temperature of core habitat within each yr
  
# Author: Erin Fedewa

# load ----
library(tidyverse)
library(ggmap)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv", skip=5)

#EBS strata data ----
sc_strata <- read_csv("./Data/crabstrata_opilio.csv")

###################################################
# data exploration ----

#Stations sampled in each year
sc_catch %>%
  group_by(AKFIN_SURVEY_YEAR) %>%
  summarise(num_stations = length(unique(GIS_STATION)))

#Plot pre-standardization data
sc_catch %>%
  filter(AKFIN_SURVEY_YEAR < 1988) %>%
  group_by(AKFIN_SURVEY_YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  distinct(GIS_STATION) %>%
ggplot() +
  geom_point(aes(x = MID_LONGITUDE, y = MID_LATITUDE), size=.5) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~AKFIN_SURVEY_YEAR)

#Earliest yrs don't survey prime snow crab habitat so 90% of habitat index
  #will be biased in these years 

###################################################

#Calculate CPUE by station for all immature snow crab 
sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX %in% 1:2,
         AKFIN_SURVEY_YEAR > 1979) %>%
  mutate(MAT_SEX = case_when((SEX == 1 & WIDTH < 95 & SHELL_CONDITION <= 2) ~ "Immature Male",
                             (SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature Female")) %>%
  filter(MAT_SEX %in% c("Immature Male", "Immature Female")) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE, GEAR_TEMPERATURE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #join to zero catch stations
  right_join(sc_strata %>%
                    filter(SURVEY_YEAR > 1979) %>%
                    distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA_SQ_NM) %>%
                    rename_all(~c("GIS_STATION", "AKFIN_SURVEY_YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
              replace_na(list(CPUE = 0)) -> cpue
#####
  
#1) Determine stations that compose average core habitat across the long-term timeseries 

#stations in 10-100 CPUE percentile range
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(AVG_CPUE = mean(CPUE)) %>%
  filter(AVG_CPUE > quantile(AVG_CPUE, 0.10)) -> perc10 #331 stations 
  
#stations in 20-100 CPUE percentile range
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(AVG_CPUE = mean(CPUE)) %>%
  filter(AVG_CPUE > quantile(AVG_CPUE, 0.20)) -> perc20 #300 stations

#stations in 50-100 CPUE percentile range
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(AVG_CPUE = mean(CPUE)) %>%
  filter(AVG_CPUE > quantile(AVG_CPUE, 0.50)) -> perc50 #187 stations
  
#Write csv for stations in 50th percentile of avg CPUE  
write.csv(perc50, file="./Data/imm_area_50perc.csv")

######

#2) Determine average bottom temperature of core habitat within each yr
cpue %>%
  group_by(AKFIN_SURVEY_YEAR) %>%
  filter(CPUE > quantile(CPUE, 0.50)) %>%
  summarise(AVG_TEMP = mean(GEAR_TEMPERATURE, na.rm = T)) -> avg_temp

#Write csv for avg bottom temp within 50th percentile home range   
write.csv(avg_temp, file="./Data/avg_temp_50perc.csv")












  
  

