# notes ----

# Determine core habitat of immature snow crab in EBS
  #1) Determine stations that compose average core habitat across the long-term timeseries 
  #2) Determine average bottom temperature of core habitat within each yr
  
# Author: Erin Fedewa

#NOTES:
#Given that stations differed pre 1980, only 1980+ dataset was used to calculate 
  #avg core habitat across timeseries
#Bottom temps need to be data corrected and imputed for final timeseries dataset

# load ----
library(tidyverse)
library(ggmap)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")
sc_catch <- read.csv("./Data/CRABHAUL_OPILIO.csv")


#EBS strata data ----
sc_strata <- read_csv("./Data/crabstrata_opilio.csv")
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")


###################################################
# data exploration ----

#Stations sampled in each year
sc_catch <- sc_catch %>%
  mutate(AKFIN_SURVEY_YEAR = lubridate::year(lubridate::parse_date_time(START_DATE, orders = "mdy"))) 

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
                    distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
                    rename_all(~c("GIS_STATION", "AKFIN_SURVEY_YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
              replace_na(list(CPUE = 0)) -> cpue
#####
  
#1) Determine stations that compose average core habitat across the long-term timeseries 

#Histogram with percentiles 
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(AVG_CPUE = mean(CPUE)) -> data
p10 <- quantile(data$AVG_CPUE, 0.10)
p20 <- quantile(data$AVG_CPUE, 0.20)
p50 <- quantile(data$AVG_CPUE, 0.50)

ggplot(data=data, aes(x=AVG_CPUE)) +
  geom_histogram(bins=15) +
  geom_vline(aes(xintercept=p10), linetype="dashed", size=1) +
  geom_vline(aes(xintercept=p20), linetype="dashed", size=1) +
  geom_vline(aes(xintercept=p50), linetype="dashed", size=1) +
  theme_bw()

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
#Lets go with the 50th percentile for defining core immature area 

#Join lat/long back in to perc50 dataset and plot
sc_strata %>%
      filter(SURVEY_YEAR == 2021) %>% #Just selecting a yr when all stations were sampled
      select(STATION_ID, LATITUDE, LONGITUDE) %>%
      dplyr::rename(GIS_STATION = STATION_ID) %>%
      right_join(perc50) -> perc50_core

# save for other analyses
write.csv(perc50_core, "output/core_stations.csv", row.names = F)


#Plot....this takes awhile to generate a stamen map FYI!!
get_stamenmap(bbox=c(-180, 54, -157, 63),
  maptype="toner-lite") %>%
  ggmap() +
  geom_point(data= perc50_core, aes(x = LONGITUDE, y = LATITUDE, size=AVG_CPUE)) +
  labs(x = "Longitude", y = "Latitude") 
  
ggsave("./Figs/core_range_map.png", width = 6, height = 4, units = 'in')
#Write csv for stations in 50th percentile of avg CPUE  
write.csv(perc50_core, file="./Data/imm_area_50perc.csv")

######

#2) Determine average bottom temperature of core habitat within each yr
cpue %>%
  group_by(AKFIN_SURVEY_YEAR) %>%
  filter(CPUE > quantile(CPUE, 0.50)) %>%
  summarise(AVG_TEMP = mean(GEAR_TEMPERATURE, na.rm = T)) -> avg_temp

#Write csv for avg bottom temp within 50th percentile home range   
write.csv(avg_temp, file="./Data/avg_temp_50perc.csv")












  
  

