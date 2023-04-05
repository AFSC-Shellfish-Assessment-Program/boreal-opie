#Double checking mean CPUE estimates produced from 
#abundance_imputation.R script lines 14-161 

#Author: Erin Fedewa
#last update: 4/4/2022

############################
## load
library(tidyverse)
library(patchwork)

## data
#EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

#Standard 375 survey stations
sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  pull(STATION_ID) ->standard

#################################
#Snow crab male 30-59 mean CPUE (no weighting)
sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4),
         SEX == 1,
         WIDTH >= 30 & WIDTH < 60,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #add zero catch stations
  right_join(sc_catch %>% 
               mutate(YEAR = floor(CRUISE/100)) %>%
               filter(HAUL_TYPE %in% c(3,4)) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)) %>%
  replace_na(list(CPUE = 0)) %>%
  mutate(log_cpue = log(CPUE +1)) %>%
  #filter for standard stations (1979 has 70 extra non-standard stations)
  filter(GIS_STATION %in% standard) -> cpue

#Calculate mean CPUE/yr across all stations
cpue %>%
  group_by(YEAR) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            mean_cpue =mean(CPUE),
            log_mean_cpue = log(mean_cpue)) -> mean_cpue

#plot
mean_cpue %>%
  ggplot(aes(YEAR, mean_cpue)) +
  geom_point() +
  geom_line() -> cpue_plot
mean_cpue %>%
  ggplot(aes(YEAR, log_mean_cpue)) +
  geom_point() +
  geom_line() -> logcpue_plot
mean_cpue %>%
  ggplot(aes(YEAR, mean_log_cpue)) +
  geom_point() +
  geom_line() -> meanlogcpue_plot
cpue_plot + logcpue_plot + meanlogcpue_plot

#Cleaning up data for imputation
#a) Remove stations that have never caught 30-59mm male snow crab 
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(station_mean = mean(CPUE)) %>% 
  filter(station_mean > 0) %>% #58 stations dropped
  #b) Remove stations < 5th CPUE percentile range 
  filter(station_mean > quantile(station_mean, 0.05)) %>% #16 additional stations dropped
  pull(GIS_STATION) -> keepers 

#Final dataset for imputation
cpue %>%
  filter(GIS_STATION %in% keepers) -> dat

