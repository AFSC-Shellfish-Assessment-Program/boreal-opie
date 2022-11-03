# notes ----

# Calculate immature snow crab abundance in EBS

#NOTE: There are no pre-1980 stratum tables on akfin so abundances were only 
  #calculated for 1980+

# Author: Erin Fedewa

# load ----
library(tidyverse)
library(ggmap)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

###################################################

# define survey year in catch
sc_catch <- sc_catch %>%
  mutate(YEAR = floor(CRUISE/100))

# rename year in strata
sc_strata <- sc_strata %>%
  rename(YEAR = SURVEY_YEAR)

# check size of stratum in each year
check <- sc_strata %>%
  group_by(YEAR, DISTRICT) %>%
  summarise(area = mean(TOTAL_AREA)) %>% 
  group_by(YEAR) %>%
  summarise(area = sum(area))

ggplot(check, aes(YEAR, area)) + 
  geom_line() +
  geom_point() # big differences!


# calculate CPUE by station for immature snow crab 
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         YEAR > 1979,
         WIDTH >= 40 & WIDTH <= 95) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 



%>%
  #join to zero catch stations
  right_join(sc_strata %>%
               filter(YEAR > 1979) %>%
               distinct(YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
               rename_all(~c("GIS_STATION", "YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
  replace_na(list(CPUE = 0)) 

# check stratum size at this point
check <- scratch %>%
  group_by(YEAR, STRATUM) %>%
  summarise(area = mean(TOTAL_AREA_SQ_NM)) %>% 
  group_by(YEAR) %>%
  summarise(area = sum(area))

ggplot(check, aes(YEAR, area)) + 
  geom_line() +
  geom_point() # same result!

# plot 1980 - 1992 samples
ggplot(filter(scratch, YEAR < 1993), aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)


%>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA_SQ_NM) %>%
  summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA_SQ_NM))) %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE = sum(ABUNDANCE)) -> imm_abundance


# calculate CPUE by station for immature snow crab 
sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         YEAR > 1979,
         WIDTH >= 40 & WIDTH <= 95) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #join to zero catch stations
  right_join(sc_strata %>%
               filter(YEAR > 1979) %>%
               distinct(YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
               rename_all(~c("GIS_STATION", "YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
  replace_na(list(CPUE = 0)) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA_SQ_NM) %>%
  summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA_SQ_NM))) %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE = sum(ABUNDANCE)) -> imm_abundance


ggplot(imm_abundance, aes(YEAR, log(ABUNDANCE))) +
  geom_point() +
  geom_line()

#Write csv as output 
write.csv(imm_abundance, file = "./Data/imm_abun.csv")
