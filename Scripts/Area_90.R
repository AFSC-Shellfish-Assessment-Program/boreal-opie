# notes ----

# Determine stations that make up 90% of immature/immature female snow crab CPUE 
# Author: Erin Fedewa

# load ----
library(tidyverse)
library(ggmap)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv", skip=5)

#EBS strata data ----
strata_ebs <- read_csv("./Data/crabstrata_opilio.csv")

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

# data  ----

#CPUE by station for all immature snow crab 
sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX %in% 1:2) %>%
  mutate(MAT_SEX = case_when((SEX == 1 & WIDTH < 95 & SHELL_CONDITION <= 2) ~ "Immature Male",
                             (SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature Female")) %>%
  filter(MAT_SEX %in% c("Immature Male", "Immature Female")) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
# compute 90% of habitat for immature snow crab
# i.e. list of stations contributing to 90% of cumulative cpue within a given year 
  group_by(AKFIN_SURVEY_YEAR) %>%
  filter(CPUE > quantile(CPUE, 0.10)) -> area90 #return stations with density > 10th percentile 
  
#Plot
area90 %>%
  filter(AKFIN_SURVEY_YEAR > 2010) %>%
  group_by(AKFIN_SURVEY_YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  distinct(GIS_STATION) %>%
  ggplot() +
  geom_point(aes(x = MID_LONGITUDE, y = MID_LATITUDE), size=.5) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~AKFIN_SURVEY_YEAR)

#right skewed data (should we be including 0 catch stations)-choose different percentile?
  #uninformative b/c caught at  alot of stations 
  
  
  
  



# function to compute 
f_area90_est <- function(x){
  x %>%
    group_by(AKFIN_SURVEY_YEAR) %>%
    arrange(-CPUE) %>% #sort by cpue (large:small)
    mutate(PROP_CPUE = CPUE/sum(CPUE),  #calculate the proportion of total cpue for each station
           cum_cpue = cumsum(PROP_CPUE)) %>%  
    filter(cum_cpue <= 0.90) %>% #T if in area90, F if not
    count() %>%
    mutate(d90 = (n + 1) * 401) %>% #add 1 station to n to push over 95%
    pull(d95)
  }

# do the estimation for immature snow crab
cpue_long %>%
  filter(!(Station %in% corner)) %>% #exclude corner stations
  group_by(Station, size_sex) %>%
  filter(n() == 32) %>% #include only stations sampled every year in 32 year timeseries
  ungroup() %>% 
  nest(-year, -size_sex) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest() %>%
  group_by(year, size_sex) %>%
  summarise(cpue = sum(num_crab) / sum(area_swept), # add a column for total cpue of each group in each year
            d95 = mean(d95))->d95 # take 'mean' just to get one value (they are all the same)

# do the estimation for immature female snow crab only
cpue_long %>%
  filter(!(Station %in% corner)) %>% #exclude corner stations
  group_by(Station, size_sex) %>%
  filter(n() == 32) %>% #include only stations sampled every year in 32 year timeseries
  ungroup() %>% 
  nest(-year, -size_sex) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest() %>%
  group_by(year, size_sex) %>%
  summarise(cpue = sum(num_crab) / sum(area_swept), # add a column for total cpue of each group in each year
            d95 = mean(d95))->d95 # take 'mean' just to get one value (they are all the same)

cbind()

write.csv(d95, file="./Output/D95_output.csv")
  
