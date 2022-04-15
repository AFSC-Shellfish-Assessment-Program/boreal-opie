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
sc_catch <- read.csv("./Data/crabhaul_opilio.csv", skip=5)

#EBS strata data ----
sc_strata <- read_csv("./Data/crabstrata_opilio.csv")

###################################################

#Calculate CPUE by station for immature snow crab 
sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX %in% 1:2,
         AKFIN_SURVEY_YEAR > 1979) %>%
  mutate(MAT_SEX = case_when((SEX == 1 & WIDTH < 95) ~ "Immature Male",
                             (SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature Female")) %>%
  filter(MAT_SEX %in% c("Immature Male", "Immature Female")) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #join to zero catch stations
  right_join(sc_strata %>%
               filter(SURVEY_YEAR > 1979) %>%
               distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA_SQ_NM) %>%
               rename_all(~c("GIS_STATION", "AKFIN_SURVEY_YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
  replace_na(list(CPUE = 0)) %>%
  #Scale to abundance by strata
  group_by(AKFIN_SURVEY_YEAR, STRATUM, TOTAL_AREA_SQ_NM) %>%
  summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA_SQ_NM))) %>%
  group_by(AKFIN_SURVEY_YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6) -> imm_abundance

#Write csv as output 
write.csv(imm_abundance, file = "./Data/imm_abun.csv")
