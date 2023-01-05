# notes ----

#1) Calculate biomass of pacific cod within imm snow crab 50th percentile home range 
#2) Calculate biomass of subarctic fish complex within imm snow crab 50th percentile home range

# Author: Erin Fedewa

#NOTES: 
#Looks like there's some missing station data prior to 1993, may be 
  #worth imputing CPUE's
#This dataset only includes 1982+, data request in for earlier yrs, and CPUE's
  #for missing stations will need to be imputed 
#Need to follow up on median CPUE approach and instead calculate station-level 
  #cod biomass...though maybe our approach for defining cod-snow crab habitat  
  #overlap needs some rethinking 

# load ----
library(tidyverse)

###################################################
# Generate groundfish master csv

#Function to import data and filter for only fish 
import <- function(filename) {
  ebs <- read_csv(filename)
  ebs %>%
    filter(SID %in% c(1:31550)) 
}

#Add all bottom trawl data files
ebs82 <- import("./Data/Groundfish Catch Data/ebs1982_1984.csv")
ebs85 <- import("./Data/Groundfish Catch Data/ebs1985_1989.csv")
ebs90 <- import("./Data/Groundfish Catch Data/ebs1990_1994.csv")
ebs95 <- import("./Data/Groundfish Catch Data/ebs1995_1999.csv")
ebs00 <- import("./Data/Groundfish Catch Data/ebs2000_2004.csv")
ebs05 <- import("./Data/Groundfish Catch Data/ebs2005_2008.csv")
ebs09 <- import("./Data/Groundfish Catch Data/ebs2009_2012.csv")
ebs13 <- import("./Data/Groundfish Catch Data/ebs2013_2016.csv")
ebs17 <- import("./Data/Groundfish Catch Data/ebs2017_2018.csv")
ebs19 <- import("./Data/Groundfish Catch Data/ebs2019.csv")
ebs21 <- import("./Data/Groundfish Catch Data/ebs2021.csv")
ebs22 <- import("./Data/Groundfish Catch Data/ebs2022.csv")

# combine datasets now and save output
bind_rows(ebs82, ebs85, ebs90, ebs95, ebs00, ebs05, ebs09, ebs13, ebs17, ebs19, ebs21, ebs22) %>%
  write_csv("./Data/groundfish_timeseries.csv")

###########################################

## Groundfish specimen data ----
gf_catch <- read.csv("./Data/groundfish_timeseries.csv")

#Core immature snow crab stations
imm_area <- read.csv("./Data/imm_area_50perc.csv")

imm_stations <- pull(imm_area, GIS_STATION)

##################
#Biomass of pacific cod within imm snow crab 50th percentile home range 

#Num of stations with catch data each yr within 187 station core region
gf_catch %>%
  filter(STATION %in% imm_stations) %>%
  group_by(YEAR) %>%
  summarise(station = length(unique(STATION))) %>%
  print(n=50)
#Looks like pre-1993 we've got missing stations 

#Plot CPUE distribution of cod by year 
gf_catch %>%
  filter(SID %in% c(21720),
         STATION %in% imm_stations) %>%
  ggplot(aes(x=WTCPUE^0.25)) +
  geom_histogram() +
  facet_wrap(~YEAR) +
  geom_vline(aes(xintercept=median(WTCPUE)), linetype="dashed", size=1) 

#Time series of median cod CPUE within immature core region
gf_catch %>%
  filter(SID %in% c(21720),
         STATION %in% imm_stations) %>%
  group_by(YEAR) %>%
  summarise(med_cod_CPUE = median(WTCPUE^0.25)) -> med_cod_cpue

ggplot(med_cod_cpue, aes(YEAR, med_cod_CPUE)) +
  geom_line() +
  geom_point()

###################################################################
#Biomass of subarctic fish complex within imm snow crab 50th percentile home range

#Plot CPUE distribution of arctic fish by year 
gf_catch %>%
  #alaska plaice, arctic cod, bering flounder, butterfly sculpin, greenland turbot
  #marbled eelpout, wattled eelpout 
  filter(SID %in% c(10285,21725,10140,21348,10115,66045,24185,10212),
         STATION %in% imm_stations) %>%
  group_by(YEAR, STATION) %>% 
  summarise(sum_CPUE = sum(WTCPUE)) %>%
  ggplot(aes(x=sum_CPUE^0.25)) +
  geom_histogram() +
  facet_wrap(~YEAR) +
  # xlim(0, 100) +
  geom_vline(aes(xintercept=median(sum_CPUE^0.25)), linetype="dashed", size=1)

#Timseries of median arctic fish CPUE within immature core region
gf_catch %>%
  filter(SID %in% c(10285,21725,10140,21348,10115,66045,24185,10212),
         STATION %in% imm_stations) %>%
  group_by(YEAR, STATION) %>% 
  summarise(sum_CPUE = sum(WTCPUE)) %>%
  group_by(YEAR) %>%
  summarise(med_arctic_CPUE = median(sum_CPUE^0.25)) -> med_arctic_cpue

ggplot(med_arctic_cpue, aes(YEAR, med_arctic_CPUE)) +
  geom_line() +
  geom_point()

#Join cod and arctic fish datasets
med_cod_cpue %>%
  left_join(med_arctic_cpue) -> final
  
#Write csv for median pcod and arctic fish CPUE within imm snow crab core habitat 
write.csv(final, file = "./Data/groundfish_med_cpue.csv")
