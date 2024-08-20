# notes ----

#1) Calculate biomass of pacific cod within imm snow crab 50th percentile home range 
#2) Calculate biomass of subarctic fish complex within imm snow crab 50th percentile home range

# Author: Erin Fedewa

# load ----
library(tidyverse)
theme_set(theme_bw())

###################################################
# Generate groundfish master csv

# #Function to import data and filter for only fish 
# import <- function(filename) {
#   ebs <- read_csv(filename)
#   ebs %>%
#     filter(SID %in% c(1:31550)) 
# }
# 
# #Add all bottom trawl data files
# ebs82 <- import("./Data/Groundfish Catch Data/ebs1982_1984.csv")
# ebs85 <- import("./Data/Groundfish Catch Data/ebs1985_1989.csv")
# ebs90 <- import("./Data/Groundfish Catch Data/ebs1990_1994.csv")
# ebs95 <- import("./Data/Groundfish Catch Data/ebs1995_1999.csv")
# ebs00 <- import("./Data/Groundfish Catch Data/ebs2000_2004.csv")
# ebs05 <- import("./Data/Groundfish Catch Data/ebs2005_2008.csv")
# ebs09 <- import("./Data/Groundfish Catch Data/ebs2009_2012.csv")
# ebs13 <- import("./Data/Groundfish Catch Data/ebs2013_2016.csv")
# ebs17 <- import("./Data/Groundfish Catch Data/ebs2017_2018.csv")
# ebs19 <- import("./Data/Groundfish Catch Data/ebs2019.csv")
# ebs21 <- import("./Data/Groundfish Catch Data/ebs2021.csv")
# ebs22 <- import("./Data/Groundfish Catch Data/ebs2022.csv")
# 
# # combine datasets now and save output
# bind_rows(ebs82, ebs85, ebs90, ebs95, ebs00, ebs05, ebs09, ebs13, ebs17, ebs19, ebs21, ebs22) %>%
#   write_csv("./Data/groundfish_timeseries.csv")

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
  summarise(station = length(unique(STATION))) -> annual_cod_n

print(annual_cod_n, n = 50)

#Looks like pre-1993 we've got missing stations 
gf_catch %>%
  filter(SID %in% c(21720),
         STATION %in% imm_stations) -> check

sum(check$WTCPUE == 0)

## make sure zero-catch stations are included

# load stratum data
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

# separate cod data
gf_catch %>%
  filter(SID %in% c(21720),
         STATION %in% imm_stations) -> cod

# join with stratum stations to include 0 catches 
stratum <- sc_strata %>%
  select(STATION_ID, SURVEY_YEAR) %>% 
  rename(STATION = STATION_ID,
         YEAR = SURVEY_YEAR) %>%
  filter(YEAR >= 1982,
         STATION %in% imm_stations)

cod <- left_join(stratum, cod) 

# replace NA catch with 0
change <- is.na(cod$WTCPUE)
cod$WTCPUE[change] <- 0


#Plot CPUE distribution of cod by year 
  ggplot(cod, aes(x=WTCPUE^0.25)) +
  geom_histogram() +
  facet_wrap(~YEAR) 

# compare distribution of 4th root and log-transformed catches
  ggplot(cod, aes(x=WTCPUE^0.25)) +
  geom_histogram(fill = "grey", color = "black")

  ggplot(cod, aes(x=log(WTCPUE+1))) +
  geom_histogram(fill = "grey", color = "black")


#Time series of mean cod CPUE within immature core region
cod %>%
  group_by(YEAR) %>%
  summarise(mean_cod_CPUE = mean(log(WTCPUE+1))) -> mean_cod_cpue

# limit to years with at least 180 stations
annual_cod_n <- annual_cod_n %>%
  filter(station >= 180)

mean_cod_cpue <- mean_cod_cpue %>%
  filter(YEAR %in% annual_cod_n$YEAR)

ggplot(mean_cod_cpue, aes(YEAR, mean_cod_CPUE)) +
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
  summarise(sum_CPUE = sum(WTCPUE))  -> arctic_cpue

arctic_cpue <- left_join(stratum, arctic_cpue)

# replace NA with 0
change <- is.na(arctic_cpue$sum_CPUE)
arctic_cpue$sum_CPUE[change] <- 0

#Time series of mean cod CPUE within immature core region
arctic_cpue %>%
  group_by(YEAR) %>%
  summarise(mean_arctic_CPUE = mean(log(sum_CPUE+1))) -> mean_arctic_cpue

# limit to years with at least 180 stations
annual_cod_n <- annual_cod_n %>%
  filter(station >= 180)

mean_arctic_cpue <- mean_arctic_cpue %>%
  filter(YEAR %in% annual_cod_n$YEAR)

ggplot(mean_arctic_cpue, aes(YEAR, mean_arctic_CPUE)) +
  geom_line() +
  geom_point()




#Join cod and arctic fish datasets
mean_cod_cpue %>%
  left_join(mean_arctic_cpue) -> final
  
#Write csv for median pcod and arctic fish CPUE within imm snow crab core habitat 
write.csv(final, file = "./Data/groundfish_mean_cpue.csv")
