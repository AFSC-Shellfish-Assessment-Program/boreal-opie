# notes ----

#1) Calculate biomass of pacific cod within imm snow crab 50th percentile home range 
#2) Calculate biomass of subarctic fish complex within imm snow crab 50th percentile home range

# Author: Erin Fedewa

# load ----
library(tidyverse)

###################################################
# Generate groundfish master csv ----

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

# combine datasets now and save output
bind_rows(ebs82, ebs85, ebs90, ebs95, ebs00, ebs05, ebs09, ebs13, ebs17, ebs19, ebs21) %>%
  write_csv("./Data/groundfish_timeseries.csv")

###########################################

## Groundfish specimen data ----
gf_catch <- read.csv("./Data/groundfish_timeseries.csv")

#Core immature snow crab stations
imm_area <- read.csv("./Data/imm_area_50perc.csv")

imm_stations <- pull(imm_area, GIS_STATION)

##################

#Biomass of pacific cod within imm snow crab 50th percentile home range 
gf_catch %>%
  filter(SID %in% c(21720, 21721, 21722),
         STATION %in% imm_stations)

#How to get from wtcpue to biomass?
#Where is strata table with areas of strata? 




#NOTE: Need to append pre-1982 data and impute missing station CPUEs, then re-run 





#Biomass of subarctic fish complex within imm snow crab 50th percentile home range

#Arctic Species- may be worthwhile to re-evaluate this complex 
gf_catch %>%
  #alaska plaice, arctic cod, bering flounder, butterfly sculpin, greenland turbot
  #marbled eelpout, wattled eelpout 
  filter(SID %in% c(10285,21725,10140,21348,10115,66045,24185,10212)) %>%
  


