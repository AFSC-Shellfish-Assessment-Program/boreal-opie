# notes ----
# Pcod/Snow Crab Spatial Overlap in EBS and NBS 
  #EBS timeseries csv is combined EBS groundfish catch data 1982-2019 only 
      #2021 groundfish data not yet ready at time of script run

##Script needs to be updated by size classes of snow crab!

# Erin Fedewa
# Last Updated 7/17/20

# load ----

library(tidyverse)
library(cowplot)
library(mgcv)
library(FNGr)
library(latticeExtra)

# data managment ----

ebs <- read_csv("./Data/ebs_timeseries.csv")
head(ebs)
nbs <- read_csv("./Data/nbs_timeseries.csv")
head(nbs)

#Calculate % overlap for EBS timeseries 
ebs %>%
  group_by(YEAR) %>% 
  mutate(TOTAL_STATIONS = n_distinct(STATION)) %>%
  filter(SID %in% c("21720", "68580")) %>%
  select(YEAR, STATION, SID, WTCPUE, TOTAL_STATIONS) %>%
  mutate(SID = case_when(SID == 68580 ~ "CRAB",
                         SID == 21720 ~ "COD")) %>%
  group_by(YEAR, STATION, SID) %>%
  pivot_wider(names_from = SID, values_from = WTCPUE) %>%
  group_by(YEAR) %>%
  # method 1 -  % of total stations that include both cod and snow crab
  # method 2 - % of positive snow crab stations that included cod
  summarise(METHOD_1 = sum((CRAB > 0 & COD > 0), na.rm = T) / mean(TOTAL_STATIONS) * 100,
            METHOD_2 = sum((CRAB > 0 & COD > 0), na.rm = T) / sum((CRAB > 0), na.rm = T) * 100) -> overlap

 
#EBS Plot 
 
  #Method 1
ggplot(aes(x = YEAR, y = METHOD_1), data = overlap) +
  geom_point() +
  geom_line() +
  labs(y = expression(atop("EBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

#Method 2
ggplot(aes(x = YEAR, y = METHOD_2), data = overlap) +
  geom_point() +
  geom_line() +
  labs(y = expression(atop("EBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  theme_bw()+
  theme(panel.grid = element_blank())
  

#Calculate % overlap for NBS timeseries
nbs %>%
  group_by(YEAR) %>% 
  mutate(TOTAL_STATIONS = n_distinct(STATIONID)) %>%
  filter(SPECIES_CODE %in% c("21720", "68580")) %>%
  select(YEAR, STATIONID, SPECIES_CODE, wCPUE, TOTAL_STATIONS) %>%
  mutate(SPECIES_CODE = case_when(SPECIES_CODE == 68580 ~ "CRAB",
                                  SPECIES_CODE == 21720 ~ "COD")) %>%
  group_by(YEAR, STATIONID, SPECIES_CODE) %>%
  pivot_wider(names_from = SPECIES_CODE, values_from = wCPUE) %>%
  group_by(YEAR) %>%
  # method 1 -  % of total stations that include both cod and snow crab
  # method 2 - % of positive snow crab stations that included cod
  summarise(METHOD_1 = sum((CRAB > 0 & COD > 0), na.rm = T) / mean(TOTAL_STATIONS) * 100,
            METHOD_2 = sum((CRAB > 0 & COD > 0), na.rm = T) / sum((CRAB > 0), na.rm = T) * 100) %>%
  add_row(YEAR = 2011:2016, METHOD_1 = NA, METHOD_2 = NA) -> nbs_overlap  #Add NA's for missing yrs 

#Let's use NBS Method 2 for draft snow crab indicator- EBS likely meaningless b/c needs to be 
  #subset for immature snow crab/large cod only 

write.csv(nbs_overlap, file = "./Output/nbs_overlap.csv")

#NBS Plot 

  #Method 1
ggplot(aes(x = YEAR, y = METHOD_1), data = nbs_overlap, na.rm = T) +
  geom_point() +
  geom_line() +
  labs(y = expression(atop("NBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

#Method 2
  ggplot(aes(x = YEAR, y = METHOD_2), data = nbs_overlap, na.rm = T) +
  geom_point(color = "#E69F00", size = 6) +
  geom_line() +
  labs(y = expression(atop("Snow crab-Pacific cod spatial overlap (%)")), x = "")+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey85", size = 0.2)) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11)) -> overlap

