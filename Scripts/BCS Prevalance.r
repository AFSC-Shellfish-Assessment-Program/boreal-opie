# notes ----]
# Calculate population prevalence of disease code 2 (i.e. visually positive BCS) 

# Author: Erin Fedewa
# last updated: 2023/1/5

# load ----
library(tidyverse)
library(ggridges)

##############################################

## EBS haul data ----
crab_ebs <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
strata_ebs <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

########################################

#compute cpue by size-sex group and disease code 2 for each station
crab_ebs %>% 
  filter(HAUL_TYPE == 3) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
         size_sex = ifelse(SEX == 1 & WIDTH_1MM < 95 | SEX == 2 & CLUTCH_SIZE == 0, "Immature",
                           ifelse(SEX == 2 & CLUTCH_SIZE >= 1 | SEX == 1 & WIDTH_1MM >= 95,  "Mature", NA)),
         bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, size_sex, bcs) %>%
  summarise(num_crab = round(sum(SAMPLING_FACTOR,na.rm = T))) %>%
  #add in data field for total population
  bind_rows(crab_ebs %>% 
              filter(HAUL_TYPE == 3) %>% 
              mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
                     bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T),
                     size_sex = "All") %>%
              group_by(YEAR, GIS_STATION, AREA_SWEPT, size_sex, bcs) %>%
              summarise(num_crab = round(sum(SAMPLING_FACTOR,na.rm = T)))) %>%
  filter(!is.na(size_sex)) %>%
  ungroup() %>%
  mutate(cpue = num_crab / AREA_SWEPT) -> catch

#abundance by size sex/BCS categories
catch %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, size_sex, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, size_sex, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  ungroup() -> abund

#Calculate BCS prevalence 
abund %>%
  group_by(YEAR,size_sex) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+(Total_abun[bcs==TRUE])))*100) -> bcs

###################################################
 
#Combined Plot 
#NOTE: 2021 has zero prev for 2 size/sex categories
bcs %>%
  filter(size_sex != "NA") %>%
  ggplot(aes(x = YEAR, y = Perc_Prevalance, group = as.factor(size_sex))) +
  geom_point(aes(colour = size_sex), size=3) +
  geom_line(aes(colour = size_sex), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  theme(legend.text=element_text(size=11)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
  theme(legend.title= element_blank()) +
  scale_x_continuous(limits = c(1988, 2022), breaks = seq(1990, 2022, 5))

#Write csv for BCS snow crab indicator
bcs %>%
  filter(size_sex != "NA") %>%
  pivot_wider(names_from = size_sex, values_from = Perc_Prevalance) %>%
  write.csv(file="./Data/bcs_prev.csv")
  

 