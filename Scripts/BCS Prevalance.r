# notes ----]
# Calculate population prevalence of disease code 2 (i.e. visually positive BCS) 

# Author: Erin Fedewa
# last updated: 2023/3/23

# load ----
library(tidyverse)
library(ggridges)

##############################################

## EBS haul data ----
crab_ebs <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
strata_ebs <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

########################################
#compute cpue by disease code 2 for mature crab only 
crab_ebs %>% 
  filter(HAUL_TYPE == 3,
         ((SEX == 2 & CLUTCH_SIZE >= 1) |
            (SEX == 1 & WIDTH_1MM >= 95))) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
         bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue 
#add zero catch stations to bcs=T/F datasets
cpue %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_ebs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_ebs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       filter(HAUL_TYPE ==3) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch
#Now our BCS categories should each contain zero catch stations so we can group by 
#bcs category to calculate abundances 

#abundance by BCS categories
catch %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  filter(!is.na(bcs)) %>%
  #calculate prevalence 
  group_by(YEAR) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+(Total_abun[bcs==TRUE])))*100) %>%
  mutate(Maturity = rep("Mature")) -> matprev

###
#compute cpue by disease code 2 for immature crab only 
crab_ebs %>% 
  filter(HAUL_TYPE == 3,
         ((SEX == 1 & WIDTH_1MM < 95 ) |
            (SEX == 2 & CLUTCH_SIZE == 0))) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
         bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue_imm
#add zero catch stations to bcs=T/F datasets
cpue_imm %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_ebs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue_imm %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_ebs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       filter(HAUL_TYPE ==3) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch_imm

#abundance by BCS categories
catch_imm %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  filter(!is.na(bcs)) %>%
  #calculate prevalence 
  group_by(YEAR) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+(Total_abun[bcs==TRUE])))*100) %>%
  mutate(Maturity = rep("Immature")) -> immprev

###
#compute cpue by disease code 2 for entire population
crab_ebs %>% 
  filter(HAUL_TYPE == 3,
         SEX %in% c(1,2)) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
         bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue_pop
#add zero catch stations to bcs=T/F datasets
cpue_pop %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_ebs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue_pop %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_ebs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       filter(HAUL_TYPE ==3) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch_pop

#abundance by BCS categories
catch_pop %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  filter(!is.na(bcs)) %>%
  #calculate prevalence 
  group_by(YEAR) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+(Total_abun[bcs==TRUE])))*100) %>%
  mutate(Maturity = rep("Population")) -> popprev

#join datasets
popprev %>%
  full_join(immprev) %>% 
  full_join(matprev) %>%
  filter(YEAR > 1988) %>%
  pivot_wider(names_from = "Maturity", values_from = "Perc_Prevalance") -> bcs

#Write csv 
write.csv(bcs, file="./Data/bcs_prev.csv")

###################################################
#Combined Plot 
bcs %>%
  pivot_longer(2:4, names_to="Maturity", values_to="Perc_prev") %>%
  ggplot(aes(x = YEAR, y = Perc_prev, group = as.factor(Maturity))) +
  geom_point(aes(colour = Maturity), size=3) +
  geom_line(aes(colour = Maturity), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  theme(legend.text=element_text(size=11)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
  theme(legend.title= element_blank()) +
  scale_x_continuous(limits = c(1988, 2022), breaks = seq(1990, 2022, 5))
  

 