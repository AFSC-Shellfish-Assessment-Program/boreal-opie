# notes ----]
# Calculate population prevalence of disease code 2 (i.e. visually positive BCS) 

# Erin Fedewa
# last updated: 2021/12/9


# load ----
library(tidyverse)
library(ggridges)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv", skip=5)

#EBS strata data ----
strata_ebs <- read_csv("./Data/crabstrata_opilio.csv")

########################################

#Snow crab haul data 
sc_catch %>%
  filter(AKFIN_SURVEY_YEAR %in% c(1988:2021)) %>%
  dplyr::select(AKFIN_SURVEY_YEAR, CRUISE, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED, NET_WIDTH, MID_LATITUDE, MID_LONGITUDE,
         HAUL_TYPE, PERFORMANCE, DISEASE_CODE) -> crab_ebs
names(crab_ebs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance", "disease_code")

## EBS strata data 
strata_ebs %>%
  filter(SURVEY_YEAR %in% c(1988:2021)) %>%
  select(STATION_ID, DISTRICT, TOWS, STRATUM, TOTAL_AREA_SQ_NM, SURVEY_YEAR) -> strata_ebs
names(strata_ebs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

################################################

#compute cpue by size-sex group and disease code 2 for each station
crab_ebs %>% 
  filter(haul_type != 17) %>% 
  mutate(size_sex = ifelse(sex == 1 & cw < 95, "Immature_male",
                        ifelse(sex == 2 & clutch >= 1, "Mature_female",
                               ifelse(sex == 1 & cw >= 95, "Mature_male",
                                       ifelse(sex == 2 & clutch == 0, "Immature_female", NA)))),
         area_swept = distance_fished * (net_width / 1000) * 0.29155335,
         bcs = ifelse(disease_code != 2 | is.na(disease_code), F, T)) %>%
  group_by(year, Station, lat, lon, area_swept, size_sex, bcs) %>%
  summarise(num_crab = round(sum(sample_factor,na.rm = T))) %>%
  bind_rows(crab_ebs %>% 
              filter(haul_type != 17) %>% 
              mutate(area_swept = distance_fished * (net_width / 1000) * 0.29155335,
                     bcs = ifelse(disease_code != 2 | is.na(disease_code), F, T),
                     size_sex = "All") %>%
              group_by(year, Station, lat, lon, area_swept, size_sex, bcs) %>%
              summarise(num_crab = round(sum(sample_factor,na.rm = T)))) %>%
  filter(!is.na(size_sex)) %>%
  ungroup() %>%
  mutate(cpue = num_crab / area_swept) -> catch

#abundance by size sex/BCS categories
catch %>%
  right_join(strata_ebs, by = c("year", "Station")) %>%
  group_by(year, stratum, size_sex, bcs) %>%
  summarise(total_area = mean(total_area),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(year, size_sex, bcs) %>%
  summarise(Total_abun=sum(abundance_mil)) %>% #sum across strata
  ungroup() -> abund
  
#Calculate BCS prevalence 
abund %>%
  group_by(year,size_sex) %>%
  summarise(Prevalance = Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+(Total_abun[bcs==TRUE]))) -> bcs

###################################################
 
#Combined Plot 
#NOTE: 2021 has zero prev for 2 size/sex categories
  bcs %>%
    ggplot(aes(x = year, y = Prevalance, group = as.factor(size_sex))) +
    geom_point(aes(colour = size_sex), size=3) +
    geom_line(aes(colour = size_sex), size=1) +
    labs(y = "Disease Prevalence", x = "") +
    theme_bw() +
    theme(legend.text=element_text(size=11)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
    theme(legend.title= element_blank()) +
    scale_x_continuous(breaks = seq(1990, 2021, 5))

#Faceted plot 
bcs %>%
  ggplot(aes(x = year, y = Prevalance)) +
  geom_point(aes(colour = size_sex), size=3) +
  geom_line(aes(colour = size_sex), size=1) +
  labs(y = "Disease Prevalence", x = "") +
  theme_bw() +
  theme(legend.text=element_text(size=11)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=12)) +
  theme(legend.title= element_blank()) +
  scale_x_continuous(breaks = seq(1990, 2021, 5)) +
  facet_wrap(~size_sex)

#Just population
bcs %>%
  filter(size_sex=="All") %>%
  ggplot(aes(x = year, y = Prevalance)) +
  geom_point(aes(colour = size_sex), size=3) +
  geom_line(aes(colour = size_sex), size=1) +
  labs(y = "Disease Prevalence", x = "") +
  theme_bw() +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
  scale_x_continuous(breaks = seq(1990, 2021, 5))

 #Write csv for snow crab indicator 
  write.csv(bcs, file="./Output/bcs_prev.csv")
  

 #################################
  
#Plot Pam's prevalence data by year
  pcr %>%
    ggplot(aes(fill=Index_site, y=Perc_inf, x=Year)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_bw() +
    labs(y = "Disease Prevalence", x = "") +
    #theme(panel.grid = element_blank()) +
    scale_fill_brewer(palette="Reds") +
    theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
    theme(axis.title.y=element_text(size=14)) +
   labs(fill="Index Site") +
    theme(legend.text=element_text(size=10)) +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5))
  
  pcr %>%
    filter(Index_site=="Index6") %>%
    ggplot(aes(y=Perc_inf, x=Year)) + 
    geom_bar(position="dodge", stat="identity", fill="blue" ) +
    theme(legend.position = "none") +
    theme_bw() +
    labs(y = "Disease Prevalence", x = "") +
    #theme(panel.grid = element_blank()) +
    theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) +
    theme(axis.title.y=element_text(size=14)) +
    labs(fill="Index Site") +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5))
  
   