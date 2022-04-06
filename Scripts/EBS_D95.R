# notes ----
# Calculate "D95" for each size sex group in EBS 1988 - 2019:
      #area of stations that make up 95% of the cumulative cpue
# Erin Fedewa
# last updated: 2020/2/5

# load ----
library(tidyverse)

#Exclude corner stations
corner <- list("QP2625","ON2625","HG2019","JI2120","IH1918",
             "GF2221","HG1918","GF2019","ON2524","PO2726",
             "IH2221","GF1918","JI2221","JI2019","JI1918",
             "HG2221","QP2726","PO2423","IH2019","PO2625",
             "QP2423","IH2120","PO2524","HG2120","GF2120",
             "QP2524")

# data ----

## EBS
catch_ebs <- read_csv("./Data/ebs_opilio_haul_data.csv")
head(catch_ebs)

# data mgmt ----

## crab data EBS
catch_ebs %>%
  filter(AKFIN_SURVEY_YEAR %in% c(1988:2019)) %>%
  select(AKFIN_SURVEY_YEAR, CRUISE, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED, NET_WIDTH, MID_LATITUDE, MID_LONGITUDE,
         HAUL_TYPE, PERFORMANCE) -> crab_ebs
names(crab_ebs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")


## compute cpue by size-sex group for each station
crab_ebs %>% 
  filter(haul_type != 17,
         #performance == 0,
         cw >= 0) %>%
  mutate(size_sex = ifelse(sex == 1 & cw >= 31 & cw <= 60.9, "male31to60",
                        ifelse(sex == 1 & cw >= 61 & cw <= 90.9, "male61to90",
                            ifelse(sex == 1 & cw >= 91 & cw <= 120.9, "male91to120",
                                  ifelse(sex == 2 & clutch >= 1, "mature_female",
                                      ifelse(sex == 2 & clutch == 0 & cw >= 31, "immature_female", NA))))),
         area_swept = distance_fished * (net_width / 1000) * 0.29155335) %>%
  group_by(year, Station, lat, lon, area_swept, size_sex) %>%
  summarise(num_crab = round(sum(sample_factor))) %>%
  filter(!is.na(area_swept)) %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  mutate(pop = sum(male31to60, male61to90, male91to120, immature_female, mature_female, na.rm = T)) %>%
  pivot_longer(c(6:12), names_to = "size_sex", values_to = "num_crab") %>%
  filter(size_sex != "NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept) %>%
  ungroup() -> cpue_long


# compute D95 by each size and sex catagory ----
# i.e. the number of stations contributing to 95% of cumulative cpue, from stations sampled in every year

# function to compute D95
f_d95_est <- function(x){
  x %>%
    arrange(-cpue) %>% #sort by cpue (large:small)
    mutate(prop_cpue = cpue/sum(cpue),  #calculate the proportion of total cpue for each station
           cum_cpue = cumsum(prop_cpue)) %>%  
    filter(cum_cpue <= 0.95) %>% #T if in d95, F if not
    count() %>%
    mutate(d95 = (n + 1) * 401) %>% #add 1 station to n to push over 95%
    pull(d95)
  }

# do the estimation
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

write.csv(d95, file="./Output/D95_output.csv")
  
