# notes ----
# Estimate missing temperature values with MICE (Multivariate Imputation by Chained Equations)
#Re-calculating avg bottom temp, cold pool extent and cold pool COD with imputed temp time series 

# authors: Mike Litzow, Erin Fedewa, Tyler Jackson
# last updated: 2020/1/9 TJ

# load ----
library(tidyverse)
library(magrittr)
library(mice)
library(Hmisc)

# define a standard error function
stderr <- function(x){sd(x) / sqrt(length(x))}

# corner station names for D95 calculations 
corner <- c("QP2625","ON2625","HG2019","JI2120","IH1918",
            "GF2221","HG1918","GF2019","ON2524","PO2726",
            "IH2221","GF1918","JI2221","JI2019","JI1918",
            "HG2221","QP2726","PO2423","IH2019","PO2625",
            "QP2423","IH2120","PO2524","HG2120","GF2120",
            "QP2524")

# raw data ----
ts <- read_csv("./Data/haul_newtimeseries.csv")
head(ts)

# data mgmt ----
ts %>%
  filter(HAUL_TYPE == 3,
         PERFORMANCE >= 0,
         SURVEY_YEAR >= 1988) %>%
  mutate(GIS_STATION = gsub("-", "", GIS_STATION)) %>%
  select(SURVEY_YEAR, GIS_STATION, GEAR_TEMPERATURE) %>%
  group_by(GIS_STATION) %>%
  mutate(TEMP_YRS = sum(!is.na(GEAR_TEMPERATURE))) %>%
  filter(TEMP_YRS >= 10) %>% #include stations with at least 10 temp observations in the 32yr timeseries 
  select(-TEMP_YRS) %>%
  pivot_wider(names_from = GIS_STATION, values_from = GEAR_TEMPERATURE) %>% 
  arrange(SURVEY_YEAR) %>%
  set_rownames(.$SURVEY_YEAR) %>%
  select(-SURVEY_YEAR) %>%
  as.matrix() -> temp_mat

# mice ----

# determine "best" predictor variables for each column
# first, get pairwise correlation coefficients
r <- rcorr(temp_mat)$r 
  r #Cross-year correlations between each station combination

# choose 20 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=20))# T for the 20 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# specify m=100 imputations, save output 
#### imp <- mice(data = temp_mat, method = "norm", m=100, predictorMatrix = pred) #Using Bayesian linear regression method
#### saveRDS(imp, "./Output/station_bottom_temp_imputations.RDS")

# read in saved imputations
imp <- readRDS("./Output/station_bottom_temp_imputations.RDS")

# generate list of completed matrices
imp_list <- list()
SURVEY_YEAR <- 1988:2019
for(i in 1:100){
  imp_list[[i]] <- cbind(SURVEY_YEAR, complete(imp, i))
  
}
head(imp_list)

# compile the completed datasets and compute a mean and standard error for estimates temps
do.call("rbind", imp_list) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "GIS_STATION", 
               values_to = "GEAR_TEMPERATURE") %>%
  nest(-SURVEY_YEAR, -GIS_STATION) %>%
  mutate(data = purrr::map(data, unlist), 
         MEAN_GT = map_dbl(data, mean),
         STD_ERR_GT = map_dbl(data, stderr)) -> temps_long

#Save per-station temp dataset 
temps_long %>%
  select(-data) ->g
write.csv(g, file="./Data/Imputed Station Temps.csv")

# take a look at standard errors, some are large...
temps_long %>%
  filter(STD_ERR_GT > 0) %>%
  pull(STD_ERR_GT) %>%
  range

# compute mean bottom temperature
temps_long %>%
  select(-data) %>%
  group_by(SURVEY_YEAR) %>%
  summarise(MEAN_BT = mean(MEAN_GT))-> avg_bt


# compute cold pool areal extent
temps_long %>%
  select(-data) %>%
  filter(MEAN_GT < 2,
         !(GIS_STATION %in% corner)) %>%
  count(SURVEY_YEAR) %>%
  mutate(CPA = n * 401) %>%
  select(-n) -> cpa

# compute cold pool COD
ts %>%
  filter(HAUL_TYPE == 3,
         PERFORMANCE >= 0,
         SURVEY_YEAR >= 1988,
         !(GIS_STATION %in% corner)) %>%
  mutate(GIS_STATION = gsub("-", "", GIS_STATION)) %>%
  select(SURVEY_YEAR, GIS_STATION, MID_LATITUDE) %>%
  right_join(temps_long, by = c("SURVEY_YEAR", "GIS_STATION")) %>%
  select(-data) %>%
  filter(MEAN_GT < 2) %>%
  group_by(SURVEY_YEAR) %>%
  summarise(CP_COD = mean(MID_LATITUDE, na.rm = T)) -> cp.cod

# combine indices and save output
avg_bt %>%
  full_join(cpa, by = "SURVEY_YEAR") %>%
  full_join(cp.cod, by = "SURVEY_YEAR") %>%
  write_csv("./Output/imputed_temperature_timeseries.csv")

#Master Timeseries Files Spatial_Env_TS.csv
  #and Spatial_Env_bySizeSex.csv updated with spatial/env indices recalculated
  #with imputed temperature timeseries 1/16/20 






















