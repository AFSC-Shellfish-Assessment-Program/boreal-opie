# male and female abundance for Fig. 1
# impute for missing stations
# using the same size / maturity designations used in borealization effects model

# load ----

library(tidyverse)
library(Hmisc)
library(patchwork)
library(moments)

theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##############################################

## data
#EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

#Standard 375 survey stations
sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  pull(STATION_ID) ->standard

## male imputation -------------------

# process catch data
# males 30-95 mm

sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4),
         SEX == 1,
         WIDTH >= 30 & WIDTH <= 95,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #add zero catch stations
  right_join(sc_catch %>% 
               mutate(YEAR = floor(CRUISE/100)) %>%
               filter(HAUL_TYPE %in% c(3,4)) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)) %>%
  replace_na(list(CPUE = 0)) %>%
  mutate(log_cpue = log(CPUE +1)) %>%
  #filter for standard stations (1979 has 70 extra non-standard stations)
  filter(GIS_STATION %in% standard) -> cpue



# Remove stations that have never caught snow crab 
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(station_mean = mean(CPUE)) %>% 
  filter(station_mean > 0) %>% #58 stations dropped
  #b) Remove stations < 5th CPUE percentile range 
  filter(station_mean > quantile(station_mean, 0.05)) %>% #16 additional stations dropped
  pull(GIS_STATION) -> keepers  

#Final dataset for imputation
cpue %>%
  filter(GIS_STATION %in% keepers) -> dat

## create bounds for imputed values for use below
# (min = 0, max = max observed)
lower_bound <- 0
upper_bound <- max(dat$log_cpue, na.rm = T) 

# get data set of excluded stations for later averaging when computing total abundance
cpue %>%
  filter(!GIS_STATION %in% keepers) -> excluded_stations

excluded_stations <- unique(excluded_stations$GIS_STATION)

## plot the data to be used for imputation

# plot means by station
plot_dat <- dat %>%
  group_by(YEAR) %>% 
  mutate(mean_log_cpue = mean(log_cpue))


ggplot(plot_dat, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")


# plot means on map
plot_map <- dat %>%
  group_by(GIS_STATION) %>% 
  mutate(mean_log_cpue = mean(log_cpue),
         latitude = mean(MID_LATITUDE),
         longitude = mean(MID_LONGITUDE))

male_abundance_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 3, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  theme(legend.position = c(0.8, 0.8))

male_abundance_map 

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 25 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

male_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

male_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "pmm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_male_30_95_stratum_drop_5.RDS")
imp_30_95_drop_5 <- readRDS("./output/abundance_imputations_male_30_95_stratum_drop_5.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_30_95_drop_5))

sum(check) # 0

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled

area <- sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  #STATION_ID %in% unique(drop_10th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA))

weights <- data.frame(STATION = colnames(complete(imp_30_95_drop_5, action = 1))) %>%
  left_join(., area)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights %>%
  group_by(DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights <- left_join(weights, summary_weights)

# set up weighted mean function
ff <- function(x) weighted.mean(x, weights$weight, na.rm = T)

# data frame of 0 cpue for excluded stations
excluded <- data.frame(STATION = str_remove_all(excluded_stations, "-"),
                       cpue = 0)


imputed_male_abundance <- imputed_male_weighted_cpue <- data.frame()

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_30_95_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
                              
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  # weighted cpue
  imputed_male_weighted_cpue <- rbind(imputed_male_weighted_cpue,
                                        data.frame(imputation = i,
                                                   year = c(1975:2019, 2021, 2022),
                                                   mean = apply(back.temp, 1, ff)))


  # loop through each year and calculate total estimated abundance
  # (sum of mean cpue*stratum areas)
  for(j in 1:nrow(back.temp)){
    # j <- 2
    
    temp.temp <- data.frame(STATION = colnames(back.temp),
                            cpue = as.vector(t(back.temp[j,])))
    
    # add in excluded stations as cpue = 0
    temp.temp <- rbind(temp.temp,
                       excluded)
    
    # join with area for each stratum
    temp.temp <- left_join(temp.temp, area)
    
    
    temp.sum <- temp.temp %>%
      group_by(DISTRICT) %>%
      summarise(mean_cpue = mean(cpue),
                total_area = mean(TOTAL_AREA),
                abundance = mean_cpue*total_area)
    
    imputed_male_abundance <- rbind(imputed_male_abundance,
                               data_frame(imputation = i,
                                          year = rownames(back.temp)[j],
                                          abundance = sum(temp.sum$abundance)))
    
  } # close j loop (years/rows)
  
print(paste("Imputation # ", i, sep = ""))
  
} # close i loop (imputation)

male_weighted_log_cpue <- data.frame(year = c(1975:2019, 2021:2022),
                                       log_mean = tapply(log(imputed_male_weighted_cpue$mean), imputed_male_weighted_cpue$year, mean, na.rm = T),
                                       sd = tapply(log(imputed_male_weighted_cpue$mean), imputed_male_weighted_cpue$year, sd, na.rm = T))

# add in NAs
xtra <- data.frame(year = 2020,
                   log_mean = NA,
                   sd = NA)

male_weighted_log_cpue <- rbind(male_weighted_log_cpue,
                         xtra) %>%
  arrange(year)  %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

# save for borealization effects analysis
write.csv(male_weighted_log_cpue, "./output/male_imputed_weighted_log_cpue.csv", row.names = F)

# summarize and plot
imputation_summary <- imputed_male_abundance %>%
  group_by(year) %>%
  summarise(abundance = mean(abundance)) %>%
  mutate(SD = tapply(imputed_male_abundance$abundance, imputed_male_abundance$year, sd))
          
# add in NAs for 2020

xtra <- data.frame(year = 2020,
                   abundance = NA,
                   SD = NA)

imputation_summary <- rbind(imputation_summary, xtra)

# and divide by 1e9 (billions of animals!)
male_imputation_plot <- imputation_summary %>%
  mutate(abundance = abundance / 1e9,
         SD = SD / 1e9)


ggplot(male_imputation_plot, aes(as.numeric(year), abundance)) +
  geom_col(fill = cb[6]) + 
  geom_errorbar(aes(ymin = abundance - 2*SD,
                    ymax = abundance + 2*SD))


# plot 1-year and 2-year differences
imputation_summary <- imputation_summary %>%
  arrange(year)

male_diff <- data.frame(year = 1976:2022,
                        diff1 = diff(imputation_summary$abundance/1e9, lag = 1, na.pad = T),
                        diff2 = c(NA,diff(imputation_summary$abundance/1e9, lag= 2, na.pad = T))) %>%
  pivot_longer(cols = -year)

ggplot(male_diff, aes(year, value, color = name)) +
  geom_point() + 
  geom_line()

## immature females ------------

sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4),
         SEX == 2,
         CLUTCH_SIZE == 0,
         WIDTH >= 30,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #add zero catch stations
  right_join(sc_catch %>% 
               mutate(YEAR = floor(CRUISE/100)) %>%
               filter(HAUL_TYPE %in% c(3,4)) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)) %>%
  replace_na(list(CPUE = 0)) %>%
  mutate(log_cpue = log(CPUE +1)) %>%
  #filter for standard stations (1979 has 70 extra non-standard stations)
  filter(GIS_STATION %in% standard) -> cpue_female

## Cleaning up data for imputation

# Remove stations that have never caught immature female snow crab 
cpue_female %>%
  group_by(GIS_STATION) %>%
  summarise(station_mean = mean(CPUE)) %>% 
  filter(station_mean > 0) %>% #58 stations dropped
  #b) Remove stations < 5th CPUE percentile range 
  filter(station_mean > quantile(station_mean, 0.05)) %>% #16 additional stations dropped
  pull(GIS_STATION) -> keepers_female

#Final dataset for imputation
cpue_female %>%
  filter(GIS_STATION %in% keepers_female) -> dat_female

## create bounds for imputed values for use below
# (min = 0, max = max observed)
lower_bound <- 0
upper_bound <- max(dat_female$log_cpue, na.rm = T) 

# get data set of excluded stations for later averaging when computing total abundance
cpue %>%
  filter(!GIS_STATION %in% keepers_female) -> excluded_stations

excluded_stations <- unique(excluded_stations$GIS_STATION)

# plot means by station
plot_dat <- dat_female %>%
  group_by(YEAR) %>% 
  mutate(mean_log_cpue = mean(log_cpue))

ggplot(plot_dat, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

# plot means on map
plot_map <- dat_female %>%
  group_by(GIS_STATION) %>% 
  mutate(mean_log_cpue = mean(log_cpue),
         latitude = mean(MID_LATITUDE),
         longitude = mean(MID_LONGITUDE))

female_drop_5_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 3, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  theme(legend.position = c(0.8, 0.8))

female_drop_5_map 


# get into matrix form for mice
dat_female <- tapply(dat_female$log_cpue, list(dat_female$YEAR, dat_female$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat_female))$r 
r #Cross-year correlations between each station combination

# choose 25 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 25 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

female_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

female_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat_female) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat_female)

imp <- mice::mice(data = dat_female, method = "pmm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_female_stratum_drop_5.RDS")
imp_female_drop_5 <- readRDS("./output/abundance_imputations_female_stratum_drop_5.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_female_drop_5))

sum(check) # 6!

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled

area <- sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  #STATION_ID %in% unique(drop_10th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA))

weights <- data.frame(STATION = colnames(complete(imp_female_drop_5, action = 1))) %>%
  left_join(., area)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights %>%
  group_by(DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights <- left_join(weights, summary_weights)


# set up weighted mean function
ff <- function(x) weighted.mean(x, weights$weight, na.rm = T)

# data frame of 0 cpue for excluded stations
excluded <- data.frame(STATION = str_remove_all(excluded_stations, "-"),
                       cpue = 0)

imputed_female_abundance <- imputed_female_cpue <- imputed_female_weighted_cpue <- data.frame()

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_female_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  # unweighted cpue
  imputed_female_cpue <- rbind(imputed_female_cpue,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             mean = apply(back.temp, 1, mean, na.rm = T)))
  
  # weighted cpue
  imputed_female_weighted_cpue <- rbind(imputed_female_weighted_cpue,
                               data.frame(imputation = i,
                                          year = c(1975:2019, 2021, 2022),
                                          mean = apply(back.temp, 1, ff)))  
  # loop through each year and calculate total estimated abundance
  # (sum of mean cpue*stratum areas)
  for(j in 1:nrow(back.temp)){
    # j <- 2
    
    temp.temp <- data.frame(STATION = colnames(back.temp),
                            cpue = as.vector(t(back.temp[j,])))
    
    # add in excluded stations as cpue = 0
    temp.temp <- rbind(temp.temp,
                       excluded)
    
    # join with area for each stratum
    temp.temp <- left_join(temp.temp, area)
    
    
    temp.sum <- temp.temp %>%
      group_by(DISTRICT) %>%
      summarise(mean_cpue = mean(cpue, na.rm = T),
                total_area = mean(TOTAL_AREA),
                abundance = mean_cpue*total_area)
    
    imputed_female_abundance <- rbind(imputed_female_abundance,
                                    data_frame(imputation = i,
                                               year = rownames(back.temp)[j],
                                               abundance = sum(temp.sum$abundance)))
    
  } # close j loop (years/rows)
  
  print(paste("Imputation # ", i, sep = ""))
  
} # close i loop (imputation)


# summarize cpue
female_log_cpue <- data.frame(year = c(1975:2019, 2021:2022),
                              log_mean = tapply(log(imputed_female_cpue$mean), imputed_female_cpue$year, mean, na.rm = T),
                              sd = tapply(log(imputed_female_cpue$mean), imputed_female_cpue$year, sd, na.rm = T))

female_weighted_log_cpue <- data.frame(year = c(1975:2019, 2021:2022),
                              log_mean = tapply(log(imputed_female_weighted_cpue$mean), imputed_female_weighted_cpue$year, mean, na.rm = T),
                              sd = tapply(log(imputed_female_weighted_cpue$mean), imputed_female_weighted_cpue$year, sd, na.rm = T))

# add in NAs
xtra <- data.frame(year = 2020,
                   log_mean = NA,
                   sd = NA)

female_log_cpue <- rbind(female_log_cpue,
                         xtra) %>%
  arrange(year) %>%
  mutate(summary = "new")

female_weighted_log_cpue <- rbind(female_weighted_log_cpue,
                         xtra) %>%
  arrange(year) %>%
  mutate(summary = "new_weighted")

# check
abundance_female <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year) %>%
  mutate( #log_mean_lag1 = lag(log_mean, 1),
         summary = "old")

check <- rbind(female_log_cpue, female_weighted_log_cpue, abundance_female) 

ggplot(check, aes(year, log_mean, color = summary)) +
  geom_point() + 
  geom_line()

cor(abundance_female$log_mean, female_weighted_log_cpue$log_mean, use = "p") 
# weighted mean = "simple" mean used in borealization effects analysis!

# summarize abundance and plot
imputation_summary <- imputed_female_abundance %>%
  group_by(year) %>%
  summarise(abundance = mean(abundance)) %>%
  mutate(SD = tapply(imputed_female_abundance$abundance, imputed_female_abundance$year, sd))

# add in NAs for 2020

xtra <- data.frame(year = 2020,
                   abundance = NA,
                   SD = NA)

imputation_summary <- rbind(imputation_summary, xtra)

# and divide by 1e9 (billions of animals!)
female_imputation_plot <- imputation_summary %>%
  mutate(abundance = abundance / 1e9,
         SD = SD / 1e9)


ggplot(female_imputation_plot, aes(as.numeric(year), abundance)) +
  geom_col(fill = cb[6]) + 
  geom_errorbar(aes(ymin = abundance - 2*SD,
                    ymax = abundance + 2*SD))


# plot 1-year and 2-year differences
imputation_summary <- imputation_summary %>%
  arrange(year)

female_diff <- data.frame(year = 1976:2022,
                        diff1 = diff(imputation_summary$abundance/1e9, lag = 1, na.pad = T),
                        diff2 = c(NA,diff(imputation_summary$abundance/1e9, lag= 2, na.pad = T))) %>%
  pivot_longer(cols = -year)

ggplot(female_diff, aes(year, value, color = name)) +
  geom_point() + 
  geom_line()

# combined plot
female_imputation_plot <- female_imputation_plot %>%
  mutate(sex = "Female")

male_imputation_plot <- male_imputation_plot %>%
  mutate(sex = "Male")

all_imputation_plot <- rbind(female_imputation_plot, male_imputation_plot)

pos_dodge = position_dodge(width = 0.3)

ggplot(all_imputation_plot, aes(as.numeric(year), abundance, color = sex)) +
  geom_point(position = pos_dodge) +
  geom_line(position = pos_dodge) + 
  geom_errorbar(aes(ymin = abundance - 2*SD,
                    ymax = abundance + 2*SD), width = 0.5, position = pos_dodge) +
  scale_color_manual(values = cb[c(2,6)])


# save data for combined Fig 1
write.csv(all_imputation_plot, "./output/imputed_male_30-95_imm_female_abundance.csv", row.names = F)

# first, back-transformed data
plot.back.trans <- data.frame(year = c(1975:2019, 2021:2022),
                              mean = tapply(log(imputed.back.trans.dat$mean), imputed.back.trans.dat$year, mean),
                              sd = tapply(log(imputed.back.trans.dat$mean), imputed.back.trans.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.back.trans <- rbind(plot.back.trans,
                         xtra)


# replace sd = 0 with NA
check <- plot.back.trans$sd == 0
plot.back.trans$sd[check] <- NA

# save
write.csv(plot.back.trans, "imputed_total_abundance.csv", row.names = F)

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))