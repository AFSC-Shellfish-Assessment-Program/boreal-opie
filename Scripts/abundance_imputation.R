# Calculate immature snow crab abundance in EBS
# Impute for missing stations

# load ----
library(tidyverse)
library(Hmisc)
library(patchwork)
library(moments)

theme_set(theme_bw())


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


#################################
#################################
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

#Calculate mean CPUE/yr across all stations
cpue %>%
  group_by(YEAR) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            mean_cpue =mean(CPUE),
            log_mean_cpue = log(mean_cpue)) -> mean_cpue

#plot
mean_cpue %>%
  ggplot(aes(YEAR, mean_cpue)) +
  geom_point() +
  geom_line() -> cpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, log_mean_cpue)) +
  geom_point() +
  geom_line() -> logcpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, mean_log_cpue)) +
  geom_point() +
  geom_line() -> meanlogcpue_plot

cpue_plot + logcpue_plot + meanlogcpue_plot

## Cleaning up data for imputation

# Remove stations that have never caught 60-95mm male snow crab 
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

male_30_95_drop_5_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
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

male_30_95_drop_5_map 

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

male_30_95_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

male_30_95_drop5_hist

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

imputed.dat <- imputed.back.trans.dat <- concentration.dat <-  data.frame()

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_30_95_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = apply(temp, 1, ff)))
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  imputed.back.trans.dat <- rbind(imputed.back.trans.dat,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             mean = apply(back.temp, 1, ff)))
  
  # concentration 
  conc <- NA
  
  for(j in 1:nrow(back.temp)){
    # j <- 1
    
    temp.temp <- as.vector(t(back.temp[j,]))
    
    temp.temp <- temp.temp[order(-temp.temp)]
    
    conc[j]  <- sum(temp.temp[1:25])/sum(temp.temp)
    
  }
  
  concentration.dat <- rbind(concentration.dat,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             concentration = conc))
  
  }

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

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# mean log-transformed data
plot.mean.log <- data.frame(year = c(1975:2019, 2021:2022),
                              mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                              sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.mean.log <- rbind(plot.mean.log,
                         xtra)


# replace sd = 0 with NA
check <- plot.mean.log$sd == 0
plot.mean.log$sd[check] <- NA

# plot
ggplot(plot.mean.log, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# plot concentration
plot.conc <- data.frame(year = c(1975:2019, 2021:2022),
                       mean = tapply(concentration.dat$concentration, concentration.dat$year, mean),
                       sd = tapply(concentration.dat$concentration, concentration.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.conc <- rbind(plot.conc,
                  xtra)


# replace sd = 0 with NA
check <- plot.conc$sd == 0
plot.conc$sd[check] <- NA

# plot
ggplot(plot.conc, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))


# combine different imputed time series and save

plot.back.trans <- plot.back.trans %>%
  rename(log_backtransformed_mean_cpue = mean,
         log_backtransformed_sd_cpue = sd)

plot.mean.log <- plot.mean.log %>%
  rename(mean_log_cpue = mean,
         sd_log_cpue = sd)

plot.conc <- plot.conc %>%
  rename(mean_concentratin = mean,
         sd_concentration = sd)

imputed_output <- left_join(plot.back.trans, plot.mean.log) %>%
  left_join(., plot.conc)


write.csv(imputed_output, "./output/male_30-95_imputed_data.csv")

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

#Calculate mean CPUE/yr across all stations
cpue_female %>%
  group_by(YEAR) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            mean_cpue =mean(CPUE),
            log_mean_cpue = log(mean_cpue)) -> mean_cpue_female

#plot
mean_cpue_female %>%
  ggplot(aes(YEAR, mean_cpue)) +
  geom_point() +
  geom_line() -> cpue_plot

mean_cpue_female %>%
  ggplot(aes(YEAR, log_mean_cpue)) +
  geom_point() +
  geom_line() -> logcpue_plot

mean_cpue_female %>%
  ggplot(aes(YEAR, mean_log_cpue)) +
  geom_point() +
  geom_line() -> meanlogcpue_plot

cpue_plot + logcpue_plot + meanlogcpue_plot

## Cleaning up data for imputation

# a) Remove stations that have never caught 60-95mm male snow crab 
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

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
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


# summarize imputation data


# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat_female, na.rm = T) 

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

imputed.dat <- imputed.back.trans.dat <- concentration.dat <-  data.frame()

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_female_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = apply(temp, 1, ff)))
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  imputed.back.trans.dat <- rbind(imputed.back.trans.dat,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             mean = apply(back.temp, 1, ff)))
  # concentration 
  conc <- NA
  
  for(j in 1:nrow(back.temp)){
    # j <- 1
    
    temp.temp <- as.vector(t(back.temp[j,]))
    
    temp.temp <- temp.temp[order(-temp.temp)]
    
    conc[j]  <- sum(temp.temp[1:25], na.rm = T)/sum(temp.temp, na.rm = T)
    
  }
  
  concentration.dat <- rbind(concentration.dat,
                             data.frame(imputation = i,
                                        year = c(1975:2019, 2021, 2022),
                                        concentration = conc))
  
  
}



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

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# mean log-transformed data
plot.mean.log <- data.frame(year = c(1975:2019, 2021:2022),
                            mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                            sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.mean.log <- rbind(plot.mean.log,
                       xtra)


# replace sd = 0 with NA
check <- plot.mean.log$sd == 0
plot.mean.log$sd[check] <- NA

# plot
ggplot(plot.mean.log, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# plot concentration
plot.conc <- data.frame(year = c(1975:2019, 2021:2022),
                        mean = tapply(concentration.dat$concentration, concentration.dat$year, mean),
                        sd = tapply(concentration.dat$concentration, concentration.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.conc <- rbind(plot.conc,
                   xtra)


# replace sd = 0 with NA
check <- plot.conc$sd == 0
plot.conc$sd[check] <- NA

# plot
ggplot(plot.conc, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))


# combine different imputed time series and save

plot.back.trans <- plot.back.trans %>%
  rename(log_backtransformed_mean_cpue = mean,
         log_backtransformed_sd_cpue = sd)

plot.mean.log <- plot.mean.log %>%
  rename(mean_log_cpue = mean,
         sd_log_cpue = sd)

plot.conc <- plot.conc %>%
  rename(mean_concentratin = mean,
         sd_concentration = sd)

imputed_output <- left_join(plot.back.trans, plot.mean.log) %>%
  left_join(., plot.conc)


write.csv(imputed_output, "./output/female_imputed_data.csv")



## combine all the plots ---------------

png("./Figs/combined_imputation_plots.png", width = 15, height = 15, units = 'in', res = 300)

ggpubr::ggarrange(male_30_59_drop_10_map, male_30_59_drop10_hist, male_30_59_stratum_drop_10,
                  male_60_95_drop_5_map, male_60_95_drop5_hist, male_60_95_stratum_drop_5,
                  female_immature_drop_10_map, female_immature_drop10_hist, female.stratum_drop_10th,
                  all_immature_drop_2_map, all_immature_drop2_hist, immature.stratum_drop_2nd,
                  ncol = 3,
                  nrow = 4,
                  labels = "auto")

dev.off()

## simple plot for talk------------------

d1 <- read.csv("./output/male3059_drop10_df.csv", row.names = 1) %>%
  dplyr::select(year, imp_log_mean, imp_sd) %>%
  dplyr::mutate(group = "30-59 mm")

xtra <- data.frame(year = 2020, imp_log_mean = NA, imp_sd = NA, group = "30-59 mm")

d1 <- rbind(d1, xtra)

d2 <- read.csv("./output/male6095_drop5_df.csv", row.names = 1) %>%
  dplyr::select(year, imp_log_mean, imp_sd) %>%
  mutate(group = "60-95 mm")
xtra <- data.frame(year = 2020, imp_log_mean = NA, imp_sd = NA, group = "60-95 mm")

d2 <- rbind(d2, xtra)

plot_dat <- rbind(d1, d2)

change <- is.na(plot_dat$imp_sd) 
plot_dat$imp_sd[change] <- 0



ggplot(plot_dat, aes(year, imp_log_mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = imp_log_mean-2*imp_sd, ymax = imp_log_mean+2*imp_sd)) +
  labs(y = "Mean ln(CPUE)") +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  theme(axis.title.x = element_blank())

ggsave("./figs/imputed_abundance_talk.png", width = 4, height = 5, units= 'in')

######################

# legal male abundance
#################################
#################################
# males 30-95 mm

sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4),
         SEX == 1,
         WIDTH >= 78) %>%
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

#Calculate mean CPUE/yr across all stations
cpue %>%
  group_by(YEAR) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            mean_cpue =mean(CPUE),
            log_mean_cpue = log(mean_cpue)) -> mean_cpue

#plot
mean_cpue %>%
  ggplot(aes(YEAR, mean_cpue)) +
  geom_point() +
  geom_line() -> cpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, log_mean_cpue)) +
  geom_point() +
  geom_line() -> logcpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, mean_log_cpue)) +
  geom_point() +
  geom_line() -> meanlogcpue_plot

cpue_plot + logcpue_plot + meanlogcpue_plot

## Cleaning up data for imputation

# Remove stations that have never caught 60-95mm male snow crab 
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

male_legal_drop_5_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
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

male_legal_drop_5_map 

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

male_legal_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

male_legal_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "pmm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_male_legal_stratum_drop_5.RDS")
imp_legal_drop_5 <- readRDS("./output/abundance_imputations_male_legal_stratum_drop_5.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_legal_drop_5))

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

weights <- data.frame(STATION = colnames(complete(imp_legal_drop_5, action = 1))) %>%
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

imputed.dat <- imputed.back.trans.dat <- concentration.dat <-  data.frame()

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_legal_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = apply(temp, 1, ff)))
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  imputed.back.trans.dat <- rbind(imputed.back.trans.dat,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             mean = apply(back.temp, 1, ff)))
  
  # concentration 
  conc <- NA
  
  for(j in 1:nrow(back.temp)){
    # j <- 1
    
    temp.temp <- as.vector(t(back.temp[j,]))
    
    temp.temp <- temp.temp[order(-temp.temp)]
    
    conc[j]  <- sum(temp.temp[1:25])/sum(temp.temp)
    
  }
  
  concentration.dat <- rbind(concentration.dat,
                             data.frame(imputation = i,
                                        year = c(1975:2019, 2021, 2022),
                                        concentration = conc))
  
}

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

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# mean log-transformed data
plot.mean.log <- data.frame(year = c(1975:2019, 2021:2022),
                            mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                            sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.mean.log <- rbind(plot.mean.log,
                       xtra)


# replace sd = 0 with NA
check <- plot.mean.log$sd == 0
plot.mean.log$sd[check] <- NA

# plot
ggplot(plot.mean.log, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# plot concentration
plot.conc <- data.frame(year = c(1975:2019, 2021:2022),
                        mean = tapply(concentration.dat$concentration, concentration.dat$year, mean),
                        sd = tapply(concentration.dat$concentration, concentration.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.conc <- rbind(plot.conc,
                   xtra)


# replace sd = 0 with NA
check <- plot.conc$sd == 0
plot.conc$sd[check] <- NA

# plot
ggplot(plot.conc, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))


# combine different imputed time series and save

plot.back.trans <- plot.back.trans %>%
  rename(log_backtransformed_mean_cpue = mean,
         log_backtransformed_sd_cpue = sd)

plot.mean.log <- plot.mean.log %>%
  rename(mean_log_cpue = mean,
         sd_log_cpue = sd)

plot.conc <- plot.conc %>%
  rename(mean_concentratin = mean,
         sd_concentration = sd)

imputed_output <- left_join(plot.back.trans, plot.mean.log) %>%
  left_join(., plot.conc)


write.csv(imputed_output, "./output/male_legal_imputed_data.csv")

ggplot(plot.back.trans, 
       aes(year, log_backtransformed_mean_cpue)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = log_backtransformed_mean_cpue - 2*log_backtransformed_sd_cpue,
                    ymax = log_backtransformed_mean_cpue + 2*log_backtransformed_sd_cpue)) +
  theme(axis.title.x = element_blank()) +
  ylab("ln(CPUE)")

ggsave("./figs/legal_male_abundance.png", width = 3.5, height = 2.5, units = 'in')


# combine with sub-legal males

# total abundance
#################################
#################################
# 
sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4)) %>%
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

#Calculate mean CPUE/yr across all stations
cpue %>%
  group_by(YEAR) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            mean_cpue =mean(CPUE),
            log_mean_cpue = log(mean_cpue)) -> mean_cpue

#plot
mean_cpue %>%
  ggplot(aes(YEAR, mean_cpue)) +
  geom_point() +
  geom_line() -> cpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, log_mean_cpue)) +
  geom_point() +
  geom_line() -> logcpue_plot

mean_cpue %>%
  ggplot(aes(YEAR, mean_log_cpue)) +
  geom_point() +
  geom_line() -> meanlogcpue_plot

cpue_plot + logcpue_plot + meanlogcpue_plot

## Cleaning up data for imputation

# Remove stations that have never caught 60-95mm male snow crab 
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

male_sub_legal_drop_5_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
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

male_sub_legal_drop_5_map 

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

all_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

all_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "pmm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/total_abundance_imputations.RDS")
imp_total_abundance_drop_5 <- readRDS("./output/total_abundance_imputations.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_total_abundance_drop_5))

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

weights <- data.frame(STATION = colnames(complete(imp_total_abundance_drop_5, action = 1))) %>%
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

imputed.dat <- imputed.back.trans.dat <- concentration.dat <-  data.frame()

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_total_abundance_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = apply(temp, 1, ff)))
  
  # and get back-transformed means
  back.temp <- exp(temp)-1
  
  imputed.back.trans.dat <- rbind(imputed.back.trans.dat,
                                  data.frame(imputation = i,
                                             year = c(1975:2019, 2021, 2022),
                                             mean = apply(back.temp, 1, ff)))
  
  # concentration 
  conc <- NA
  
  for(j in 1:nrow(back.temp)){
    # j <- 1
    
    temp.temp <- as.vector(t(back.temp[j,]))
    
    temp.temp <- temp.temp[order(-temp.temp)]
    
    conc[j]  <- sum(temp.temp[1:25])/sum(temp.temp)
    
  }
  
  concentration.dat <- rbind(concentration.dat,
                             data.frame(imputation = i,
                                        year = c(1975:2019, 2021, 2022),
                                        concentration = conc))
  
}

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

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# mean log-transformed data
plot.mean.log <- data.frame(year = c(1975:2019, 2021:2022),
                            mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                            sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.mean.log <- rbind(plot.mean.log,
                       xtra)


# replace sd = 0 with NA
check <- plot.mean.log$sd == 0
plot.mean.log$sd[check] <- NA

# plot
ggplot(plot.mean.log, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))

# plot concentration
plot.conc <- data.frame(year = c(1975:2019, 2021:2022),
                        mean = tapply(concentration.dat$concentration, concentration.dat$year, mean),
                        sd = tapply(concentration.dat$concentration, concentration.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.conc <- rbind(plot.conc,
                   xtra)


# replace sd = 0 with NA
check <- plot.conc$sd == 0
plot.conc$sd[check] <- NA

# plot
ggplot(plot.conc, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))


# combine different imputed time series and save

plot.back.trans <- plot.back.trans %>%
  rename(log_backtransformed_mean_cpue = mean,
         log_backtransformed_sd_cpue = sd)

plot.mean.log <- plot.mean.log %>%
  rename(mean_log_cpue = mean,
         sd_log_cpue = sd)

plot.conc <- plot.conc %>%
  rename(mean_concentratin = mean,
         sd_concentration = sd)

imputed_output <- left_join(plot.back.trans, plot.mean.log) %>%
  left_join(., plot.conc)


write.csv(imputed_output, "./output/male_total_abundance_imputed_data.csv")

ggplot(plot.back.trans, 
       aes(year, log_backtransformed_mean_cpue)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = log_backtransformed_mean_cpue - 2*log_backtransformed_sd_cpue,
                    ymax = log_backtransformed_mean_cpue + 2*log_backtransformed_sd_cpue)) +
  theme(axis.title.x = element_blank()) +
  ylab("ln(CPUE)")

ggsave("./figs/total_abundance_male_abundance.png", width = 3.5, height = 2.5, units = 'in')

