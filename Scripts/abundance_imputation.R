# notes ----

# Calculate immature snow crab abundance in EBS

#NOTE: There are no pre-1980 stratum tables on akfin so abundances were only 
  #calculated for 1980+

# Author: Erin Fedewa

# load ----
library(tidyverse)
library(ggmap)
library(Hmisc)

theme_set(theme_bw())


##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

###################################################

# evaluate available hauls
haul <- read.csv("./Data/haul.csv")

# restrict to haul_type = 3 and plot
haul <- haul %>%
  filter(HAUL_TYPE == 3) 

ggplot(haul, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# load core stations from Erin's analysis
core <- read.csv("./output/core_stations.csv")

# combine with haul
core <- core %>%
  select(GIS_STATION) %>%
  left_join(., haul)

# plot these
ggplot(core, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# plot these
ggplot(core, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# define survey year in catch
sc_catch <- sc_catch %>%
  mutate(YEAR = floor(CRUISE/100))

# rename year in strata
sc_strata <- sc_strata %>%
  rename(YEAR = SURVEY_YEAR)

## now attempt for the entire stratum---------------------

# start over with loading up data

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

# evaluate available hauls
haul <- read.csv("./Data/haul.csv")

# restrict to haul_type = 3 and plot
haul <- haul %>%
  filter(HAUL_TYPE == 3) 

ggplot(haul, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# combine stratum with haul
use <- unique(sc_strata$STATION_ID)
  
stratum <- haul %>%
  filter(GIS_STATION %in% use) %>%
  rename(YEAR = SURVEY_YEAR)

# plot these
ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

## get strata areas for weighted mean abundance --------------

area <- sc_strata %>%
  filter(SURVEY_YEAR == 2022)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(count = n(),
            area = mean(TOTAL_AREA))

## all immature for the snow crab stratum --------------------
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX %in% 1:2) %>%
  mutate(MAT_SEX = case_when((SEX == 1 & WIDTH >= 30 & WIDTH <= 95 & SHELL_CONDITION <= 2) ~ "Immature Male",
                             (SEX == 2 & WIDTH >= 30 & CLUTCH_SIZE == 0) ~ "Immature Female")) %>%
  filter(MAT_SEX %in% c("Immature Male", "Immature Female")) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 


# join with stratum stations to include 0 catches 
stratum <- haul %>%
  filter(GIS_STATION %in% use) %>%
  rename(YEAR = SURVEY_YEAR)

stratum <- stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(stratum$CPUE)

sum(change)

stratum$CPUE[change] <- 0

# check samples per year
check <- stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
stratum$log_cpue <- log(stratum$CPUE + 1)

# plot histogram of mean cpue per station
check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# plot zeros

zero <- check %>%
  filter(mean_log_cpue == 0)

ggplot(zero, aes(longitude, latitude)) +
  geom_point() +
  xlim(-180, -155) +
  ylim(54, 62)

# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# View(drop_0)


# # drop catches below 1
# 
# drop_1 <- check %>%
#   filter(mean_log_cpue >= 1)
# 
# ggplot(drop_1, aes(mean_log_cpue)) + 
#   geom_histogram(bins = 30, fill = "grey", color = "black")
# 
# ggplot(drop_1, aes(longitude, latitude, color = mean_log_cpue)) +
#   geom_point(size = 4, shape = 15) +
#   scale_colour_gradient(
#     low = "purple",
#     high = "red",
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "colour"
#   )

# drop catches < 5th quantile

drop_5th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.05))

ggplot(drop_5th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

all_immature_drop5_map <- ggplot(drop_5th, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour") +
  theme(axis.title = element_blank(),
        legend.position = c(0.8,0.8))

all_immature_drop5_map

# remove mean_log_cpue < 5th quantile from annual catch data

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_5th$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

hist(r[pred])

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_5.RDS")
imp_immature_drop_5 <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_5.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_immature_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_5th$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_drop_5th <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum > 5th percentile")

immature.stratum_drop_5th


## total immature, stratum, drop 0 only ----------------------------

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_0$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_0.RDS")
imp <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_0.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 13


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_0$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_drop_0 <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum, drop 0 catch")

immature.stratum_drop_0




## total immature, stratum, drop > 2nd percentile---------------

# drop catches < 2nd quantile

drop_2nd <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.02))

ggplot(drop_2nd, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_2nd, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 2nd quantile from annual catch data

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_2nd$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_2nd.RDS")
imp_all <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_2nd.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp_all))

sum(check) # 0!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_all, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_2nd$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_drop_2nd <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum > 2nd percentile")

immature.stratum_drop_2nd

## male 30-59, whole stratum, drop 10th quantile -------------

scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         WIDTH >= 30 & WIDTH < 60,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 


# join with stratum stations to include 0 catches 
stratum <- stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(stratum$CPUE)

sum(change)

stratum$CPUE[change] <- 0

# check samples per year
check <- stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
stratum$log_cpue <- log(stratum$CPUE + 1)

# plot histogram of mean cpue per station

check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# plot zeros

zero <- check %>%
  filter(mean_log_cpue == 0)

ggplot(zero, aes(longitude, latitude)) +
  geom_point() +
  xlim(-180, -155) +
  ylim(54, 62)

# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# View(drop_0)

# drop catches < 10th quantile

drop_10th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.1))

ggplot(drop_10th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_10th, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 10th quantile from annual catch data

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_10th$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_male_30-59_stratum_drop_10th.RDS")
imp_30_59_drop_10 <- readRDS("./output/abundance_imputations_male_30-59_stratum_drop_10th.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp_30_59_drop_10))

sum(check) # 4!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_30_59_drop_10, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_10th$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

male_30_59_stratum_drop_10th <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Male 30-59 mm, stratum > 10th percentile")

male_30_59_stratum_drop_10th


## male 60-95, whole stratum, drop 5th quantile --------------------
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         WIDTH >= 60 & WIDTH <= 95,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 


# join with stratum stations to include 0 catches 
stratum <- stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(stratum$CPUE)

sum(change)

stratum$CPUE[change] <- 0

# check samples per year
check <- stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
stratum$log_cpue <- log(stratum$CPUE + 1)

# plot histogram of mean cpue per station

check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# plot zeros

zero <- check %>%
  filter(mean_log_cpue == 0)

ggplot(zero, aes(longitude, latitude)) +
  geom_point() +
  xlim(-180, -155) +
  ylim(54, 62)

# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# drop catches < 5th quantile

drop_5th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.05))

ggplot(drop_5th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_5th, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 2nd quantile from annual catch data

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_5th$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_male_60_95_stratum_drop_5.RDS")
imp_60_95_drop_5 <- readRDS("./output/abundance_imputations_male_60_95_stratum_drop_5.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_60_95_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_5th$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

male_60_95_stratum_drop_5th <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Male 60-95, stratum > 5th percentile")

male_60_95_stratum_drop_5th






## female immature > 30 mm ---------------------------
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 2,
         CLUTCH_SIZE == 0,
         WIDTH >= 30,
         SHELL_CONDITION == 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 

# join with stratum stations to include 0 catches 
stratum <- haul %>%
  filter(GIS_STATION %in% use) %>%
  rename(YEAR = SURVEY_YEAR)

stratum <- stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(stratum$CPUE)

sum(change)

stratum$CPUE[change] <- 0

# check samples per year
check <- stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
stratum$log_cpue <- log(stratum$CPUE + 1)

# plot histogram of mean cpue per station

check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# plot zeros

zero <- check %>%
  filter(mean_log_cpue == 0)

ggplot(zero, aes(longitude, latitude)) +
  geom_point() +
  xlim(-180, -155) +
  ylim(54, 62)

# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# drop 10th quantile
drop_10th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.1))

ggplot(drop_10th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_10th, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 10th quantile from annual catch data

dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_10th$GIS_STATION))


# remove mean_log_cpue = 0 from annual catch data

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

hist(r[pred])

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_female_immature_stratum_drop_10.RDS")
imp_female_drop_10 <- readRDS("./output/abundance_imputations_female_immature_stratum_drop_10.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp_female_drop_10))

sum(check) # 4!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_female_drop_10, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_10th$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

female_stratum_drop_10th <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Female > 30 mm, stratum > 10th percentile")

female_stratum_drop_10th


###################

## look at the influence of corner stations

# start over with loading up data

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

# evaluate available hauls
haul <- read.csv("./Data/haul.csv")

# restrict to haul_type = 3 and plot
haul <- haul %>%
  filter(HAUL_TYPE == 3) 

ggplot(haul, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# get station names in stratum
use <- unique(sc_strata$STATION_ID)

# corner station names are 6 characters long, regular are 4 characters long
stations <- data.frame(station = use,
                       length = str_length(use))

corner <- stations %>%
  filter(length > 5)

regular <- stations %>%
  filter(length < 5)

regular_stratum <- haul %>%
  filter(GIS_STATION %in% regular$station) %>%
  rename(YEAR = SURVEY_YEAR)

corner_stratum <- haul %>%
  filter(GIS_STATION %in% corner$station) %>%
  rename(YEAR = SURVEY_YEAR)

# plot these
ggplot(regular_stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

ggplot(corner_stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

## impute all immature, all non-0 mean cpue stations, regular stratum (no corners)--------------

scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX %in% 1:2) %>%
  mutate(MAT_SEX = case_when((SEX == 1 & WIDTH >= 30 & WIDTH <= 95 & SHELL_CONDITION <= 2) ~ "Immature Male",
                             (SEX == 2 & WIDTH >= 30 & CLUTCH_SIZE == 0) ~ "Immature Female")) %>%
  filter(MAT_SEX %in% c("Immature Male", "Immature Female")) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 


# join with stratum stations to include 0 catches 
regular_stratum <- regular_stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(regular_stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(regular_stratum$CPUE)

sum(change)

regular_stratum$CPUE[change] <- 0

# check samples per year
check <- regular_stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
regular_stratum$log_cpue <- log(regular_stratum$CPUE + 1)

# plot histogram of mean cpue per station

check <- regular_stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)


ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# View(drop_0)

dat <- regular_stratum

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_corners.RDS")
imp <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_corners.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- regular_stratum %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_no_corner <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum, no corner")

immature.stratum_no_corner

## total immature, drop corners, > 2nd quantile--------------

# drop catches < 2nd quantile

drop_2nd <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.02))

ggplot(drop_2nd, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_2nd, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 2nd quantile from annual catch data

dat <- regular_stratum %>%
  filter(regular_stratum$GIS_STATION %in% unique(drop_2nd$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_corners_drop_2nd.RDS")
imp <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_corners_drop_2nd.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_2nd$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_drop_2nd <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum > 2nd percentile")

immature.stratum_drop_2nd

## male 30-60 mm, stratum, no corners, > 2nd quantile--------------
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         WIDTH >= 30 & WIDTH < 60,
         SHELL_CONDITION <= 2) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 

regular_stratum <- haul %>%
  filter(GIS_STATION %in% regular$station) %>%
  rename(YEAR = SURVEY_YEAR)

# join with stratum stations to include 0 catches 
regular_stratum <- regular_stratum %>%
  select(GIS_STATION, YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  left_join(., scratch)

ggplot(regular_stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(regular_stratum$CPUE)

sum(change)

regular_stratum$CPUE[change] <- 0

# check samples per year
check <- regular_stratum %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
regular_stratum$log_cpue <- log(regular_stratum$CPUE + 1)

# plot histogram of mean cpue per station

check <- regular_stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(check, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )

# remove stations with mean cpue = 0

drop_0 <- check %>%
  filter(mean_log_cpue > 0)


ggplot(drop_0, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# looks at histogram for non-zero catches
ggplot(drop_0, aes(mean_log_cpue)) +
  geom_histogram(bins = 50, fill = "grey", color = "black")

# View(drop_0)

dat <- regular_stratum

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# drop catches < 5th quantile

drop_5th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.05))

ggplot(drop_5th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(drop_5th, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 4, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )


# remove mean_log_cpue < 2nd quantile from annual catch data

dat <- regular_stratum %>%
  filter(regular_stratum$GIS_STATION %in% unique(drop_5th$GIS_STATION))

# check stations per year
count <- dat %>%
  group_by(YEAR) %>%
  summarise(count = n())

count


# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations_male_30-60__stratum_drop_corners_drop_5th.RDS")
imp <- readRDS("./output/abundance_imputations_male_30-60_stratum_drop_corners_drop_5th.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = rowMeans(temp, na.rm = T)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# summarize raw (non-imputed data) to plot

raw.dat <- stratum %>%
  filter(stratum$GIS_STATION %in% unique(drop_2nd$GIS_STATION)) %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

immature.stratum_drop_2nd <- ggplot(plot.dat, aes(year, log_mean)) +
  geom_line(data = raw.dat, aes(year, log_mean), color = "red") +
  geom_point(data = raw.dat, aes(year, log_mean), color = "red") + 
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum > 2nd percentile")

immature.stratum_drop_2nd


%>%
  #join to zero catch stations
  right_join(sc_strata %>%
               distinct(YEAR, STATION_ID, AREA_SWEPT) %>%
               rename_all(~c("GIS_STATION", "YEAR",
                            "TOTAL_AREA_SQ_NM"))) %>%
  replace_na(list(CPUE = 0)) 

# check stratum size at this point
check <- scratch %>%
  group_by(YEAR, STRATUM) %>%
  summarise(area = mean(TOTAL_AREA_SQ_NM)) %>% 
  group_by(YEAR) %>%
  summarise(area = sum(area))

ggplot(check, aes(YEAR, area)) + 
  geom_line() +
  geom_point() # same result!

# plot 1980 - 1992 samples
ggplot(filter(scratch, YEAR < 1993), aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)


%>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA_SQ_NM) %>%
  summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA_SQ_NM))) %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE = sum(ABUNDANCE)) -> imm_abundance


# calculate CPUE by station for immature snow crab 
sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         YEAR > 1979,
         WIDTH >= 40 & WIDTH <= 95) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT,MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #join to zero catch stations
  right_join(sc_strata %>%
               filter(YEAR > 1979) %>%
               distinct(YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
               rename_all(~c("GIS_STATION", "YEAR",
                             "STRATUM", "TOTAL_AREA_SQ_NM"))) %>%
  replace_na(list(CPUE = 0)) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA_SQ_NM) %>%
  summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA_SQ_NM))) %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE = sum(ABUNDANCE)) -> imm_abundance


ggplot(imm_abundance, aes(YEAR, log(ABUNDANCE))) +
  geom_point() +
  geom_line()

#Write csv as output 
write.csv(imm_abundance, file = "./Data/imm_abun.csv")
