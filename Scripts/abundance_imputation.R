# Calculate immature snow crab abundance in EBS
# Impute for missing stations


# load ----
library(tidyverse)
library(Hmisc)

theme_set(theme_bw())


##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

# evaluate available hauls
haul <- read.csv("./Data/haul.csv")

# restrict to haul_type = 3 and plot
haul <- haul %>%
  filter(HAUL_TYPE %in% 3:4) 

ggplot(haul, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# define survey year in catch
sc_catch <- sc_catch %>%
  mutate(YEAR = floor(CRUISE/100))

# rename year in strata
sc_strata <- sc_strata %>%
  rename(YEAR = SURVEY_YEAR)

# combine stratum with haul
use <- unique(sc_strata$STATION_ID)
  
stratum <- haul %>%
  filter(GIS_STATION %in% use) %>%
  rename(YEAR = SURVEY_YEAR)

# plot these
ggplot(stratum, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

## all immature for the snow crab stratum, drop <2nd quantile --------------------
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

# remove stations with mean cpue = 0
check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

# drop catches < 2nd quantile
drop_2nd <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.02))

ggplot(drop_2nd, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

all_immature_drop_2_map <-  ggplot(drop_2nd, aes(longitude, latitude, color = mean_log_cpue)) +
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

all_immature_drop_2_map

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

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

all_immature_drop2_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

all_immature_drop2_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_all_immature_stratum_drop_2nd.RDS")
imp_immature_drop_2nd <- readRDS("./output/abundance_imputations_all_immature_stratum_drop_2nd.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_immature_drop_2nd))

sum(check) # 0!

# also create df to save mean annual temp and sampling day for each imputed temperature data set
imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)
lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled
area <- sc_strata %>%
  filter(YEAR == 2022) %>%
         #STATION_ID %in% unique(drop_0$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA))

weights <- data.frame(STATION = colnames(complete(imp_immature_drop_2nd, action = 1))) %>%
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

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_immature_drop_2nd, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
  
  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  log_mean = apply(temp, 1, ff)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# replace sd = 0 with NA
check <- plot.dat$sd == 0
plot.dat$sd[check] <- NA

# summarize raw (non-imputed data) to plot
# first, get weights
raw_weights <- sc_strata %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(YEAR, STATION, DISTRICT, TOTAL_AREA)

weights1 <- distinct(data.frame(YEAR = sc_catch$YEAR, STATION = sc_catch$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION, "-")) %>%
  left_join(., raw_weights)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights1 %>%
  group_by(YEAR, DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights2 <- left_join(weights1, summary_weights)

weights2 %>%
  select(STATION, YEAR, DISTRICT, weight) %>%
  distinct() %>%
  as_tibble() %>%
  rename(GIS_STATION = STATION) -> new_wts

stratum %>%
  mutate(GIS_STATION = str_remove_all(GIS_STATION, "-")) -> stratum # remove dash between GIS station to match new_wts

raw.dat <- left_join(stratum, new_wts) 

# replace NA weights with 1
change <- is.na(raw.dat$weight)
raw.dat$weight[change] <- 1

plot.raw <- raw.dat %>%
  mutate(weighted_log_cpue = log_cpue*weight) %>%
  group_by(YEAR) %>%
  summarise(sum_weight = sum(weight),
            sum_weighted_log_cpue = sum(weighted_log_cpue))%>%
  mutate(log_mean = sum_weighted_log_cpue / sum_weight) %>%
  rename(year = YEAR) 

unweighted.raw <- raw.dat %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

# Plot
immature.stratum_drop_2nd <- ggplot(plot.dat, aes(year, log_mean, color = "black")) +
  geom_line(data = plot.raw, aes(year, log_mean, color = "red")) +
  geom_point(data = plot.raw, aes(year, log_mean, color = "red")) + 
  geom_line(data = unweighted.raw, aes(year, log_mean, color = "blue")) +
  geom_point(data = unweighted.raw, aes(year, log_mean, color = "blue")) + 
  geom_line() +
  scale_colour_manual(name = element_blank(), values = c("black", "red", "blue"), labels = c("Imputed", "Unweighted", "Weighted"))+
  geom_point(size = 2) +
  theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Total immature, stratum > 2nd percentile") + 
  theme(axis.title.x = element_blank())

immature.stratum_drop_2nd

# Create csv of imputed cpue, imputed sd, raw weighted cpue, raw unweighted cpue, and # stations sampled by year
allimm_drop2_df <- cbind(data.frame(year = plot.dat$year, imp_log_mean = plot.dat$log_mean, imp_sd = plot.dat$sd),
                         plot.raw %>% select(log_mean) %>% rename(wtd_log_mean = log_mean),
                         unweighted.raw %>% select(log_mean) %>% rename(unwtd_log_mean = log_mean),
                         count %>% select(count) %>% rename(n_stations = count))

write.csv(allimm_drop2_df, "./output/allimm_drop2_df.csv")

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

# remove stations with mean cpue = 0
drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

# drop catches < 10th quantile
drop_10th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.1))

ggplot(drop_10th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

male_30_59_drop_10_map <- ggplot(drop_10th, aes(longitude, latitude, color = mean_log_cpue)) +
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

male_30_59_drop_10_map 

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

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

male_30_59_drop10_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

male_30_59_drop10_hist

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

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled
area <- sc_strata %>%
  filter(YEAR == 2022) %>%
         #STATION_ID %in% unique(drop_10th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA))

weights <- data.frame(STATION = colnames(complete(imp_30_59_drop_10, action = 1))) %>%
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
                                  log_mean = apply(temp, 1, ff)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# replace sd = 0 with NA
check <- plot.dat$sd == 0
plot.dat$sd[check] <- NA

# summarize raw (non-imputed data) to plot
# first, get weights
raw_weights <- sc_strata %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(YEAR, STATION, DISTRICT, TOTAL_AREA)

weights1 <- distinct(data.frame(YEAR = sc_catch$YEAR, STATION = sc_catch$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION, "-")) %>%
  left_join(., raw_weights)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights1 %>%
  group_by(YEAR, DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights2 <- left_join(weights1, summary_weights)

weights2 %>%
  select(STATION, YEAR, DISTRICT, weight) %>%
  distinct() %>%
  as_tibble() %>%
  rename(GIS_STATION = STATION) -> new_wts

stratum %>%
  mutate(GIS_STATION = str_remove_all(GIS_STATION, "-")) -> stratum # remove dash between GIS station to match new_wts

raw.dat <- left_join(stratum, new_wts) 

# replace NA weights with 1
change <- is.na(raw.dat$weight)
raw.dat$weight[change] <- 1

plot.raw <- raw.dat %>%
  mutate(weighted_log_cpue = log_cpue*weight) %>%
  group_by(YEAR) %>%
  summarise(sum_weight = sum(weight),
            sum_weighted_log_cpue = sum(weighted_log_cpue))%>%
  mutate(log_mean = sum_weighted_log_cpue / sum_weight) %>%
  rename(year = YEAR) 

unweighted.raw <- raw.dat %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

# Plot
male_30_59_stratum_drop_10 <- ggplot(plot.dat, aes(year, log_mean, color = "black")) +
  geom_line(data = plot.raw, aes(year, log_mean, color = "red")) +
  geom_point(data = plot.raw, aes(year, log_mean, color = "red")) + 
  geom_line(data = unweighted.raw, aes(year, log_mean, color = "blue")) +
  geom_point(data = unweighted.raw, aes(year, log_mean, color = "blue")) + 
  geom_line() +
  scale_colour_manual(name = element_blank(), values = c("black", "red", "blue"), labels = c("Imputed", "Unweighted", "Weighted"))+
  geom_point(size = 2) +
  theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Male 30-59 mm, stratum > 10th percentile") + 
  theme(axis.title.x = element_blank())

male_30_59_stratum_drop_10

# Create csv of imputed cpue, imputed sd, raw weighted cpue, raw unweighted cpue, and # stations sampled by year
male3059_drop10_df <- cbind(data.frame(year = plot.dat$year, imp_log_mean = plot.dat$log_mean, imp_sd = plot.dat$sd),
                         plot.raw %>% select(log_mean) %>% rename(wtd_log_mean = log_mean),
                         unweighted.raw %>% select(log_mean) %>% rename(unwtd_log_mean = log_mean),
                         count %>% select(count) %>% rename(n_stations = count))

write.csv(male3059_drop10_df , "./output/male3059_drop10_df.csv")


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

# remove stations with mean cpue = 0
drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

# drop catches < 5th quantile
drop_5th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.05))

ggplot(drop_5th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

male_60_95_drop_5_map <- ggplot(drop_5th, aes(longitude, latitude, color = mean_log_cpue)) +
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

male_60_95_drop_5_map

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

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

male_60_95_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

male_60_95_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_male_60_95_stratum_drop_5.RDS")
imp_60_95_drop_5 <- readRDS("./output/abundance_imputations_male_60_95_stratum_drop_5.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp_60_95_drop_5))

sum(check) # 0


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled
area <- sc_strata %>%
  filter(YEAR == 2022) %>%
         #STATION_ID %in% unique(drop_5th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area <- area %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA))

weights <- data.frame(STATION = colnames(complete(imp_60_95_drop_5, action = 1))) %>%
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
                                  log_mean = apply(temp, 1, ff)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# replace sd = 0 with NA
check <- plot.dat$sd == 0
plot.dat$sd[check] <- NA

# summarize raw (non-imputed data) to plot
# first, get weights
raw_weights <- sc_strata %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(YEAR, STATION, DISTRICT, TOTAL_AREA)

weights1 <- distinct(data.frame(YEAR = sc_catch$YEAR, STATION = sc_catch$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION, "-")) %>%
  left_join(., raw_weights)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights1 %>%
  group_by(YEAR, DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights2 <- left_join(weights1, summary_weights)

weights2 %>%
  select(STATION, YEAR, DISTRICT, weight) %>%
  distinct() %>%
  as_tibble() %>%
  rename(GIS_STATION = STATION) -> new_wts

stratum %>%
  mutate(GIS_STATION = str_remove_all(GIS_STATION, "-")) -> stratum # remove dash between GIS station to match new_wts

raw.dat <- left_join(stratum, new_wts) 

# replace NA weights with 1
change <- is.na(raw.dat$weight)
raw.dat$weight[change] <- 1

plot.raw <- raw.dat %>%
  mutate(weighted_log_cpue = log_cpue*weight) %>%
  group_by(YEAR) %>%
  summarise(sum_weight = sum(weight),
            sum_weighted_log_cpue = sum(weighted_log_cpue))%>%
  mutate(log_mean = sum_weighted_log_cpue / sum_weight) %>%
  rename(year = YEAR) 

unweighted.raw <- raw.dat %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

# Plot
male_60_95_stratum_drop_5 <- ggplot(plot.dat, aes(year, log_mean, color = "black")) +
  geom_line(data = plot.raw, aes(year, log_mean, color = "red")) +
  geom_point(data = plot.raw, aes(year, log_mean, color = "red")) + 
  geom_line(data = unweighted.raw, aes(year, log_mean, color = "blue")) +
  geom_point(data = unweighted.raw, aes(year, log_mean, color = "blue")) + 
  geom_line() +
  scale_colour_manual(name = element_blank(), values = c("black", "red", "blue"), labels = c("Imputed", "Unweighted", "Weighted"))+
  geom_point(size = 2) +
  theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Male 60-95 mm, stratum > 5th percentile") + 
  theme(axis.title.x = element_blank())

male_60_95_stratum_drop_5

# Create csv of imputed cpue, imputed sd, raw weighted cpue, raw unweighted cpue, and # stations sampled by year
male6095_drop5_df <- cbind(data.frame(year = plot.dat$year, imp_log_mean = plot.dat$log_mean, imp_sd = plot.dat$sd),
                            plot.raw %>% select(log_mean) %>% rename(wtd_log_mean = log_mean),
                            unweighted.raw %>% select(log_mean) %>% rename(unwtd_log_mean = log_mean),
                            count %>% select(count) %>% rename(n_stations = count))

write.csv(male6095_drop5_df , "./output/male6095_drop5_df.csv")


## female immature > 30 mm ---------------------------
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 2,
         CLUTCH_SIZE == 0,
         WIDTH >= 30,
         SHELL_CONDITION <= 2) %>%
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

# remove stations with mean cpue = 0
check <- stratum %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue),
            latitude = mean(MID_LATITUDE),
            longitude = mean(MID_LONGITUDE))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

drop_0 <- check %>%
  filter(mean_log_cpue > 0)

ggplot(drop_0, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

# drop 10th quantile
drop_10th <- drop_0 %>%
  filter(mean_log_cpue > quantile(drop_0$mean_log_cpue, 0.1))

ggplot(drop_10th, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")

female_immature_drop_10_map <- ggplot(drop_10th, aes(longitude, latitude, color = mean_log_cpue)) +
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

female_immature_drop_10_map

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

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

female_immature_drop10_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

female_immature_drop10_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)
# 
imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/abundance_imputations_female_immature_stratum_drop_10.RDS")
imp_female_drop_10 <- readRDS("./output/abundance_imputations_female_immature_stratum_drop_10.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp_female_drop_10))

sum(check) # 2!


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) 

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled
area_full <- sc_strata %>%
  filter(YEAR == 2022) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area_full <- area_full %>%
  group_by(DISTRICT) %>%
  summarise(full_count = n())

area_10 <- sc_strata %>%
  filter(YEAR == 2022)%>%
         #STATION_ID %in% unique(drop_10th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

sum_area_10 <- area_10 %>%
  group_by(DISTRICT) %>%
  summarise(total_count = n(),
            total_area = mean(TOTAL_AREA)) %>%
  left_join(., sum_area_full) %>%
  mutate(area_10 = (total_count / full_count) * total_area) %>%
  select(DISTRICT, area_10) %>%
  rename(weight = area_10)

weights <- data.frame(STATION = colnames(complete(imp_female_drop_10, action = 1))) %>%
  left_join(., area_10) %>%
  select(-TOTAL_AREA) %>%
  left_join(.,sum_area_10)

# # get proportion of each stratum that is included in the imputation
# summary_weights <- weights %>%
#   group_by(DISTRICT) %>%
#   summarise(count_sampled = n()) %>%
#   left_join(., sum_area) %>%
#   mutate(weight = (count_sampled / total_count) * total_area) %>%
#   select(DISTRICT, weight)

# weights <- left_join(weights, summary_weights)

# set up weighted mean function
ff <- function(x) weighted.mean(x, weights$weight, na.rm = T)

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
                                  log_mean = apply(temp, 1, ff)))
}


# View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))

# replace sd = 0 with NA
check <- plot.dat$sd == 0
plot.dat$sd[check] <- NA

# summarize raw (non-imputed data) to plot
# first, get weights
raw_weights <- sc_strata %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(YEAR, STATION, DISTRICT, TOTAL_AREA)

weights1 <- distinct(data.frame(YEAR = sc_catch$YEAR, STATION = sc_catch$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION, "-")) %>%
  left_join(., raw_weights)

# get proportion of each stratum that is included in the imputation
summary_weights <- weights1 %>%
  group_by(YEAR, DISTRICT) %>%
  summarise(count_sampled = n()) %>%
  left_join(., sum_area) %>%
  mutate(weight = (count_sampled / total_count) * total_area) %>%
  select(DISTRICT, weight)

weights2 <- left_join(weights1, summary_weights)

weights2 %>%
  select(STATION, YEAR, DISTRICT, weight) %>%
  distinct() %>%
  as_tibble() %>%
  rename(GIS_STATION = STATION) -> new_wts

stratum %>%
  mutate(GIS_STATION = str_remove_all(GIS_STATION, "-")) -> stratum # remove dash between GIS station to match new_wts

raw.dat <- left_join(stratum, new_wts) 

# replace NA weights with 1
change <- is.na(raw.dat$weight)
raw.dat$weight[change] <- 1

plot.raw <- raw.dat %>%
  mutate(weighted_log_cpue = log_cpue*weight) %>%
  group_by(YEAR) %>%
  summarise(sum_weight = sum(weight),
            sum_weighted_log_cpue = sum(weighted_log_cpue))%>%
  mutate(log_mean = sum_weighted_log_cpue / sum_weight) %>%
  rename(year = YEAR) 

unweighted.raw <- raw.dat %>%
  group_by(YEAR) %>%
  rename(year = YEAR) %>%
  summarise(log_mean = mean(log_cpue))

# Plot
female.stratum_drop_10th <- ggplot(plot.dat, aes(year, log_mean, color = "black")) +
  geom_line(data = plot.raw, aes(year, log_mean, color = "red")) +
  geom_point(data = plot.raw, aes(year, log_mean, color = "red")) + 
  geom_line(data = unweighted.raw, aes(year, log_mean, color = "blue")) +
  geom_point(data = unweighted.raw, aes(year, log_mean, color = "blue")) + 
  geom_line() +
  scale_colour_manual(name = element_blank(), values = c("black", "red", "blue"), labels = c("Imputed", "Unweighted", "Weighted"))+
  geom_point(size = 2) +
  theme(legend.position = "bottom")+
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd)) +
  ggtitle("Female immature, stratum > 10th percentile") + 
  theme(axis.title.x = element_blank())

female.stratum_drop_10th

# Create csv of imputed cpue, imputed sd, raw weighted cpue, raw unweighted cpue, and # stations sampled by year
female_drop10_df <- cbind(data.frame(year = plot.dat$year, imp_log_mean = plot.dat$log_mean, imp_sd = plot.dat$sd),
                            plot.raw %>% select(log_mean) %>% rename(wtd_log_mean = log_mean),
                            unweighted.raw %>% select(log_mean) %>% rename(unwtd_log_mean = log_mean),
                            count %>% select(count) %>% rename(n_stations = count))

write.csv(female_drop10_df , "./output/female_drop10_df.csv")

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
