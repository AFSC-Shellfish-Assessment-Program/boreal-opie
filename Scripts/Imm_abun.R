# notes ----

# Calculate immature snow crab abundance in EBS

#NOTE: There are no pre-1980 stratum tables on akfin so abundances were only 
  #calculated for 1980+

# Author: Erin Fedewa

# load ----
library(tidyverse)
library(ggmap)
library(Hmisc)

##############################################

## EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")
sc_catch <- read.csv("./Data/CRABHAUL_OPILIO.csv")


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

# # stations sampled at least 40 times
# count.stations <- haul %>%
#   filter(HAUL_TYPE == 3) %>%
#   group_by(GIS_STATION) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count))
# 
# View(count.stations)

# # load core stations from Erin's analysis
# core <- read.csv("./output/core_stations.csv")
# 
# # combine with haul
# core <- core %>%
#   select(GIS_STATION) %>%
#   left_join(., haul)
# 
# # plot these
# ggplot(core, aes(MID_LONGITUDE, MID_LATITUDE)) +
#   geom_point() +
#   facet_wrap(~SURVEY_YEAR)


# # limit to 33+ (regular stations!)
# 
# use.stations <- count.stations %>%
#   filter(count >= 33)

## restrict to stratum stations

# get count of station samples in stratum
count <- sc_strata %>%
  group_by(STATION_ID) %>%
  summarise(count = n()) %>%
  arrange(count)

View(count) # all good - 33 or more

use <- unique(count$STATION_ID)

core <- haul %>% 
  filter(GIS_STATION %in% use)

# plot these
ggplot(core, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# look at histogram of mean log catch by year to select stations



# define survey year in catch
sc_catch <- sc_catch %>%
  mutate(YEAR = floor(CRUISE/100))

# rename year in strata
sc_strata <- sc_strata %>%
  rename(YEAR = SURVEY_YEAR)

# # check size of stratum in each year
# check <- sc_strata %>%
#   group_by(YEAR, DISTRICT) %>%
#   summarise(area = mean(TOTAL_AREA)) %>% 
#   group_by(YEAR) %>%
#   summarise(area = sum(area))
# 
# ggplot(check, aes(YEAR, area)) + 
#   geom_line() +
#   geom_point() # big differences!


# calculate CPUE by station for immature male snow crab 
scratch <- sc_catch %>%
  filter(HAUL_TYPE == 3, 
         SEX == 1,
         WIDTH >= 45 & WIDTH <= 95) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) 

# join with core stations to include 0 catches 
core <- core %>%
  select(GIS_STATION, SURVEY_YEAR, MID_LONGITUDE, MID_LATITUDE) %>%
  rename(YEAR = SURVEY_YEAR) %>%
  left_join(., scratch)

ggplot(core, aes(MID_LONGITUDE, MID_LATITUDE)) +
  geom_point() +
  facet_wrap(~YEAR)

# replace NA CPUE with 0
change <- is.na(core$CPUE)
core$CPUE[change] <- 0

# check samples per year
check <- core %>%
  group_by(YEAR) %>%
  summarise(count = n())

check

# log transform CPUE
core$log_cpue <- log(core$CPUE + 1)

# plot histogram of mean cpue per station

check <- core %>%
  group_by(GIS_STATION) %>%
  summarise(mean_log_cpue = mean(log_cpue))

ggplot(check, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  

# get into matrix form for mice
dat <- tapply(core$log_cpue, list(core$YEAR, core$GIS_STATION), mean)

# examine correlations
r <- rcorr(as.matrix(dat))$r 
r #Cross-year correlations between each station combination

hist(r)
# choose 25 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=15))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) 

saveRDS(imp, "./output/abundance_imputations.RDS")
imp <- readRDS("./output/abundance_imputations.RDS")


# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check) # 0


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.dat <- data.frame()

# this is clunky but should work!

# create bounds for imputed values 
# (min = 0, max = 0.5 sd above observed)

lower_bound <- 0
upper_bound <- max(dat, na.rm = T) + 0.5*sd(dat, na.rm = T)

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
                                  log_mean = rowMeans(temp)))
}


View(imputed.dat)

# not using dplyr because sd using summarize cannot handle sd = 0
plot.dat <- data.frame(year = c(1975:2019, 2021:2022),
                       log_mean = tapply(imputed.dat$log_mean, imputed.dat$year, mean),
                       sd = tapply(imputed.dat$log_mean, imputed.dat$year, sd))
  
ggplot(plot.dat, aes(year, log_mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = log_mean - 2*sd,
                    ymax = log_mean + 2*sd))


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
