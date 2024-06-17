library(tidyverse)
library(lubridate)
library(mice)
library(Hmisc)

dat <- read.csv("./data/HAUL_NEWTIMESERIES.csv")

head(dat)

# load immature opilio core habitat
imm_area <- read.csv("./Data/imm_area_50perc.csv")

# plot what we have
plot.dat <- dat %>%
  select(SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE) %>%
  rename(LATITUDE = MID_LATITUDE,
         LONGITUDE = MID_LONGITUDE)

# add imm_area
add_area <- imm_area %>%
  select(GIS_STATION) %>%
  mutate(core = "core")

plot.dat <- left_join(plot.dat, add_area)

change <- is.na(plot.dat$core)
plot.dat$core[change] <- "non-core"

ggplot(filter(plot.dat, SURVEY_YEAR <= 1982), aes(LONGITUDE, LATITUDE, color = core)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

# different names in different years??
check79 <- plot.dat %>%
  filter(SURVEY_YEAR == 1979,
         core == "non-core")

check80 <- plot.dat %>%
  filter(SURVEY_YEAR == 1980,
         core == "core")

whats.this <- intersect(unique(check79$GIS_STATION), unique(check80$GIS_STATION)) # none

# let's use a list of core stations sampled in 1977 to generate our bottom temp index
use <- plot.dat %>%
  filter(SURVEY_YEAR == 1975,
          core == "core")
  
use.stations <- use$GIS_STATION

plot.dat <- plot.dat %>%
  filter(GIS_STATION %in% use.stations)

ggplot(filter(plot.dat, SURVEY_YEAR <= 1982), aes(LONGITUDE, LATITUDE, color = core)) +
  geom_point() +
  facet_wrap(~SURVEY_YEAR)

check <- plot.dat %>%
  group_by(SURVEY_YEAR) %>%
  summarise(missing = length(use.stations)-n())

check

# extract gear temp, date, year, gis_station

dat <- dat %>%
  mutate(julian=yday(parse_date_time(START_TIME, "d-m-y", "US/Alaska"))) %>%
  select(julian, GEAR_TEMPERATURE, SURVEY_YEAR, GIS_STATION) %>%
  rename(bottom.temp = GEAR_TEMPERATURE, year = SURVEY_YEAR, station = GIS_STATION) %>%
  filter(station %in% use.stations)



# look at these stations sampled more than once a year 

check <- dat %>%
  group_by(year, station) %>%
  dplyr::summarize(count = n()) %>%
  filter(count > 1)

check # none!


## impute missing values --------------------------

# here is our workflow going forwards:

# 1) make data into tow df with row = year, columns = station, and data = either day or temp

# 2) impute each 100 times

# 3) use mean temps / dates from imputed data frames in a model to estimate year and sst:station effects
# (for core opilio area only!)


# set up for multiple imputation

# first, clean up dat
dat.julian <- dat %>%
  select(julian, year, station) %>%
  pivot_wider(names_from = station, values_from = julian) %>%
  arrange(year) %>%
  select(-year)

# check number missing
f <- function(x) sum(is.na(x))

check <- apply(dat.julian, 1, f)
check

# examine correlations
r <- rcorr(as.matrix(dat.julian))$r 
r #Cross-year correlations between each station combination


# choose 30 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=30))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

colnames(pred) <- colnames(dat.julian) <- str_remove_all(colnames(pred), "-")

imp <- mice(data = dat.julian, method = "norm.predict", m=100)#, pred = pred) #Using Bayesian linear regression method
saveRDS(imp, "./output/station_julian_day_imputations.RDS")
imp <- readRDS("./output/station_julian_day_imputations.RDS")

str(imp$imp)

View(complete(imp))

# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check)
check <- which(is.na(complete(imp)))
check
# this looks great....


# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.day <- data.frame()

# this is clunky but should work!

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  imputed.day <- rbind(imputed.day,
                         data.frame(imputation = i,
                                    year = c(1975:2019, 2021, 2022),
                                    mean.day = rowMeans(temp)))
}
  
imputed.day <- imputed.day %>%
  pivot_wider(names_from = imputation,
              values_from = mean.day)

# all are identical!

plot.dat <- data.frame(year = c(1975:2019, 2021, 2022),
                       imputed.mean.day = rowMeans(imputed.day[,2:101]))
  
add.dat <- dat %>%
  group_by(year) %>%
  dplyr::summarize(observed.mean.day = mean(julian, na.rm = T))


plot.dat <- left_join(plot.dat, add.dat) %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value, color = name)) +
  geom_line() +
  geom_point()

## now impute temperature ----------------------
# first, clean up dat
dat.temp <- dat %>%
  select(bottom.temp, year, station) %>%
  pivot_wider(names_from = station, values_from = bottom.temp) %>%
  arrange(year) %>%
  select(-year)

# check number missing
f <- function(x) sum(is.na(x))

check <- apply(dat.julian, 1, f)
check

colnames(dat.temp) <- str_remove_all(colnames(dat.temp), "-")

imp <- mice(data = dat.temp, method = "norm.predict", m=100)#, pred = pred) #Using Bayesian linear regression method
saveRDS(imp, "./output/station_bottom_temp_imputations.RDS")
imp <- readRDS("./output/station_bottom_temp_imputations.RDS")

imputed.temp <- data.frame()

# this is clunky but should work!

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  imputed.temp <- rbind(imputed.temp,
                       data.frame(imputation = i,
                                  year = c(1975:2019, 2021, 2022),
                                  mean.temp = rowMeans(temp)))
}

imputed.temp <- imputed.temp %>%
  pivot_wider(names_from = imputation,
              values_from = mean.temp)

# all are identical!

plot.dat <- data.frame(year = c(1975:2019, 2021, 2022),
                       imputed.mean.temp = rowMeans(imputed.temp[,2:101]))

add.temp <- dat %>%
  group_by(year) %>%
  dplyr::summarize(observed.mean.temp = mean(bottom.temp, na.rm = T))


plot.dat <- left_join(plot.dat, add.temp) %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value, color = name)) +
  geom_line() +
  geom_point()


imputed.dat <- data.frame(year_fac = as.factor(c(1975:2019, 2021, 2022)),
                          mean.day = rowMeans(imputed.day[,2:101]),
                       mean.temp = rowMeans(imputed.temp[,2:101]))


# get station-level data for GAM
station.imputed.temp <- data.frame()


for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i) %>%
    mutate(year = c(1975:2019, 2021, 2022)) %>%
    pivot_longer(cols = -year) 
  
  station.imputed.temp <- rbind(station.imputed.temp,
                        temp)
}

station.imputed.temp <- station.imputed.temp %>%
  group_by(year, name) %>%
  summarise(bottom.temp = mean(value)) %>%
  rename(station = name)

# now day
imp <- readRDS("./output/station_julian_day_imputations.RDS")


station.imputed.day <- data.frame()



for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i) %>%
    mutate(year = c(1975:2019, 2021, 2022)) %>%
    pivot_longer(cols = -year) 
  
  station.imputed.day <- rbind(station.imputed.day,
                                temp)
}

station.imputed.day <- station.imputed.day %>%
  group_by(year, name) %>%
  summarise(day = mean(value)) %>%
  rename(station = name)

gam.dat <- left_join(station.imputed.day, station.imputed.temp) %>%
  mutate(year_fac = as.factor(year),
         station_fac = as.factor(station))


mod <- mgcv::gamm(bottom.temp ~ year_fac + s(day, k = 5), random = list(station=  ~ 1), data = gam.dat)
summary(mod$gam)
plot(mod$gam)

newdat <- data.frame(year_fac = as.factor(c(1975:2019, 2021, 2022)),
                     day = mean(gam.dat$day))

annual.temp <- predict(mod$gam, newdata = newdat)


estimated.temp <- data.frame(year = c(1975:2022),
                   bottom.temp = c(annual.temp[1:45], NA, annual.temp[46:47]))

ggplot(estimated.temp, aes(year, bottom.temp)) +
  geom_point() + 
  geom_line()

write.csv(estimated.temp, "./Data/date_corrected_bottom_temp.csv", row.names = F)

