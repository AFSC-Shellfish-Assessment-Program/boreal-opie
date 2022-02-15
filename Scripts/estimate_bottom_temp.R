library(tidyverse)
library(lubridate)
library(mice)

dat <- read.csv("./data/haul.csv")

head(dat)

# extract gear temp, date, year, gis_station

dat <- dat %>%
  mutate(julian=yday(parse_date_time(START_TIME, "d-m-y", "US/Alaska"))) %>%
  select(julian, GEAR_TEMPERATURE, SURVEY_YEAR, GIS_STATION) %>%
  rename(bottom.temp = GEAR_TEMPERATURE, year = SURVEY_YEAR, station = GIS_STATION)

# check for  which stations are sampled each year
ff <- function(x) sum(!is.na(x))

check <- tapply(dat$julian, list(dat$year, dat$station), ff)
View(check)

# quite a few are sampled only once!

check <- dat %>%
  group_by(station) %>%
  summarise(count = ff(julian)) %>% 
  arrange(count)

check            

ggplot(check, aes(count)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

# so we can cut off above 20...
# but it also looks like we have some that were sampled > once in a year (> 46 samples in 46-year time series)
# might make the multiple imputation more complicated (i.e., if we want to use a year x station table for imputation!)

# first, cut off all the n < 20 sets

keep <- check$station[check$count > 20]

dat <- dat %>%
  filter(station %in% keep)

# look at these stations sampled more than once a year 

check <- dat %>%
  group_by(year, station) %>%
  summarize(count = n()) %>%
  filter(count > 1)

View(check)

# make sure there are a max of 2 tows per station/year combo
max(check$count) # yes

# aha - I think these are the BBRKC retows!

# easiest/best not to use these for temperature
# bad to throw away information, but these are out of the norm re julian:station combinations,
# so might be too complicated to try to analyze (in addition to the difficulties in mice)

# need to go through the situations where 2 tows are made and use the first

# set up an index for the ones we want to drop
dat$station_year = paste(dat$station, dat$year, sep = "_")

check$station_year = paste(check$station, check$year, sep = "_")

# now pick these out of dat

doubles <- dat %>%
  filter(station_year %in% check$station_year)

# confirm this has the double sample events

check.double <- doubles %>%
  group_by(station_year) %>%
  summarise(count = n())

unique(check.double$count) # good...

# now select the maximum date for each station-year
# these are the ones to drop from dat

drop.check <- doubles %>%
  group_by(station_year) %>%
  summarise(julian = max(julian)) %>% # picks the second day in each station-year combo
  mutate(station_year_julian = paste(station_year, julian, sep = "_")) 

# and drop these from dat
dat$station_year_julian <- paste(dat$station_year, dat$julian, sep = "_")
  
dat <- dat %>%
  filter(!station_year_julian %in% drop.check$station_year_julian)

# check
c.dat <- dat %>%
  group_by(year, station) %>%
  summarise(count = n())

unique(c.dat$count) # looks good

## input missing values --------------------------

# here is our workflow going forwards:

# 1) make data into df with row = year, columns = station-day and station-temp

# 2) impute 100 times

# 3) combine these 100 imputed data frames back into full versions of the data (including julian day)

# 4) feed this list of imputed data frames into a brms model to estimate year and sst:station effects


# set up for multiple imputation

# first, clean up dat
dat <- dat %>%
  select(julian, bottom.temp, year, station)

stations <- unique(dat$station) 

impute.this <- data.frame(year = c(1975:2019, 2021))

for(i in 1:length(stations)){
  
  # i <- 1
  
  station.dat <- dat %>%
    filter(station == stations[i]) %>%
    select(-station)
  
  names(station.dat)[1] <- paste("j", i, sep = "")
  names(station.dat)[2] <- paste("t", i, sep = "")
  
  
  impute.this <- left_join(impute.this, station.dat)
  
}

View(impute.this)

# check station that is returning NA for imputed date!
check <- which(stations == "GF2120")

index.check <- grep("248", names(impute.this))
impute.this[,index.check]

# remove year!
impute.this <- impute.this %>%
  select(-year)

# examine correlations
cors <- cor(impute.this, use = "p")
View(cors)

# limit to columns with missing cases!
ff <- function(x) sum(!is.na(x))

check <- apply(impute.this, 2, ff)

max(check) # that's right!

keep <- check[check < 46]

keep

missing.cors <- cors[,colnames(cors) %in% names(keep)]

# figure out how many examples would be included for each correlation cutoff

# remove cor = 1
change <- missing.cors == 1 

missing.cors[change] <- NA

ff <- function(x) sum(abs(x) > 0.65, na.rm = T)

count <- apply(missing.cors, 2, ff)

range(count)

# big ranges!

# let's figure out the mincor that will produce 15 predictors for each
plyr::round_any(18.129, 0.01, f = floor)
ff <- function(x) plyr::round_any(sort(abs(x), decreasing = T)[20], 0.001, f = floor)

mincors <- apply(cors, 2, ff)


# impute 100 times - using 15 predictors for each
pred <- quickpred(impute.this, mincor = as.vector(mincors))

table(rowSums(pred)) 
# perfect - at least 19 predictors for each, and vast majority are in 15-25 range
# so we have some buffer for dropped predictors in logged events

# and impute
imp <- mice(impute.this, m = 100, pred = pred)

head(imp$loggedEvents)
unique(imp$loggedEvents$meth)
sum(imp$loggedEvents$meth == "pmm")

# so we're dealing with 1 collinear variable and all the rest are pmm

str_count((imp$loggedEvents$out[2]), pattern = ",")

dropped <- data.frame(dropped = 1 + str_count((imp$loggedEvents$out), pattern = ","))

str(dropped)

ggplot(dropped, aes(dropped)) +
  geom_histogram(bins = 6, fill = "grey", color = "black")

# so not lots of variables dropped - enough remaining for prediction! 
# will need to look into the pmm warning (haven't found much immediately), but proceeding for now

# now need to put imp into correct form for analysis

str(imp$imp)

View(complete(imp))

# are there NAs in complete(imp)?

check <- is.na(complete(imp))

sum(check)
check <- which(is.na(complete(imp)))
complete
# this looks great....

# save!
saveRDS(imp, "./output/imputed_julian_temp")

# need to replace index values from "stations" with station name, then put into long form 
# with year, station, julian, and temperature columns

# then add to list of data frames

imputed.data <- list()

# also create df to save mean annual temp and sampling day for each imputed temperature data set

imputed.means <- data.frame()

# this is clunky but should work!

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  
  # split colname
  split1 <- str_sub(colnames(temp), 1L, 1L) 
  
  # change to "julian" and "temperature"
  change <- split1 == "j"
  split1[change] <- "julian"
  
  change <- split1 == "t"
  split1[change] <- "temperature" 
  
  # combine with stations
  colnames(temp) <- str_c(split1, stations, sep = "_")
  
  # split into temp and day then recombine
  t1 <- grep("temp", colnames(temp))
  temperature.dat <- temp[,t1]  
  
  t2 <- grep("julian", colnames(temp))
  julian.dat <- temp[,t2]   
  
  temperature.dat$year <- julian.dat$year <- c(1975:2019, 2021)
  
  #cool - pivot both longer, split out leading text (julian/temperature) and join by year / station!
  
  temperature.dat <- temperature.dat %>%
    pivot_longer(cols = -year) %>%
    mutate(name = str_sub(name, 13L, -1L)) %>%
    rename(temperature = value)
  
  julian.dat <- julian.dat %>%
    pivot_longer(cols = -year) %>%
    mutate(name = str_sub(name, 8L, -1L)) %>%
    rename(julian = value)
  
  find <- is.na(julian.dat$julian)
  sum(find)
  julian.dat[find,]
  #join into single df
  imputed <- as.data.frame(left_join(temperature.dat, julian.dat))
  
  ff <- function(x) sum(!is.na())
  
  temp.imputed <- imputed %>%
    select(year, temperature) %>%
    na.omit() %>%
    group_by(year) %>%
    summarise(mean.temp = mean(temperature),
              n.temp = n())
  
  temp.day <- imputed %>%
    select(year, julian) %>%
    na.omit() %>%
    group_by(year) %>%
    summarise(mean.day= mean(julian),
              n.day = n())
  
  temp.temp <- left_join(temp.imputed, temp.day)
  
  temp.temp$imputation = i
    
  # add to summary of each imputation
  imputed.means <- rbind(imputed.means, temp.temp)
    
  # and add to list of dfs for analysis in brms
  imputed.data[[i]] <- imputed 
  
  }
  
check <- is.na(imputed.means)
sum(check)

plot.means <- imputed
ggplot(imputed.means, aes())
  
  # examine imputation
plot.imp <- imputed.means %>%
  group_by(year) %>%
  summarise(mean.temp = mean(mean.temp),
            mean.day = mean(mean.day))
# can't use summarise for sd b/c years with no NA will return NA (sd = 0)

plot.imp$sd.temp <- tapply(imputed.means$mean.temp, imputed.means$year, sd)
plot.imp$sd.day <- tapply(imputed.means$mean.day, imputed.means$year, sd)

ggplot(plot.imp, aes(year, mean.temp)) +
  geom_point() +
  geom_errorbar(aes(ymax = mean.temp + 2*sd.temp, ymin = mean.temp - 2*sd.temp))

ggplot(plot.imp, aes(year, mean.day)) +
  geom_point() +
  geom_errorbar(aes(ymax = mean.day + 2*sd.day, ymin = mean.day - 2*sd.day))

ggplot(plot.imp, aes(mean.day, mean.temp)) +
  geom_point()

# so the only thing that seems unclear is the 6 incidences when day cannot be estimated by multiple imputation

# pass along to brms - seems that the effect will just be a reduction in # of df passed to brms
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./scripts/stan_utils.R")

theme_set(theme_bw())

## brms: setup ---------------------------------------------

## Define model formulas

brms_formula <-  bf(temperature ~ as.factor(year) + s(julian, k = 3, by = name))

## fit --------------------------------------
bottom_temp_brm <- brm_multiple(brms_formula,
                    data = imputed.data,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10), 
                    verbose = T)



codR_dfa_brm  <- add_criterion(codR_dfa_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_dfa_brm, file = "output/codR_dfa_brm.rds")

codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
check_hmc_diagnostics(codR_dfa_brm$fit)
neff_lowest(codR_dfa_brm$fit)
rhat_highest(codR_dfa_brm$fit)
summary(codR_dfa_brm)
bayes_R2(codR_dfa_brm)
plot(codR_dfa_brm$criteria$loo, "k")
plot(conditional_effects(codR_dfa_brm), ask = FALSE)
y <- dfa$model
yrep_codR_dfa_brm  <- fitted(codR_dfa_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_dfa_brm[sample(nrow(yrep_codR_dfa_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("codR_dfa_brm")
pdf("./figs/trace_codR_dfa_brm.pdf", width = 6, height = 4)
trace_plot(codR_dfa_brm$fit)
dev.off()