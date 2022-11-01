# process directed fishery bycatch data for predicting 2020 survey abundance

library(tidyverse)

dat <- read.csv("./Data/raw_df.csv", row.names = 1)

str(dat)

# get julian day, month and year

dat <- dat %>%
  mutate(julian = lubridate::yday(lubridate::parse_date_time(dat$sampdate, orders = "mdy")),
         month = lubridate::month(lubridate::parse_date_time(dat$sampdate, orders = "mdy")),
         year = lubridate::year(lubridate::parse_date_time(dat$sampdate, orders = "mdy")))

unique(dat$spn)
length(unique(dat$spn))

unique(dat$year)

dat <- dat %>%
  mutate(trip_pot = paste(trip, spn, sep = "_"))
length(unique(dat$trip_pot))

summary <- dat %>%
  select(year, month, trip_pot) %>%
  group_by(trip_pot) %>% 
  summarise(month = mean(month),
            year = mean(year)) %>%
  group_by(year, month) %>%
  summarise(count = n()) 

View(summary)

head(dat)

summary <- dat %>%
  filter(sex == 2) %>%
  group_by(year, month, clutch) %>%
  summarise(count = n()) 

# split into male and female immature
male <- dat %>%
  filter(sex == 1, size < 95) 

ggplot(male, aes(size)) +
  geom_histogram()

# check maturity codes for females
check <- dat %>%
  filter(sex == 2) %>%
  group_by(clutch, maturity) %>%
  summarise(count = n())

View(check)

female <- dat %>%
  filter(sex == 2)

ggplot(female, aes(size)) + 
  geom_histogram(bins = 30) +
  facet_wrap(~maturity, scales = "free_y") 

# looks like they're all immature!


# catch per pot 

catch <- male %>%
  group_by(trip_pot) %>%
  summarise(catch_95 = n())

# get a df with one entry for each pot
unique_pots <- data.frame(trip_pot = unique(dat$trip_pot))

# join to dat (one entry per pot)
temp.dat <- dat %>% 
  group_by(trip_pot) %>% 
  summarise(julian = mean(julian),
            month = mean(month),
            year = mean(year),
            soaktime = mean(soaktime))

unique_pots <- left_join(unique_pots, temp.dat) 
unique_pots <- left_join(unique_pots, catch)

unique_pots <- unique_pots %>%
  mutate(january_year = if_else(month > 8, year+1, year))

# check for NA
check <- is.na(unique_pots$catch_95)
sum(check)

# change to 0!
unique_pots$catch_95[check] <- 0

unique_pots$log_catch <- log(unique_pots$catch_95+1)
min(unique_pots$log_catch)

ggplot(unique_pots, aes(julian, log_catch)) +
  geom_point() +
  facet_wrap(~january_year) +
  geom_smooth(method = "lm")

summary <- unique_pots %>%
  mutate(month_fac = as.factor(as.character(month))) %>%
  group_by(january_year, month_fac) %>%
  summarize(count = n())

order <- data.frame(month_fac = levels(summary$month_fac),
                    order = c(3,1,2,4:8))

summary <- left_join(summary, order)

summary$month_fac <- reorder(summary$month_fac, summary$order)

ggplot(summary, aes(month_fac, count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~january_year)

# keep months with at least 50 pots
keep <- summary %>%
  filter(count >= 50)

catch <- unique_pots %>%
  mutate(month_fac = as.factor(as.character(month))) %>%
  group_by(january_year, month_fac) %>%
  summarise(mean_catch = mean(catch_95))

keep_catch <- left_join(keep, catch)

# fill in NAs

all_months <- data.frame(january_year = 1995:2022,
                         month_fac = rep(1:12, each = 28)) %>%
  filter(month_fac %in% c(1:6, 11:12)) %>%
  mutate(month_fac = as.factor(as.character(month_fac)))


all_months <- left_join(all_months, order)

all_months$month_fac <- reorder(all_months$month_fac, all_months$order)                         

all_months <- left_join(all_months, keep_catch)

ggplot(all_months, aes(january_year, mean_catch)) +
  geom_line() +
  geom_point() +
  facet_wrap(~month_fac)

all_months <- all_months %>%
  filter(month_fac %in% c(1:4))

jan_apr_catch <- unique_pots %>%
  filter(month %in% 1:4) %>%
  group_by(january_year) %>%
  summarise(mean_catch = mean(catch_95))

ggplot(jan_apr_catch, aes(january_year, mean_catch)) +
  geom_line() +
  geom_point() 
