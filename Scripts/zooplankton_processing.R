# create an index of annual Calanus and Pseudocalanus abundance
# controlling for spatial-temporal differences in sampling among years

library(tidyverse)
library(mgcv)

#Last update: 1/13/2023 (modified dataset from D. Kimmel)

theme_set(theme_bw())

# load data
dat <- read.csv("./Data/Calanus_Pseudo_Combined.csv")

##############################
# Data wrangling
dat <- dat %>%
  mutate(julian = 
           lubridate::yday(lubridate::parse_date_time(paste(DAY, MONTH, YEAR, sep = "-"), order = "dmy")),
         year_fac = as.factor(YEAR)) %>%
        mutate(abundance = exp(log10_EST_NUM_PERM3)) %>%
         rename(log_abundance = log10_EST_NUM_PERM3) %>%
          rename_with(tolower)

#Mean date sampled 
dat %>%
  group_by(year) %>%
  summarise(mean_date = mean(julian, na.rm=T),
            min_date = min(julian, na.rm=T),
            max_date = max(julian, na.rm=T))
#plot date sampled 
dat %>%
  ggplot(aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year)  

#Plot stations sampled 
dat %>%
  group_by(year, lat, lon) %>%
  distinct(station_name) %>%
  ggplot() +
  geom_point(aes(x =lon, y = lat), size=.5) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~year) #very few stations prior to 2005 

#plot # stations sampled
dat %>%
  group_by(year) %>%
  summarise(num_stations = n_distinct(station_name, na.rm=T)) %>%
  ggplot(aes(x=year, y=num_stations, group=1)) +
  geom_point() +
  geom_line() 

#plot psuedocalanus mean abundance: full timeseries
dat %>%
  filter(taxa == "Pseudocalanus spp.") %>%
  group_by(year) %>%
  summarise(mean_pseudo = mean(log_abundance, na.rm=T))%>%
  ggplot(aes(x=year, y=mean_pseudo)) +
  geom_point() +
  geom_line()

#plot calanus mean abundance: full timeseries 
dat %>%
  filter(taxa == "Calanus glacialis") %>%
  group_by(year) %>%
  summarise(mean_calanus = mean(log_abundance, na.rm=T))%>%
  ggplot(aes(x=year, y=mean_calanus)) +
  geom_point() +
  geom_line()

#Ideally I think we'd correct for a subset of stations sampled 
  #consistently across years but it looks like station id/lat/long
  #is never the same across years. Lets filter for stations 2005+ and run a GAM

########
# fit model to pseudocalanus to control for lat,long and julian day and 
# estimate an annual abundance

pseudo.dat <- dat %>%
  filter(taxa == "Pseudocalanus spp.", 
         year > 2004)

pseudo.mod <- gam(log_abundance ~ te(lon, lat) + s(julian, k = 5) + year_fac, 
                  data = pseudo.dat)

summary(pseudo.mod) 
gam.check(pseudo.mod) #Yikes...going to proceed for now but really need to look
  #at some different model structures 

#Back transform and extract year coefficient (2005 is our intercept)
c(coef(pseudo.mod)[1], coef(pseudo.mod)[1] + coef(pseudo.mod)[2:16]) -> est

year <- data.frame(year = c(2005:2019, 2021))
cbind(est,year) -> dat.2
as_tibble(dat.2) %>%
  rename(Pseudocalanus = est) -> pseudo.dat

#Plot 
ggplot(pseudo.dat, aes(year, Pseudocalanus)) +
  geom_point() +
  geom_line()

#########
# fit model to calanus to control for lat,long and julian day and 
# estimate an annual abundance

calanus.dat <- dat %>%
  filter(taxa == "Calanus glacialis", 
         year > 2004)

calanus.mod <- gam(log_abundance ~ te(lon, lat) + s(julian, k = 5) + year_fac, 
                  data = calanus.dat)

summary(calanus.mod) 
gam.check(calanus.mod) #much better than psuedocalanus model diagnostics 

#Back transform and extract year coefficient (2005 is our intercept)
c(coef(calanus.mod)[1], coef(calanus.mod)[1] + coef(calanus.mod)[2:16]) -> est.cal

year <- data.frame(year = c(2005:2019, 2021))
cbind(est.cal,year) -> dat.3
as_tibble(dat.3) %>%
  rename(Calanus_glacialis = est.cal) -> cal.dat

#Plot 
ggplot(cal.dat, aes(year, Calanus_glacialis)) +
  geom_point() +
  geom_line()

#Join datasets
pseudo.dat %>%
  full_join(cal.dat) %>%
  pivot_longer(c(1,3), names_to = "taxa", values_to = "log_abundance") -> zoo_dat

# and save
write.csv(zoo_dat, "./Data/summarized_zooplankton.csv", row.names = F)

