# create an index of annual Calanus and Pseudocalanus abundance
# controlling for spatial-temporal differences in sampling among years

library(tidyverse)
library(mgcv)

theme_set(theme_bw())

# load data
dat <- read.csv("./Data/Calanus_Pseudo_Combined.csv", row.names = 1)

##############################
# Data wrangling
dat <- dat %>%
  mutate(julian = 
           lubridate::yday(lubridate::parse_date_time(paste(DAY, MONTH, YEAR, sep = "-"), order = "dmy")),
         year_fac = as.factor(YEAR)) %>%
        mutate(abundance = exp(log10_EST_NUM_PERM3)) %>%
         rename(lat = LAT,
                long = LON,
                log_abundance = log10_EST_NUM_PERM3)

#Mean date sampled 
dat %>%
  group_by(YEAR) %>%
  summarise(mean_date = mean(julian, na.rm=T),
            min_date = min(julian, na.rm=T),
            max_date = max(julian, na.rm=T))
#plot date sampled 
dat.bs %>%
  ggplot(aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year) # big differences 

#Plot stations sampled 
dat %>%
  group_by(YEAR, lat, long) %>%
  distinct(STATION_NAME) %>%
  ggplot() +
  geom_point(aes(x =long, y = lat), size=.5) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~YEAR) #big differences 

#plot # stations sampled
dat %>%
  group_by(YEAR) %>%
  summarise(num_stations = n_distinct(STATION_NAME, na.rm=T)) %>%
  ggplot(aes(x=YEAR, y=num_stations)) +
  geom_point() +
  geom_line() #oh boy....

#plot psuedocalanus abundance
dat %>%
  filter(TAXA == "Pseudocalanus spp.") %>%
  group_by(YEAR) %>%
  summarise(total_pseudo = sum(abundance, na.rm=T))%>%
  mutate(total_log_pseudo = log(total_pseudo)) %>%
  ggplot(aes(x=YEAR, y=total_log_pseudo)) +
  geom_point() +
  geom_line()

#plot calanus abundance
dat %>%
  filter(TAXA == "Calanus glacialis") %>%
  group_by(YEAR) %>%
  summarise(total_cal = sum(abundance, na.rm=T))%>%
  mutate(total_log_cal = log(total_cal)) %>%
  ggplot(aes(x=YEAR, y=total_log_cal)) +
  geom_point() +
  geom_line()

#Number of stations is a big issue.....

#####
# fit model to pseudocalanus to control for lat,long and julian day and 
# estimate an annual abundance

pseudo.dat <- dat %>%
  filter(TAXA == "Pseudocalanus spp.")

pseudo.mod <- gam(log_abundance ~ te(long, lat) + s(julian, k = 5) + year_fac, 
                  data = pseudo.dat)

summary(pseudo.mod) 
gam.check(pseudo.mod) #Yikes...we really shouldn't be date correcting with these diagnostics

#Pausing here.....

#Plot 
ggplot(pred.dat, aes(year, log_abundance)) +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  facet_wrap(~group, scales = "free_y", ncol = 1)

# and save
write.csv(pred.dat, "./Data/summarized_zooplankton.csv", row.names = F)

