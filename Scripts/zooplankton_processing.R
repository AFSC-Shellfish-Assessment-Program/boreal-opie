# create an index of annual Calanus and Pseudocalanus abundance
# controlling for spatial-temporal differences in sampling among years

library(tidyverse)
library(mgcv)

theme_set(theme_bw())

# load data
dat <- read.csv("./Data/Calanus_Pseudo_Combined.csv", row.names = 1)

# add a Julian day column and turn year into a factor

dat <- dat %>%
  mutate(julian = 
           lubridate::yday(lubridate::parse_date_time(paste(DAY, MONTH, YEAR, sep = "-"), order = "dmy")),
         year_fac = as.factor(YEAR)) %>%
         rename(lat = LAT,
                long = LON,
                log_abundance = log10_EST_NUM_PERM3)


# fit model to each spp. group to control for lat,long and julian day and 
# estimate an annual abundance

pseudo.dat <- dat %>%
  filter(TAXA == "Pseudocalanus spp.")

pseudo.mod <- gam(log_abundance ~ te(long, lat) + s(julian, k = 5) + year_fac, 
                  data = pseudo.dat)

summary(pseudo.mod)
# plot(pseudo.mod)


##

calanus.dat <- dat %>%
  filter(TAXA == "Calanus glacialis")

calanus.mod <- gam(log_abundance ~ te(long, lat) + s(julian, k = 5) + year_fac, 
                  data = calanus.dat)

summary(calanus.mod)
# plot(calanus.mod)

# now create new data frame with average lat, long, and day and predict abundance for each year

new.pseudo.dat <- data.frame(julian = mean(pseudo.dat$julian),
                      lat = mean(pseudo.dat$lat),
                      long = mean(pseudo.dat$long),
                      year_fac = unique(pseudo.dat$year_fac))

pseudo.pred <- predict(pseudo.mod, newdata = new.pseudo.dat, se = T)


new.calanus.dat <- data.frame(julian = mean(calanus.dat$julian),
                             lat = mean(calanus.dat$lat),
                             long = mean(calanus.dat$long),
                             year_fac = unique(calanus.dat$year_fac))

calanus.pred <- predict(calanus.mod, newdata = new.calanus.dat, se = T)


pred.dat <- rbind(data.frame(year = new.pseudo.dat$year_fac,
                             log_abundance = pseudo.pred$fit,
                             LCI = pseudo.pred$fit - 1.96 * pseudo.pred$se.fit,
                             UCI = pseudo.pred$fit + 1.96 * pseudo.pred$se.fit,
                             group = "Pseudocalanus"),
                  data.frame(year = new.calanus.dat$year_fac,
                             log_abundance = calanus.pred$fit,
                             LCI = calanus.pred$fit - 1.96 * calanus.pred$se.fit,
                             UCI = calanus.pred$fit + 1.96 * calanus.pred$se.fit,
                             group = "Calanus glacialis")) %>%
  mutate(year = as.numeric(as.character(year)))

ggplot(pred.dat, aes(year, log_abundance)) +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  facet_wrap(~group, scales = "free_y", ncol = 1)

# and save
write.csv(pred.dat, "./Data/summarized_zooplankton.csv", row.names = F)

