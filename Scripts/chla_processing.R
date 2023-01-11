# Goal: create an index of annual chla size fraction and total chla, 
# controlling for spatial-temporal differences in sampling among years

#Data source includes summer/fall depth integrated (top 50m) chla biomass  
  #from the Bering Sea fall BASIS survey 2003-2018, filtered by BSIERP regions 3&6
#NOTES: chl-a size fraction ratio is avg of >10um/avg total chla

library(tidyverse)
library(mgcv)
library(lubridate)
library(raster)
library(sf)
library(rgdal)

#Author: Erin Fedewa
#Last update: 1/10/2023

# load data
dat <- read.csv("./Data/chla_raw.csv")

#load BSIERP region shapefiles 
bsregions <- readOGR("./Data/BSIERP shapefiles", layer="BSIERP_regions_2012")

###############################################
#Spatially subset chla data for BSIERP regions 3 & 6:
dat %>%
  separate(gear_time, c("date", "time"), sep= "\\s", extra = "drop", fill = "right") %>%
  mutate(date = mdy(date),
        julian = yday(date),
        year = as.factor(year)) -> dat.new

#Assign coordinates to x and y variables
x <- dat.new$gear_long
y <- dat.new$gear_lat

#Define coordinate system
d <- data.frame(lon=x, lat=y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4269")

#Get coordinate system string
CRSargs <-bsregions@proj4string

#convert lat-lon file to this projection
CRS.new <- CRS(CRSargs@projargs)

#Convert data to map projection of shapefile
d.new <- spTransform(d, CRS.new)

#Assign each row of data to a specific polygon found in the polygons file
q <- over(d.new, as(bsregions, "SpatialPolygons"))

#Append this variable to chla data set
dat.new$bsregion <- NULL
dat.new$bsregion <- as.character(q)

#Filter by the SE Bering Sea regions 3 and 6
dat.new %>%
  filter(bsregion %in% c(3,6)) -> dat.bs

####################################
#Spatiotemporal data exploration:

#Stations sampled each yr 
dat.bs %>%
  group_by(year) %>%
  summarise(num_stations = length(unique(station_ID)))
#No data for 2013/2017/2019 

#Plot stations sampled 
dat.bs %>%
  group_by(year, gear_lat, gear_long) %>%
  distinct(station_ID) %>%
  ggplot() +
  geom_point(aes(x = gear_long, y = gear_lat), size=.5) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~year)

#Mean date sampled 
dat.bs %>%
  group_by(year) %>%
  summarise(mean_date = mean(julian, na.rm=T),
            min_date = min(julian, na.rm=T),
            max_date = max(julian, na.rm=T)) 
#plot
dat.bs %>%
  ggplot(aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year) # big differences 

#Avg chla ratio by year (no date correction)
dat.bs %>%
  group_by(year) %>%
  summarise(total_ratio = mean(avg_ratio, na.rm=T))%>%
  ggplot(aes(year, total_ratio, group=1)) +
  geom_point() +
  geom_line() 

#Avg total chla by year (no date correction)
dat.bs %>%
  group_by(year) %>%
  summarise(total_chla = mean(avg_totalchla, na.rm=T))%>%
  ggplot(aes(year, total_chla, group=1)) +
  geom_point() +
  geom_line() 

###############################################
# fit model to chl-a ratio to control for lat,long and julian day and 
# estimate annual chl-a size fraction 

chla.mod <- gam(avg_ratio ~ te(gear_long, gear_lat) + s(julian, k = 5) + year, 
                  data = dat.bs)

summary(chla.mod) 
plot(chla.mod)
gam.check(chla.mod) #not great...lots of positive residuals

#4th root transform response and re-run
dat.bs %>%
  filter(avg_ratio > 0) -> dat.new.fourth #Can't have zero observations here 
chla.mod.2 <- gam(avg_ratio^0.25 ~ te(gear_long, gear_lat) + s(julian, k = 5) + year, 
                data = dat.new.fourth)

summary(chla.mod.2)
plot(chla.mod.2)
gam.check(chla.mod.2) #much better

#Back transform and extract year coefficient (2003 is our intercept)
c(coef(chla.mod.2)[1], coef(chla.mod.2)[1] + coef(chla.mod.2)[2:14])^4 -> est

year <- data.frame(year = c(2003:2012, 2014:2016, 2018))
cbind(est,year) -> dat
as_tibble(dat) %>%
  rename(avg_chla_ratio = est) -> ratio.dat

#######################################
# fit model to total chl-a to control for lat,long and julian day and 
# estimate annual avg chl-a concentration (fall bloom size)

dat.bs %>%
  filter(avg_totalchla > 0) -> dat.new.total #Can't have zero observations here 
chla.mod.3 <- gam(avg_totalchla^0.25 ~ te(gear_long, gear_lat) + s(julian, k = 5) + year, 
                  data = dat.new.total)

summary(chla.mod.3)
plot(chla.mod.3)
gam.check(chla.mod.3) 

#Back transform and extract year coefficient (2003 is our intercept)
c(coef(chla.mod.3)[1], coef(chla.mod.3)[1] + coef(chla.mod.3)[2:14])^4 -> est

year <- data.frame(year = c(2003:2012, 2014:2016, 2018))
cbind(est,year) -> dat
as_tibble(dat) %>%
  rename(total_avg_chla = est) %>%
  full_join(ratio.dat) -> chla.final

#plot total chla
chla.final %>%
  ggplot(aes(year, total_avg_chla )) +
  geom_point() + 
  geom_line() #Trends are similar to raw data plot 

# Save product
write.csv(chla.final, "./Data/summarized_chla.csv", row.names=FALSE)



