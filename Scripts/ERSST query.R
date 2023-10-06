    # produce time series for annual, winter and spring/summer mean SST for E Bering
    
    library(tidyverse)
    library(ncdf4)
    library(zoo)
    library(maps)
    library(mapdata)
    library(chron)
    library(fields)
    library(oce)
    
    # set palettes
    new.col <- oceColorsPalette(64)
    cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    # set theme
    theme_set(theme_bw())
    
    
    ## load and process ERSST ------------------------
    
    # download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2022-03-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(66)][(180):1:(204)]", "~temp")
    # download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2022-12-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(66)][(180):1:(204)]", "~temp")
    
    # load and process SST data
    # nc <- nc_open("~temp")
    
    nc <- nc_open("./Data/nceiErsstv5_b990_926e_15da.nc")
    
    # full time series
    nc <- nc_open("./Data/nceiErsstv5_845b_f996_87bc.nc")
    # process
    
    ncvar_get(nc, "time")   # seconds since 1-1-1970
    raw <- ncvar_get(nc, "time")
    h <- raw/(24*60*60)
    d <- dates(h, origin = c(1,1,1970))
    m <- months(d)
    yr <- years(d)
    
    x <- ncvar_get(nc, "longitude")
    y <- ncvar_get(nc, "latitude")
    
    SST <- ncvar_get(nc, "sst", verbose = F)
    
    SST <- aperm(SST, 3:1)
    
    SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y, length(x))
    lon <- rep(x, each = length(y))
    dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
    
    # plot to check
    
    temp.mean <- colMeans(SST, na.rm=T)
    z <- t(matrix(temp.mean,length(y)))
    image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x, y, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")
    
    # extract area used in SST attribution!
    ebs.x <- c(183, 183, 203, 203, 191) 
    ebs.y <- c(53, 65, 65, 57.5, 53)
    
    polygon(ebs.x, ebs.y, border = "red", lwd = 2)


# define cells within polygon and plot to check

ebs.sst <- as.data.frame(SST)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(180, 205), ylim = c(50, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# create time series

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

# calculate annual weighted mean temp
monthly.sst <- apply(ebs.sst, 1, weighted.cell.mean)
annual.sst <- tapply(monthly.sst, yr, mean)

plot.sst <- data.frame(year = 1854:2022,
                       sst = annual.sst)

# save for plotting later
write.csv(plot.sst, "./output/annual_sst.csv", row.names = F)


ggplot(plot.sst, aes(year, annual.sst)) +
  geom_point(size = 1) +
  geom_line() +
  geom_smooth(se = F) +
  scale_x_continuous(breaks = seq(1860, 2020, 20)) +
  theme(axis.title.x = element_blank()) +
  ylab("Mean SST (°C)")

ggsave("./figs/annual_sst_ts.png", width = 4, height = 3, units = "in")


# create a function to compute monthly anomalies
monthly.anomalies <- function(x) tapply(x, m, mean) 

# short year and month vectors for 1950-2021
yr.1950.2022 <- yr[yr %in% 1950:2022]
m.1950.2022 <- m[yr %in% 1950:2022]

# and define winter year
winter.year <- if_else(m.1950.2022  %in% c("Nov", "Dec"), as.numeric(as.character(yr.1950.2022)) +1, as.numeric(as.character(yr.1950.2022)))
winter.year <- winter.year[m.1950.2022  %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# process
  # first, calculate monthly anomalies
  mu <- apply(ebs.sst, 2, monthly.anomalies)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location
  mu <- rbind(mu, mu[1:3,]) # add trailing months (Jan-Mar 2022)

  anom <- ebs.sst - mu   # compute matrix of anomalies
  
  # calculate weighted monthly means
  monthly.anom <- apply(anom, 1, weighted.cell.mean)

  
      
  # make data frame for plotting
  obs.plot <- data.frame(monthly.anom = monthly.anom,
                         date = lubridate::parse_date_time(x = paste(yr,as.numeric(m),"01"),orders="ymd",tz="America/Anchorage"))

  ggplot(obs.plot, aes(date, monthly.anom)) +
    geom_line(lwd = 0.02) +
    ylab("Monthly anomaly") +
    theme(axis.title.x = element_blank())
  
  # looks good

  
  ## now create time series of annual means
  ## unsmoothed, and two- and three-year running means
  
  time.series <- anomaly.time.series <- data.frame()
 
  # calculate monthly mean temp weighted by area for 1950-2022
  monthly.sst <- apply(ebs.sst[yr %in% 1950:2022,], 1, weighted.cell.mean) 
  
  # calculate annual means
  temp.annual <- tapply(monthly.sst, yr[yr %in% 1950:2022], mean) 
  
  # drop trailing year (2022 is incomplete)
  trailing <- max(names(temp.annual))
  temp.annual[names(temp.annual) %in% trailing] <- NA
  
  temp.2yr <- rollmean(temp.annual, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
  
  temp.3yr <- rollmean(temp.annual, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry

  # calculate winter means
  temp.winter.monthly.sst <- monthly.sst[m.1950.2022 %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

  temp.winter <- tapply(temp.winter.monthly.sst, winter.year, mean)
  
  # change leading year to NA b/c incomplete
  leading <- min(names(temp.winter))

  temp.winter[names(temp.winter) %in% leading] <- NA
  
  temp.winter.2yr <- rollmean(temp.winter, 2, fill = NA, align = "left")
  
  temp.winter.3yr <- rollmean(temp.winter, 3, fill = NA, align = "center")
  
  
  # combine into data frame of time series
  time.series <- rbind(time.series,
                            data.frame(year = 1950:2022,
                                       annual.unsmoothed = temp.annual[names(temp.annual) %in% 1950:2022],
                                       annual.two.yr.running.mean = temp.2yr[names(temp.annual) %in% 1950:2022],
                                       annual.three.yr.running.mean = temp.3yr[names(temp.annual) %in% 1950:2022],
                                       winter.unsmoothed = temp.winter[names(temp.winter) %in% 1950:2022],
                                       winter.two.yr.running.mean = temp.winter.2yr[names(temp.winter) %in% 1950:2022],
                                       winter.three.yr.running.mean = temp.winter.3yr[names(temp.winter) %in% 1950:2022]))

  ## now calculate the data as anomalies wrt 1950-1999
  # calculate annual anomalies
  annual.climatology.mean <- mean(temp.annual[names(temp.annual) %in% 1950:1999])
    
  annual.climatology.sd <- sd(temp.annual[names(temp.annual) %in% 1950:1999])
  
  temp.annual.anom <- (temp.annual - annual.climatology.mean) / annual.climatology.sd
  
  temp.anom.2yr <- rollmean(temp.annual.anom, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
  
  temp.anom.3yr <- rollmean(temp.annual.anom, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
  
  # calculate winter anomalies
  winter.climatology.mean <- mean(temp.winter[names(temp.winter) %in% 1950:1999], na.rm = T)
  
  winter.climatology.sd <- sd(temp.winter[names(temp.winter) %in% 1950:1999], na.rm = T)
  
  temp.winter.anom <- (temp.winter - winter.climatology.mean) / winter.climatology.sd
  
  temp.winter.anom.2yr <- rollmean(temp.winter.anom, 2, fill = NA, align = "left")
  
  temp.winter.anom.3yr <- rollmean(temp.winter.anom, 3, fill = NA, align = "center")
  
  # combine into data frame of time series by region
  anomaly.time.series <- rbind(anomaly.time.series,
                            data.frame(year = 1950:2021,
                                       annual.anomaly.unsmoothed = temp.annual.anom[names(temp.annual.anom) %in% 1950:2021],
                                       annual.anomaly.two.yr.running.mean = temp.anom.2yr[names(temp.annual.anom) %in% 1950:2021],
                                       annual.anomaly.three.yr.running.mean = temp.anom.3yr[names(temp.annual.anom) %in% 1950:2021],
                                       winter.anomaly.unsmoothed = temp.winter.anom[names(temp.winter.anom) %in% 1950:2021],
                                       winter.anomaly.two.yr.running.mean = temp.winter.anom.2yr[names(temp.winter.anom) %in% 1950:2021],
                                       winter.anomaly.three.yr.running.mean = temp.winter.anom.3yr[names(temp.winter.anom) %in% 1950:2021]))
  
    



# plot

plot.dat <- time.series %>%
  select(year,
         annual.unsmoothed,
         winter.unsmoothed) %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scale = "free_y", ncol = 1) 
  
  
# and save 
write.csv(time.series, "./Data/ersst_time_series.csv", row.names = F)
write.csv(anomaly.time.series, "./Data/ersst_anomaly_time_series.csv", row.names = F)