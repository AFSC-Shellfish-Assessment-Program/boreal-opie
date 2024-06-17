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
    
    # nc <- nc_open("./Data/nceiErsstv5_b990_926e_15da.nc")
    
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

tiff("./figs/Extended_Data_Fig_1.tiff", width = 4, height = 3, units = "in", res = 300)

ggplot(plot.sst, aes(year, annual.sst)) +
  geom_point(size = 1) +
  geom_line() +
  geom_smooth(se = F) +
  scale_x_continuous(breaks = seq(1860, 2020, 20)) +
  theme(axis.title.x = element_blank()) +
  ylab("Mean SST (°C)")

dev.off()


