    # produce ice cover time series for Bering (55-64N, 180-160W)
    
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
    
    
    ## load and process ------------------------
    
    # note that there are separate time series for 1950-1978 and 1979-present
    
    nc1 <- nc_open("./Data/ERA5_ice_1950-1978.nc")
    
    # process
    
    ncvar_get(nc1, "time")   # hours since 1-1-1900
    raw <- ncvar_get(nc1, "time")
    h <- raw/24
    d1 <- dates(h, origin = c(1,1,1900))
    m1 <- months(d1)
    yr1 <- years(d1)
    
    x1 <- ncvar_get(nc1, "longitude")
    y1 <- ncvar_get(nc1, "latitude")
    
    ice1 <- ncvar_get(nc1, "siconc", verbose = F)
    dim(ice1) # 87 long, 37 lat, 203 months
    
    # reverse lat for plotting
    ice1 <- ice1[,37:1,]
    
    # reverse y too
    y1 <- rev(y1)
    
    ice1 <- aperm(ice1, 3:1)
    
    ice1 <- matrix(ice1, nrow=dim(ice1)[1], ncol=prod(dim(ice1)[2:3]))
    
    # plot to check
    
    ice.mean <- colMeans(ice1, na.rm=T)
    z <- t(matrix(ice.mean,length(y1)))
    image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x1, y1, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  
    
    # now the second time series
    
    nc2 <- nc_open("./Data/ERA5_ice_1979-2022.nc")
    
    # process
    
    ncvar_get(nc2, "time")   # hours since 1-1-1900
    raw <- ncvar_get(nc2, "time")
    h <- raw/24
    d2 <- dates(h, origin = c(1,1,1900))
    m2 <- months(d2)
    yr2 <- years(d2)
    
    x2 <- ncvar_get(nc2, "longitude")
    y2 <- ncvar_get(nc2, "latitude")
    
    # expver <-  ncvar_get(nc2, "expver", verbose = F)
    # expver # 1 and 5??
    
   
    ice2 <- ncvar_get(nc2, "siconc", verbose = F)
    dim(ice2) # 87 long, 37 lat, 2 expver, 203 months
    
    # expver1 - this is ERA5
    
    # ice2 <- ice2[,,1,]
    
    # reverse lat for plotting
    ice2 <- ice2[,37:1,]
    
    # reverse y too
    y2 <- rev(y2)
    
    ice2 <- aperm(ice2, 3:1)
    
    ice2 <- matrix(ice2, nrow=dim(ice2)[1], ncol=prod(dim(ice2)[2:3]))
    
    
    # plot to check
    
    ice.mean <- colMeans(ice2, na.rm=T)
    z <- t(matrix(ice.mean,length(y2)))
    image.plot(x2,y2,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x2, y2, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  
    
    
    # check dimensions
    identical(x1, x2)
    identical(y1,y2)
    

    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y1, length(x1))
    lon <- rep(x1, each = length(y1))
    
    ice <- rbind(ice1, ice2)
    
    # drop E of 165 and N of 63
    drop <- lon > -165 | lat > 63
    ice[,drop] <- NA
    
    # plot to check
    ice.mean <- colMeans(ice, na.rm=T)
    z <- t(matrix(ice.mean,length(y1)))
    image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x1, y1, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3") # perfecto
    
    dimnames(ice) <- list(as.character(c(d1, d2)), paste("N", lat, "E", lon, sep=""))
    
    f <- function(x) colMeans(x, na.rm = T)
    
    m <- c(as.character(m1), as.character(m2))
    
    yr <- c(as.numeric(as.character(yr1)), as.numeric(as.character(yr2)))
    
    means <- data.frame(month = m,
                        year = as.numeric(as.character(yr)),
                        ice = rowMeans(ice, na.rm = T)) 
    
    
    ggplot(means, aes(year, ice, color = month)) +
      geom_line()
    
    # drop Oct - Dec
    means <- means %>%
      filter(!month %in% c("Oct", "Nov", "Dec"))
    
    
    # pivot wider
    means <- means %>% 
      pivot_wider(values_from = ice, names_from = month) %>%
      filter(year %in% 1953:2021) # fixed values through 1952, 2022 incomplete!
    
    means[,2:5] <- apply(means[,2:5], 2, scale)
    
    plot <- means %>%
      pivot_longer(cols = -year)
    
    ggplot(plot, aes(year, value, color = name)) +
      geom_line()
    
    
    # generate Jan-Feb and Mar-Apr means
    means$JanFeb_ice <- apply(means[,2:3], 1, mean)
    means$MarApr_ice <- apply(means[,4:5], 1, mean)
    
    # clean up
    means <- means %>%
      select(year, JanFeb_ice, MarApr_ice)
    
    # save ERA5
    write.csv(means, "./Data/ice.csv", row.names = F)
    
    plot <- means %>%
      pivot_longer(cols = -year)
    
    ggplot(plot, aes(year, value, color = name)) +
      geom_line() +
      geom_point()

    
# compare with ORAS5
oras <- read.csv("./Data/ORAS5_summarized_ice.csv")
oras$data = "ORAS5"
          
means$data <- "ERA5"

oras <- oras %>%
  pivot_longer(cols = c(-year, -data))


means <- means %>%
  pivot_longer(cols = c(-year, -data))


compare <- rbind(oras, filter(means, year %in% oras$year))

compare <- compare %>%
  pivot_wider(names_from = c(data, name), values_from = value) %>%
  mutate(era = if_else(year < 1979, "1972-1978", "1979-2018"))

ggplot(compare, aes(ORAS5_JanFeb_ice, ERA5_JanFeb_ice)) +
  geom_point() +
  facet_wrap(~era)

ggplot(compare, aes(ORAS5_MarApr_ice, ERA5_MarApr_ice)) +
  geom_point() +
  facet_wrap(~era)

compare <- rbind(oras, filter(means, year %in% oras$year))
          
ggplot(compare, aes(year, value, color = data)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, ncol = 1)
         
# pretty similar!!


