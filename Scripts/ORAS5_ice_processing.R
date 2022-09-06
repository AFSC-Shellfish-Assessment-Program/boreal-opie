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


# load file

# URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_eed0_4cbd_14b1.nc?ileadfra[(1972-01-15):1:(2018-12-15T00:00:00Z)][(55.5):1:(62.5)][(180):1:(195)]"


nc <- nc_open("./Data/hawaii_soest_eed0_4cbd_14b1_79fd_149d_f84a.nc")

nc


x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

ncvar_get(nc, "time")   # second since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

ice <- ncvar_get(nc, "ileadfra", verbose = F)
dim(ice) 


# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
ice <- aperm(ice, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
ice <- matrix(ice, nrow=dim(ice)[1], ncol=prod(dim(ice)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(ice) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# check for values < 0!
check <- ice < 0
sum(check, na.rm = T) 0!


ice.mean <- colMeans(ice, na.rm=T)
z <- t(matrix(ice.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x, y, z, add=T)
map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  

# save Jan-Feb and Mar-Apr means for comparison with ERA5

f <- function(x) colMeans(x, na.rm = T)

m <- as.character(m)

yr <- as.numeric(as.character(yr))

means <- data.frame(month = m,
                    year = yr,
                    ice = rowMeans(ice, na.rm = T)) 


ggplot(means, aes(year, ice, color = month)) +
  geom_line()

# keep Jan - Apr
means <- means %>%
  filter(month %in% c("Jan", "Feb", "Mar", "Apr"))


# pivot wider
means <- means %>% 
  pivot_wider(values_from = ice, names_from = month)

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

plot <- means %>%
  pivot_longer(cols = -year)

ggplot(plot, aes(year, value, color = name)) +
  geom_line() +
  geom_point()


# save 
write.csv(means, "./Data/ORAS5_summarized_ice.csv", row.names = F)
