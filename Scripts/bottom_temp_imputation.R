# estimate missing bottom temp values for alternate borealization effects models
# using ice time series to predict bottom temp

library(tidyverse)
library(corrplot)

theme_set(theme_bw())

# Updated 1/10/23 with 2022 data by EJF


## data processing -------------------------------------

#sea ice
d1 <- read.csv("./Data/ice.csv")

d1 <- d1 %>%
  rename(`Jan-Feb ice cover` = JanFeb_ice,
         `Mar-Apr ice cover` = MarApr_ice) %>%
  filter(year >= 1972) %>%
  pivot_longer(cols = -year)

#bottom temp
d2 <- read.csv("./Data/date_corrected_bottom_temp.csv")

d2 <- d2 %>%
  rename(value = bottom.temp) %>%
  mutate(name = 'Bottom temperature')


dat <- rbind(d1, d2) %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  dplyr::select(-year)

# clean up dat names for imputation
names(dat) <- c("Jan.Feb.ice", "Mar.Apr.ice", "Bottom.temp")

# examine correlations
r <- cor(as.matrix(dat), use = "p")
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 25 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE


imp <- mice::mice(data = dat, method = "pmm", m=100, predictorMatrix = pred)

saveRDS(imp, "./output/bottom_temp_imputation.RDS")
imp <- readRDS("./output/bottom_temp_imputation.RDS")

check <- is.na(imp)
sum(check)


# summarize (mean and SD) imputed dat

imputed.dat <-  data.frame()


for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp, action = i)
  

  imputed.dat <- rbind(imputed.dat,
                       data.frame(imputation = i,
                                  year = c(1972:2022),
                                  Bottom_temp = temp$Bottom.temp))
  
}



summary.dat <- data.frame(year = 1972:2022,
                          Bottom_temp = tapply(imputed.dat$Bottom_temp, imputed.dat$year, mean, na.rm = T),
                          SD = tapply(imputed.dat$Bottom_temp, imputed.dat$year, sd, na.rm = T))


# replace SD = 0 with NA
change <- summary.dat$SD == 0
summary.dat$SD[change] <- NA

ggplot(summary.dat, aes(year, Bottom_temp)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = Bottom_temp - 2*SD,
                    ymax = Bottom_temp + 2*SD))

# and save
write.csv(summary.dat, "./output/imputed_bottom_temp.csv", row.names = F)
  